# Performance Pitfall: uproot AwkwardExtensionArray with `pd.concat`

## Summary

Reading a ROOT file with `uproot` using `library="pd"` returns string branches
(`Char[N]/C`) as `AwkwardExtensionArray` objects.  Passing these columns
through `pd.concat` triggers a Python-level element-by-element fallback that is
**~8000× slower** than the equivalent operation on a regular pandas object array.

Observed impact in a real analysis: **74 minutes** per file instead of **~50 seconds**.

---

## Root cause

uproot interprets fixed-length C-string branches (`Char[32]/C`, `Char[64]/C`,
etc.) using its `AsStrings()` reader.  In pandas mode (`library="pd"`) this
produces an `AwkwardExtensionArray` — a custom `pandas.api.extensions.ExtensionArray`
subtype from the `awkward-array` library.

Pandas dispatches `pd.concat` for `ExtensionArray` columns through
`ExtensionArray._concat_same_type()`.  When two `AwkwardExtensionArray` instances
are concatenated, this method iterates the arrays at Python speed rather than
delegating to a NumPy kernel.  For N rows the cost is O(N) Python object creations
and function calls.

Benchmark (300 000-row chunk, two string columns):

| operation | time |
|-----------|------|
| `pd.concat([awk_df, awk_df])` | **87 s** |
| `pd.concat([obj_df, obj_df])` (after fix) | **11 ms** |
| `np.asarray(series)` conversion | 300 ms (one-time, per column per chunk) |

---

## How to detect this in your own code

### Sign 1 — column dtype

```python
import uproot

f = uproot.open("file.root")
df = f["Tree"].arrays(["some_string_branch"], library="pd")
print(df["some_string_branch"].dtype)   # prints: awkward
print(type(df["some_string_branch"].iloc[0]))  # prints: <class 'str'>
```

The dtype name contains `"awkward"`.  The elements are already Python `str`
objects (uproot already decoded the null-terminated C-string), but they are
wrapped in an awkward container that pandas cannot efficiently handle.

### Sign 2 — sudden slowdown on the second loop iteration

If you are iterating a large tree in chunks and buffering a "leftover" partial
event from one chunk to the next with `pd.concat`, the **first chunk is fast**
and every subsequent chunk is ~100–10 000× slower:

```
chunk 0:  0.55 s   (no leftover to concat)
chunk 1: 89.74 s   (concat awk_leftover + awk_new_chunk)
chunk 2: 90.24 s   (same)
...
```

### Sign 3 — profiling shows `pd.concat` dominates

```python
import cProfile
cProfile.run("pd.concat([df1, df2], ignore_index=True)")
# Shows millions of calls into awkward internals
```

### Quick diagnostic

```python
def has_awkward_cols(df):
    return [c for c in df.columns if "awkward" in str(df[c].dtype)]

chunk = next(iter(tree.iterate(step_size=100_000, library="pd")))
print(has_awkward_cols(chunk))   # non-empty → you are affected
```

---

## The fix

Convert `AwkwardExtensionArray` columns to regular numpy object arrays
**immediately after reading**, before any `pd.concat`, filter, or groupby
that touches those columns.

```python
# In the chunked-read loop, add these lines right after receiving each chunk:
for chunk in tree.iterate(columns, step_size=chunk_size, library="pd"):
    # Fix: materialise awkward string columns into plain numpy object arrays.
    # Without this, pd.concat([leftover, chunk]) takes ~90 s per chunk
    # instead of <1 ms because pandas has no fast path for AwkwardExtensionArray.
    for col in ["detType", "particle"]:      # replace with your string column names
        if "awkward" in str(chunk[col].dtype):
            chunk[col] = np.asarray(chunk[col])
    ...
```

`np.asarray()` materialises the elements into a regular `numpy` object array of
Python `str` objects in ~300 ms per 300k rows — a one-time cost that pays for
itself on the very first concat.

### Alternative: read as numpy and convert manually

```python
arrays = tree.arrays(string_cols, library="np")   # returns numpy fixed-length bytes
df[col] = arrays[col].astype(str)                  # fast vectorised conversion
```

### Alternative: avoid string columns entirely

Map the string values to integer codes once at the start of the script and work
with integers throughout.  Integers have no AwkwardExtensionArray overhead:

```python
DET_CODE  = {name: i for i, name in enumerate(["DriftGas", "PlasticScint", ...])}
PART_CODE = {"e-": 0, "e+": 1}

chunk["det_id"]  = np.asarray(chunk["detType"]).astype(object)
chunk["det_id"]  = chunk["det_id"].map(DET_CODE).fillna(-1).astype(int)
```

---

## Secondary bug: NaN from left-merge corrupts boolean columns

When a `pd.merge(..., how="left")` joins an EventTree (all events) with a
hit-summary derived from a subset of HitTree chunks, events that appear in
the EventTree but not in the hit-summary get `NaN` for all hit-derived bool
columns.

In pandas, a `bool` column that receives even one `NaN` is silently promoted
to `float64` (numpy booleans cannot represent NaN).  Any subsequent use of
`.values` for numpy fancy indexing then fails with:

```
IndexError: arrays used as indices must be of integer (or boolean) type
```

### Fix

Call `.fillna(False).astype(bool)` on every hit-derived boolean column before
using it as a mask or index:

```python
# Instead of:
mask = merged.get("double_trig", pd.Series(False, index=merged.index))
arr[mask.values]   # crashes when NaN present

# Do:
mask = (merged.get("double_trig", pd.Series(False, index=merged.index))
        .fillna(False).astype(bool))
arr[mask.values]   # safe
```

---

## Affected versions

Reproduced with:

- `uproot` 5.x (uses `awkward` 2.x backend for `AsStrings()`)
- `pandas` 1.x / 2.x
- `awkward` 2.x

The issue is present any time `uproot` returns a column as `AwkwardExtensionArray`
and that column is subsequently passed through `pd.concat`.  String-valued ROOT
branches (`Char[N]/C`, `TString`, `std::string`) are the most common trigger.
Numeric branches are unaffected.

---

## Checklist for uproot + pandas code

- [ ] Print `df.dtypes` for the first chunk and look for `awkward` dtype entries
- [ ] If any exist, add `chunk[col] = np.asarray(chunk[col])` before any `pd.concat`
- [ ] After a `left`-merge, treat all hit-derived bool columns as potentially `float64`; use `.fillna(False).astype(bool)` before boolean indexing
- [ ] Profile with `cProfile` or manual `time.perf_counter()` checkpoints around each `pd.concat` call if unexpected slowdowns appear
