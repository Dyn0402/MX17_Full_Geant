# The E0 (0⁺ → 0⁺) pair branch in ⁴He — what Alberto is asking, in plain terms

**For:** Dylan (intuition-building) · **Re:** question 2
**Date:** 2026-06-14 · **Status:** explainer draft — physics first, decisions later

The goal here is just to understand the question clearly before we decide
whether/how to put it in the generator. No jargon dumps; where a term is
unavoidable I define it.

---

## 1. The level scheme in one picture

The states involved sit just above the `n + ³He` and `p + ³H` thresholds in
⁴He (Alberto's figure):

```
   E* (MeV)
   23.30  2⁻
   21.84  2⁻
   21.21  0⁻   <-- n + ³He threshold is at 20.58
   20.21  0⁺   <-- the "subthreshold" 0⁺ state
   ─────────────────────────────────────
    0.00  0⁺   ground state of ⁴He
```

Both the 20.21 MeV state and the ⁴He ground state are **0⁺** — spin 0, positive
parity. The transition between them is what Alberto means by "the E0
transition." "Subthreshold" just means the 0⁺ level lies *below* the 20.58 MeV
`n + ³He` threshold, so when a slow neutron is captured the system is fed mainly
into the states at/above threshold and only leaks down into the 0⁺ via its tail
— hence "small branch."

## 2. Why 0⁺ → 0⁺ is special: no real photon allowed

A photon carries away at least one unit of angular momentum (it has no spin-0
mode). To emit a single real γ you need the initial and final spins to differ
by at least 1, **or** at least allow a 0→0 with one unit — and a `0 → 0`
transition can do *neither*. The rule is hard: **a 0⁺ → 0⁺ transition cannot
emit a single real gamma ray at all.** (This is the same selection rule that
makes the famous 0⁺→0⁺ transitions elsewhere "γ-forbidden.")

So how does the nucleus get rid of its 20.21 MeV? Through processes that don't
require a real photon:

- **Internal conversion** — hand the energy to an atomic electron. Negligible
  here (helium, almost no bound electrons, and a gas).
- **Internal pair creation (E0 pair emission)** — if the available energy
  exceeds `2mₑc² = 1.022 MeV` (here it's 20.21 MeV, hugely above), the nucleus
  can create an `e⁺e⁻` pair *directly*, via a virtual photon. **This is the
  only practical de-excitation path for the 0⁺ state, and it produces exactly
  the `e⁺e⁻` pairs our experiment detects.**

This is the "E0" in the question: **E** = electric, **0** = monopole
(zero units of angular momentum). It is a pure pair (+ conversion) transition —
there is no γ line to go with it.

## 3. How this differs from the IPC we already simulate

We currently generate `e⁺e⁻` pairs from the **1⁻ continuum** near 20.58 MeV
("E1", electric dipole). That is *internal pair conversion of an allowed γ
transition*: the 1⁻ → 0⁺ decay can emit a real photon, and a small fraction of
the time that photon is "internally converted" into a pair instead. Two things
are different for E0:

| | E1 (what we have) | E0 (the missing branch) |
|---|---|---|
| Parent state | 1⁻ continuum, ~20.58 MeV | 0⁺ subthreshold, ~20.21 MeV |
| Real γ competing? | **yes** — pairs are a small IPC tail of a γ line | **no** — pairs are the whole transition |
| Transition energy | 20.58 MeV (our `transition_energy_MeV`) | **20.21 MeV** (different) |
| Pair angular / energy sharing | dipole (E1) pattern | **monopole (E0) pattern — different** |

The practical upshots for us:

1. **Different transition energy** (20.21 vs 20.58 MeV) → slightly different
   total pair energy and boost.
2. **Different pair kinematics.** The opening-angle and the energy-sharing
   between e⁺ and e⁻ follow the E0 monopole distribution, which is *not* the
   same shape as the E1/IPC distribution our generator samples. Since our whole
   measurement is an opening-angle / invariant-mass analysis, the *shape*
   matters, not just the rate.
3. **No accompanying γ**, which subtly changes the in-detector background
   picture for these events relative to an E1 transition of the same energy.

## 4. Why Alberto flags it as non-negligible

His estimate is that the 0⁺ E0 pair yield is *of the same order* as the E1 pair
yield from the 1⁻ continuum — i.e. **comparable to a background component we
already include**, even though both are tiny next to the dominant (n,p)t
channel. If that holds, leaving E0 out means our simulated `e⁺e⁻` background
(and therefore the shape we fit the X17 signal against) is missing a piece of
the same size as a piece we kept. That is a modelling completeness issue, not a
rate-table issue.

## 5. What "including it" would mean for the sim (not deciding yet)

For when we come back to this — the generator already has the machinery
(`X17PrimaryGenerator::GenerateIPC` builds a pair from a parent of a given mass
and transition energy). Adding E0 would be a third pair sub-generator mixed in
by a branching weight, needing three physics inputs we'd have to pin down:

1. the **E0 branching weight** relative to the E1/IPC component (Alberto's
   "same order of magnitude" — we'd want a number);
2. the **transition energy** (20.21 MeV);
3. the **E0 pair energy-sharing / opening-angle distribution** (the monopole
   pair-conversion shape — this is the part that needs a real reference, e.g. a
   pair internal-conversion coefficient / spectrum for E0).

None of that is hard to wire in once we agree on the inputs. The honest gap
right now is physics inputs (especially #1 and #3), not code.

---

## Questions to resolve with Alberto (once the physics is clear)

- Where does the "same order of magnitude as E1" estimate come from — a
  measured/evaluated E0 strength for this 0⁺ state, or a theoretical estimate?
  That sets the branching weight #1.
- Do we have a reference for the E0 pair energy-sharing and angular
  distribution at ~20 MeV, or should we use the standard monopole
  pair-conversion formula? That sets #3.
- Is the 0⁻ (21.21 MeV) state worth the same treatment, or is it negligible?

Once these are answered, the generator change is small and we can re-run the
IPC background pool with the E0 component included and see how much the
opening-angle background shape moves.
