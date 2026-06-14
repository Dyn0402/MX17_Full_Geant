# Run guide: neutron-path event displays (lxplus)

The per-step trajectory dump is now implemented in the simulation
(`--trajdump`), and the plotting script is written and smoke-tested locally.
What remains needs the actual sim, which only runs on lxplus (Geant4 is not
installed on the workstation). Steps below.

## 1. Build (lxplus, with the campaign's Geant4 + ROOT)

```bash
source scripts/setup_lxplus.sh      # same env as the neutron campaign
cmake -S . -B build && cmake --build build -j
```

> Note: the C++ trajectory dump (`SimConfig::trajDump`,
> `SteppingAction::DumpTrajectoryStep`, the `--trajdump` flag) was written
> here but **not compile-tested locally** — no Geant4 on the workstation.
> Watch the first build for typos.

## 2. Generate trajectories at a few energies

`--trajdump N` writes `<out>_traj[_t<tid>].csv` for the first N events. Run
single-thread (`-t 1`) and pin a narrow energy window per run so each file is
~mono-energetic. ~15–20 events each is plenty; thermal neutrons take many
elastic steps, so keep N small.

```bash
N=15
for tag in thermal kev mev; do
  case $tag in
    thermal) EMIN=0.02   EMAX=0.05   ;;   # ~25 meV
    kev)     EMIN=1e3    EMAX=3e3    ;;   # ~1 keV
    mev)     EMIN=1e6    EMAX=3e6    ;;   # ~1 MeV
  esac
  ./build/mx17_full_sim -t 1 -n 4000 --trajdump $N \
     --neutron data/fluxEAR2-Ph3_in_different_units.root data/lamda2DvsEn_EAR2.root \
     --emin $EMIN --emax $EMAX -o traj_$tag
done
```

(`-n` is larger than `N` because not every beam neutron reaches the gas; the
dump stops after the first `N` events regardless.) Copy the three
`traj_*_traj*.csv` back to the workstation.

## 3. Make the figures (locally)

```bash
PY=/home/dylan/PycharmProjects/nTof_x17/venv/bin/python
$PY scripts/make_neutron_event_display.py \
    traj_thermal_traj*.csv traj_kev_traj*.csv traj_mev_traj*.csv \
    -o docs/neutron_histories/figs/event_display.pdf --max-events 3
```

Each panel draws one neutron history over the real capsule cross-section,
coloured by **fraction of incident KE remaining** (so the slowing-down shows
the same way at every energy), with charged secondaries and the terminal
vertex marked. Title shows incident energy, step count and fate
(`neutronInelastic` = (n,p)t, `nCapture` = radiative).

## 4. Optional: the E-at-capture vs E-incident scatter

The same CSVs (or, with higher stats, a normal `--neutron` run reading the
terminal `cap_*`/`neutron_E_eV` branches) give the quantitative companion:
KE at the terminal step vs incident energy, showing how much moderation
precedes the excitation. Ping me with the CSVs and I'll add it +
fold both figures into `neutron_histories_note.md`.
