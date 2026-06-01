#!/usr/bin/env bash
# setup_lxplus.sh
# Sources Geant4 11.2 + ROOT 6.30 from CVMFS on lxplus9 (AlmaLinux 9).
# Run this ONCE in your login shell before building or submitting jobs.
#
# Usage:
#   source scripts/setup_lxplus.sh
#
# If you later want a different G4 version, check available builds:
#   ls /cvmfs/geant4.cern.ch/geant4/
#   ls /cvmfs/sft.cern.ch/lcg/releases/Geant4/

# ---- Sanity check: are we on lxplus? ----
if [[ ! -d /cvmfs/geant4.cern.ch ]]; then
    echo "ERROR: /cvmfs/geant4.cern.ch not found. Are you on lxplus with CVMFS mounted?"
    return 1 2>/dev/null || exit 1
fi

# ---- Geant4 11.2.p02 on AlmaLinux9 / gcc13 ----
# This is the latest stable build available on lxplus9 as of 2025.
# Falls back to 11.1 if 11.2 isn't present yet on your node.
G4_BASE="/cvmfs/geant4.cern.ch/geant4"
G4_VER=""

for try_ver in "11.2.p02" "11.2.p01" "11.2" "11.1.p02" "11.1"; do
    candidate="${G4_BASE}/${try_ver}"
    if [[ -d "$candidate" ]]; then
        G4_VER="$try_ver"
        G4_DIR="$candidate"
        break
    fi
done

if [[ -z "$G4_VER" ]]; then
    echo "ERROR: Could not find any Geant4 installation under $G4_BASE"
    echo "Available:"
    ls "$G4_BASE" 2>/dev/null || echo "(none)"
    return 1 2>/dev/null || exit 1
fi

# Pick the right architecture string (try gcc13 for lxplus9, fall back to gcc11)
G4_ARCH=""
for try_arch in \
    "x86_64-el9-gcc13-optdeb-MT" \
    "x86_64-centos9-gcc13-optdeb-MT" \
    "x86_64-el9-gcc11-optdeb-MT" \
    "x86_64-centos9-gcc11-optdeb-MT"; do
    if [[ -d "${G4_DIR}/${try_arch}" ]]; then
        G4_ARCH="$try_arch"
        break
    fi
done

if [[ -z "$G4_ARCH" ]]; then
    echo "ERROR: No matching architecture build found in ${G4_DIR}/"
    echo "Available builds:"
    ls "${G4_DIR}/" 2>/dev/null
    return 1 2>/dev/null || exit 1
fi

G4_INSTALL="${G4_DIR}/${G4_ARCH}"
echo "==> Sourcing Geant4 ${G4_VER} [${G4_ARCH}]"
source "${G4_INSTALL}/bin/geant4.sh"
source "${G4_INSTALL}/CMake-setup.sh" 2>/dev/null || true

# ---- ROOT 6.30 via LCG (compatible with gcc13 / el9) ----
# LCG_106 ships ROOT 6.30.x with a gcc13 el9 build.
LCG_ROOT=""
for try_lcg in "LCG_106" "LCG_105" "LCG_104"; do
    candidate="/cvmfs/sft.cern.ch/lcg/views/${try_lcg}/x86_64-el9-gcc13-opt/setup.sh"
    if [[ -f "$candidate" ]]; then
        LCG_ROOT="$candidate"
        LCG_TAG="$try_lcg"
        break
    fi
done

if [[ -n "$LCG_ROOT" ]]; then
    echo "==> Sourcing ROOT via ${LCG_TAG}"
    source "$LCG_ROOT"
else
    echo "WARNING: Could not find LCG ROOT setup. Will build without ROOT (CSV output only)."
    echo "Check: ls /cvmfs/sft.cern.ch/lcg/views/"
fi

# ---- CMake (LCG ships a recent one) ----
CMAKE_BIN="/cvmfs/sft.cern.ch/lcg/contrib/CMake/3.26.2/Linux-x86_64/bin"
if [[ -d "$CMAKE_BIN" ]]; then
    export PATH="${CMAKE_BIN}:$PATH"
fi

echo ""
echo "==> Environment ready:"
echo "    Geant4  : $G4_VER  ($G4_INSTALL)"
which cmake && cmake --version | head -1
echo ""
echo "==> Now run: bash scripts/build.sh"
