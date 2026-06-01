#!/usr/bin/env bash
# build.sh
# Builds the mm_sim executable.
# Must run AFTER sourcing setup_lxplus.sh (or setup_local.sh).
#
# Usage:
#   bash scripts/build.sh [clean]
#
# Options:
#   clean  -- remove build dir first

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="${SRC_DIR}/build"

if [[ "$1" == "clean" ]]; then
    echo "==> Cleaning build dir: $BUILD_DIR"
    rm -rf "$BUILD_DIR"
fi

echo "==> Source dir : $SRC_DIR"
echo "==> Build dir  : $BUILD_DIR"

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Configure
cmake "$SRC_DIR" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${SRC_DIR}/install"

# Build (use all available cores, capped at 8 to be polite on lxplus login nodes)
NCPU=$(nproc 2>/dev/null || echo 4)
NCPU=$(( NCPU > 8 ? 8 : NCPU ))
echo "==> Building with $NCPU cores..."
make -j"$NCPU"

echo ""
echo "==> Build successful!"
echo "    Executable: ${BUILD_DIR}/mx17_full_sim"
echo ""
echo "==> Quick test (1000 X17 pair events, ArCF4):"
echo "    ${BUILD_DIR}/mx17_full_sim -n 1000 -o /tmp/test_x17"
