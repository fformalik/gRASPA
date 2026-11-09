#!/bin/bash
# Build the vanilla (no ML) gRASPA on Perlmutter from src_clean into <repo>/build-vanilla.
# Applies Step 6 (experimental::filesystem) to the copy only, so src_clean stays clean in git.
# cuda/12.4 no longer exists on Perlmutter -> nvidia/24.5 + cudatoolkit/12.4.
# Usage: Cluster-Setup/NERSC_Installation_Nov2025/build_vanilla_perlmutter.sh   (from anywhere)
set -e
REPO=$(cd "$(dirname "$0")/../.." && pwd)
BUILD=$REPO/build-vanilla

module load PrgEnv-nvidia
module load nvidia/24.5
module load cudatoolkit/12.4

rm -rf "$BUILD"; mkdir -p "$BUILD"
cp "$REPO"/src_clean/*.h "$REPO"/src_clean/*.cuh "$REPO"/src_clean/*.cu "$REPO"/src_clean/*.cpp "$BUILD"/
cd "$BUILD"
sed -i "s/std::filesystem/std::experimental::filesystem/g" *
sed -i "s/<filesystem>/<experimental\/filesystem>/g" *

CXX="nvc++"
LINKFLAG="-lstdc++fs -L/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/cuda/lib64 -L/usr/lib64/"
NVCFLAG="-O3 -std=c++20 -Minline -mp -target=gpu -gpu=cc80 -cuda"

$CXX $NVCFLAG -c axpy.cu $LINKFLAG
$CXX $NVCFLAG -c main.cpp $LINKFLAG
$CXX $NVCFLAG -c read_data.cpp $LINKFLAG
$CXX $NVCFLAG -c data_struct.cpp $LINKFLAG
$CXX $NVCFLAG -c VDW_Coulomb.cu $LINKFLAG
$CXX $NVCFLAG main.o read_data.o axpy.o data_struct.o VDW_Coulomb.o -o nvc_main.x $LINKFLAG

# record which commit this binary was built from
{ git -C "$REPO" log -1 --format='%H %s'; git -C "$REPO" status --porcelain src_clean; } > BUILD_INFO
echo "built $BUILD/nvc_main.x from $(head -1 BUILD_INFO)"
