# PetIBM - utilities
---

A collection of scripts and C++ codes to pre- and post-process the numerical solution of a PetIBM (0.2) simulation.


## Dependencies (last tested)
---

* PETSc (3.8.1)


## Contents
---

* `petibm-vorticity2d` (compute the vorticity field from the 2D velocity field)
* `petibm-vorticity3d` (compute the vorticity field from the 3D velocity field)
* `petibm-interpolation2d` (interpolate the 2D field values from one grid to another)
* `petibm-interpolation3d` (interpolate the 3D field values from one grid to another)


## Installation
---

```
PETIBM_UTILITIES_DIR=<directory>
cd $PETIBM_UTILITIES_DIR
mkdir build
cd build
PETSC_DIR=<petsc-dir>
PETSC_ARCH=<petsc-arch>
$PETIBM_UTILITIES_DIR/configure \
    --prefix=$PETIBM_UTILITIES_DIR/install \
    CXX=$PETSC_DIR/$PETSC_ARCH/bin/mpicxx \
    CXXFLAGS="-O3 -Wall -Wno-deprecated -std=c++11" \
    --with-petsc-dir=$PETSC_DIR \
    --with-petsc-arch=$PETSC_ARCH
make all
make install
export PATH=$PETIBM_UTILITIES_DIR/install/bin:$PATH
```
