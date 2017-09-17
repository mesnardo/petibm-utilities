# PetIBM - utilities
---

A collection of scripts and C++ codes to pre- and post-process the numerical solution of a PetIBM (0.2) simulation.


## Dependencies
---

* PETSc (3.7.4)


## Contents
---

* `petibm-vorticity2d` (compute the 2D vorticity field from the 2D velocity field)


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
