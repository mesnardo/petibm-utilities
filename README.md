# petibm-utilities
---

Collection of applications for processing the numerical solution from [PetIBM](https://github.com/barbagroup/PetIBM).

Tested with PetIBM 0.2 and 0.3.


## Dependencies (last tested)
---

* [PETSc](https://www.mcs.anl.gov/petsc/) (3.8.1)


## Installation
---

Clone the GitHub repository:
```
git clone https://github.com/mesnardo/petibm-utilities.git
```

Configure, build, and install `petibm-utilities`:
```
cd petibm-utilities
PETIBM_UTILS_DIR=$PWD
mkdir build
cd build
$PETIBM_UTILS_DIR/configure \
    --prefix=$PETIBM_UTILS_DIR/install \
    CXX=mpicxx \
    CXXFLAGS="-O3 -Wall -Wno-deprecated -std=c++11" \
    --with-petsc-dir=<petsc-dir> \
    --with-petsc-arch=<petsc-arch> \
    --with-petibm=0.2  # or --with-petibm=0.3
make all
make install
export PATH=$PETIBM_UTILS_DIR/install/bin:$PATH
```

## Applications
---

* `petibm-vorticity2d` (compute and write the vorticity field from the 2D velocity field)
* `petibm-vorticity3d` (compute and write the vorticity field from the 3D velocity field)
* `petibm-interpolate` (interpolate field values from one grid to another)
* `petibm-convert` (write a field originally written in HDF5 format to PETSc binary format, or the other way around)
* `petibm-crop` (write a field and its grid by keeping values in a given sub-domain)

For each application, use the command-line `<application-name> -help intro` to print the command-line options.
