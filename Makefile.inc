# Options common to all Makefiles

PETSC_DIR = $(HOME)/src/petsc/3.7.4
PETSC_ARCH = darwin-gnu-opt

include $(PETSC_DIR)/lib/petsc/conf/variables
include $(PETSC_DIR)/lib/petsc/conf/rules

CXX = g++-6
CXXFLAGS = -O3 -std=c++11 -Wall -Wno-deprecated
