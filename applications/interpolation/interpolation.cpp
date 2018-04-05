/*! Interpolates a PetIBM field A from grid A to grid B.
 * \file interpolation.cpp
 */

#include <string>
#include <sys/stat.h>

#include <petscsys.h>

#include "petibm-utilities/field.h"
#include "petibm-utilities/grid.h"
#include "petibm-utilities/misc.h"

#ifndef DIMENSIONS
#define DIMENSIONS 2
#endif


int main(int argc, char **argv)
{
	PetscErrorCode ierr;
	PetibmField fieldA, fieldB;
	PetibmFieldCtx fieldACtx, fieldBCtx;
	PetibmGrid gridA, gridB;
	PetibmGridCtx gridACtx, gridBCtx;
	const PetscInt dim = DIMENSIONS;
	std::string outdir;

	ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);
	
	ierr = PetibmOptionsInsertFile(nullptr); CHKERRQ(ierr);
	ierr = PetibmGetDirectory(
		&outdir, "-output_directory", PETSC_TRUE); CHKERRQ(ierr);

	// Create and read the grid A
	ierr = PetibmGridGetOptions("gridA_", &gridACtx); CHKERRQ(ierr);
	ierr = PetibmGridCtxPrintf("Grid A", gridACtx); CHKERRQ(ierr);
	ierr = PetibmGridInitialize(gridACtx, gridA); CHKERRQ(ierr);
#ifdef PETIBM_0_3
	ierr = PetibmGridHDF5Read(gridACtx.path, gridACtx.name, gridA); CHKERRQ(ierr);
#elif PETIBM_0_2
	ierr = PetibmGridHDF5Read(gridACtx.path, gridA); CHKERRQ(ierr);
#endif
	ierr = PetibmGridSetBoundaryPoints(
		gridACtx.starts, gridACtx.ends, gridA); CHKERRQ(ierr);
	// Create and read the field A
	ierr = PetibmFieldGetOptions("fieldA_", &fieldACtx); CHKERRQ(ierr);
	ierr = PetibmFieldCtxPrintf("Field A", fieldACtx); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(fieldACtx, gridA, fieldA); CHKERRQ(ierr);
	ierr = PetibmFieldHDF5Read(
		fieldACtx.path, fieldACtx.name, fieldA); CHKERRQ(ierr);
	ierr = PetibmFieldSetBoundaryPoints(
		fieldACtx.bc_value, fieldA); CHKERRQ(ierr);

	// Create and read the grid B
	ierr = PetibmGridGetOptions("gridB_", &gridBCtx); CHKERRQ(ierr);
	ierr = PetibmGridCtxPrintf("Grid B", gridBCtx); CHKERRQ(ierr);
	Vec coords[dim];
	ierr = VecCreateSeq(PETSC_COMM_SELF, gridBCtx.nx, coords); CHKERRQ(ierr);
#ifdef PETIBM_0_3
	ierr = PetibmGridlineHDF5Read(
		gridBCtx.path, gridBCtx.name, "x", coords[0]); CHKERRQ(ierr);
#elif PETIBM_0_2
	ierr = PetibmGridlineHDF5Read(gridBCtx.path, "x", coords[0]); CHKERRQ(ierr);
#endif
	ierr = VecCreateSeq(PETSC_COMM_SELF, gridBCtx.ny, coords+1); CHKERRQ(ierr);
#ifdef PETIBM_0_3
	ierr = PetibmGridlineHDF5Read(
		gridBCtx.path, gridBCtx.name, "y", coords[1]); CHKERRQ(ierr);
#elif PETIBM_0_2
	ierr = PetibmGridlineHDF5Read(gridBCtx.path, "y", coords[1]); CHKERRQ(ierr);
#endif
	if (dim == 3)
	{
		ierr = VecCreateSeq(PETSC_COMM_SELF, gridBCtx.nz, coords+2); CHKERRQ(ierr);
#ifdef PETIBM_0_3
		ierr = PetibmGridlineHDF5Read(
			gridBCtx.path, gridBCtx.name, "z", coords[2]); CHKERRQ(ierr);
#elif PETIBM_0_2
		ierr = PetibmGridlineHDF5Read(gridBCtx.path, "z", coords[2]); CHKERRQ(ierr);
#endif
	}
	ierr = PetibmGridInitialize(gridA, coords, gridB); CHKERRQ(ierr);
#ifdef PETIBM_0_3
	ierr = PetibmGridHDF5Read(
		gridBCtx.path, gridBCtx.name, gridB); CHKERRQ(ierr);
#elif PETIBM_0_2
	ierr = PetibmGridHDF5Read(gridBCtx.path, gridB); CHKERRQ(ierr);
#endif
	ierr = PetibmGridSetBoundaryPoints(
		gridBCtx.starts, gridBCtx.ends, gridB); CHKERRQ(ierr);
	// Create the field B
	ierr = PetibmFieldGetOptions("fieldB_", &fieldBCtx); CHKERRQ(ierr);
	ierr = PetibmFieldCtxPrintf("Field B", fieldBCtx); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(fieldBCtx, gridB, fieldB); CHKERRQ(ierr);
	ierr = PetibmFieldSetBoundaryPoints(
		fieldBCtx.bc_value, fieldB); CHKERRQ(ierr);

	// Interpolate field A (defined on grid A) onto field B (defined on grid B)
	ierr = PetibmFieldInterpolate(gridA, fieldA, gridB, fieldB); CHKERRQ(ierr);

	ierr = PetibmFieldHDF5Write(
		fieldBCtx.path, fieldBCtx.name, fieldB); CHKERRQ(ierr);

	ierr = PetibmFieldDestroy(fieldA); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(fieldB); CHKERRQ(ierr);
	ierr = PetibmGridDestroy(gridA); CHKERRQ(ierr);
	ierr = PetibmGridDestroy(gridB); CHKERRQ(ierr);
	ierr = PetscFinalize(); CHKERRQ(ierr);

	return 0;
} // main
