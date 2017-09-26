/*! Interpolates a PetIBM field from one grid to another.
 * \file interpolation.cpp
 */

#include <string>

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

	ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);
	
	// create fieldA
	ierr = PetibmGridGetOptions("fieldA_", &gridACtx); CHKERRQ(ierr);
	ierr = PetibmGridInitialize(gridACtx, gridA); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) gridA.x.da, nullptr, "-gridA_x_dmda_view"); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) gridA.y.da, nullptr, "-gridA_y_dmda_view"); CHKERRQ(ierr);
	if (dim == 3)
	{
		ierr = PetscObjectViewFromOptions(
			(PetscObject) gridA.z.da, nullptr, "-gridA_z_dmda_view"); CHKERRQ(ierr);
	}
	ierr = PetibmGridHDF5Read(gridACtx.path, gridA); CHKERRQ(ierr);
	ierr = PetibmGridSetBoundaryPoints(
		gridACtx.starts, gridACtx.ends, gridA); CHKERRQ(ierr);
	ierr = PetibmFieldGetOptions("fieldA_", &fieldACtx); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(fieldACtx, gridA, fieldA); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) fieldA.da, nullptr, "-fieldA_dmda_view"); CHKERRQ(ierr);
	ierr = PetibmFieldHDF5Read(
		fieldACtx.path, fieldACtx.name, fieldA); CHKERRQ(ierr);
	ierr = PetibmFieldSetBoundaryPoints(
		fieldACtx.bc_value, fieldA); CHKERRQ(ierr);

	// create fieldB
	ierr = PetibmGridGetOptions("fieldB_", &gridBCtx); CHKERRQ(ierr);
	Vec coords[dim];
	ierr = VecCreateSeq(PETSC_COMM_SELF, gridBCtx.nx, coords); CHKERRQ(ierr);
	ierr = PetibmGridlineHDF5Read(gridBCtx.path, "x", coords[0]); CHKERRQ(ierr);
	ierr = VecCreateSeq(PETSC_COMM_SELF, gridBCtx.ny, coords+1); CHKERRQ(ierr);
	ierr = PetibmGridlineHDF5Read(gridBCtx.path, "y", coords[1]); CHKERRQ(ierr);
	if (dim == 3)
	{
		ierr = VecCreateSeq(PETSC_COMM_SELF, gridBCtx.nz, coords+2); CHKERRQ(ierr);
		ierr = PetibmGridlineHDF5Read(gridBCtx.path, "z", coords[2]); CHKERRQ(ierr);
	}
	ierr = PetibmGridInitialize(gridA, coords, gridB); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) gridB.x.da, nullptr, "-gridB_x_dmda_view"); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) gridB.y.da, nullptr, "-gridB_y_dmda_view"); CHKERRQ(ierr);
	if (dim == 3)
	{
		ierr = PetscObjectViewFromOptions(
			(PetscObject) gridB.z.da, nullptr, "-gridB_z_dmda_view"); CHKERRQ(ierr);
	}
	ierr = PetibmGridHDF5Read(gridBCtx.path, gridB); CHKERRQ(ierr);
	ierr = PetibmGridSetBoundaryPoints(
		gridBCtx.starts, gridBCtx.ends, gridB); CHKERRQ(ierr);
	ierr = PetibmFieldGetOptions("fieldB_", &fieldBCtx); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(fieldBCtx, gridB, fieldB); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) fieldB.da, nullptr, "-fieldB_dmda_view"); CHKERRQ(ierr);
	ierr = PetibmFieldSetBoundaryPoints(
		fieldBCtx.bc_value, fieldB); CHKERRQ(ierr);

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
