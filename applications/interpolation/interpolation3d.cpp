/*! Interpolates the 3D PetIBM solution from one grid to another.
 * \file interpolation2d.cpp
 */

#include <string>

#include <petscsys.h>

#include "petibm-utilities/field.h"
#include "petibm-utilities/grid.h"
#include "petibm-utilities/misc.h"


int main(int argc, char **argv)
{
	PetscErrorCode ierr;
	DMBoundaryType bType_x,
	               bType_y,
	               bType_z;

	ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);

	PetibmField fieldA;
	PetibmFieldCtx fieldACtx;
	PetibmGridCtx gridACtx;
	ierr = PetibmGridGetOptions("gridA_", &gridACtx); CHKERRQ(ierr);
	ierr = PetibmFieldGetOptions("fieldA_", &fieldACtx); CHKERRQ(ierr);
	bType_x = (gridACtx.periodic_x) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	bType_y = (gridACtx.periodic_y) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	bType_z = (gridACtx.periodic_z) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
	                    bType_x, bType_y, bType_z,
	                    DMDA_STENCIL_BOX,
	                    gridACtx.nx, gridACtx.ny, gridACtx.nz,
	                    PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
	                    1, 1,
	                    nullptr, nullptr, nullptr,
	                    &fieldA.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) fieldA.da, nullptr, "-fieldA_dmda_view"); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(fieldA); CHKERRQ(ierr);
	ierr = PetibmGridInitialize(gridACtx, fieldA.grid); CHKERRQ(ierr);
	ierr = PetibmGridHDF5Read(gridACtx.path, fieldA.grid); CHKERRQ(ierr);
	ierr = PetibmFieldHDF5Read(
		fieldACtx.path, fieldACtx.name, fieldA); CHKERRQ(ierr);
	ierr = PetibmFieldExternalGhostPointsSet(
		fieldA, fieldACtx.bc_value); CHKERRQ(ierr);

	PetibmField fieldB;
	PetibmFieldCtx fieldBCtx;
	PetibmGridCtx gridBCtx;
	ierr = PetibmGridGetOptions("gridB_", &gridBCtx); CHKERRQ(ierr);
	ierr = PetibmFieldGetOptions("fieldB_", &fieldBCtx); CHKERRQ(ierr);
	bType_x = (gridBCtx.periodic_x) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	bType_y = (gridBCtx.periodic_y) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	bType_z = (gridBCtx.periodic_z) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
	                    bType_x, bType_y, bType_z,
	                    DMDA_STENCIL_BOX,
	                    gridBCtx.nx, gridBCtx.ny, gridBCtx.nz,
	                    PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
	                    1, 1,
	                    nullptr, nullptr, nullptr,
	                    &fieldB.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) fieldB.da, nullptr, "-fieldB_dmda_view"); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(fieldB); CHKERRQ(ierr);
	ierr = PetibmGridInitialize(gridBCtx, fieldB.grid); CHKERRQ(ierr);
	ierr = PetibmGridHDF5Read(gridBCtx.path, fieldB.grid); CHKERRQ(ierr);
	ierr = PetibmFieldExternalGhostPointsSet(
		fieldB, fieldBCtx.bc_value); CHKERRQ(ierr);

	ierr = PetibmFieldInterpolate3D(
		fieldA, fieldB, fieldBCtx.bc_value); CHKERRQ(ierr);

	ierr = PetibmFieldHDF5Write(
		fieldBCtx.path, fieldBCtx.name, fieldB); CHKERRQ(ierr);

	ierr = PetibmFieldDestroy(fieldA); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(fieldB); CHKERRQ(ierr);
	ierr = PetscFinalize(); CHKERRQ(ierr);

	return 0;
} // main
