/*! Interpolates the 2D PetIBM solution from one grid to another.
 * \file interpolate2d.cpp
 */

#include <string>

#include <petscsys.h>

#include "petibm-utilities/mesh.h"
#include "petibm-utilities/misc.h"
#include "petibm-utilities/field.h"


int main(int argc, char **argv)
{
	PetscErrorCode ierr;
	DMBoundaryType bType_x,
	               bType_y;

	ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);

	PetibmMeshInfo meshAInfo;
	PetibmFieldInfo fieldAInfo;
	ierr = PetibmMeshInfoGetOptions("fieldA_", &meshAInfo); CHKERRQ(ierr);
	ierr = PetibmFieldInfoGetOptions("fieldA_", &fieldAInfo); CHKERRQ(ierr);

	PetibmMeshInfo meshBInfo;
	PetibmFieldInfo fieldBInfo;
	ierr = PetibmMeshInfoGetOptions("fieldB_", &meshBInfo); CHKERRQ(ierr);
	ierr = PetibmFieldInfoGetOptions("fieldB_", &fieldBInfo); CHKERRQ(ierr);

	PetibmField fieldA;
	bType_x = (meshAInfo.periodic_x) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	bType_y = (meshAInfo.periodic_y) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	ierr = DMDACreate2d(PETSC_COMM_WORLD,
	                    bType_x, bType_y,
	                    DMDA_STENCIL_BOX,
	                    meshAInfo.nx, meshAInfo.ny,
	                    PETSC_DECIDE, PETSC_DECIDE, 1, 1, nullptr, nullptr,
	                    &fieldA.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) fieldA.da, nullptr, "-fieldA_dmda_view"); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(fieldA); CHKERRQ(ierr);
	ierr = PetibmFieldReadGrid(meshAInfo.path, fieldA); CHKERRQ(ierr);
	ierr = PetibmFieldReadValues(
		fieldAInfo.path, fieldAInfo.name, fieldA); CHKERRQ(ierr);
	ierr = PetibmFieldExternalGhostPointsSet(
		fieldA, fieldAInfo.bcValue); CHKERRQ(ierr);

	PetibmField fieldB;
	bType_x = (meshBInfo.periodic_x) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	bType_y = (meshBInfo.periodic_y) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	ierr = DMDACreate2d(PETSC_COMM_WORLD,
	                    bType_x, bType_y,
	                    DMDA_STENCIL_BOX,
	                    meshBInfo.nx, meshBInfo.ny,
	                    PETSC_DECIDE, PETSC_DECIDE, 1, 1, nullptr, nullptr,
	                    &fieldB.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) fieldB.da, nullptr, "-fieldB_dmda_view"); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(fieldB); CHKERRQ(ierr);
	ierr = PetibmFieldReadGrid(meshBInfo.path, fieldB); CHKERRQ(ierr);
	ierr = PetibmFieldExternalGhostPointsSet(
		fieldB, fieldBInfo.bcValue); CHKERRQ(ierr);

	ierr = PetibmFieldInterpolate2D(fieldA, fieldB); CHKERRQ(ierr);

	ierr = PetibmFieldWriteValues(
		fieldBInfo.path, fieldBInfo.name, fieldB); CHKERRQ(ierr);

	ierr = PetibmFieldDestroy(fieldA); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(fieldB); CHKERRQ(ierr);
	ierr = PetscFinalize(); CHKERRQ(ierr);

	return 0;
} // main
