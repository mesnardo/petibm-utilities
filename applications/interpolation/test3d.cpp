/*! Tests 3D interpolation.
 * \file test3d.cpp
 */

#include <cmath>
#include <string>

#include <petscsys.h>

#include "petibm-utilities/mesh.h"
#include "petibm-utilities/misc.h"
#include "petibm-utilities/field.h"


int main(int argc, char **argv)
{
	PetscErrorCode ierr;
	DMBoundaryType bType_x,
	               bType_y,
	               bType_z;
	PetscReal xstart = 0.0, xend = 1.0,
	          ystart = 1.0, yend = 2.0,
	          zstart = 2.0, zend = 3.0,
	          h,
	          gammaA = 1.01, gammaB = 1.02;
	PetscInt i, j, k, ite, nt = 100;
	PetscBool found;

	ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);

	ierr = PetscOptionsGetInt(
		nullptr, nullptr, "-nt", &nt, &found); CHKERRQ(ierr);

	PetibmMeshInfo meshAInfo;
	meshAInfo.nx = 8;
	meshAInfo.ny = 8;
	meshAInfo.nz = 8;

	PetibmMeshInfo meshBInfo;
	meshBInfo.nx = 4;
	meshBInfo.ny = 4;
	meshBInfo.nz = 4;

	PetibmField fieldA;
	bType_x = (meshAInfo.periodic_x) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	bType_y = (meshAInfo.periodic_y) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	bType_z = (meshAInfo.periodic_z) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
	                    bType_x, bType_y, bType_z,
	                    DMDA_STENCIL_BOX,
	                    meshAInfo.nx, meshAInfo.ny, meshAInfo.nz,
	                    PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
	                    1, 1, nullptr, nullptr, nullptr,
	                    &fieldA.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) fieldA.da, nullptr, "-fieldA_dmda_view"); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(fieldA); CHKERRQ(ierr);
	ierr = VecSet(fieldA.global, 1.2345);
	ierr = DMGlobalToLocalBegin(
		fieldA.da, fieldA.global, INSERT_VALUES, fieldA.local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(
		fieldA.da, fieldA.global, INSERT_VALUES, fieldA.local); CHKERRQ(ierr);
	ierr = PetibmFieldExternalGhostPointsSet(fieldA, 1.2345); CHKERRQ(ierr);

	PetscReal *xA, *yA, *zA;
	ierr = VecGetArray(fieldA.x, &xA); CHKERRQ(ierr);
	ierr = VecGetArray(fieldA.y, &yA); CHKERRQ(ierr);
	ierr = VecGetArray(fieldA.z, &zA); CHKERRQ(ierr);
	h = (xend - xstart) * (gammaA - 1.0) / (std::pow(gammaA, meshAInfo.nx) - 1.0);
	xA[0] = xstart + 0.5 * h;
	for (i=0; i<meshAInfo.nx-1; i++)
		xA[i+1] = xA[i] + h * std::pow(gammaA, i);
	h = (yend - ystart) * (gammaA - 1.0) / (std::pow(gammaA, meshAInfo.ny) - 1.0);
	yA[0] = ystart + 0.5 * h;
	for (j=0; j<meshAInfo.ny-1; j++)
		yA[j+1] = yA[j] + h * std::pow(gammaA, j);
	h = (zend - zstart) * (gammaA - 1.0) / (std::pow(gammaA, meshAInfo.nz) - 1.0);
	zA[0] = zstart + 0.5 * h;
	for (k=0; k<meshAInfo.nz-1; k++)
		zA[k+1] = zA[k] + h * std::pow(gammaA, k);
	ierr = VecRestoreArray(fieldA.x, &xA); CHKERRQ(ierr);
	ierr = VecRestoreArray(fieldA.y, &yA); CHKERRQ(ierr);
	ierr = VecRestoreArray(fieldA.z, &zA); CHKERRQ(ierr);

	PetibmField fieldB;
	bType_x = (meshBInfo.periodic_x) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	bType_y = (meshBInfo.periodic_y) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	bType_z = (meshBInfo.periodic_z) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
	                    bType_x, bType_y, bType_z,
	                    DMDA_STENCIL_BOX,
	                    meshBInfo.nx, meshBInfo.ny, meshBInfo.nz,
	                    PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
	                    1, 1, nullptr, nullptr, nullptr,
	                    &fieldB.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) fieldB.da, nullptr, "-fieldB_dmda_view"); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(fieldB); CHKERRQ(ierr);
	ierr = PetibmFieldExternalGhostPointsSet(fieldB, 1.2345); CHKERRQ(ierr);

	PetscReal *xB, *yB, *zB;
	ierr = VecGetArray(fieldB.x, &xB); CHKERRQ(ierr);
	ierr = VecGetArray(fieldB.y, &yB); CHKERRQ(ierr);
	ierr = VecGetArray(fieldB.z, &zB); CHKERRQ(ierr);
	h = (xend - xstart) * (gammaB - 1.0) / (std::pow(gammaB, meshBInfo.nx) - 1.0);
	xB[0] = xstart + 0.5 * h;
	for (i=0; i<meshBInfo.nx-1; i++)
		xB[i+1] = xB[i] + h * std::pow(gammaB, i);
	h = (yend - ystart) * (gammaB - 1.0) / (std::pow(gammaB, meshBInfo.ny) - 1.0);
	yB[0] = ystart + 0.5 * h;
	for (j=0; j<meshBInfo.ny-1; j++)
		yB[j+1] = yB[j] + h * std::pow(gammaB, j);
	h = (zend - zstart) * (gammaB - 1.0) / (std::pow(gammaB, meshBInfo.nz) - 1.0);
	zB[0] = zstart + 0.5 * h;
	for (k=0; k<meshBInfo.nz-1; k++)
		zB[k+1] = zB[k] + h * std::pow(gammaB, k);
	ierr = VecRestoreArray(fieldB.x, &xB); CHKERRQ(ierr);
	ierr = VecRestoreArray(fieldB.y, &yB); CHKERRQ(ierr);
	ierr = VecRestoreArray(fieldB.z, &zB); CHKERRQ(ierr);

	for (ite=0; ite<nt; ite++)
	{
		if (ite % 4 == 0 or ite % 4 == 1)
		{
			ierr = PetibmFieldInterpolate3D(fieldA, fieldB); CHKERRQ(ierr);
		}
		else
		{
			ierr = PetibmFieldInterpolate3D(fieldB, fieldA); CHKERRQ(ierr);
		}
	}

	PetscMPIInt rank;
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
	if (rank == 0)
	{
		ierr = PetibmFieldWriteGrid("gridA.h5", fieldA); CHKERRQ(ierr);
		ierr = PetibmFieldWriteGrid("gridB.h5", fieldB); CHKERRQ(ierr);
	}
	ierr = PetibmFieldWriteValues("fieldA.h5", "phi", fieldA); CHKERRQ(ierr);
	ierr = PetibmFieldWriteValues("fieldB.h5", "phi", fieldB); CHKERRQ(ierr);

	ierr = PetibmFieldDestroy(fieldA); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(fieldB); CHKERRQ(ierr);
	ierr = PetscFinalize(); CHKERRQ(ierr);

	return 0;
} // main
