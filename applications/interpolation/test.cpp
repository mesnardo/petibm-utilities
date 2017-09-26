/*! Tests interpolation of field values from one grid to another.
 * \file test.cpp
 */

#include <cmath>
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
	PetscReal *xA, *yA, *zA, *xB, *yB, *zB;
	Vec coords[3];
	DMDALocalInfo info;
	PetscInt i, j, k;
	PetscReal starts[3] = {0.0, 1.0, 2.0},
	          ends[3] = {1.0, 2.0, 4.0},
	          h;
	PetscBool found, fine2coarse = PETSC_FALSE;
	const PetscInt dim = DIMENSIONS;

	ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);

	ierr = PetscOptionsGetBool(nullptr, nullptr, "-fine_to_coarse",
	                           &fine2coarse, &found); CHKERRQ(ierr);

	gridACtx.nx = (fine2coarse) ? 8 : 6;
	gridACtx.ny = (fine2coarse) ? 7 : 5;
	if (dim == 3)
		gridACtx.nz = (fine2coarse) ? 6 : 4;
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
	
	h = (ends[0] - starts[0]) / gridACtx.nx;
	ierr = DMDAGetLocalInfo(gridA.x.da, &info); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(gridA.x.da, gridA.x.coords, &xA); CHKERRQ(ierr);
	for (i=info.xs; i<info.xs+info.xm; i++)
		xA[i] = starts[0] + (0.5 + i) * h;
	ierr = DMDAVecRestoreArray(gridA.x.da, gridA.x.coords, &xA); CHKERRQ(ierr);
	h = (ends[1] - starts[1]) / gridACtx.ny;
	ierr = DMDAGetLocalInfo(gridA.y.da, &info); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(gridA.y.da, gridA.y.coords, &yA); CHKERRQ(ierr);
	for (j=info.xs; j<info.xs+info.xm; j++)
		yA[j] = starts[1] + (0.5 + j) * h;
	ierr = DMDAVecRestoreArray(gridA.y.da, gridA.y.coords, &yA); CHKERRQ(ierr);
	if (dim == 3)
	{
		h = (ends[2] - starts[2]) / gridACtx.nz;
		ierr = DMDAGetLocalInfo(gridA.z.da, &info); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(gridA.z.da, gridA.z.coords, &zA); CHKERRQ(ierr);
		for (k=info.xs; k<info.xs+info.xm; k++)
			zA[k] = starts[2] + (0.5 + k) * h;
		ierr = DMDAVecRestoreArray(gridA.z.da, gridA.z.coords, &zA); CHKERRQ(ierr);
	}
	ierr = PetibmGridSetBoundaryPoints(starts, ends, gridA); CHKERRQ(ierr);
	ierr = PetibmGridHDF5Write("gridA.h5", gridA); CHKERRQ(ierr);

	fieldACtx.bc_value = 1.2345;
	ierr = PetibmFieldInitialize(fieldACtx, gridA, fieldA); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) fieldA.da, nullptr, "-fieldA_dmda_view"); CHKERRQ(ierr);
	ierr = VecSet(fieldA.global, 1.2345);
	ierr = PetibmFieldSetBoundaryPoints(1.2345, fieldA); CHKERRQ(ierr);

	gridBCtx.nx = (fine2coarse) ? 6 : 8;
	gridBCtx.ny = (fine2coarse) ? 5 : 7;
	if (dim == 3)
		gridBCtx.nz = (fine2coarse) ? 4 : 6;
	ierr = VecCreateSeq(PETSC_COMM_SELF, gridBCtx.nx, coords); CHKERRQ(ierr);
	h = (ends[0] - starts[0]) / gridBCtx.nx;
	ierr = VecGetArray(coords[0], &xB); CHKERRQ(ierr);
	for (i=0; i<gridBCtx.nx; i++)
		xB[i] = starts[0] + (0.5 + i) * h;
	ierr = VecRestoreArray(coords[0], &xB); CHKERRQ(ierr);
	ierr = VecCreateSeq(PETSC_COMM_SELF, gridBCtx.ny, coords+1); CHKERRQ(ierr);
	h = (ends[1] - starts[1]) / gridBCtx.ny;
	ierr = VecGetArray(coords[1], &yB); CHKERRQ(ierr);
	for (j=0; j<gridBCtx.ny; j++)
		yB[j] = starts[1] + (0.5 + j) * h;
	ierr = VecRestoreArray(coords[1], &yB); CHKERRQ(ierr);
	if (dim == 3)
	{
		ierr = VecCreateSeq(PETSC_COMM_SELF, gridBCtx.nz, coords+2); CHKERRQ(ierr);
		h = (ends[2] - starts[2]) / gridBCtx.nz;
		ierr = VecGetArray(coords[2], &zB); CHKERRQ(ierr);
		for (k=0; k<gridBCtx.nz; k++)
			zB[k] = starts[2] + (0.5 + k) * h;
		ierr = VecRestoreArray(coords[2], &zB); CHKERRQ(ierr);
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
	
	h = (ends[0] - starts[0]) / gridBCtx.nx;
	ierr = DMDAGetLocalInfo(gridB.x.da, &info); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(gridB.x.da, gridB.x.coords, &xB); CHKERRQ(ierr);
	for (i=info.xs; i<info.xs+info.xm; i++)
		xB[i] = starts[0] + (0.5 + i) * h;
	ierr = DMDAVecRestoreArray(gridB.x.da, gridB.x.coords, &xB); CHKERRQ(ierr);
	h = (ends[1] - starts[1]) / gridBCtx.ny;
	ierr = DMDAGetLocalInfo(gridB.y.da, &info); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(gridB.y.da, gridB.y.coords, &yB); CHKERRQ(ierr);
	for (j=info.xs; j<info.xs+info.xm; j++)
		yB[j] = starts[1] + (0.5 + j) * h;
	ierr = DMDAVecRestoreArray(gridB.y.da, gridB.y.coords, &yB); CHKERRQ(ierr);
	if (dim == 3)
	{
		h = (ends[2] - starts[2]) / gridBCtx.nz;
		ierr = DMDAGetLocalInfo(gridB.z.da, &info); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(gridB.z.da, gridB.z.coords, &zB); CHKERRQ(ierr);
		for (k=info.xs; k<info.xs+info.xm; k++)
			zB[k] = starts[2] + (0.5 + k) * h;
		ierr = DMDAVecRestoreArray(gridB.z.da, gridB.z.coords, &zB); CHKERRQ(ierr);
	}
	ierr = PetibmGridSetBoundaryPoints(starts, ends, gridB); CHKERRQ(ierr);
	ierr = PetibmGridHDF5Write("gridB.h5", gridB); CHKERRQ(ierr);

	ierr = PetibmFieldInitialize(fieldBCtx, gridB, fieldB); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) fieldB.da, nullptr, "-fieldB_dmda_view"); CHKERRQ(ierr);
	ierr = PetibmFieldSetBoundaryPoints(1.2345, fieldB); CHKERRQ(ierr);

	ierr = PetibmFieldInterpolate(
		gridA, fieldA, gridB, fieldB); CHKERRQ(ierr);

	ierr = PetibmFieldHDF5Write("fieldA.h5", "phi", fieldA); CHKERRQ(ierr);
	ierr = PetibmFieldHDF5Write("fieldB.h5", "phi", fieldB); CHKERRQ(ierr);

	ierr = PetibmFieldDestroy(fieldA); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(fieldB); CHKERRQ(ierr);
	ierr = VecDestroy(coords); CHKERRQ(ierr);
	ierr = VecDestroy(coords+1); CHKERRQ(ierr);
	if (dim == 3)
	{
		ierr = VecDestroy(coords+2); CHKERRQ(ierr);
	}
	ierr = PetscFinalize(); CHKERRQ(ierr);

	return 0;
} // main
