/*! Crops the solution in a given inner box.
 * \file crop.cpp
 */

#include <string>
#include <cstring>
#include <sys/stat.h>

#include <petscsys.h>

#include "petibm-utilities/field.h"
#include "petibm-utilities/grid.h"
#include "petibm-utilities/misc.h"


struct AppCtx
{
	PetscReal x_start, x_end,
	          y_start, y_end,
	          z_start, z_end;
}; // AppCtx


PetscErrorCode AppGetOptions(const char prefix[], AppCtx *ctx)
{
	PetscErrorCode ierr;
	PetscBool found;

	PetscFunctionBeginUser;

	ierr = PetscOptionsGetReal(
		nullptr, prefix, "-x_start", &ctx->x_start, &found); CHKERRQ(ierr);
	// get ending point in the x-direction
	ierr = PetscOptionsGetReal(
		nullptr, prefix, "-x_end", &ctx->x_end, &found); CHKERRQ(ierr);
	// get starting point in the y-direction
	ierr = PetscOptionsGetReal(
		nullptr, prefix, "-y_start", &ctx->y_start, &found); CHKERRQ(ierr);
	// get ending point in the y-direction
	ierr = PetscOptionsGetReal(
		nullptr, prefix, "-y_end", &ctx->y_end, &found); CHKERRQ(ierr);
	// get starting point in the z-direction
	ierr = PetscOptionsGetReal(
		nullptr, prefix, "-z_start", &ctx->z_start, &found); CHKERRQ(ierr);
	// get ending point in the z-direction
	ierr = PetscOptionsGetReal(
		nullptr, prefix, "-z_end", &ctx->z_end, &found); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // AppGetOptions


PetscErrorCode PetibmGridlineGetBoundingIndices(
	const PetibmGridline line, const PetscReal x_start, const PetscReal x_end,
	PetscInt &idx_start, PetscInt &idx_end)
{
	PetscErrorCode ierr;
	PetscReal *arr;
	DMDALocalInfo info;

	PetscFunctionBeginUser;

	ierr = DMDAGetLocalInfo(line.da, &info); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(line.da, line.coords, &arr); CHKERRQ(ierr);
	PetscInt i = info.xs;
	while (arr[i] < x_start and i < info.xs+info.xm)
		i++;
	idx_start = i;
	while (arr[i] < x_end and i < info.xs+info.xm)
		i++;
	idx_end = i-1;
	ierr = DMDAVecRestoreArray(line.da, line.coords, &arr); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmGridlineGetBoundingIndices


PetscErrorCode PetibmGridlineCrop(
	const PetibmGridline lineA, const PetscReal start, const PetscReal end,
	PetibmGridline &lineB)
{
	PetscErrorCode ierr;
	const PetscInt *lxA;
	PetscInt *lxB;
	PetscReal *lineA_arr, *lineB_arr;
	PetscInt m;
	PetscMPIInt rank;
	DMDALocalInfo info;

	PetscFunctionBeginUser;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

	ierr = DMDAGetOwnershipRanges(
		lineA.da, &lxA, nullptr, nullptr); CHKERRQ(ierr);
	ierr = DMDAGetInfo(lineA.da,
	                   nullptr,
	                   nullptr, nullptr, nullptr,
	                   &m, nullptr, nullptr,
	                   nullptr, nullptr,
	                   nullptr, nullptr, nullptr,
	                   nullptr); CHKERRQ(ierr);
	ierr = PetscMalloc(m*sizeof(*lxB), &lxB); CHKERRQ(ierr);
	ierr = PetscMemcpy(lxB, lxA, m*sizeof(*lxB)); CHKERRQ(ierr);

	// Count number of point on gridline between starting and ending points
	ierr = DMDAVecGetArray(lineA.da, lineA.coords, &lineA_arr); CHKERRQ(ierr);
	ierr = DMDAGetLocalInfo(lineA.da, &info); CHKERRQ(ierr);
	PetscInt xm = info.xm;
	for (PetscInt i=info.xs; i<info.xs+info.xm; i++)
	{
		if (!(lineA_arr[i] >= start and lineA_arr[i] <= end))
		{
			if (m == 1)
				lxB[0]--;
			else
				lxB[rank]--;
			xm--;
		}
	}
	ierr = MPI_Allgather(
		lxB+rank, 1, MPIU_INT, lxB, 1, MPIU_INT, PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lineA.da, lineA.coords, &lineA_arr); CHKERRQ(ierr);
	// Create the DMDA object for the sub-gridline
	if (m == 1)
	{
		ierr = DMDACreate1d(PETSC_COMM_SELF,
	                      DM_BOUNDARY_GHOSTED, xm, 1, 1, lxB,
	                      &lineB.da); CHKERRQ(ierr);
	}
	else
	{
		PetscInt M = 0;
		ierr = MPI_Allreduce(
			&xm, &M, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
		ierr = DMDACreate1d(PETSC_COMM_WORLD,
	                      DM_BOUNDARY_GHOSTED, M, 1, 1, lxB,
	                      &lineB.da); CHKERRQ(ierr);	
	}
	ierr = PetscFree(lxB); CHKERRQ(ierr);
	ierr = DMSetFromOptions(lineB.da); CHKERRQ(ierr);
	ierr = DMSetUp(lineB.da); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(lineB.da, &lineB.coords); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(lineB.da, &lineB.local); CHKERRQ(ierr);

	// Fill coordinates of the sub-gridline
	ierr = DMDAVecGetArray(lineA.da, lineA.coords, &lineA_arr); CHKERRQ(ierr);
	ierr = VecGetArray(lineB.local, &lineB_arr); CHKERRQ(ierr);
	PetscInt idx = 1;
	for (PetscInt i=info.xs; i<info.xs+info.xm; i++)
	{
		if (lineA_arr[i] >= start and lineA_arr[i] <= end)
		{
			lineB_arr[idx] = lineA_arr[i];
			idx++;
		}
	}
	ierr = DMDAVecRestoreArray(lineA.da, lineA.coords, &lineA_arr); CHKERRQ(ierr);
	ierr = VecRestoreArray(lineB.local, &lineB_arr); CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(
		lineB.da, lineB.local, INSERT_VALUES, lineB.coords); CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(
		lineB.da, lineB.local, INSERT_VALUES, lineB.coords); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmGridlineCrop


PetscErrorCode PetibmGridCrop(
	const PetibmGrid gridA, const AppCtx ctx, PetibmGrid &gridB)
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	// Crop the gridline along the x-direction
	ierr = PetibmGridlineCrop(
		gridA.x, ctx.x_start, ctx.x_end, gridB.x); CHKERRQ(ierr);
	// Crop the gridline along the y-direction
	ierr = PetibmGridlineCrop(
		gridA.y, ctx.y_start, ctx.y_end, gridB.y); CHKERRQ(ierr);
	if (gridA.dim == 3)
	{
		gridB.dim = 3;
		// Crop the gridline along the z-direction
		ierr = PetibmGridlineCrop(
			gridA.z, ctx.z_start, ctx.z_end, gridB.z); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
} // PetibmGridCrop


PetscErrorCode PetibmFieldCrop2d(
	const PetibmField fieldA,
	const PetscInt I_start, const PetscInt I_end,
	const PetscInt J_start, const PetscInt J_end,
	PetibmField &fieldB)
{
	PetscErrorCode ierr;
	PetscReal **arrA, **arrB;
	DMDALocalInfo info;
	PetscInt i, j;

	PetscFunctionBeginUser;

	ierr = DMDAVecGetArray(fieldA.da, fieldA.global, &arrA); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fieldB.da, fieldB.global, &arrB); CHKERRQ(ierr);
	ierr = DMDAGetLocalInfo(fieldB.da, &info); CHKERRQ(ierr);
	PetscInt I=I_start, J=J_start;
	for (j=info.ys; j<info.ym; j++)
	{
		for (i=info.xs; i<info.xm; i++)
		{
			arrB[j][i] = arrA[J][I];
			I++;
		}
		J++;
	}
	ierr = DMDAVecRestoreArray(fieldA.da, fieldA.global, &arrA); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fieldB.da, fieldB.global, &arrB); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmFieldCrop2d


PetscErrorCode PetibmFieldCrop3d(
	const PetibmField fieldA,
	const PetscInt I_start, const PetscInt I_end,
	const PetscInt J_start, const PetscInt J_end,
	const PetscInt K_start, const PetscInt K_end,
	PetibmField &fieldB)
{
	PetscErrorCode ierr;
	PetscReal ***arrA, ***arrB;
	DMDALocalInfo info;
	PetscInt i, j, k;

	PetscFunctionBeginUser;

	ierr = DMDAVecGetArray(fieldA.da, fieldA.global, &arrA); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fieldB.da, fieldB.global, &arrB); CHKERRQ(ierr);
	ierr = DMDAGetLocalInfo(fieldB.da, &info); CHKERRQ(ierr);
	PetscInt I=I_start, J=J_start, K=K_start;
	for (k=info.zs; k<info.zm; k++)
	{
		for (j=info.ys; j<info.ym; j++)
		{
			for (i=info.xs; i<info.xm; i++)
			{
				arrB[k][j][i] = arrA[K][J][I];
				I++;
			}
			J++;
		}
		K++;
	}
	ierr = DMDAVecRestoreArray(fieldA.da, fieldA.global, &arrA); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fieldB.da, fieldB.global, &arrB); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmFieldCrop3d


PetscErrorCode PetibmFieldCrop(
	const PetibmGrid gridA, const PetibmField fieldA,
	const AppCtx ctx, PetibmField &fieldB)
{
	PetscErrorCode ierr;
	PetscInt I_start, I_end,
	         J_start, J_end,
	         K_start, K_end;

	PetscFunctionBeginUser;

	ierr = PetibmGridlineGetBoundingIndices(
		gridA.x, ctx.x_start, ctx.x_end, I_start, I_end); CHKERRQ(ierr);
	ierr = PetibmGridlineGetBoundingIndices(
		gridA.y, ctx.y_start, ctx.y_end, J_start, J_end); CHKERRQ(ierr);
	if (gridA.dim == 3)
	{
		ierr = PetibmGridlineGetBoundingIndices(
			gridA.z, ctx.z_start, ctx.z_end, K_start, K_end); CHKERRQ(ierr);
		ierr = PetibmFieldCrop3d(fieldA,
		                         I_start, I_end, J_start, J_end, K_start, K_end,
		                         fieldB); CHKERRQ(ierr);
	}
	else
	{
		ierr = PetibmFieldCrop2d(fieldA,
		                         I_start, I_end, J_start, J_end,
		                         fieldB); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
} // PetibmFieldCrop


int main(int argc, char **argv)
{
	PetscErrorCode ierr;
	PetibmField fieldA, fieldB;
	PetibmFieldCtx fieldACtx, fieldBCtx;
	PetibmGrid gridA, gridB;
	PetibmGridCtx gridACtx, gridBCtx;
	AppCtx ctx;
	std::string outdir, datadir, filepath;

	ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);
	
	ierr = PetibmOptionsInsertFile(nullptr); CHKERRQ(ierr);
	ierr = PetibmGetDirectory(
		&datadir, "-data_directory", "solution"); CHKERRQ(ierr);
	ierr = PetibmGetDirectory(
		&outdir, "-output_directory", "output", PETSC_TRUE); CHKERRQ(ierr);
	ierr = AppGetOptions(nullptr, &ctx); CHKERRQ(ierr);

	// Create and read the grid A
	ierr = PetibmGridGetOptions("gridA_", &gridACtx); CHKERRQ(ierr);
	ierr = PetibmGridCtxPrintf("Grid A", gridACtx); CHKERRQ(ierr);
	ierr = PetibmGridInitialize(gridACtx, gridA); CHKERRQ(ierr);
	ierr = PetibmGridHDF5Read(
		gridACtx.path, gridACtx.name, gridA); CHKERRQ(ierr);

	// Crop grid A to get sub grid B and write
	std::strcpy(gridBCtx.name, gridACtx.name);
	ierr = PetibmGridCrop(gridA, ctx, gridB); CHKERRQ(ierr);
	filepath = outdir + "/grid.h5";
	ierr = PetibmGridHDF5Write(filepath, gridBCtx.name, gridB); CHKERRQ(ierr);

	// Create and read the field A
	ierr = PetibmFieldGetOptions("fieldA_", &fieldACtx); CHKERRQ(ierr);
	ierr = PetibmFieldCtxPrintf("Field A", fieldACtx); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(fieldACtx, gridA, fieldA); CHKERRQ(ierr);
	ierr = PetibmFieldHDF5Read(
		fieldACtx.path, fieldACtx.name, fieldA); CHKERRQ(ierr);

	// Create field B and crop field A and write field B
	std::strcpy(fieldBCtx.name, fieldACtx.name);
	ierr = PetibmFieldInitialize(fieldBCtx, gridB, fieldB); CHKERRQ(ierr);
	ierr = PetibmFieldCrop(gridA, fieldA, ctx, fieldB); CHKERRQ(ierr);
	filepath = outdir + "/0000010.h5";
	ierr = PetibmFieldHDF5Write(filepath, fieldBCtx.name, fieldB); CHKERRQ(ierr);

	// Clean workspace
	ierr = PetibmFieldDestroy(fieldA); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(fieldB); CHKERRQ(ierr);
	ierr = PetibmGridDestroy(gridA); CHKERRQ(ierr);
	ierr = PetibmGridDestroy(gridB); CHKERRQ(ierr);
	ierr = PetscFinalize(); CHKERRQ(ierr);

	return 0;
} // main
