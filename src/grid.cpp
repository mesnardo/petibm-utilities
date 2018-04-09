/*! Implementation of the functions related to the grid.
 * \file grid.cpp
 */

#include <fstream>

#include <petscviewerhdf5.h>

#include "petibm-utilities/grid.h"
#include "petibm-utilities/misc.h"


/*! Gets options from command-line or config file.
 *
 * \param prefix String to prepend the name of the options.
 * \param ctx The PetibmGridCtx structure to fill.
 */
PetscErrorCode PetibmGridGetOptions(
	const char prefix[], PetibmGridCtx *ctx)
{
	PetscErrorCode ierr;
	PetscBool found;

	PetscFunctionBeginUser;

	// get path of the file with the gridline stations
	ierr = PetscOptionsGetString(nullptr, prefix, "-path", ctx->path,
	                             sizeof(ctx->path), &found); CHKERRQ(ierr);
	// get name of the grid
	ierr = PetscOptionsGetString(nullptr, prefix, "-name", ctx->name,
	                             sizeof(ctx->name), &found); CHKERRQ(ierr);
	// get number of points in the x-direction
	ierr = PetscOptionsGetInt(
		nullptr, prefix, "-nx", &ctx->nx, &found); CHKERRQ(ierr);
	// get number of points in the y-direction
	ierr = PetscOptionsGetInt(
		nullptr, prefix, "-ny", &ctx->ny, &found); CHKERRQ(ierr);
	// get number of points in the z-direction
	ierr = PetscOptionsGetInt(
		nullptr, prefix, "-nz", &ctx->nz, &found); CHKERRQ(ierr);
	// get starting point in the x-direction
	ierr = PetscOptionsGetReal(
		nullptr, prefix, "-x_start", &ctx->starts[0], &found); CHKERRQ(ierr);
	// get ending point in the x-direction
	ierr = PetscOptionsGetReal(
		nullptr, prefix, "-x_end", &ctx->ends[0], &found); CHKERRQ(ierr);
	// get starting point in the y-direction
	ierr = PetscOptionsGetReal(
		nullptr, prefix, "-y_start", &ctx->starts[1], &found); CHKERRQ(ierr);
	// get ending point in the y-direction
	ierr = PetscOptionsGetReal(
		nullptr, prefix, "-y_end", &ctx->ends[1], &found); CHKERRQ(ierr);
	// get starting point in the z-direction
	ierr = PetscOptionsGetReal(
		nullptr, prefix, "-z_start", &ctx->starts[2], &found); CHKERRQ(ierr);
	// get ending point in the z-direction
	ierr = PetscOptionsGetReal(
		nullptr, prefix, "-z_end", &ctx->ends[2], &found); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmGridGetOptions


PetscErrorCode PetibmGridCtxPrintf(
	const std::string name, const PetibmGridCtx ctx)
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	ierr = PetscPrintf(PETSC_COMM_WORLD, "+ %s:\n", name.c_str()); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\t- path: %s\n",
	                   ctx.path); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\t- size: %d x %d",
	                   ctx.nx, ctx.ny); CHKERRQ(ierr);
	if (ctx.nz > 0)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " x %d\n", ctx.nz); CHKERRQ(ierr);
	}
	else
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "\n"); CHKERRQ(ierr);
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\t- x: [%g, %g]\n",
	                   ctx.starts[0], ctx.ends[0]); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\t- y: [%g, %g]\n",
	                   ctx.starts[1], ctx.ends[1]); CHKERRQ(ierr);
	if (ctx.nz > 0)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "\t- z: [%g, %g]\n",
		                   ctx.starts[2], ctx.ends[2]); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
} // PetibmGridCtxPrintf


/*! Initializes the grid based on the context.
 *
 * Creates a 1D DMDA object for each direction and creates the local global
 * vectors associated with each DMDA.
 * The domain is decomposed along the y-direction only.
 *
 * \param ctx The grid context.
 * \param grid The grid to initialize (passed by reference).
 */
PetscErrorCode PetibmGridInitialize(
	const PetibmGridCtx ctx, PetibmGrid &grid)
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	// gridline in x-direction
	ierr = DMDACreate1d(PETSC_COMM_SELF,
	                    DM_BOUNDARY_GHOSTED, ctx.nx, 1, 1, nullptr,
	                    &grid.x.da); CHKERRQ(ierr);
	ierr = DMSetFromOptions(grid.x.da); CHKERRQ(ierr);
	ierr = DMSetUp(grid.x.da); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(grid.x.da, &grid.x.coords); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(grid.x.da, &grid.x.local); CHKERRQ(ierr);
	// gridline in y-direction
	ierr = DMDACreate1d(PETSC_COMM_WORLD,
	                    DM_BOUNDARY_GHOSTED, ctx.ny, 1, 1, nullptr,
	                    &grid.y.da); CHKERRQ(ierr);
	ierr = DMSetFromOptions(grid.y.da); CHKERRQ(ierr);
	ierr = DMSetUp(grid.y.da); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(grid.y.da, &grid.y.coords); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(grid.y.da, &grid.y.local); CHKERRQ(ierr);
	if (ctx.nz > 0)
	{
		grid.dim = 3;
		// gridline in z-direction
		ierr = DMDACreate1d(PETSC_COMM_SELF,
		                    DM_BOUNDARY_GHOSTED, ctx.nz, 1, 1, nullptr,
		                    &grid.z.da); CHKERRQ(ierr);
		ierr = DMSetFromOptions(grid.z.da); CHKERRQ(ierr);
		ierr = DMSetUp(grid.z.da); CHKERRQ(ierr);
		ierr = DMCreateGlobalVector(grid.z.da, &grid.z.coords); CHKERRQ(ierr);
		ierr = DMCreateLocalVector(grid.z.da, &grid.z.local); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
} // PetibmGridInitialize


/*! Initializes a grid, providing the coordinates, based on a reference grid.
 *
 * Creates a 1D DMDA objects for each direction and creates the local and global
 * vectors associated with each DMDA.
 * The coordinates are provided as an array of sequential vectors, each vectors
 * containing the stations along a gridline.
 * The decomposition followed the decomposition of a reference grid so that
 * physically closed points are located on the same process.
 *
 * \param other The grid used as a reference.
 * \param coords The stations along each direction.
 * \param grid The grid to initialize (passed by reference).
 */
PetscErrorCode PetibmGridInitialize(
	const PetibmGrid other, const Vec coords[], PetibmGrid &grid)
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	ierr = PetibmGridlineInitialize(other.x, coords[0], grid.x); CHKERRQ(ierr);
	ierr = PetibmGridlineInitialize(other.y, coords[1], grid.y); CHKERRQ(ierr);
	if (other.dim == 3)
	{
		grid.dim = 3;
		ierr = PetibmGridlineInitialize(other.z, coords[2], grid.z); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
} // PetibmGridInitialize


/*! Initializes a gridline, providing the coordinates, based on a reference
 * gridline.
 *
 * Creates a 1D DMDA object for the direction and creates the local and global
 * vectors associated with it.
 * The coordinates are provided as a sequential vector containing the stations
 * along the gridline.
 * The decomposition followed the decomposition of a reference gridline so that
 * physically closed points are located on the same process.
 *
 * \param other The gridline used as a reference.
 * \param coords The stations along the direction.
 * \param line The gridline to initialize (passed by reference).
 */
PetscErrorCode PetibmGridlineInitialize(
	PetibmGridline other, const Vec coords, PetibmGridline &line)
{
	PetscErrorCode ierr;
	PetscInt *lx;
	PetscInt M, m;
	DMBoundaryType b;
	DMDALocalInfo info;
	PetscMPIInt rank, size;

	PetscFunctionBeginUser;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
	ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

	// get number of coordinates
	ierr = VecGetSize(coords, &M); CHKERRQ(ierr);

	// get info from other (model) 1D DMDA
	ierr = DMDAGetInfo(other.da,
	                   nullptr,
	                   nullptr, nullptr, nullptr,
	                   &m, nullptr, nullptr,
	                   nullptr, nullptr,
	                   &b, nullptr, nullptr,
	                   nullptr); CHKERRQ(ierr);
	ierr = PetscMalloc(m*sizeof(*lx), &lx); CHKERRQ(ierr);
	if (m == 1)
	{
		lx[0] = M;
		ierr = DMDACreate1d(
			PETSC_COMM_SELF, b, M, 1, 1, nullptr, &line.da); CHKERRQ(ierr);
		ierr = DMSetFromOptions(line.da); CHKERRQ(ierr);
		ierr = DMSetUp(line.da); CHKERRQ(ierr);
	}
	else if (m == size)
	{	
		ierr = DMDAGetLocalInfo(other.da, &info); CHKERRQ(ierr);
		PetscReal *arr;
		ierr = PetibmGridlineGlobalToLocal(other); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(other.da, other.local, &arr); CHKERRQ(ierr);
		PetscReal start = arr[info.gxs],
		          end = arr[info.xs+info.xm-1];
		if (info.xs+info.xm == info.mx)
			end = arr[info.xs+info.xm];
		ierr = PetibmGetNumPoints1D(coords, start, end, lx+rank); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(other.da, other.coords, &arr); CHKERRQ(ierr);
		ierr = MPI_Allgather(MPI_IN_PLACE, 1, MPIU_INT,
		                     lx, 1, MPIU_INT, PETSC_COMM_WORLD); CHKERRQ(ierr);
		ierr = DMDACreate1d(
			PETSC_COMM_WORLD, b, M, 1, 1, lx, &line.da); CHKERRQ(ierr);
		ierr = DMSetFromOptions(line.da); CHKERRQ(ierr);
		ierr = DMSetUp(line.da); CHKERRQ(ierr);
	}
	else
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,
		        "Function only supports with 1 or comm's size");
	ierr = PetscFree(lx); CHKERRQ(ierr);

	ierr = DMCreateGlobalVector(line.da, &line.coords); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(line.da, &line.local); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmGridlineInitialize


/*! Sets the stations at external boundary points for all directions.
 *
 * \param starts List of starting points.
 * \param ends List of ending points.
 * \param grid The grid to modify (passed by reference).
 */
PetscErrorCode PetibmGridSetBoundaryPoints(
	const PetscReal starts[], const PetscReal ends[],PetibmGrid &grid)
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	ierr = PetibmGridlineSetBoundaryPoints(
		starts[0], ends[0], grid.x); CHKERRQ(ierr);
	ierr = PetibmGridlineSetBoundaryPoints(
		starts[1], ends[1], grid.y); CHKERRQ(ierr);
	if (grid.dim == 3)
	{
		ierr = PetibmGridlineSetBoundaryPoints(
			starts[2], ends[2], grid.z); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
} // PetibmGridSetBoundaryPoints


/*! Sets the starting end ending points for a gridline.
 *
 * \param start The starting point.
 * \param end The ending point.
 * \param line The gridline to modify (passed by reference).
 */
PetscErrorCode PetibmGridlineSetBoundaryPoints(
	const PetscReal start, const PetscReal end, PetibmGridline &line)
{
	PetscErrorCode ierr;
	DMDALocalInfo info;
	PetscReal *arr;

	PetscFunctionBeginUser;

	ierr = DMDAGetLocalInfo(line.da, &info); CHKERRQ(ierr);
	ierr = PetibmGridlineGlobalToLocal(line); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(line.da, line.local, &arr); CHKERRQ(ierr);
	if (info.xs == 0)
		arr[-1] = start;
	if (info.xs+info.xm == info.mx)
		arr[info.mx] = end;
	ierr = DMDAVecRestoreArray(line.da, line.local, &arr); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmGridlineSetBoundaryPoints


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
	const PetibmGrid gridA, const PetibmGridCtx ctx, PetibmGrid &gridB)
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	// Crop the gridline along the x-direction
	ierr = PetibmGridlineCrop(
		gridA.x, ctx.starts[0], ctx.ends[0], gridB.x); CHKERRQ(ierr);
	// Crop the gridline along the y-direction
	ierr = PetibmGridlineCrop(
		gridA.y, ctx.starts[1], ctx.ends[1], gridB.y); CHKERRQ(ierr);
	if (gridA.dim == 3)
	{
		gridB.dim = 3;
		// Crop the gridline along the z-direction
		ierr = PetibmGridlineCrop(
			gridA.z, ctx.starts[2], ctx.ends[2], gridB.z); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
} // PetibmGridCrop


/*! Inserts global values into local vectors for the grid.
 *
 * \param grid The grid to work on (passed by reference).
 */
PetscErrorCode PetibmGridGlobalToLocal(PetibmGrid &grid)
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	ierr = PetibmGridlineGlobalToLocal(grid.x); CHKERRQ(ierr);
	ierr = PetibmGridlineGlobalToLocal(grid.y); CHKERRQ(ierr);
	if (grid.dim == 3)
	{
		ierr = PetibmGridlineGlobalToLocal(grid.z); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
} // PetibmGridGlobalToLocal


/*! Inserts global values into local vector for the gridline.
 *
 * \param line The gridline to work on (passed by reference).
 */
PetscErrorCode PetibmGridlineGlobalToLocal(PetibmGridline &line)
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	ierr = DMGlobalToLocalBegin(
		line.da, line.coords, INSERT_VALUES, line.local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(
		line.da, line.coords, INSERT_VALUES, line.local); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmGridlineGlobalToLocal


/*! Destroys a PetibmGrid structure.
 *
 * \param grid The PetibmGrid structure to destroy (passed by reference).
 */
PetscErrorCode PetibmGridDestroy(PetibmGrid &grid)
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	ierr = PetibmGridlineDestroy(grid.x); CHKERRQ(ierr);
	ierr = PetibmGridlineDestroy(grid.x); CHKERRQ(ierr);
	if (grid.dim == 3)
	{
		ierr = PetibmGridlineDestroy(grid.z); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
} // PetibmGridDestroy


/*! Destroys a PetibmGridline structure.
 *
 * \param line The PetibmGridline structure to destroy (passed by reference).
 */
PetscErrorCode PetibmGridlineDestroy(PetibmGridline &line)
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	ierr = VecDestroy(&line.coords); CHKERRQ(ierr);
	ierr = VecDestroy(&line.local);
	ierr = DMDestroy(&line.da); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmGridlineDestroy


/*! Reads the gridlines from file.
 *
 * The gridlines should be stored in HDF5 format in the same file.
 *
 * \param filepath Path of the input file.
 * \param varname Name of variable (group name in the HDF5).
 * \param grid The PetibmGrid object to fill (passed by reference).
 */
PetscErrorCode PetibmGridHDF5Read(
	const std::string filepath, const std::string varname, PetibmGrid &grid)
{
	PetscErrorCode ierr;
	PetscViewer viewer;

	PetscFunctionBeginUser;

	ierr = PetscViewerHDF5Open(
		PETSC_COMM_SELF, filepath.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
#ifndef PETIBM_0_2
	ierr = PetscViewerHDF5PushGroup(viewer, varname.c_str()); CHKERRQ(ierr);
#endif
	// read x-stations
	ierr = PetscObjectSetName((PetscObject) grid.x.coords, "x"); CHKERRQ(ierr);
	ierr = VecLoad(grid.x.coords, viewer); CHKERRQ(ierr);
	// read y-stations
	ierr = PetscObjectSetName((PetscObject) grid.y.coords, "y"); CHKERRQ(ierr);
	ierr = VecLoad(grid.y.coords, viewer); CHKERRQ(ierr);
	if (grid.dim == 3)
	{
		// read z-stations
		ierr = PetscObjectSetName((PetscObject) grid.z.coords, "z"); CHKERRQ(ierr);
		ierr = VecLoad(grid.z.coords, viewer); CHKERRQ(ierr);
	}

	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmGridHDF5Read


/*! Reads the gridline stations from a file.
 *
 * The stations along the gridline should be stored in HDF5 format.
 *
 * \param filepath The path of the input file.
 * \param varname The name of the variable.
 * \param name The name of the gridline (the direction).
 * \param line The sequential vector to fill (passed by reference).
 */
PetscErrorCode PetibmGridlineHDF5Read(
	const std::string filepath, const std::string varname,
	const std::string name, Vec &line)
{
	PetscErrorCode ierr;
	PetscViewer viewer;

	PetscFunctionBeginUser;

	ierr = PetscViewerHDF5Open(
		PETSC_COMM_SELF, filepath.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
#ifndef PETIBM_0_2
	ierr = PetscViewerHDF5PushGroup(viewer, varname.c_str()); CHKERRQ(ierr);
#endif
	ierr = PetscObjectSetName(
		(PetscObject) line, name.c_str()); CHKERRQ(ierr);
	ierr = VecLoad(line, viewer); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmGridlineHDF5Read


/*! Writes the gridlines into file in HDF5 format.
 *
 * \param filepath Path of the output file.
 * \param varname Name of the grid.
 * \param grid The PetibmGrid structure containing the gridlines.
 */
PetscErrorCode PetibmGridHDF5Write(
	const std::string filepath, const std::string varname, const PetibmGrid grid)
{
	PetscErrorCode ierr;
	PetscViewer viewer;
	PetscFileMode mode;

	PetscFunctionBeginUser;

	std::ifstream infile(filepath.c_str());
	if (infile.good())
		mode = FILE_MODE_APPEND;
	else
		mode = FILE_MODE_WRITE;

	ierr = PetscViewerHDF5Open(
		PETSC_COMM_SELF, filepath.c_str(), mode, &viewer); CHKERRQ(ierr);
#ifndef PETIBM_0_2
	ierr = PetscViewerHDF5PushGroup(viewer, varname.c_str()); CHKERRQ(ierr);
#endif
	// write x-stations
	ierr = PetscObjectSetName((PetscObject) grid.x.coords, "x"); CHKERRQ(ierr);
	ierr = VecView(grid.x.coords, viewer); CHKERRQ(ierr);
	// write y-stations
	ierr = PetscObjectSetName((PetscObject) grid.y.coords, "y"); CHKERRQ(ierr);
	ierr = VecView(grid.y.coords, viewer); CHKERRQ(ierr);
	if (grid.dim == 3)
	{
		// write z-stations
		ierr = PetscObjectSetName((PetscObject) grid.z.coords, "z"); CHKERRQ(ierr);
		ierr = VecView(grid.z.coords, viewer); CHKERRQ(ierr);
	}

	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmGridHDF5Write


PetscErrorCode PetibmGridCreateSeq(PetibmGridCtx ctx, PetibmGrid &grid)
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	ierr = VecCreateSeq(
		PETSC_COMM_SELF, ctx.nx, &grid.x.coords); CHKERRQ(ierr);
	ierr = VecCreateSeq(
		PETSC_COMM_SELF, ctx.ny, &grid.y.coords); CHKERRQ(ierr);
	if (ctx.nz > 0)
	{
		grid.dim = 3;
		ierr = VecCreateSeq(
			PETSC_COMM_SELF, ctx.nz, &grid.z.coords); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
} // PetibmGridCreateSeq
