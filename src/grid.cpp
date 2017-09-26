/*! Implementation of the functions for the PetibmGrid structure.
 * \file grid.cpp
 */

#include <petscviewerhdf5.h>

#include "petibm-utilities/grid.h"


/*! Gets options from command-line or config file.
 *
 * \param prefix String to prepend the name of the options.
 * \param ctx The PetibmGridCtx structure to fill.
 */
PetscErrorCode PetibmGridGetOptions(
	const char prefix[], PetibmGridCtx *ctx)
{
	PetscErrorCode ierr;
	char path[PETSC_MAX_PATH_LEN];
	PetscBool found;

	PetscFunctionBeginUser;

	ierr = PetscOptionsGetString(nullptr, prefix, "-config_file",
	                             path, sizeof(path), &found); CHKERRQ(ierr);
	if (found)
	{
		ierr = PetscOptionsInsertFile(
			PETSC_COMM_WORLD, nullptr, path, PETSC_FALSE); CHKERRQ(ierr);
	}
	ierr = PetscOptionsGetString(nullptr, prefix, "-grid_path", ctx->path,
	                             sizeof(ctx->path), &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(
		nullptr, prefix, "-nx", &ctx->nx, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(
		nullptr, prefix, "-ny", &ctx->ny, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(
		nullptr, prefix, "-nz", &ctx->nz, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(nullptr, prefix, "-periodic_x",
	                           &ctx->periodic_x, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(nullptr, prefix, "-periodic_y",
	                           &ctx->periodic_y, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(nullptr, prefix, "-periodic_z",
	                           &ctx->periodic_z, &found); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmGridGetOptions


/*! Initializes a PetibmGrid structure.
 *
 * \param grid The PetibmGrid structure to initialize (passed by reference).
 */
PetscErrorCode PetibmGridInitialize(const PetibmGridCtx ctx, PetibmGrid &grid)
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	ierr = VecCreateSeq(PETSC_COMM_SELF, ctx.nx, &grid.x); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, ctx.ny, &grid.y); CHKERRQ(ierr);
  if (ctx.nz > 0)
  {
  	grid.dim = 3;
  	ierr = VecCreateSeq(PETSC_COMM_SELF, ctx.nz, &grid.z); CHKERRQ(ierr);
  }

	PetscFunctionReturn(0);
} // PetibmGridInitialize


/*! Destroys a PetibmGrid structure.
 *
 * \param grid The PetibmGrid structure to destroy (passed by reference).
 */
PetscErrorCode PetibmGridDestroy(PetibmGrid &grid)
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	ierr = VecDestroy(&grid.x); CHKERRQ(ierr);
	ierr = VecDestroy(&grid.y); CHKERRQ(ierr);
	if (grid.dim == 3)
	{
		ierr = VecDestroy(&grid.z); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
} // PetibmGridDestroy


/*! Reads the gridlines from file.
 *
 * \param filepath Path of the input file.
 * \param grid The PetibmGrid object to fill (passed by reference).
 */
PetscErrorCode PetibmGridHDF5Read(
	const std::string filepath, PetibmGrid &grid)
{
	PetscErrorCode ierr;
	PetscViewer viewer;

	PetscFunctionBeginUser;

	ierr = PetscViewerHDF5Open(
		PETSC_COMM_SELF, filepath.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
	// read x-stations
	ierr = PetscObjectSetName((PetscObject) grid.x, "x"); CHKERRQ(ierr);
	ierr = VecLoad(grid.x, viewer); CHKERRQ(ierr);
	// read y-stations
	ierr = PetscObjectSetName((PetscObject) grid.y, "y"); CHKERRQ(ierr);
	ierr = VecLoad(grid.y, viewer); CHKERRQ(ierr);
	if (grid.dim == 3)
	{
		// read z-stations
		ierr = PetscObjectSetName((PetscObject) grid.z, "z"); CHKERRQ(ierr);
		ierr = VecLoad(grid.z, viewer); CHKERRQ(ierr);
	}

	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmGridHDF5Read


/*! Writes the gridlines into file.
 *
 * \param filepath Path of the output file.
 * \param grid The PetibmGrid object containing the gridlines.
 */
PetscErrorCode PetibmGridHDF5Write(
	const std::string filepath, const PetibmGrid grid)
{
	PetscErrorCode ierr;
	PetscViewer viewer;
	PetscMPIInt rank;

	PetscFunctionBeginUser;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

	if (rank == 0)
	{
		ierr = PetscViewerHDF5Open(
			PETSC_COMM_SELF, filepath.c_str(), FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
		// write x-stations
		ierr = PetscObjectSetName((PetscObject) grid.x, "x"); CHKERRQ(ierr);
		ierr = VecView(grid.x, viewer); CHKERRQ(ierr);
		// write y-stations
		ierr = PetscObjectSetName((PetscObject) grid.y, "y"); CHKERRQ(ierr);
		ierr = VecView(grid.y, viewer); CHKERRQ(ierr);
		if (grid.dim == 3)
		{
			// write z-stations
			ierr = PetscObjectSetName((PetscObject) grid.z, "z"); CHKERRQ(ierr);
			ierr = VecView(grid.z, viewer); CHKERRQ(ierr);
		}

		ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
} // PetibmGridHDF5Write
