/*! Definition of the structure and functions related to the grid.
 * \file grid.h
 */

#pragma once

#include <string>

#include <petscsys.h>
#include <petscvec.h>


typedef struct
{
	PetscInt dim = 2;
	Vec x, y, z;
} PetibmGrid;

typedef struct
{
	char path[PETSC_MAX_PATH_LEN];
	PetscInt nx = 0,
	         ny = 0,
	         nz = 0;
	PetscBool periodic_x = PETSC_FALSE,
	          periodic_y = PETSC_FALSE,
	          periodic_z = PETSC_FALSE;
} PetibmGridCtx;

/*! Gets options from command-line or config file.
 *
 * \param prefix String to prepend the name of the options.
 * \param grid The PetibmGrid object to fill.
 */
PetscErrorCode PetibmGridGetOptions(
	const char prefix[], PetibmGridCtx *ctx);

/*! Initializes a PetibmGrid structure.
 *
 * \param grid The PetibmGrid structure to initialize (passed by reference).
 */
PetscErrorCode PetibmGridInitialize(const PetibmGridCtx ctx, PetibmGrid &grid);

/*! Destroys a PetibmGrid structure.
 *
 * \param grid The PetibmGrid structure to destroy (passed by reference).
 */
PetscErrorCode PetibmGridDestroy(PetibmGrid &grid);

/*! Reads the gridlines from file.
 *
 * \param filepath Path of the input file.
 * \param grid The PetibmGrid object to fill (passed by reference).
 */
PetscErrorCode PetibmGridHDF5Read(
	const std::string filepath, PetibmGrid &grid);

/*! Writes the gridlines into file.
 *
 * \param filepath Path of the output file.
 * \param grid The PetibmGrid object containing the gridlines.
 */
PetscErrorCode PetibmGridHDF5Write(
	const std::string filepath, const PetibmGrid grid);
