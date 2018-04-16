/*! Definition of miscellaneous functions.
 * \file misc.h
 */

#pragma once

#include <string>
#include <vector>

#include <petscsys.h>
#include <petscvec.h>

#include "petibm-utilities/grid.h"


/*! Gets a directory from the command-line arguments.
 *
 * \param directory [out] Directory (passed by pointer).
 * \param key [in] The key option one is seeking.
 * \param defval [in] Default value to use if not found.
 * \param create [in] Create the directory if PETSC_TRUE (default if PETSC_FALSE).
 */
PetscErrorCode PetibmGetDirectory(
	std::string *directory, const char key[],
	const std::string defval, const PetscBool create=PETSC_FALSE);


/*! Gets a file path from the command-line arguments.
 *
 * \param fielpath [out] File path (passed by pointer).
 * \param key [in] The key option one is seeking.
 * \param defval [in] Default value to use if not found.
 */
PetscErrorCode PetibmGetFilePath(
	std::string *filepath, const char key[], const std::string defval);


/*! Loads options from configuration file.
 *
 * \param prefix [in] Prefix of "-config_file".
 */
PetscErrorCode PetibmOptionsInsertFile(const char prefix[]);


/*! Counts the number of points within provided boundaries.
 *
 * \param x [in] The sequential vector containing the points.
 * \param start [in] The starting boundary.
 * \param end [in] The ending boundary.
 * \param num [out] The number of points (passed by pointer).
 */
PetscErrorCode PetibmGetNumPoints1D(
	const Vec x, const PetscReal start, const PetscReal end, PetscInt *num);


/*! Get the global index of the inferior neighbor from a gridline for each
 * station along another gridline.
 *
 * \param lineB [in] The gridline where each station will be associated with a neighbor.
 * \param lineA [in] The reference gridline where to find neighbors.
 * \param Iv [out] Standard vector that will contains the indices (passed by reference).
 */
PetscErrorCode PetibmGetNeighbors1D(
	PetibmGridline lineB, PetibmGridline lineA, std::vector<PetscInt> &Iv);


/*! Helper function to print on-process vector in a sequential manner.
 *
 * \param v [in] The vector to print.
 */
PetscErrorCode PetibmVecView(const Vec v);


/*! Helper function to print on-process array of integers in a sequential manner.
 *
 * \param N [in] The number of elements.
 * \param v [in] The array of integers to print.
 */
PetscErrorCode PetibmIntView(const PetscInt N, const PetscInt v[]);
