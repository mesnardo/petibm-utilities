/*! Definition of miscellaneous functions.
 * \file misc.h
 */

#pragma once

#include <string>
#include <vector>

#include <petscsys.h>
#include <petscvec.h>

#include "petibm-utilities/grid.h"


/*! Gets the directory from the command-line.
 *
 * If the command-line parameter `-directory <path>` is not found, the present
 * working directory is used.
 *
 * \param directory Directory of the numerical solution (passed by pointer).
 */
PetscErrorCode PetibmGetDirectory(
	std::string *directory, const char key[], const PetscBool create=PETSC_FALSE);


PetscErrorCode PetibmGetFilePath(std::string *filepath, const char key[]);


/*! Loads options from configuration file.
 *
 * \param prefix Prefix of "-config_file".
 */
PetscErrorCode PetibmOptionsInsertFile(const char prefix[]);


/*! Counts the number of points within provided boundaries.
 *
 * \param x The sequential vector containing the points.
 * \param start The starting boundary.
 * \param end The ending boundary.
 * \param num The number of points (passed by pointer).
 */
PetscErrorCode PetibmGetNumPoints1D(
	const Vec x, const PetscReal start, const PetscReal end, PetscInt *num);


/*! Get the global index of the inferior neighbor from a gridline for each
 * station along another gridline.
 *
 * \param lineB The gridline where each station will be associated with a neighbor.
 * \param lineA The reference gridline where to find neighbors.
 * \param Iv Standard vector that will contains the indices (passed by reference).
 */
PetscErrorCode PetibmGetNeighbors1D(
	PetibmGridline lineB, PetibmGridline lineA, std::vector<PetscInt> &Iv);


/*! Helper function to print on-process vector in a sequential manner.
 *
 * \param v The vector to print.
 */
PetscErrorCode PetibmVecView(const Vec v);


/*! Helper function to print on-process array of integers in a sequential manner.
 *
 * \param N The number of elements.
 * \param v The array of integers to print.
 */
PetscErrorCode PetibmIntView(const PetscInt N, const PetscInt v[]);
