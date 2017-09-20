/*! Definition of miscellaneous functions.
 * \file misc.hpp
 */

#pragma once

#include <string>

#include <petscsys.h>
#include <petscvec.h>


/*! Gets the directory from the command-line.
 *
 * \param directory Directory of the numerical solution.
 */
PetscErrorCode PetibmGetDirectory(std::string *directory);

/*! Finds index of closest point on left side.
 *
 * \param x_i Point for which we look the neighbor.
 * \param x Vec that contains all neighbors.
 * \param index Index of the closest neighbor.
 * \param found Set to PETSC_TRUE is closest neighbor found.
 */
PetscErrorCode PetibmGetNeighborIndex1D(
	const PetscReal x_i, const Vec x, PetscInt *index, PetscBool *found);
