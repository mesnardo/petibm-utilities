/*! Definition of miscellaneous functions.
 * \file misc.hpp
 */

#pragma once

#include <string>

#include <petscsys.h>


/*! Gets the directory from the command-line.
 *
 * \param directory Directory of the numerical solution.
 */
PetscErrorCode PetibmGetDirectory(std::string *directory);
