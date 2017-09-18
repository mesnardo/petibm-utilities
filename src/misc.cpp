/*! Implementation of miscellaneous functions.
 * \file misc.cpp
 */

#include "petibm-utilities/misc.hpp"


/*! Gets the directory from the command-line.
 *
 * \param directory Directory of the numerical solution.
 */
PetscErrorCode PetibmGetDirectory(std::string *directory)
{
	PetscErrorCode ierr;
	char dir[PETSC_MAX_PATH_LEN];
	PetscBool found;

	PetscFunctionBeginUser;

	ierr = PetscOptionsGetString(
		nullptr, nullptr, "-directory", dir, sizeof(dir), &found); CHKERRQ(ierr);
	*directory = (!found) ? "." : dir;

	PetscFunctionReturn(0);
} // PetibmGetDirectory
