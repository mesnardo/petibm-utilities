/*! Implementation of miscellaneous functions.
 * \file misc.cpp
 */

#include "petibm-utilities/misc.h"


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


/*! Finds index of closest point on left side.
 *
 * \param x_i Point for which we look the neighbor.
 * \param x Vec that contains all neighbors.
 * \param index Index of the closest neighbor.
 * \param found Set to PETSC_TRUE is closest neighbor found.
 */
PetscErrorCode PetibmGetNeighborIndex1D(
	const PetscReal x_i, const Vec x, PetscInt *index, PetscBool *found)
{
	PetscErrorCode ierr;
	PetscReal *x_a;
	PetscInt i, i_start = *index;
	PetscInt n;

	PetscFunctionBeginUser;

	ierr = VecGetLocalSize(x, &n); CHKERRQ(ierr);
	ierr = VecGetArray(x, &x_a); CHKERRQ(ierr);
	for (i=i_start; i<n-1; i++)
	{
		if (x_a[i] <= x_i and x_i < x_a[i+1])
		{
			*index = i;
			*found = PETSC_TRUE;
			break;
		}
	}
	ierr = VecRestoreArray(x, &x_a); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmGetNeighborIndex1D
