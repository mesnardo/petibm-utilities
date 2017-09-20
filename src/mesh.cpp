/*! Definition of the structure and functions related to the mesh.
 * \file mesh.cpp
 */

#include <petscsys.h>

#include "petibm-utilities/mesh.h"


/*! Gets options from command-line or config file.
 *
 * \param prefix String to prepend the name of the options.
 * \param info The PetibmMeshInfo to fill.
 */
PetscErrorCode PetibmMeshInfoGetOptions(
	const char prefix[], PetibmMeshInfo *info)
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
	ierr = PetscOptionsGetString(nullptr, prefix, "-grid_path", info->path,
	                             sizeof(info->path), &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(
		nullptr, prefix, "-nx", &info->nx, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(
		nullptr, prefix, "-ny", &info->ny, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(
		nullptr, prefix, "-nz", &info->nz, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(nullptr, prefix, "-periodic_x",
	                           &info->periodic_x, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(nullptr, prefix, "-periodic_y",
	                           &info->periodic_y, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(nullptr, prefix, "-periodic_z",
	                           &info->periodic_z, &found); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmMeshInfoGetOptions
