/*! Definition of the structure and functions related to the mesh.
 * \file mesh.hpp
 */

#pragma once


typedef struct
{
	char path[PETSC_MAX_PATH_LEN];
  PetscInt nx = 0,
           ny = 0,
           nz = 0;
  PetscBool periodic_x = PETSC_FALSE,
            periodic_y = PETSC_FALSE,
            periodic_z = PETSC_FALSE;
} PetibmMeshInfo;

/*! Gets options from command-line or config file.
 *
 * \param prefix String to prepend the name of the options.
 * \param info The MeshInfo to fill.
 */
PetscErrorCode PetibmMeshInfoGetOptions(
	const char prefix[], PetibmMeshInfo *info);
