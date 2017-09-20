/*! Definition of the structure PetibmField.
 * \file field.hpp
 */

#pragma once

#include <petscsys.h>
#include <petscdmda.h>


typedef struct
{
  DM da;
  Vec x, y, z;
  Vec global, local;
} PetibmField;

typedef struct
{
	char path[PETSC_MAX_PATH_LEN];
	char name[PETSC_MAX_PATH_LEN];
	PetscReal bcValue = 0.0;
} PetibmFieldInfo;

/*! Initializes a PetibmField structure.
 *
 * Creates the vectors based on the DMDA object.
 *
 * \param field The PetibmField structure to initialize (passed by reference).
 */
PetscErrorCode PetibmFieldInitialize(PetibmField &field);

/*! Destroys the PETSc objects of a PetibmField structure.
 *
 * \param field The PetibmField structure (passed by reference).
 */
PetscErrorCode PetibmFieldDestroy(PetibmField &field);

/*! Reads the field values from a HDF5 file.
 *
 * \param filePath Path of the file to read.
 * \param name The name of the field.
 * \param field The PetibmField structure (passed by reference).
 */
PetscErrorCode PetibmFieldReadValues(
	std::string filePath, std::string name, PetibmField &field);

/*! Reads a HDF5 grid file for a given field.
 *
 * \param filePath Path of the grid file to read.
 * \param field The PetibmField structure (passed by reference).
 */
PetscErrorCode PetibmFieldReadGrid(std::string filePath, PetibmField &field);

/*! Writes a grid into a HDF5 file.
 *
 * \param filePath Path of the file to write into.
 * \param field PetibmField structure that contains the grid.
 */
PetscErrorCode PetibmFieldWriteGrid(std::string filePath, PetibmField field);

/*! Writes the field values into a HDF5 file.
 *
 * \param filePath Path of the file to write into.
 * \param name Name of the field.
 * \param field PetibmField structure.
 */
PetscErrorCode PetibmFieldWriteValues(
	std::string filePath, std::string name, PetibmField field);

/*! Interpolates 2D field values from one mesh to another.
 *
 * \param fieldA PetibmField to be interpolated.
 * \param fieldB PetibmField on which the solution is interpolated.
 */
PetscErrorCode PetibmFieldInterpolate2D(
	const PetibmField fieldA, PetibmField &fieldB);

/*! Interpolates 3D field values from one mesh to another.
 *
 * \param fieldA PetibmField to be interpolated.
 * \param fieldB PetibmField on which the solution is interpolated.
 */
PetscErrorCode PetibmFieldInterpolate3D(
	const PetibmField fieldA, PetibmField &fieldB);

/*! Sets the value at external ghost points.
 *
 * \param field The PetibmField to update.
 * \param value Value at the external ghost points.
 */
PetscErrorCode  PetibmFieldExternalGhostPointsSet(
	PetibmField field, PetscReal value);

/*! Gets options from command-line or config file.
 *
 * \param prefix String to prepend to options.
 * \param info The PetibmFieldInfo structure to fill.
 */
PetscErrorCode PetibmFieldInfoGetOptions(
	const char prefix[], PetibmFieldInfo *info);
