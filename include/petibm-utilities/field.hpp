/*! Definition of the structure Field.
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
} Field;

/*! Initializes a Field structure.
 *
 * Creates the vectors based on the DMDA object.
 *
 * \param field The Field structure to initialize (passed by reference).
 */
PetscErrorCode FieldInitialize(Field &field);

/*! Destroys the PETSc objects of a Field structure.
 *
 * \param field The Field structure (passed by reference).
 */
PetscErrorCode FieldDestroy(Field &field);

/*! Reads the field values from a HDF5 file.
 *
 * \param filePath Path of the file to read.
 * \param name The name of the field.
 * \param field The field structure (passed by reference).
 */
PetscErrorCode FieldReadValues(
	std::string filePath, std::string name, Field &field);

/*! Reads a HDF5 grid file for a given field.
 *
 * \param filePath Path of the grid file to read.
 * \param field The field structure (passed by reference).
 */
PetscErrorCode FieldReadGrid(std::string filePath, Field &field);

/*! Writes a grid into a HDF5 file.
 *
 * \param filePath Path of the file to write into.
 * \param field Field structure that contains the grid.
 */
PetscErrorCode FieldWriteGrid(std::string filePath, Field field);

/*! Writes the field values into a HDF5 file.
 *
 * \param filePath Path of the file to write into.
 * \param name Name of the field.
 * \param field Field structure.
 */
PetscErrorCode FieldWriteValues(
	std::string filePath, std::string name, Field field);
