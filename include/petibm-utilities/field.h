/*! Definition of the structure PetibmField.
 * \file field.hpp
 */

#pragma once

#include <petscsys.h>
#include <petscdmda.h>

#include "petibm-utilities/grid.h"


typedef struct
{
	PetibmGrid grid;
	DM da;
	Vec global, local;
} PetibmField;

typedef struct
{
	char path[PETSC_MAX_PATH_LEN];
	char name[PETSC_MAX_PATH_LEN];
	PetscReal bc_value = 0.0;
} PetibmFieldCtx;

/*! Gets options from command-line or config file.
 *
 * \param prefix String to prepend to options.
 * \param ctx The PetibmFieldCtx structure to fill.
 */
PetscErrorCode PetibmFieldGetOptions(
	const char prefix[], PetibmFieldCtx *ctx);

/*! Initializes a PetibmField structure.
 *.
 * \param field The PetibmField structure to initialize (passed by reference).
 */
PetscErrorCode PetibmFieldInitialize(PetibmField &field);

/*! Destroys the PETSc objects of a PetibmField structure.
 *
 * \param field The PetibmField structure (passed by reference).
 */
PetscErrorCode PetibmFieldDestroy(PetibmField &field);

/*! Reads the field values from file.
 *
 * \param filepath Path of the input file.
 * \param name The name of the field.
 * \param field The PetibmField structure (passed by reference).
 */
PetscErrorCode PetibmFieldHDF5Read(
	const std::string filepath, const std::string name, PetibmField &field);

/*! Writes the field values into file.
 *
 * \param filepath Path of the output file.
 * \param name Name of the field.
 * \param field PetibmField structure.
 */
PetscErrorCode PetibmFieldHDF5Write(
	const std::string filepath, const std::string name, const PetibmField field);

/*! Interpolates 2D field values from one mesh to another.
 *
 * \param fieldA PetibmField to be interpolated.
 * \param fieldB PetibmField on which the solution is interpolated.
 * \param bc_value Value to use near boundary when no neighbor is found.
 */
PetscErrorCode PetibmFieldInterpolate2D(
	const PetibmField fieldA, PetibmField &fieldB, const PetscReal bc_value);

/*! Interpolates 3D field values from one mesh to another.
 *
 * \param fieldA PetibmField to be interpolated.
 * \param fieldB PetibmField on which the solution is interpolated.
 * \param bc_value Value to use near boundary when no neighbor is found.
 */
PetscErrorCode PetibmFieldInterpolate3D(
	const PetibmField fieldA, PetibmField &fieldB, const PetscReal bc_value);

/*! Sets the value at external ghost points.
 *
 * \param field The PetibmField to update.
 * \param value Value at the external ghost points.
 */
PetscErrorCode  PetibmFieldExternalGhostPointsSet(
	PetibmField field, PetscReal value);
