/*! Definition of the structure PetibmField and related functions.
 * \file field.h
 */

#pragma once

#include <petscsys.h>
#include <petscdmda.h>

#include "petibm-utilities/grid.h"


/*! Structure holding the field (decomposition and vectors).
 */
struct PetibmField
{
	DM da;  /// the decomposition of the field
	Vec global,  /// parallel vector containing the field values
	    local;  /// sequential ghosted vector containing the field values on process.
}; // PetibmField


/*! Structure holding information about the field.
 */
struct PetibmFieldCtx
{
	char path[PETSC_MAX_PATH_LEN];  /// path of the file containing field values
	char name[PETSC_MAX_PATH_LEN];  /// name of the field
	PetscReal bc_value = 0.0;  /// value to set at external boundary points
	PetscBool periodic_x = PETSC_FALSE,  /// is field periodic in the x-direction?
	          periodic_y = PETSC_FALSE,  /// is field periodic in the y-direction?
	          periodic_z = PETSC_FALSE;  /// is field periodic in the z-direction?
}; // PetibmFieldCtx


/*! Gets options from command-line or config file.
 *
 * \param prefix String to prepend to options.
 * \param ctx The PetibmFieldCtx structure to fill (passed by pointer).
 */
PetscErrorCode PetibmFieldGetOptions(
	const char prefix[], PetibmFieldCtx *ctx);


/*! Initializes a PetibmField structure based on the grid.
 *
 * Creates the DMDA object and local and global vectors associated with it.
 * Creates the vectors based on the DMDA object.
 * The decomposition of the field follows the decomposition of the provided grid.
 *
 * \param ctx The context.
 * \param grid The grid used a reference for domain decomposition of the field.
 * \param field The field to initialize (passed by reference).
 */
PetscErrorCode PetibmFieldInitialize(
	const PetibmFieldCtx ctx, const PetibmGrid grid, PetibmField &field);


/*! Initializes a PetibmField structure based on a given DMDA object.
 *
 * Creates the vectors based on the provided DMDA object.
 *
 * \param da The DMDA object.
 * \param field The field to initialize (passed by reference).
 */
PetscErrorCode PetibmFieldInitialize(const DM da, PetibmField &field);


/*! Sets the value at external boundary points.
 *
 * \param value The value on the external boundaries.
 * \param field The field to modify (passed by reference).
 */
PetscErrorCode PetibmFieldSetBoundaryPoints(
	const PetscReal value, PetibmField &field);


/*! Inserts values from global vector into local vector.
 *
 * \param field The field to work on (passed by reference).
 */
PetscErrorCode PetibmFieldGlobalToLocal(PetibmField &field);


/*! Destroys the PETSc objects of a PetibmField structure.
 *
 * \param field The PetibmField structure (passed by reference).
 */
PetscErrorCode PetibmFieldDestroy(PetibmField &field);


/*! Reads the field values stored in given format from file.
 *
 * \param filepath Path of the input file.
 * \param name The name of the field.
 * \param viewerType PETSc viewer type.
 * \param field The PetibmField structure (passed by reference).
 */
PetscErrorCode PetibmFieldRead(const std::string filepath,
                               const std::string name,
                               const PetscViewerType viewerType,
                               PetibmField &field);


/*! Reads the field values stored in HDF5 format from file.
 *
 * \param filepath Path of the input file.
 * \param name The name of the field.
 * \param field The PetibmField structure (passed by reference).
 */
PetscErrorCode PetibmFieldHDF5Read(
	const std::string filepath, const std::string name, PetibmField &field);


/*! Reads the field values stored in binary format from file.
 *
 * \param filepath Path of the input file.
 * \param name The name of the field.
 * \param field The PetibmField structure (passed by reference).
 */
PetscErrorCode PetibmFieldBinaryRead(
	const std::string filepath, const std::string name, PetibmField &field);


/*! Writes the field values into file in HDF5 format.
 *
 * \param filepath Path of the output file.
 * \param name Name of the field.
 * \param field PetibmField structure.
 */
PetscErrorCode PetibmFieldHDF5Write(
	const std::string filepath, const std::string name, const PetibmField field);


/*! Writes the field values into file in binary format.
 *
 * \param filepath Path of the output file.
 * \param name Name of the field.
 * \param field PetibmField structure.
 */
PetscErrorCode PetibmFieldBinaryWrite(
	const std::string filepath, const std::string name, const PetibmField field);


/*! Interpolates field A associated with grid A onto grid B.
 *
 * \param gridA The grid to interpolate from.
 * \param fieldA The field to interpolate.
 * \param gridB The grid to interpolate on.
 * \param fieldB The resulting interpolated field (passed by reference).
 */
PetscErrorCode PetibmFieldInterpolate(
	PetibmGrid gridA, PetibmField fieldA, PetibmGrid gridB, PetibmField &fieldB);


/*! Interpolates a 2D field A associated with grid A onto grid B.
 *
 * Performs bi-linear interpolation.
 *
 * \param gridA The grid to interpolate from.
 * \param fieldA The field to interpolate.
 * \param gridB The grid to interpolate on.
 * \param fieldB The resulting interpolated field (passed by reference).
 */
PetscErrorCode PetibmFieldInterpolate2D(
	PetibmGrid gridA, PetibmField fieldA, PetibmGrid gridB, PetibmField &fieldB);


/*! Interpolates a 3D field A associated with grid A onto grid B.
 *
 * Performs tri-linear interpolation.
 *
 * \param gridA The grid to interpolate from.
 * \param fieldA The field to interpolate.
 * \param gridB The grid to interpolate on.
 * \param fieldB The resulting interpolated field (passed by reference).
 */
PetscErrorCode PetibmFieldInterpolate3D(
	PetibmGrid gridA, PetibmField fieldA, PetibmGrid gridB, PetibmField &fieldB);
