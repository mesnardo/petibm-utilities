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


/*! Gets options from command-line or configuration file.
 *
 * \param prefix [in] String to prepend to options.
 * \param ctx [out] The PetibmFieldCtx structure to fill (passed by pointer).
 */
PetscErrorCode PetibmFieldGetOptions(
	const char prefix[], PetibmFieldCtx *ctx);


/*! Prints field parameters.
 *
 * \param name [in] Name of the field.
 * \param ctx [in] Parameters of the field.
 */
PetscErrorCode PetibmFieldCtxPrintf(
	const std::string name, const PetibmFieldCtx ctx);


/*! Initializes a PetibmField structure based on the grid.
 *
 * Creates the DMDA object and local and global vectors associated with it.
 * Creates the vectors based on the DMDA object.
 * The decomposition of the field follows the decomposition of the provided grid.
 *
 * \param ctx [in] The context.
 * \param grid [in] The grid used a reference for domain decomposition of the field.
 * \param field [out] The field to initialize (passed by reference).
 */
PetscErrorCode PetibmFieldInitialize(
	const PetibmFieldCtx ctx, const PetibmGrid grid, PetibmField &field);


/*! Initializes a PetibmField structure based on a given DMDA object.
 *
 * Creates the vectors based on the provided DMDA object.
 *
 * \param da [in] The DMDA object.
 * \param field [out] The field to initialize (passed by reference).
 */
PetscErrorCode PetibmFieldInitialize(const DM da, PetibmField &field);


/*! Sets the value at external boundary points.
 *
 * \param value [in] The value on the external boundaries.
 * \param field [out] The field to modify (passed by reference).
 */
PetscErrorCode PetibmFieldSetBoundaryPoints(
	const PetscReal value, PetibmField &field);


/*! Crops a 2D PetibmField given the starting and ending indices and returns
 *  the sub-field.
 *
 * \param fieldA [in] The PetibmField to crop.
 * \param I_start [in] Starting index in the x-direction.
 * \param I_end [in] Ending index in the x-direction.
 * \param J_start [in] Starting index in the y-direction.
 * \param J_end [in] Ending index in the y-direction.
 * \param fieldB [out] The resulting sub PetibmField (passed by reference).
 */
PetscErrorCode PetibmFieldCrop2d(
	const PetibmField fieldA,
	const PetscInt I_start, const PetscInt I_end,
	const PetscInt J_start, const PetscInt J_end,
	PetibmField &fieldB);


/*! Crops a 3D PetibmField given the starting and ending indices and returns
 *  the sub-field.
 *
 * \param fieldA [in] The PetibmField to crop.
 * \param I_start [in] Starting index in the x-direction.
 * \param I_end [in] Ending index in the x-direction.
 * \param J_start [in] Starting index in the y-direction.
 * \param J_end [in] Ending index in the y-direction.
 * \param K_start [in] Starting index in the z-direction.
 * \param K_end [in] Ending index in the z-direction.
 * \param fieldB [out] The resulting sub PetibmField (passed by reference).
 */
PetscErrorCode PetibmFieldCrop3d(
	const PetibmField fieldA,
	const PetscInt I_start, const PetscInt I_end,
	const PetscInt J_start, const PetscInt J_end,
	const PetscInt K_start, const PetscInt K_end,
	PetibmField &fieldB);


/*! Crops a PetibmField and returns the sub-field.
 *
 * \param gridA [in] The PetibmGrid of the PetibmField to crop.
 * \param fieldA [in] The PetibmField to crop.
 * \param ctx [in] The context with parameters to crop.
 * \param fieldB [out] The resulting sub PetibmField (passed by reference).
 */
PetscErrorCode PetibmFieldCrop(
	const PetibmGrid gridA, const PetibmField fieldA,
	const PetibmGridCtx ctx, PetibmField &fieldB);


/*! Inserts values from global vector into local vector.
 *
 * \param field [out] The field to work on (passed by reference).
 */
PetscErrorCode PetibmFieldGlobalToLocal(PetibmField &field);


/*! Destroys the PETSc objects of a PetibmField structure.
 *
 * \param field [out] The PetibmField structure (passed by reference).
 */
PetscErrorCode PetibmFieldDestroy(PetibmField &field);


/*! Reads the field values stored in given format from file.
 *
 * \param comm [in] MPI communicator.
 * \param filepath [in] Path of the input file.
 * \param name [in] The name of the field.
 * \param viewerType [in] PETSc viewer type.
 * \param field [out] The PetibmField structure (passed by reference).
 */
PetscErrorCode PetibmFieldRead(const MPI_Comm comm,
                               const std::string filepath,
                               const std::string name,
                               const PetscViewerType viewerType,
                               PetibmField &field);


/*! Reads the field values stored in HDF5 format from file.
 *
 * \param comm [in] MPI communicator.
 * \param filepath [in] Path of the input file.
 * \param name [in] The name of the field.
 * \param field [out] The PetibmField structure (passed by reference).
 */
PetscErrorCode PetibmFieldHDF5Read(
	const MPI_Comm comm, const std::string filepath, const std::string name,
	PetibmField &field);


/*! Reads the field values stored in binary format from file.
 *
 * \param comm [in] MPI communicator.
 * \param filepath [in] Path of the input file.
 * \param name [in] The name of the field.
 * \param field [out] The PetibmField structure (passed by reference).
 */
PetscErrorCode PetibmFieldBinaryRead(
	const MPI_Comm comm, const std::string filepath, const std::string name,
	PetibmField &field);


/*! Writes the field values to file in given format.
 *
 * \param comm [in] MPI communicator.
 * \param filepath [in] Path of the output file.
 * \param name [in] The name of the field.
 * \param viewerType [in] PETSc viewer type.
 * \param field [in] The PetibmField structure.
 */
PetscErrorCode PetibmFieldWrite(const MPI_Comm comm,
                                const std::string filepath,
                                const std::string name,
                                const PetscViewerType viewerType,
                                const PetibmField field);


/*! Writes the field values into file in HDF5 format.
 *
 * \param comm [in] MPI communicator.
 * \param filepath [in] Path of the output file.
 * \param name [in] Name of the field.
 * \param field [in] PetibmField structure.
 */
PetscErrorCode PetibmFieldHDF5Write(
	const MPI_Comm comm, const std::string filepath, const std::string name,
	const PetibmField field);


/*! Writes the field values into file in binary format.
 *
 * \param comm [in] MPI communicator.
 * \param filepath [in] Path of the output file.
 * \param name [in] Name of the field.
 * \param field [in] PetibmField structure.
 */
PetscErrorCode PetibmFieldBinaryWrite(
	const MPI_Comm comm, const std::string filepath, const std::string name,
	const PetibmField field);


/*! Interpolates field A associated with grid A onto grid B.
 *
 * \param gridA [in] The grid to interpolate from.
 * \param fieldA [in] The field to interpolate.
 * \param gridB [in] The grid to interpolate on.
 * \param fieldB [out] The resulting interpolated field (passed by reference).
 */
PetscErrorCode PetibmFieldInterpolate(
	PetibmGrid gridA, PetibmField fieldA, PetibmGrid gridB, PetibmField &fieldB);


/*! Interpolates a 2D field A associated with grid A onto grid B.
 *
 * Performs bi-linear interpolation.
 *
 * \param gridA [in] The grid to interpolate from.
 * \param fieldA [in] The field to interpolate.
 * \param gridB [in] The grid to interpolate on.
 * \param fieldB [out] The resulting interpolated field (passed by reference).
 */
PetscErrorCode PetibmFieldInterpolate2d(
	PetibmGrid gridA, PetibmField fieldA, PetibmGrid gridB, PetibmField &fieldB);


/*! Interpolates a 3D field A associated with grid A onto grid B.
 *
 * Performs tri-linear interpolation.
 *
 * \param gridA [in] The grid to interpolate from.
 * \param fieldA [in] The field to interpolate.
 * \param gridB [in] The grid to interpolate on.
 * \param fieldB [out] The resulting interpolated field (passed by reference).
 */
PetscErrorCode PetibmFieldInterpolate3d(
	PetibmGrid gridA, PetibmField fieldA, PetibmGrid gridB, PetibmField &fieldB);


/*! Creates a 2D DMDA object for a PetibmField.
 *
 * \param gridCtx [in] The parameters of the grid.
 * \param fieldCtx [in] The parameters of the field.
 * \param da [out] The PETSc DMDA object (passed by reference).
 */
PetscErrorCode PetibmFieldDMDACreate2d(
	const PetibmGridCtx gridCtx, const PetibmFieldCtx fieldCtx, DM &da);


/*! Creates a 3D DMDA object for a PetibmField.
 *
 * \param gridCtx [in] The parameters of the grid.
 * \param fieldCtx [in] The parameters of the field.
 * \param da [out] The PETSc DMDA object (passed by reference).
 */
PetscErrorCode PetibmFieldDMDACreate3d(
	const PetibmGridCtx gridCtx, const PetibmFieldCtx fieldCtx, DM &da);


/*! Creates a DMDA object for a PetibmField.
 *
 * \param gridCtx [in] The parameters of the grid.
 * \param fieldCtx [in] The parameters of the field.
 * \param da [out] The PETSc DMDA object (passed by reference).
 */
PetscErrorCode PetibmFieldDMDACreate(
	const PetibmGridCtx gridCtx, const PetibmFieldCtx fieldCtx, DM &da);


/*! Creates a 2D DMDA object for a PetibmField based on a model 2D DMDA object.
 *
 * The resulting DMDA object follows the same domain decomposition
 * than the model DMDA object.
 *
 * \param name [in] Name of the PetibmField.
 * \param da_in [in] The model 2D DMDA object.
 * \param da [out] The resulting DMDA object (passed by reference).
 */
PetscErrorCode PetibmFieldDMDACreate2d(
	const std::string name, const DM da_in, DM &da);


/*! Creates a 3D DMDA object for a PetibmField based on a model 3D DMDA object.
 *
 * The resulting DMDA object follows the same domain decomposition
 * than the model DMDA object.
 *
 * \param name [in] Name of the PetibmField.
 * \param da_in [in] The model 3D DMDA object.
 * \param da [out] The resulting DMDA object (passed by reference).
 */
PetscErrorCode PetibmFieldDMDACreate3d(
	const std::string name, const DM da_in, DM &da);
