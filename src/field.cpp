/*! Implementation of the functions for the PetibmField structure.
 * \file field.cpp
 */

#include <petscviewerhdf5.h>

#include "petibm-utilities/field.hpp"


/*! Initializes a PetibmField structure.
 *
 * Creates the vectors based on the DMDA object.
 *
 * \param field The PetibmField structure to initialize (passed by reference).
 */
PetscErrorCode PetibmFieldInitialize(PetibmField &field)
{
  PetscErrorCode ierr;
  DMDALocalInfo info;

  PetscFunctionBeginUser;

  ierr = DMCreateGlobalVector(field.da, &field.global); CHKERRQ(ierr);
  ierr = DMCreateLocalVector(field.da, &field.local);
  ierr = DMDAGetLocalInfo(field.da, &info);
  ierr = VecCreateSeq(PETSC_COMM_SELF, info.mx, &field.x); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, info.my, &field.y); CHKERRQ(ierr);
  if (info.dim == 3)
  {
  	ierr = VecCreateSeq(PETSC_COMM_SELF, info.mz, &field.z); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
} // PetibmFieldInitialize


/*! Destroys the PETSc objects of a PetibmField structure.
 *
 * \param field The PetibmField structure (passed by reference).
 */
PetscErrorCode PetibmFieldDestroy(PetibmField &field)
{
	PetscErrorCode ierr;
	DMDALocalInfo info;

	PetscFunctionBeginUser;

	ierr = DMDAGetLocalInfo(field.da, &info);
	ierr = VecDestroy(&field.x); CHKERRQ(ierr);
	ierr = VecDestroy(&field.y); CHKERRQ(ierr);
	if (info.dim == 3)
	{
		ierr = VecDestroy(&field.z); CHKERRQ(ierr);
	}
	ierr = VecDestroy(&field.global); CHKERRQ(ierr);
	ierr = VecDestroy(&field.local); CHKERRQ(ierr);
	ierr = DMDestroy(&field.da); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmFieldDestroy


/*! Reads the field values from a HDF5 file.
 *
 * \param filePath Path of the file to read.
 * \param name The name of the field.
 * \param field The PetibmField structure (passed by reference).
 */
PetscErrorCode PetibmFieldReadValues(
	std::string filePath, std::string name, PetibmField &field)
{
	PetscErrorCode ierr;
	PetscViewer viewer;

	PetscFunctionBeginUser;

	ierr = PetscObjectSetName((PetscObject) field.global, name.c_str()); CHKERRQ(ierr);
	ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
	ierr = PetscViewerSetType(viewer, PETSCVIEWERHDF5); CHKERRQ(ierr);
	ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ); CHKERRQ(ierr);
	ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
	ierr = VecLoad(field.global, viewer); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmFieldReadValues


/*! Reads a HDF5 grid file for a given field.
 *
 * \param filePath Path of the grid file to read.
 * \param field The PetibmField structure (passed by reference).
 */
PetscErrorCode PetibmFieldReadGrid(std::string filePath, PetibmField &field)
{
	PetscErrorCode ierr;
	PetscViewer viewer;
	DMDALocalInfo info;

	PetscFunctionBeginUser;

	ierr = DMDAGetLocalInfo(field.da, &info); CHKERRQ(ierr);

	ierr = PetscViewerHDF5Open(
		PETSC_COMM_SELF, filePath.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
	// read x-stations
	ierr = PetscObjectSetName((PetscObject) field.x, "x"); CHKERRQ(ierr);
	ierr = VecLoad(field.x, viewer); CHKERRQ(ierr);
	// read y-stations
	ierr = PetscObjectSetName((PetscObject) field.y, "y"); CHKERRQ(ierr);
	ierr = VecLoad(field.y, viewer); CHKERRQ(ierr);
	if (info.dim == 3)
	{
		// read z-stations
		ierr = PetscObjectSetName((PetscObject) field.z, "z"); CHKERRQ(ierr);
		ierr = VecLoad(field.z, viewer); CHKERRQ(ierr);
	}

	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmFieldReadGrid


/*! Writes a grid into a HDF5 file.
 *
 * \param filePath Path of the file to write into.
 * \param field PetibmField structure that contains the grid.
 */
PetscErrorCode PetibmFieldWriteGrid(std::string filePath, PetibmField field)
{
	PetscErrorCode ierr;
	PetscViewer viewer;
	DMDALocalInfo info;

	PetscFunctionBeginUser;

	ierr = DMDAGetLocalInfo(field.da, &info); CHKERRQ(ierr);

	ierr = PetscViewerHDF5Open(
		PETSC_COMM_SELF, filePath.c_str(), FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
	// read x-stations
	ierr = PetscObjectSetName((PetscObject) field.x, "x"); CHKERRQ(ierr);
	ierr = VecView(field.x, viewer); CHKERRQ(ierr);
	// read y-stations
	ierr = PetscObjectSetName((PetscObject) field.y, "y"); CHKERRQ(ierr);
	ierr = VecView(field.y, viewer); CHKERRQ(ierr);
	if (info.dim == 3)
	{
		// read z-stations
		ierr = PetscObjectSetName((PetscObject) field.z, "z"); CHKERRQ(ierr);
		ierr = VecView(field.z, viewer); CHKERRQ(ierr);
	}

	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmFieldWriteGrid


/*! Writes the field values into a HDF5 file.
 *
 * \param filePath Path of the file to write into.
 * \param name Name of the field.
 * \param field PetibmField structure.
 */
PetscErrorCode PetibmFieldWriteValues(
	std::string filePath, std::string name, PetibmField field)
{
	PetscErrorCode ierr;
	PetscViewer viewer;

	PetscFunctionBeginUser;

	ierr = PetscObjectSetName(
		(PetscObject) field.global, name.c_str()); CHKERRQ(ierr);
	ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
	ierr = PetscViewerSetType(viewer, PETSCVIEWERHDF5); CHKERRQ(ierr);
	ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
	ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
	ierr = VecView(field.global, viewer); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmFieldWriteValues
