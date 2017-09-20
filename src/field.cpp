/*! Implementation of the functions for the PetibmField structure.
 * \file field.cpp
 */

#include <petscviewerhdf5.h>

#include "petibm-utilities/field.h"
#include "petibm-utilities/misc.h"


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


/*! Interpolates 2D field values from one mesh to another.
 *
 * \param fieldA PetibmField to be interpolated.
 * \param fieldB PetibmField on which the solution is interpolated.
 */
PetscErrorCode PetibmFieldInterpolate2D(
	const PetibmField fieldA, PetibmField &fieldB)
{
	PetscErrorCode ierr;
	DMDALocalInfo info;
	PetscInt i, j, I = 0, J = 0;
	PetscBool found_x = PETSC_FALSE, found_y = PETSC_FALSE;
	PetscReal *xA, *yA, *xB, *yB;
	PetscReal **vA, **vB;
	PetscReal v1, v2;

	PetscFunctionBeginUser;

	ierr = DMGlobalToLocalBegin(
		fieldA.da, fieldA.global, INSERT_VALUES, fieldA.local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(
		fieldA.da, fieldA.global, INSERT_VALUES, fieldA.local); CHKERRQ(ierr);

	ierr = DMDAGetLocalInfo(fieldB.da, &info); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(fieldA.da, fieldA.local, &vA); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fieldB.da, fieldB.global, &vB); CHKERRQ(ierr);
	ierr = VecGetArray(fieldA.x, &xA); CHKERRQ(ierr);
	ierr = VecGetArray(fieldA.y, &yA); CHKERRQ(ierr);
	ierr = VecGetArray(fieldB.x, &xB); CHKERRQ(ierr);
	ierr = VecGetArray(fieldB.y, &yB); CHKERRQ(ierr);
	for (j=info.ys; j<info.ys+info.ym; j++)
	{
		ierr = PetibmGetNeighborIndex1D(
			yB[j], fieldA.y, &J, &found_y); CHKERRQ(ierr);
		I = 0;
		for (i=info.xs; i<info.xs+info.xm;i++)
		{
			ierr = PetibmGetNeighborIndex1D(
				xB[i], fieldA.x, &I, &found_x); CHKERRQ(ierr);
			v1 = ((xA[I+1] - xB[i]) / (xA[I+1] - xA[I]) * vA[J][I] +
			      (xB[i] - xA[I]) / (xA[I+1] - xA[I]) * vA[J][I+1]);
			v2 = ((xA[I+1] - xB[i]) / (xA[I+1] - xA[I]) * vA[J+1][I] +
			      (xB[i] - xA[I]) / (xA[I+1] - xA[I]) * vA[J+1][I+1]);
			vB[j][i] = ((yA[J+1] - yB[j]) / (yA[J+1] - yA[J]) * v1 +
			            (yB[j] - yA[J]) / (yA[J+1] - yA[J]) * v2);
		}
	}
	ierr = VecRestoreArray(fieldA.x, &xA); CHKERRQ(ierr);
	ierr = VecRestoreArray(fieldA.y, &yA); CHKERRQ(ierr);
	ierr = VecRestoreArray(fieldB.x, &xB); CHKERRQ(ierr);
	ierr = VecRestoreArray(fieldB.y, &yB); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fieldA.da, fieldA.local, &vA); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fieldB.da, fieldB.global, &vB); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmFieldInterpolate2D


/*! Interpolates 3D field values from one mesh to another.
 *
 * \param fieldA PetibmField to be interpolated.
 * \param fieldB PetibmField on which the solution is interpolated.
 */
PetscErrorCode PetibmFieldInterpolate3D(
	const PetibmField fieldA, PetibmField &fieldB)
{
	PetscErrorCode ierr;
	DMDALocalInfo info;
	PetscInt i, j, k, I = 0, J = 0, K=0;
	PetscBool found_x = PETSC_FALSE, found_y = PETSC_FALSE, found_z = PETSC_FALSE;
	PetscReal *xA, *yA, *zA, *xB, *yB, *zB;
	PetscReal ***vA, ***vB;
	PetscReal v1, v2, v3, v4, v12, v34;

	PetscFunctionBeginUser;

	ierr = DMGlobalToLocalBegin(
		fieldA.da, fieldA.global, INSERT_VALUES, fieldA.local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(
		fieldA.da, fieldA.global, INSERT_VALUES, fieldA.local); CHKERRQ(ierr);

	ierr = DMDAGetLocalInfo(fieldB.da, &info); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(fieldA.da, fieldA.local, &vA); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fieldB.da, fieldB.global, &vB); CHKERRQ(ierr);
	ierr = VecGetArray(fieldA.x, &xA); CHKERRQ(ierr);
	ierr = VecGetArray(fieldA.y, &yA); CHKERRQ(ierr);
	ierr = VecGetArray(fieldA.z, &zA); CHKERRQ(ierr);
	ierr = VecGetArray(fieldB.x, &xB); CHKERRQ(ierr);
	ierr = VecGetArray(fieldB.y, &yB); CHKERRQ(ierr);
	ierr = VecGetArray(fieldB.z, &zB); CHKERRQ(ierr);
	for (k=info.zs; k<info.zs+info.zm; k++)
	{
		ierr = PetibmGetNeighborIndex1D(
			zB[k], fieldA.z, &K, &found_z); CHKERRQ(ierr);
		J = 0;
		for (j=info.ys; j<info.ys+info.ym; j++)
		{
			ierr = PetibmGetNeighborIndex1D(
				yB[j], fieldA.y, &J, &found_y); CHKERRQ(ierr);
			I = 0;
			for (i=info.xs; i<info.xs+info.xm;i++)
			{
				ierr = PetibmGetNeighborIndex1D(
					xB[i], fieldA.x, &I, &found_x); CHKERRQ(ierr);
				v1 = ((xA[I+1] - xB[i]) / (xA[I+1] - xA[I]) * vA[K][J][I] +
				      (xB[i] - xA[I]) / (xA[I+1] - xA[I]) * vA[K][J][I+1]);
				v2 = ((xA[I+1] - xB[i]) / (xA[I+1] - xA[I]) * vA[K][J+1][I] +
				      (xB[i] - xA[I]) / (xA[I+1] - xA[I]) * vA[K][J+1][I+1]);
				v3 = ((xA[I+1] - xB[i]) / (xA[I+1] - xA[I]) * vA[K+1][J][I] +
				      (xB[i] - xA[I]) / (xA[I+1] - xA[I]) * vA[K+1][J][I+1]);
				v4 = ((xA[I+1] - xB[i]) / (xA[I+1] - xA[I]) * vA[K+1][J+1][I] +
				      (xB[i] - xA[I]) / (xA[I+1] - xA[I]) * vA[K+1][J+1][I+1]);
				v12 = ((yA[J+1] - yB[j]) / (yA[J+1] - yA[J]) * v1 +
				       (yB[j] - yA[J]) / (yA[J+1] - yA[J]) * v2);
				v34 = ((yA[J+1] - yB[j]) / (yA[J+1] - yA[J]) * v3 +
				       (yB[j] - yA[J]) / (yA[J+1] - yA[J]) * v4);
				vB[k][j][i] = ((zA[K+1] - zB[k]) / (zA[K+1] - zA[K]) * v12 +
				              (zB[k] - zA[K]) / (zA[K+1] - zA[K]) * v34);
			}
		}
	}
	ierr = VecRestoreArray(fieldA.x, &xA); CHKERRQ(ierr);
	ierr = VecRestoreArray(fieldA.y, &yA); CHKERRQ(ierr);
	ierr = VecRestoreArray(fieldA.z, &zA); CHKERRQ(ierr);
	ierr = VecRestoreArray(fieldB.x, &xB); CHKERRQ(ierr);
	ierr = VecRestoreArray(fieldB.y, &yB); CHKERRQ(ierr);
	ierr = VecRestoreArray(fieldB.z, &zB); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fieldA.da, fieldA.local, &vA); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fieldB.da, fieldB.global, &vB); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmFieldInterpolate3D


/*! Sets the value at external ghost points.
 *
 * \param field The PetibmField to update.
 * \param value Value at the external ghost points.
 */
PetscErrorCode  PetibmFieldExternalGhostPointsSet(
	PetibmField field, PetscReal value)
{
	PetscErrorCode ierr;
	DMDALocalInfo info;
	PetscInt i, j, k;

	PetscFunctionBeginUser;

	ierr = DMGlobalToLocalBegin(
		field.da, field.global, INSERT_VALUES, field.local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(
		field.da, field.global, INSERT_VALUES, field.local); CHKERRQ(ierr);

	ierr = DMDAGetLocalInfo(field.da, &info); CHKERRQ(ierr);
	if (info.dim == 2)
	{
		PetscReal **arr;
		ierr = DMDAVecGetArray(field.da, field.local, &arr); CHKERRQ(ierr);
		for (i=info.xs; i<info.xs+info.xm; i++)
			if (i == 0 or i == info.mx-1)
				for (j=info.gys; j<info.gys+info.gym; j++)
					arr[j][i] = value;
		for (j=info.ys; j<info.ys+info.ym; j++)
			if (j == 0 or j == info.my-1)
				for (i=info.gxs; i<info.gxs+info.gxm; i++)
					arr[j][i] = value;
		ierr = DMDAVecRestoreArray(field.da, field.local, &arr); CHKERRQ(ierr);
	}
	else if (info.dim == 3)
	{
		PetscReal ***arr;
		ierr = DMDAVecGetArray(field.da, field.local, &arr); CHKERRQ(ierr);
		for (i=info.xs; i<info.xs+info.xm; i++)
			if (i == 0 or i == info.mx-1)
				for (k=info.gzs; k<info.gzs+info.gzm; k++)
					for (j=info.gys; j<info.gys+info.gym; j++)
						arr[k][j][i] = value;
		for (j=info.ys; j<info.ys+info.ym; j++)
			if (j == 0 or j == info.my-1)
				for (k=info.gzs; k<info.gzs+info.gzm; k++)
					for (i=info.gxs; i<info.gxs+info.gxm; i++)
						arr[k][j][i] = value;
		for (k=info.zs; k<info.zs+info.zm; k++)
			if (k == 0 or k == info.mz-1)
				for (j=info.gys; j<info.gys+info.gym; j++)
					for (i=info.gxs; i<info.gxs+info.gxm; i++)
						arr[k][j][i] = value;
		ierr = DMDAVecRestoreArray(field.da, field.local, &arr); CHKERRQ(ierr);
	}
	else
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,
		        "Function only supports 2D or 3D fields");

	PetscFunctionReturn(0);
} // PetibmFieldExternalGhostPointsSet


/*! Gets options from command-line or config file.
 *
 * \param prefix String to prepend to options.
 * \param info The PetibmFieldInfo structure to fill.
 */
PetscErrorCode PetibmFieldInfoGetOptions(
	const char prefix[], PetibmFieldInfo *info)
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
	ierr = PetscOptionsGetString(nullptr, prefix, "-path", info->path,
	                             sizeof(info->path), &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetString(nullptr, prefix, "-name", info->name,
	                             sizeof(info->name), &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(
		nullptr, prefix, "-bc_value", &info->bcValue, &found); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmFieldInfoGetOptions
