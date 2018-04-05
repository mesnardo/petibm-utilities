/*! Implementation of the functions for the PetibmField structure.
 * \file field.cpp
 */

#include <vector>
#include <cstring>
#include <fstream>

#include <petscviewerhdf5.h>

#include "petibm-utilities/field.h"
#include "petibm-utilities/misc.h"


/*! Gets options from command-line or config file.
 *
 * \param prefix String to prepend to options.
 * \param ctx The PetibmFieldCtx structure to fill (passed by pointer).
 */
PetscErrorCode PetibmFieldGetOptions(
	const char prefix[], PetibmFieldCtx *ctx)
{
	PetscErrorCode ierr;
	PetscBool found;

	PetscFunctionBeginUser;

	// get path of file containing field values
	ierr = PetscOptionsGetString(nullptr, prefix, "-path", ctx->path,
	                             sizeof(ctx->path), &found); CHKERRQ(ierr);
	// get name of the field
	ierr = PetscOptionsGetString(nullptr, prefix, "-name", ctx->name,
	                             sizeof(ctx->name), &found); CHKERRQ(ierr);
	// get external boundary value
	ierr = PetscOptionsGetReal(
		nullptr, prefix, "-bc_value", &ctx->bc_value, &found); CHKERRQ(ierr);
	// check if periodic in the x-direction
	ierr = PetscOptionsGetBool(nullptr, prefix, "-periodic_x",
	                           &ctx->periodic_x, &found); CHKERRQ(ierr);
	// check if periodic in the y-direction
	ierr = PetscOptionsGetBool(nullptr, prefix, "-periodic_y",
	                           &ctx->periodic_y, &found); CHKERRQ(ierr);
	// check if periodic in the z-direction
	ierr = PetscOptionsGetBool(nullptr, prefix, "-periodic_z",
	                           &ctx->periodic_z, &found); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmFieldGetOptions


PetscErrorCode PetibmFieldCtxPrintf(
	const std::string name, const PetibmFieldCtx ctx)
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	ierr = PetscPrintf(PETSC_COMM_WORLD, "+ %s:\n", name.c_str()); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\t- name: %s\n",
	                   ctx.name); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\t- path: %s\n",
	                   ctx.path); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\t- boundary value: %f\n",
	                   ctx.bc_value); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmFieldCtxPrintf


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
	const PetibmFieldCtx ctx, const PetibmGrid grid, PetibmField &field)
{
  PetscErrorCode ierr;
  PetscInt M, N, P, m, n, p;
  DMBoundaryType bx, by, bz;
  const PetscInt *lx, *ly, *lz;

  PetscFunctionBeginUser;

  // get type of boundary conditions
  bx = (ctx.periodic_x) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
  by = (ctx.periodic_y) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
  bz = (ctx.periodic_z) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
  // get info from gridline in x-direction
	ierr = DMDAGetInfo(grid.x.da,
	                   nullptr,
	                   &M, nullptr, nullptr,
	                   &m, nullptr, nullptr,
	                   nullptr, nullptr,
	                   nullptr, nullptr, nullptr,
	                   nullptr); CHKERRQ(ierr);
	ierr = DMDAGetOwnershipRanges(
		grid.x.da, &lx, nullptr, nullptr); CHKERRQ(ierr);
	// get info from gridline in y-direction
	ierr = DMDAGetInfo(grid.y.da,
	                   nullptr,
	                   &N, nullptr, nullptr,
	                   &n, nullptr, nullptr,
	                   nullptr, nullptr,
	                   nullptr, nullptr, nullptr,
	                   nullptr); CHKERRQ(ierr);
	ierr = DMDAGetOwnershipRanges(
		grid.y.da, &ly, nullptr, nullptr); CHKERRQ(ierr);
	if (grid.dim == 3)
	{
		// get info from gridline in z-direction
		ierr = DMDAGetInfo(grid.z.da,
		                   nullptr,
		                   &P, nullptr, nullptr,
		                   &p, nullptr, nullptr,
		                   nullptr, nullptr,
		                   nullptr, nullptr, nullptr,
		                   nullptr); CHKERRQ(ierr);
		ierr = DMDAGetOwnershipRanges(
			grid.z.da, &lz, nullptr, nullptr); CHKERRQ(ierr);
		// create 3D DMDA
		ierr = DMDACreate3d(PETSC_COMM_WORLD,
		                    bx, by, bz,
		                    DMDA_STENCIL_BOX,
		                    M, N, P,
		                    m, n, p,
		                    1, 1,
		                    lx, ly, lz,
		                    &field.da); CHKERRQ(ierr);
	}
	else
	{
		// create 2D DMDA
		ierr = DMDACreate2d(PETSC_COMM_WORLD,
		                    bx, by,
		                    DMDA_STENCIL_BOX,
		                    M, N,
		                    m, n,
		                    1, 1,
		                    lx, ly,
		                    &field.da); CHKERRQ(ierr);
	}
	ierr = DMSetFromOptions(field.da); CHKERRQ(ierr);
	ierr = DMSetUp(field.da); CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(field.da, &field.global); CHKERRQ(ierr);
  ierr = DMCreateLocalVector(field.da, &field.local); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // PetibmFieldInitialize


/*! Initializes a PetibmField structure based on a given DMDA object.
 *
 * Creates the vectors based on the provided DMDA object.
 *
 * \param da The DMDA object.
 * \param field The field to initialize (passed by reference).
 */
PetscErrorCode PetibmFieldInitialize(const DM da, PetibmField &field)
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	field.da = da;
	ierr = DMCreateGlobalVector(da, &field.global); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(da, &field.local); CHKERRQ(ierr);	

	PetscFunctionReturn(0);
} // PetibmFieldInitialize


/*! Sets the value at external boundary points.
 *
 * \param value The value on the external boundaries.
 * \param field The field to modify (passed by reference).
 */
PetscErrorCode PetibmFieldSetBoundaryPoints(
	const PetscReal value, PetibmField &field)
{
	PetscErrorCode ierr;
	DMDALocalInfo info;
	PetscInt i, j, k;

	PetscFunctionBeginUser;

	ierr = DMDAGetLocalInfo(field.da, &info); CHKERRQ(ierr);
	if (info.dim == 2)
	{
		PetscReal **arr;
		ierr = DMDAVecGetArray(field.da, field.local, &arr); CHKERRQ(ierr);
		if (info.bx != DM_BOUNDARY_PERIODIC)
		{
			if (info.xs == 0)
				for (j=info.gys; j<info.gys+info.gym; j++)
					arr[j][-1] = value;
			if (info.xs+info.xm == info.mx)
				for (j=info.gys; j<info.gys+info.gym; j++)
					arr[j][info.mx] = value;
		}
		if (info.by != DM_BOUNDARY_PERIODIC)
		{
			if (info.ys == 0)
				for (i=info.gxs; i<info.gxs+info.gxm; i++)
					arr[-1][i] = value;
			if (info.ys+info.ym == info.my)
				for (i=info.gxs; i<info.gxs+info.gxm; i++)
					arr[info.my][i] = value;
		}
		ierr = DMDAVecRestoreArray(field.da, field.local, &arr); CHKERRQ(ierr);
	}
	else if (info.dim == 3)
	{
		PetscReal ***arr;
		ierr = DMDAVecGetArray(field.da, field.local, &arr); CHKERRQ(ierr);
		if (info.bx != DM_BOUNDARY_PERIODIC)
		{
			if (info.xs == 0)
				for (k=info.gzs; k<info.gzs+info.gzm; k++)
					for (j=info.gys; j<info.gys+info.gym; j++)
						arr[k][j][-1] = value;
			if (info.xs+info.xm == info.mx)
				for (k=info.gzs; k<info.gzs+info.gzm; k++)
					for (j=info.gys; j<info.gys+info.gym; j++)
						arr[k][j][info.mx] = value;
		}
		if (info.by != DM_BOUNDARY_PERIODIC)
		{
			if (info.ys == 0)
				for (k=info.gzs; k<info.gzs+info.gzm; k++)
					for (i=info.gxs; i<info.gxs+info.gxm; i++)
						arr[k][-1][i] = value;
			if (info.ys+info.ym == info.my)
				for (k=info.gzs; k<info.gzs+info.gzm; k++)
					for (i=info.gxs; i<info.gxs+info.gxm; i++)
						arr[k][info.my][i] = value;
		}
		if (info.bz != DM_BOUNDARY_PERIODIC)
		{
			if (info.zs == 0)
				for (j=info.gys; j<info.gys+info.gym; j++)
					for (i=info.gxs; i<info.gxs+info.gxm; i++)
						arr[-1][j][i] = value;
			if (info.zs+info.zm == info.mz)
				for (j=info.gys; j<info.gys+info.gym; j++)
					for (i=info.gxs; i<info.gxs+info.gxm; i++)
						arr[info.mz][j][i] = value;
		}
		ierr = DMDAVecRestoreArray(field.da, field.local, &arr); CHKERRQ(ierr);
	}
	else
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,
		        "Function only supports 2D or 3D fields");

	PetscFunctionReturn(0);
} // PetibmFieldSetBoundaryPoints


/*! Inserts values from global vector into local vector.
 *
 * \param field The field to work on (passed by reference).
 */
PetscErrorCode PetibmFieldGlobalToLocal(PetibmField &field)
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	ierr = DMGlobalToLocalBegin(
		field.da, field.global, INSERT_VALUES, field.local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(
		field.da, field.global, INSERT_VALUES, field.local); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmFieldGlobalToLocal


/*! Destroys the PETSc objects of a PetibmField structure.
 *
 * \param field The PetibmField structure (passed by reference).
 */
PetscErrorCode PetibmFieldDestroy(PetibmField &field)
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	ierr = VecDestroy(&field.global); CHKERRQ(ierr);
	ierr = VecDestroy(&field.local); CHKERRQ(ierr);
	ierr = DMDestroy(&field.da); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmFieldDestroy


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
                               PetibmField &field)
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	if (std::strcmp(viewerType, PETSCVIEWERBINARY) == 0)
	{
		ierr = PetibmFieldBinaryRead(filepath, name, field); CHKERRQ(ierr);
	}
	else if (std::strcmp(viewerType, PETSCVIEWERHDF5) == 0)
	{
		ierr = PetibmFieldHDF5Read(filepath, name, field); CHKERRQ(ierr);
	}
	else
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,
		        "Function only supports Binary and HDF5 viewers");

	PetscFunctionReturn(0);
} // PetibmFieldRead


/*! Reads the field values stored in HDF5 format from file.
 *
 * \param filepath Path of the input file.
 * \param name The name of the field.
 * \param field The PetibmField structure (passed by reference).
 */
PetscErrorCode PetibmFieldHDF5Read(
	const std::string filepath, const std::string name, PetibmField &field)
{
	PetscErrorCode ierr;
	PetscViewer viewer;

	PetscFunctionBeginUser;

	ierr = PetscObjectSetName(
		(PetscObject) field.global, name.c_str()); CHKERRQ(ierr);
	ierr = PetscViewerHDF5Open(
		PETSC_COMM_SELF, filepath.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
	ierr = VecLoad(field.global, viewer); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmFieldHDF5Read


/*! Reads the field values stored in binary format from file.
 *
 * \param filepath Path of the input file.
 * \param name The name of the field.
 * \param field The PetibmField structure (passed by reference).
 */
PetscErrorCode PetibmFieldBinaryRead(
	const std::string filepath, const std::string name, PetibmField &field)
{
	PetscErrorCode ierr;
	PetscViewer viewer;

	PetscFunctionBeginUser;

	ierr = PetscObjectSetName(
		(PetscObject) field.global, name.c_str()); CHKERRQ(ierr);
	ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
	ierr = PetscViewerSetType(viewer, PETSCVIEWERBINARY); CHKERRQ(ierr);
	ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ); CHKERRQ(ierr);
	ierr = PetscViewerFileSetName(viewer, filepath.c_str()); CHKERRQ(ierr);
	ierr = VecLoad(field.global, viewer); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmFieldBinaryRead


/*! Writes the field values into file in HDF5 format.
 *
 * \param filepath Path of the output file.
 * \param name Name of the field.
 * \param field PetibmField structure.
 */
PetscErrorCode PetibmFieldHDF5Write(
	const std::string filepath, const std::string name, const PetibmField field)
{
	PetscErrorCode ierr;
	PetscViewer viewer;
	PetscFileMode mode;

	PetscFunctionBeginUser;

	std::ifstream infile(filepath.c_str());
	if (infile.good())
		mode = FILE_MODE_APPEND;
	else
		mode = FILE_MODE_WRITE;

	ierr = PetscObjectSetName(
		(PetscObject) field.global, name.c_str()); CHKERRQ(ierr);
	ierr = PetscViewerHDF5Open(
		PETSC_COMM_SELF, filepath.c_str(), mode, &viewer); CHKERRQ(ierr);
	ierr = VecView(field.global, viewer); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmFieldHDF5Write


/*! Writes the field values into file in binary format.
 *
 * \param filepath Path of the output file.
 * \param name Name of the field.
 * \param field PetibmField structure.
 */
PetscErrorCode PetibmFieldBinaryWrite(
	const std::string filepath, const std::string name, const PetibmField field)
{
	PetscErrorCode ierr;
	PetscViewer viewer;

	PetscFunctionBeginUser;

	ierr = PetscObjectSetName(
		(PetscObject) field.global, name.c_str()); CHKERRQ(ierr);
	ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
	ierr = PetscViewerSetType(viewer, PETSCVIEWERBINARY); CHKERRQ(ierr);
	ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
	ierr = PetscViewerFileSetName(viewer, filepath.c_str()); CHKERRQ(ierr);
	ierr = VecView(field.global, viewer); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmFieldBinaryWrite


/*! Interpolates field A associated with grid A onto grid B.
 *
 * \param gridA The grid to interpolate from.
 * \param fieldA The field to interpolate.
 * \param gridB The grid to interpolate on.
 * \param fieldB The resulting interpolated field (passed by reference).
 */
PetscErrorCode PetibmFieldInterpolate(
	PetibmGrid gridA, PetibmField fieldA, PetibmGrid gridB, PetibmField &fieldB)
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	if (gridA.dim == 2)
	{
		ierr = PetibmFieldInterpolate2d(
			gridA, fieldA, gridB, fieldB); CHKERRQ(ierr);
	}
	else if (gridA.dim == 3)
	{
		ierr = PetibmFieldInterpolate3d(
			gridA, fieldA, gridB, fieldB); CHKERRQ(ierr);
	}
	else
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,
		        "Function only supports 2D or 3D fields");

	PetscFunctionReturn(0);
} // PetibmFieldInterpolate


/*! Interpolates a 2D field A associated with grid A onto grid B.
 *
 * Performs bi-linear interpolation.
 *
 * \param gridA The grid to interpolate from.
 * \param fieldA The field to interpolate.
 * \param gridB The grid to interpolate on.
 * \param fieldB The resulting interpolated field (passed by reference).
 */
PetscErrorCode PetibmFieldInterpolate2d(
	PetibmGrid gridA, PetibmField fieldA, PetibmGrid gridB, PetibmField &fieldB)
{
	PetscErrorCode ierr;
	DMDALocalInfo info;
	PetscInt i, j, I, J;
	PetscReal *xA, *yA, *xB, *yB;
	PetscReal **vA, **vB;
	PetscReal v1, v2;
	std::vector<PetscInt> Iv, Jv;

	PetscFunctionBeginUser;

	ierr = PetibmFieldGlobalToLocal(fieldA); CHKERRQ(ierr);
	ierr = PetibmGridGlobalToLocal(gridA); CHKERRQ(ierr);
	ierr = PetibmGridGlobalToLocal(gridB); CHKERRQ(ierr);

	ierr = PetibmGetNeighbors1D(gridB.x, gridA.x, Iv); CHKERRQ(ierr);
	ierr = PetibmGetNeighbors1D(gridB.y, gridA.y, Jv); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(fieldA.da, fieldA.local, &vA); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fieldB.da, fieldB.global, &vB); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(gridA.x.da, gridA.x.local, &xA); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(gridA.y.da, gridA.y.local, &yA); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(gridB.x.da, gridB.x.local, &xB); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(gridB.y.da, gridB.y.local, &yB); CHKERRQ(ierr);
	ierr = DMDAGetLocalInfo(fieldB.da, &info); CHKERRQ(ierr);
	for (j=info.ys; j<info.ys+info.ym; j++)
	{
		J = Jv[j-info.ys];
		for (i=info.xs; i<info.xs+info.xm;i++)
		{
			I = Iv[i-info.xs];
			v1 = ((xA[I+1] - xB[i]) / (xA[I+1] - xA[I]) * vA[J][I] +
			      (xB[i] - xA[I]) / (xA[I+1] - xA[I]) * vA[J][I+1]);
			v2 = ((xA[I+1] - xB[i]) / (xA[I+1] - xA[I]) * vA[J+1][I] +
			      (xB[i] - xA[I]) / (xA[I+1] - xA[I]) * vA[J+1][I+1]);
			vB[j][i] = ((yA[J+1] - yB[j]) / (yA[J+1] - yA[J]) * v1 +
			            (yB[j] - yA[J]) / (yA[J+1] - yA[J]) * v2);
		}
	}
	ierr = DMDAVecRestoreArray(gridA.x.da, gridA.x.local, &xA); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(gridA.y.da, gridA.y.local, &yA); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(gridB.x.da, gridB.x.local, &xB); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(gridB.y.da, gridB.y.local, &yB); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fieldA.da, fieldA.local, &vA); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fieldB.da, fieldB.global, &vB); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmFieldInterpolate2d


/*! Interpolates a 3D field A associated with grid A onto grid B.
 *
 * Performs tri-linear interpolation.
 *
 * \param gridA The grid to interpolate from.
 * \param fieldA The field to interpolate.
 * \param gridB The grid to interpolate on.
 * \param fieldB The resulting interpolated field (passed by reference).
 */
PetscErrorCode PetibmFieldInterpolate3d(
	PetibmGrid gridA, PetibmField fieldA, PetibmGrid gridB, PetibmField &fieldB)
{
	PetscErrorCode ierr;
	DMDALocalInfo info;
	PetscInt i, j, k, I, J, K;
	PetscReal *xA, *yA, *zA, *xB, *yB, *zB;
	PetscReal ***vA, ***vB;
	PetscReal v1, v2, v3, v4, v12, v34;
	std::vector<PetscInt> Iv, Jv, Kv;

	PetscFunctionBeginUser;

	ierr = PetibmFieldGlobalToLocal(fieldA); CHKERRQ(ierr);
	ierr = PetibmGridGlobalToLocal(gridA); CHKERRQ(ierr);
	ierr = PetibmGridGlobalToLocal(gridB); CHKERRQ(ierr);

	ierr = PetibmGetNeighbors1D(gridB.x, gridA.x, Iv); CHKERRQ(ierr);
	ierr = PetibmGetNeighbors1D(gridB.y, gridA.y, Jv); CHKERRQ(ierr);
	ierr = PetibmGetNeighbors1D(gridB.z, gridA.z, Kv); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(fieldA.da, fieldA.local, &vA); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fieldB.da, fieldB.global, &vB); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(gridA.x.da, gridA.x.local, &xA); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(gridA.y.da, gridA.y.local, &yA); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(gridA.z.da, gridA.z.local, &zA); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(gridB.x.da, gridB.x.local, &xB); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(gridB.y.da, gridB.y.local, &yB); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(gridB.z.da, gridB.z.local, &zB); CHKERRQ(ierr);
	ierr = DMDAGetLocalInfo(fieldB.da, &info); CHKERRQ(ierr);
	for (k=info.zs; k<info.zs+info.zm; k++)
	{
		K = Kv[k-info.zs];
		for (j=info.ys; j<info.ys+info.ym; j++)
		{
			J = Jv[j-info.ys];
			for (i=info.xs; i<info.xs+info.xm;i++)
			{
				I = Iv[i-info.xs];
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
	ierr = DMDAVecRestoreArray(gridA.x.da, gridA.x.local, &xA); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(gridA.y.da, gridA.y.local, &yA); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(gridA.z.da, gridA.z.local, &zA); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(gridB.x.da, gridB.x.local, &xB); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(gridB.y.da, gridB.y.local, &yB); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(gridB.z.da, gridB.z.local, &zA); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fieldA.da, fieldA.local, &vA); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fieldB.da, fieldB.global, &vB); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmFieldInterpolate3d


PetscErrorCode PetibmFieldDMDACreate2d(
	const PetibmGridCtx gridCtx, const PetibmFieldCtx fieldCtx, DM &da)
{
	PetscErrorCode ierr;
	DMBoundaryType bType_x, bType_y;

	PetscFunctionBeginUser;

	bType_x = (fieldCtx.periodic_x) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	bType_y = (fieldCtx.periodic_y) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	ierr = DMDACreate2d(PETSC_COMM_WORLD,
	                    bType_x, bType_y,
	                    DMDA_STENCIL_STAR,
	                    gridCtx.nx, gridCtx.ny,
	                    PETSC_DECIDE, PETSC_DECIDE,
	                    1, 1,
	                    nullptr, nullptr,
	                    &da); CHKERRQ(ierr);
	ierr = DMSetFromOptions(da); CHKERRQ(ierr);
	ierr = DMSetUp(da); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmFieldDMDACreate2d


PetscErrorCode PetibmFieldDMDACreate3d(
	const PetibmGridCtx gridCtx, const PetibmFieldCtx fieldCtx, DM &da)
{
	PetscErrorCode ierr;
	DMBoundaryType bType_x, bType_y, bType_z;

	PetscFunctionBeginUser;

	bType_x = (fieldCtx.periodic_x) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	bType_y = (fieldCtx.periodic_y) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	bType_z = (fieldCtx.periodic_z) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
	                    bType_x, bType_y, bType_z,
	                    DMDA_STENCIL_STAR,
	                    gridCtx.nx, gridCtx.ny, gridCtx.nz,
	                    PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
	                    1, 1,
	                    nullptr, nullptr, nullptr,
	                    &da); CHKERRQ(ierr);
	ierr = DMSetFromOptions(da); CHKERRQ(ierr);
	ierr = DMSetUp(da); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmFieldDMDACreate3d


PetscErrorCode PetibmFieldDMDACreate2d(
	const std::string name, const DM da_in, DM &da)
{
	PetscErrorCode ierr;
	const PetscInt *plx, *ply;
	PetscInt *lx, *ly;
	PetscInt M, N, m, n;
	DMBoundaryType bType_x, bType_y;

	PetscFunctionBeginUser;

	ierr = DMDAGetOwnershipRanges(da_in, &plx, &ply, nullptr); CHKERRQ(ierr);
	ierr = DMDAGetInfo(da_in,
	                   nullptr,
	                   &M, &N, nullptr,
	                   &m, &n, nullptr,
	                   nullptr, nullptr,
	                   &bType_x, &bType_y, nullptr,
	                   nullptr); CHKERRQ(ierr);
	// create DMDA and vector for velocity in x-direction
	ierr = PetscMalloc(m*sizeof(*lx), &lx); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*ly), &ly); CHKERRQ(ierr);
	ierr = PetscMemcpy(lx, plx, m*sizeof(*lx)); CHKERRQ(ierr);
	ierr = PetscMemcpy(ly, ply, n*sizeof(*ly)); CHKERRQ(ierr);
	if (name == "ux" and bType_x != DM_BOUNDARY_PERIODIC)
	{
		lx[m-1]--;
		M--;
	}
	if (name == "uy" and bType_y != DM_BOUNDARY_PERIODIC)
	{
		ly[n-1]--;
		N--;
	}
	if (name == "wz")
	{
		lx[m-1]--;
		ly[n-1]--;
		M--;
		N--;
	}
	ierr = DMDACreate2d(PETSC_COMM_WORLD,
	                    bType_x, bType_y,
	                    DMDA_STENCIL_BOX,
	                    M, N, m, n, 1, 1, lx, ly,
	                    &da); CHKERRQ(ierr);
	ierr = DMSetFromOptions(da); CHKERRQ(ierr);
	ierr = DMSetUp(da); CHKERRQ(ierr);
	ierr = PetscFree(lx); CHKERRQ(ierr);
	ierr = PetscFree(ly); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmFieldDMDACreate2d


PetscErrorCode PetibmFieldDMDACreate3d(
	const std::string name, const DM da_in, DM &da)
{
	PetscErrorCode ierr;
	const PetscInt *plx, *ply, *plz;
	PetscInt *lx, *ly, *lz;
	PetscInt M, N, P, m, n, p;
	DMBoundaryType bType_x, bType_y, bType_z;

	PetscFunctionBeginUser;

	ierr = DMDAGetOwnershipRanges(da_in, &plx, &ply, &plz); CHKERRQ(ierr);
	ierr = DMDAGetInfo(da_in,
	                   nullptr,
	                   &M, &N, &P,
	                   &m, &n, &p,
	                   nullptr, nullptr,
	                   &bType_x, &bType_y, &bType_z,
	                   nullptr); CHKERRQ(ierr);
	// create DMDA and vector for velocity in x-direction
	ierr = PetscMalloc(m*sizeof(*lx), &lx); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*ly), &ly); CHKERRQ(ierr);
	ierr = PetscMalloc(p*sizeof(*lz), &lz); CHKERRQ(ierr);
	ierr = PetscMemcpy(lx, plx, m*sizeof(*lx)); CHKERRQ(ierr);
	ierr = PetscMemcpy(ly, ply, n*sizeof(*ly)); CHKERRQ(ierr);
	ierr = PetscMemcpy(lz, plz, p*sizeof(*lz)); CHKERRQ(ierr);
	if (name == "ux" and bType_x != DM_BOUNDARY_PERIODIC)
	{
		lx[m-1]--;
		M--;
	}
	if (name == "uy" and bType_y != DM_BOUNDARY_PERIODIC)
	{
		ly[n-1]--;
		N--;
	}
	if (name == "uz" and bType_z != DM_BOUNDARY_PERIODIC)
	{
		lz[p-1]--;
		P--;
	}
	if (name == "wx")
	{
		ly[n-1]--;
		lz[p-1]--;
		N--;
		P--;
	}
	if (name == "wz")
	{
		lx[m-1]--;
		ly[n-1]--;
		M--;
		N--;
	}
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
	                    bType_x, bType_y, bType_z,
	                    DMDA_STENCIL_BOX,
	                    M, N, P, m, n, p, 1, 1, lx, ly, lz,
	                    &da); CHKERRQ(ierr);
	ierr = DMSetFromOptions(da); CHKERRQ(ierr);
	ierr = DMSetUp(da); CHKERRQ(ierr);
	ierr = PetscFree(lx); CHKERRQ(ierr);
	ierr = PetscFree(ly); CHKERRQ(ierr);
	ierr = PetscFree(lz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmFieldDMDACreate3d
