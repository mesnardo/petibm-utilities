/*! Converts the numerical solution from one format to another.
 * \file convert.cpp
 */

#include <petscsys.h>
#include <petscdmda.h>

#include "petibm-utilities/field.h"
#include "petibm-utilities/grid.h"
#include "petibm-utilities/misc.h"


static char help[] = "petibm-convert (0.1.0)\n\n" \
"Convert a field solution from/to PETSc binary format to/from HDF5 format.\n\n" \
"Usage: petibm-convert [arg]\n" \
"Options and arguments:\n" \
"  -config_file <path>\tInsert options and arguments from a given file\n" \
"  -source <path>\tPath of the field to convert\n" \
"  -destination <path>\tPath of the converted field to write\n" \
"  -hdf52binary <bool>\tConvert from HDF5 to PETSc binary (otherwise binary -> HDF5)\n" \
"  -name <string>\tName of the field\n" \
"  -nx <int>\tNumber of points in the x-direction\n" \
"  -ny <int>\tNumber of points in the y-direction\n" \
"  -nz <int>\tNumber of points in the z-direction\n" \
"  -periodic_x\tUse periodic boundary conditions in the x-direction\n" \
"  -periodic_y\tUse periodic boundary conditions in the y-direction\n" \
"  -periodic_z\tUse periodic boundary conditions in the z-direction\n" \
"\n"
;


struct AppCtx
{
	char source[PETSC_MAX_PATH_LEN];
	char destination[PETSC_MAX_PATH_LEN];
	PetscBool hdf52binary = PETSC_FALSE;
}; // AppCtx


PetscErrorCode AppGetOptions(const char prefix[], AppCtx *ctx)
{
	PetscErrorCode ierr;
	PetscBool found;

	PetscFunctionBeginUser;

	ierr = PetscOptionsGetString(
		nullptr, prefix, "-source", ctx->source,
		sizeof(ctx->source), &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetString(
		nullptr, prefix, "-destination", ctx->destination,
		sizeof(ctx->destination), &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(
		nullptr, prefix, "-hdf52binary", &ctx->hdf52binary, &found); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // AppGetOptions


int main(int argc, char **argv)
{
	PetscErrorCode ierr;
	DM da;
	PetibmGridCtx gridCtx;
	PetibmFieldCtx fieldCtx;
	AppCtx appCtx;
	PetibmField field;

	ierr = PetscInitialize(&argc, &argv, nullptr, help); CHKERRQ(ierr);

	// parse command-line options
	ierr = PetibmOptionsInsertFile(nullptr); CHKERRQ(ierr);
	ierr = PetibmGridGetOptions(nullptr, &gridCtx); CHKERRQ(ierr);
	ierr = PetibmFieldGetOptions(nullptr, &fieldCtx); CHKERRQ(ierr);
	ierr = AppGetOptions(nullptr, &appCtx); CHKERRQ(ierr);

	// create DMDA object
	ierr = PetibmFieldDMDACreate(gridCtx, fieldCtx, da); CHKERRQ(ierr);

	// initialize, read, and write
	ierr = PetibmFieldInitialize(da, field); CHKERRQ(ierr);
	std::string directory;
	ierr = PetibmGetParentDirectory(appCtx.destination, directory); CHKERRQ(ierr);
	ierr = PetibmCreateDirectory(directory); CHKERRQ(ierr);
	if (appCtx.hdf52binary)
	{
		ierr = PetibmFieldHDF5Read(
			PETSC_COMM_WORLD, appCtx.source, fieldCtx.name, field); CHKERRQ(ierr);
		ierr = PetibmFieldBinaryWrite(
			PETSC_COMM_WORLD, appCtx.destination, fieldCtx.name, field); CHKERRQ(ierr);
	}
	else
	{
		ierr = PetibmFieldBinaryRead(
			PETSC_COMM_WORLD, appCtx.source, fieldCtx.name, field); CHKERRQ(ierr);
		ierr = PetibmFieldHDF5Write(
			PETSC_COMM_WORLD, appCtx.destination, fieldCtx.name, field); CHKERRQ(ierr);
	}

	ierr = PetibmFieldDestroy(field); CHKERRQ(ierr);
	ierr = PetscFinalize(); CHKERRQ(ierr);
	return 0;
} // main
