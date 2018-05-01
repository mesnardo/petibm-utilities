/*! Interpolates a PetIBM field A from grid A to grid B.
 * \file interpolate.cpp
 */

#include <string>
#include <sys/stat.h>

#include <petscsys.h>

#include "petibm-utilities/field.h"
#include "petibm-utilities/grid.h"
#include "petibm-utilities/misc.h"


static char help[] = "petibm-interpolate (0.1.0)\n\n" \
"Interpolate the field solution from one grid to another.\n\n" \
"Usage: petibm-interpolate [arg]\n" \
"Options and arguments:\n" \
"  -config_file <path>\tInsert options and arguments from a given file\n" \
"  -gridA_name <string>\tName of the grid (PetIBM-0.3)\n" \
"  -gridA_path <string>\tPath of the grid file\n" \
"  -gridA_nx <int>\tNumber of cells in the x-direction\n" \
"  -gridA_ny <int>\tNumber of cells in the y-direction\n" \
"  -gridA_nz <int>\tNumber of cells in the z-direction\n" \
"  -gridA_x_start <float>\tBottom-left corner x-coordinate of the domain\n" \
"  -gridA_y_start <float>\tBottom-left corner y-coordinate of the domain\n" \
"  -gridA_z_start <float>\tBottom-left corner z-coordinate of the domain\n" \
"  -gridA_x_end <float>\tTop-right corner x-coordinate of the domain\n" \
"  -gridA_y_end <float>\tTop-right corner y-coordinate of the domain\n" \
"  -gridA_z_end <float>\tTop-right corner z-coordinate of the domain\n" \
"  -fieldA_name <string>\tName of the field to interpolate\n" \
"  -fieldA_path <path>\tPath of the field to interpolate\n" \
"  -fieldA_periodic_x\tUse periodic boundary conditions in the x-direction\n" \
"  -fieldA_periodic_y\tUse periodic boundary conditions in the y-direction\n" \
"  -fieldA_periodic_z\tUse periodic boundary conditions in the z-direction\n" \
"  -fieldA_bc_value <float>\tValue of the field at boundaries\n" \
"  -gridB_name <string>\tName of the grid (PetIBM-0.3)\n" \
"  -gridB_path <string>\tPath of the grid file\n" \
"  -gridB_nx <int>\tNumber of cells in the x-direction\n" \
"  -gridB_ny <int>\tNumber of cells in the y-direction\n" \
"  -gridB_nz <int>\tNumber of cells in the z-direction\n" \
"  -gridB_x_start <float>\tBottom-left corner x-coordinate of the domain\n" \
"  -gridB_y_start <float>\tBottom-left corner y-coordinate of the domain\n" \
"  -gridB_z_start <float>\tBottom-left corner z-coordinate of the domain\n" \
"  -gridB_x_end <float>\tTop-right corner x-coordinate of the domain\n" \
"  -gridB_y_end <float>\tTop-right corner y-coordinate of the domain\n" \
"  -gridB_z_end <float>\tTop-right corner z-coordinate of the domain\n" \
"  -fieldB_name <string>\tName of the interpolated field\n" \
"  -fieldB_path <path>\tPath of the output file for the interpolated field\n" \
"  -fieldB_bc_value <float>\tValue of the field at boundaries\n" \
"  -input_binary\tRead the PetIBM field written in PETSc binary format (PetIBM-0.2)\n" \
"  -output_binary\tWrite the PetIBM field in PETSc binary format (PetIBM-0.2)\n" \
"\n"
;


int main(int argc, char **argv)
{
	PetscErrorCode ierr;
	PetibmField fieldA, fieldB;
	PetibmFieldCtx fieldACtx, fieldBCtx;
	PetibmGrid gridA, gridB;
	PetibmGridCtx gridACtx, gridBCtx;
	Vec coords[3];
	std::string outdir;

	ierr = PetscInitialize(&argc, &argv, nullptr, help); CHKERRQ(ierr);
	
	ierr = PetibmOptionsInsertFile(nullptr); CHKERRQ(ierr);
#ifdef PETIBM_0_2
	// get the input and output formats
	PetscBool input_binary = PETSC_FALSE,
	          output_binary = PETSC_FALSE,
	          found = PETSC_FALSE;
	ierr = PetscOptionsGetBool(
		nullptr, nullptr, "-input_binary", &input_binary, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(
		nullptr, nullptr, "-output_binary", &output_binary, &found); CHKERRQ(ierr);
	PetscViewerType inViewerType = PETSCVIEWERHDF5,
	                outViewerType = PETSCVIEWERHDF5;
	if (input_binary == PETSC_TRUE)
		inViewerType = PETSCVIEWERBINARY;
	if (output_binary == PETSC_TRUE)
		outViewerType = PETSCVIEWERBINARY;
#endif

	// Create and read the grid A
	ierr = PetibmGridGetOptions("gridA_", &gridACtx); CHKERRQ(ierr);
	ierr = PetibmGridInitialize(gridACtx, gridA); CHKERRQ(ierr);
	ierr = PetibmGridHDF5Read(
		PETSC_COMM_SELF, gridACtx.path, gridACtx.name, gridA); CHKERRQ(ierr);
	ierr = PetibmGridSetBoundaryPoints(
		gridACtx.starts, gridACtx.ends, gridA); CHKERRQ(ierr);
	// Create and read the field A
	ierr = PetibmFieldGetOptions("fieldA_", &fieldACtx); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(fieldACtx, gridA, fieldA); CHKERRQ(ierr);
#ifdef PETIBM_0_2
	ierr = PetibmFieldRead(PETSC_COMM_WORLD, fieldACtx.path, fieldACtx.name,
	                       inViewerType, fieldA); CHKERRQ(ierr);
#else
	ierr = PetibmFieldHDF5Read(
		PETSC_COMM_WORLD, fieldACtx.path, fieldACtx.name, fieldA); CHKERRQ(ierr);
#endif
	ierr = PetibmFieldSetBoundaryPoints(
		fieldACtx.bc_value, fieldA); CHKERRQ(ierr);

	// Create and read the grid B
	ierr = PetibmGridGetOptions("gridB_", &gridBCtx); CHKERRQ(ierr);
	ierr = VecCreateSeq(PETSC_COMM_SELF, gridBCtx.nx, coords); CHKERRQ(ierr);
	ierr = PetibmGridlineHDF5Read(
		PETSC_COMM_SELF, gridBCtx.path, gridBCtx.name, "x", coords[0]); CHKERRQ(ierr);
	ierr = VecCreateSeq(PETSC_COMM_SELF, gridBCtx.ny, coords+1); CHKERRQ(ierr);
	ierr = PetibmGridlineHDF5Read(
		PETSC_COMM_SELF, gridBCtx.path, gridBCtx.name, "y", coords[1]); CHKERRQ(ierr);
	if (gridBCtx.nz > 0)
	{
		ierr = VecCreateSeq(PETSC_COMM_SELF, gridBCtx.nz, coords+2); CHKERRQ(ierr);
		ierr = PetibmGridlineHDF5Read(
			PETSC_COMM_SELF, gridBCtx.path, gridBCtx.name, "z", coords[2]); CHKERRQ(ierr);
	}
	ierr = PetibmGridInitialize(gridA, coords, gridB); CHKERRQ(ierr);
	ierr = PetibmGridHDF5Read(
		PETSC_COMM_SELF, gridBCtx.path, gridBCtx.name, gridB); CHKERRQ(ierr);
	ierr = PetibmGridSetBoundaryPoints(
		gridBCtx.starts, gridBCtx.ends, gridB); CHKERRQ(ierr);
	// Create the field B
	ierr = PetibmFieldGetOptions("fieldB_", &fieldBCtx); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(fieldBCtx, gridB, fieldB); CHKERRQ(ierr);
	ierr = PetibmFieldSetBoundaryPoints(
		fieldBCtx.bc_value, fieldB); CHKERRQ(ierr);

	// Interpolate field A (defined on grid A) onto field B (defined on grid B)
	ierr = PetibmFieldInterpolate(gridA, fieldA, gridB, fieldB); CHKERRQ(ierr);

	std::string directory;
	ierr = PetibmGetParentDirectory(fieldBCtx.path, directory); CHKERRQ(ierr);
	ierr = PetibmCreateDirectory(directory); CHKERRQ(ierr);
#ifdef PETIBM_0_2
	ierr = PetibmFieldWrite(PETSC_COMM_WORLD, fieldBCtx.path, fieldBCtx.name,
	                        outViewerType, fieldB); CHKERRQ(ierr);
#else
	ierr = PetibmFieldHDF5Write(
		PETSC_COMM_WORLD, fieldBCtx.path, fieldBCtx.name, fieldB); CHKERRQ(ierr);
#endif

	ierr = PetibmFieldDestroy(fieldA); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(fieldB); CHKERRQ(ierr);
	ierr = PetibmGridDestroy(gridA); CHKERRQ(ierr);
	ierr = PetibmGridDestroy(gridB); CHKERRQ(ierr);
	for (PetscInt d=0; d<gridA.dim; d++)
	{
		ierr = VecDestroy(coords+d); CHKERRQ(ierr);
	}
	ierr = PetscFinalize(); CHKERRQ(ierr);

	return 0;
} // main
