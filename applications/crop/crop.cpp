/*! Crops the solution in a given inner box.
 * \file crop.cpp
 */

#include <string>
#include <cstring>
#include <sys/stat.h>
#include <iomanip>

#include <petscsys.h>

#include "petibm-utilities/field.h"
#include "petibm-utilities/grid.h"
#include "petibm-utilities/misc.h"
#include "petibm-utilities/timestep.h"


static char help[] = "petibm-crop (0.1.0)\n\n" \
"Crop the numerical solution in a given sub-domain.\n\n" \
"Usage: petibm-crop [arg]\n" \
"Options and arguments:\n" \
"  -config_file <path>\tInsert options and arguments from a given file\n" \
"  -data_directory <path>\tData directory [default='solution']\n" \
"  -output_directory <path>\tOutput directory [default='output']\n" \
"  -gridA_name <string>\tName of the grid (PetIBM-0.3)\n" \
"  -gridA_path <string>\tPath of the grid file\n" \
"  -gridA_nx <int>\tNumber of cells in the x-direction\n" \
"  -gridA_ny <int>\tNumber of cells in the y-direction\n" \
"  -gridA_nz <int>\tNumber of cells in the z-direction\n" \
"  -fieldA_name <string>\tName of the field\n" \
"  -input_binary\tRead the PetIBM field written in PETSc binary format (PetIBM-0.2)\n" \
"  -output_binary\tWrite the PetIBM field in PETSc binary format (PetIBM-0.2)\n" \
"  -x_start <float>\tBottom-left corner x-coordinate of the sub-domain\n" \
"  -y_start <float>\tBottom-left corner y-coordinate of the sub-domain\n" \
"  -z_start <float>\tBottom-left corner z-coordinate of the sub-domain\n" \
"  -x_end <float>\tTop-right corner x-coordinate of the sub-domain\n" \
"  -y_end <float>\tTop-right corner y-coordinate of the sub-domain\n" \
"  -z_end <float>\tTop-right corner z-coordinate of the sub-domain\n" \
"  -nstart <int>\tStarting time-step\n" \
"  -nend <int>\tEnding time-step\n" \
"  -nstep <int>\tTime-step increment\n" \
"\n"
;


int main(int argc, char **argv)
{
	PetscErrorCode ierr;
	PetibmField fieldA, fieldB;
	PetibmFieldCtx fieldACtx, fieldBCtx;
	PetibmGrid gridA, gridB;
	PetibmGridCtx gridACtx, gridBCtx, cropCtx;
	PetibmTimeStepCtx stepCtx;
	std::string outdir, datadir, filepath;

	ierr = PetscInitialize(&argc, &argv, nullptr, help); CHKERRQ(ierr);
	
	ierr = PetibmOptionsInsertFile(nullptr); CHKERRQ(ierr);
	ierr = PetibmGetDirectory(
		&datadir, "-data_directory", "solution"); CHKERRQ(ierr);
	ierr = PetibmGetDirectory(
		&outdir, "-output_directory", "output", PETSC_TRUE); CHKERRQ(ierr);
	ierr = PetibmTimeStepGetOptions(nullptr, &stepCtx); CHKERRQ(ierr);
	ierr = PetibmGridGetOptions(nullptr, &cropCtx); CHKERRQ(ierr);
#ifdef PETIBM_0_2
	// get the input and output formats
	PetscBool input_binary = PETSC_FALSE,
	          output_binary = PETSC_FALSE,
	          found = PETSC_FALSE;
	ierr = PetscOptionsGetBool(
		nullptr, nullptr, "-input_binary", &input_binary, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(
		nullptr, nullptr, "-output_binary", &output_binary, &found); CHKERRQ(ierr);
	std::string inExt = ".h5",
	            outExt = ".h5";
	PetscViewerType inViewerType = PETSCVIEWERHDF5,
	                outViewerType = PETSCVIEWERHDF5;
	if (input_binary == PETSC_TRUE)
	{
		inExt = ".dat"; inViewerType = PETSCVIEWERBINARY;
	}
	if (output_binary == PETSC_TRUE)
	{
		outExt = ".dat"; outViewerType = PETSCVIEWERBINARY;
	}
#endif

	// Create and read the grid A
	ierr = PetibmGridGetOptions("gridA_", &gridACtx); CHKERRQ(ierr);
	ierr = PetibmGridInitialize(gridACtx, gridA); CHKERRQ(ierr);
	ierr = PetibmGridHDF5Read(
		PETSC_COMM_WORLD, gridACtx.path, gridACtx.name, gridA); CHKERRQ(ierr);

	// Crop grid A to get sub grid B and write grid B
	std::strcpy(gridBCtx.name, gridACtx.name);
	ierr = PetibmGridCrop(gridA, cropCtx, gridB); CHKERRQ(ierr);
#ifndef PETIBM_0_2
	filepath = outdir + "/grid.h5";
#else
	std::string griddir = outdir + "/grids";
	ierr = PetibmCreateDirectory(griddir); CHKERRQ(ierr);
	filepath = griddir + "/" + gridBCtx.name + ".h5";
#endif
	ierr = PetibmGridHDF5Write(
		PETSC_COMM_WORLD, filepath, gridBCtx.name, gridB); CHKERRQ(ierr);

	// Create field A
	ierr = PetibmFieldGetOptions("fieldA_", &fieldACtx); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(fieldACtx, gridA, fieldA); CHKERRQ(ierr);

	// Create field B
	std::strcpy(fieldBCtx.name, fieldACtx.name);
	ierr = PetibmFieldInitialize(fieldBCtx, gridB, fieldB); CHKERRQ(ierr);

	for (PetscInt step=stepCtx.start; step<=stepCtx.end; step+=stepCtx.step)
	{
		ierr = PetscPrintf(
			PETSC_COMM_WORLD, "[time-step %d]\n", step); CHKERRQ(ierr);
#ifndef PETIBM_0_2
		// Read field A from file
		std::stringstream ss;
		ss << std::setfill('0') << std::setw(7) << step << ".h5";
		filepath = datadir + "/" + ss.str();
		ierr = PetibmFieldHDF5Read(
			PETSC_COMM_WORLD, filepath, fieldACtx.name, fieldA); CHKERRQ(ierr);
#else
		// Read field A from file
		std::stringstream ss;
		ss << std::setfill('0') << std::setw(7) << step;
		filepath = datadir + "/" + ss.str() + "/" + fieldACtx.name + inExt;
		ierr = PetibmFieldRead(PETSC_COMM_WORLD, filepath, fieldACtx.name,
		                       inViewerType, fieldA); CHKERRQ(ierr);
#endif
		// Crop field A to fill field B
		ierr = PetibmFieldCrop(gridA, fieldA, cropCtx, fieldB); CHKERRQ(ierr);
#ifndef PETIBM_0_2
		filepath = outdir + "/" + ss.str();
		ierr = PetibmFieldHDF5Write(
			PETSC_COMM_WORLD, filepath, fieldBCtx.name, fieldB); CHKERRQ(ierr);
#else
		std::string folder = outdir + "/" + ss.str();
		ierr = PetibmCreateDirectory(folder); CHKERRQ(ierr);
		filepath = folder + "/" + fieldBCtx.name + outExt;
		ierr = PetibmFieldWrite(PETSC_COMM_WORLD, filepath, fieldBCtx.name,
		                        outViewerType, fieldB); CHKERRQ(ierr);
#endif
	}

	// Clean workspace
	ierr = PetibmFieldDestroy(fieldA); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(fieldB); CHKERRQ(ierr);
	ierr = PetibmGridDestroy(gridA); CHKERRQ(ierr);
	ierr = PetibmGridDestroy(gridB); CHKERRQ(ierr);
	ierr = PetscFinalize(); CHKERRQ(ierr);

	return 0;
} // main
