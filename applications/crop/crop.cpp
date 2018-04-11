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


int main(int argc, char **argv)
{
	PetscErrorCode ierr;
	PetibmField fieldA, fieldB;
	PetibmFieldCtx fieldACtx, fieldBCtx;
	PetibmGrid gridA, gridB;
	PetibmGridCtx gridACtx, gridBCtx, cropCtx;
	PetibmTimeStepCtx stepCtx;
	std::string outdir, datadir, filepath;

	ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);
	
	ierr = PetibmOptionsInsertFile(nullptr); CHKERRQ(ierr);
	ierr = PetibmGetDirectory(
		&datadir, "-data_directory", "solution"); CHKERRQ(ierr);
	ierr = PetibmGetDirectory(
		&outdir, "-output_directory", "output", PETSC_TRUE); CHKERRQ(ierr);
	ierr = PetibmTimeStepGetOptions(nullptr, &stepCtx); CHKERRQ(ierr);
	ierr = PetibmGridGetOptions(nullptr, &cropCtx); CHKERRQ(ierr);

	// Create and read the grid A
	ierr = PetibmGridGetOptions("gridA_", &gridACtx); CHKERRQ(ierr);
	ierr = PetibmGridCtxPrintf("Grid A", gridACtx); CHKERRQ(ierr);
	ierr = PetibmGridInitialize(gridACtx, gridA); CHKERRQ(ierr);
	ierr = PetibmGridHDF5Read(
		gridACtx.path, gridACtx.name, gridA); CHKERRQ(ierr);

	// Crop grid A to get sub grid B and write grid B
	std::strcpy(gridBCtx.name, gridACtx.name);
	ierr = PetibmGridCrop(gridA, cropCtx, gridB); CHKERRQ(ierr);
#ifndef PETIBM_0_2
	filepath = outdir + "/grid.h5";
#else
	std::string griddir = outdir + "/grids";
	mkdir(griddir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	filepath = griddir + "/" + gridBCtx.name + ".h5";
#endif
	ierr = PetibmGridHDF5Write(filepath, gridBCtx.name, gridB); CHKERRQ(ierr);

	// Create field A
	ierr = PetibmFieldGetOptions("fieldA_", &fieldACtx); CHKERRQ(ierr);
	ierr = PetibmFieldCtxPrintf("Field A", fieldACtx); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(fieldACtx, gridA, fieldA); CHKERRQ(ierr);

	// Create field B
	std::strcpy(fieldBCtx.name, fieldACtx.name);
	ierr = PetibmFieldInitialize(fieldBCtx, gridB, fieldB); CHKERRQ(ierr);

	for (PetscInt step=stepCtx.start; step<=stepCtx.end; step+=stepCtx.step)
	{
		ierr = PetscPrintf(
			PETSC_COMM_WORLD, "[time-step %d] Cropping...\n", step); CHKERRQ(ierr);
#ifndef PETIBM_0_2
		// Read field A from file
		std::stringstream ss;
		ss << std::setfill('0') << std::setw(7) << step << ".h5";
		filepath = datadir + "/" + ss.str();
#else
		// Read field A from file
		std::stringstream ss;
		ss << std::setfill('0') << std::setw(7) << step;
		filepath = datadir + "/" + ss.str() + "/" + fieldACtx.name + ".h5";
#endif
		ierr = PetibmFieldHDF5Read(
			filepath, fieldACtx.name, fieldA); CHKERRQ(ierr);
		// Crop field A to fill field B
		ierr = PetibmFieldCrop(gridA, fieldA, cropCtx, fieldB); CHKERRQ(ierr);
#ifndef PETIBM_0_2
		filepath = outdir + "/" + ss.str();
#else
		std::string folder = outdir + "/" + ss.str();
		mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		filepath = folder + "/" + fieldBCtx.name + ".h5";
#endif
		ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);
		ierr = PetibmFieldHDF5Write(filepath, fieldBCtx.name, fieldB); CHKERRQ(ierr);
	}

	// Clean workspace
	ierr = PetibmFieldDestroy(fieldA); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(fieldB); CHKERRQ(ierr);
	ierr = PetibmGridDestroy(gridA); CHKERRQ(ierr);
	ierr = PetibmGridDestroy(gridB); CHKERRQ(ierr);
	ierr = PetscFinalize(); CHKERRQ(ierr);

	return 0;
} // main
