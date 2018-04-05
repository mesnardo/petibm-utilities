/*! Computes the vorticity field from the 2D velocity vector field.
 * \file vorticity2d.cpp
 */

#include <iomanip>
#include <iostream>
#include <sys/stat.h>

#include <petscsys.h>
#include <petscdmda.h>

#include "petibm-utilities/field.h"
#include "petibm-utilities/grid.h"
#include "petibm-utilities/misc.h"
#include "petibm-utilities/timestep.h"
#include "petibm-utilities/vorticity.h"


int main(int argc, char **argv)
{
	PetscErrorCode ierr;
	std::string directory, datadir, outdir, griddir, gridpath;
	PetibmGrid grid, gridux, griduy, gridwz;
	PetibmGridCtx gridCtx, griduxCtx, griduyCtx, gridwzCtx;
	PetibmField ux, uy, wz;
	PetibmFieldCtx fieldCtx;
	PetibmTimeStepCtx stepCtx;
	DM da;
	PetscInt ite;
	PetscMPIInt rank;
	PetscBool found = PETSC_FALSE,
	          binary_format = PETSC_FALSE;

	ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

	// parse command-line options
	ierr = PetibmOptionsInsertFile(nullptr); CHKERRQ(ierr);
	directory = ".";
	ierr = PetibmGetDirectory(&directory, "-directory"); CHKERRQ(ierr);
	ierr = PetibmGetDirectory(&datadir, "-data_directory"); CHKERRQ(ierr);
	ierr = PetibmGetDirectory(&outdir, "-output_directory", PETSC_TRUE); CHKERRQ(ierr);
	ierr = PetibmGridGetOptions(nullptr, &gridCtx); CHKERRQ(ierr);
	ierr = PetibmFieldGetOptions(nullptr, &fieldCtx); CHKERRQ(ierr);
	ierr = PetibmTimeStepGetOptions(nullptr, &stepCtx); CHKERRQ(ierr);
#ifdef PETIBM_0_3
	gridpath = directory+"/grid.h5";
	ierr = PetibmGetFilePath(&gridpath, "-grid_path"); CHKERRQ(ierr);
#elif PETIBM_0_2
	griddir = directory+"/grids";
	ierr = PetibmGetDirectory(&griddir, "-grid_directory"); CHKERRQ(ierr);
#endif
	ierr = PetscOptionsGetBool(
		nullptr, nullptr, "-binary_format", &binary_format, &found); CHKERRQ(ierr);

	// read cell-centered gridline stations
	ierr = PetibmGridCreateSeq(gridCtx, grid); CHKERRQ(ierr);
#ifdef PETIBM_0_3
	ierr = PetibmGridHDF5Read(gridpath, "p", grid); CHKERRQ(ierr);
#elif PETIBM_0_2
	ierr = PetibmGridHDF5Read(griddir+"/cell-centered.h5", grid); CHKERRQ(ierr);
#endif
	// read staggered gridline stations for x-velocity
	griduxCtx.nx = (fieldCtx.periodic_x) ? gridCtx.nx : gridCtx.nx-1;
	griduxCtx.ny = gridCtx.ny;
	ierr = PetibmGridCreateSeq(griduxCtx, gridux); CHKERRQ(ierr);
#ifdef PETIBM_0_3
	ierr = PetibmGridHDF5Read(gridpath, "u", gridux); CHKERRQ(ierr);
#elif PETIBM_0_2
	ierr = PetibmGridHDF5Read(griddir+"/staggered-x.h5", gridux); CHKERRQ(ierr);
#endif
	// read staggered gridline stations for y-velocity
	griduyCtx.nx = gridCtx.nx;
	griduyCtx.ny = (fieldCtx.periodic_y) ? gridCtx.ny : gridCtx.ny-1;
	ierr = PetibmGridCreateSeq(griduyCtx, griduy); CHKERRQ(ierr);
#ifdef PETIBM_0_3
	ierr = PetibmGridHDF5Read(gridpath, "v", griduy); CHKERRQ(ierr);
#elif PETIBM_0_2
	ierr = PetibmGridHDF5Read(griddir+"/staggered-y.h5", griduy); CHKERRQ(ierr);
#endif
	// create grid for z-vorticity
	gridwzCtx.nx = gridCtx.nx - 1;
	gridwzCtx.ny = gridCtx.ny - 1;
	ierr = PetibmGridCreateSeq(gridwzCtx, gridwz); CHKERRQ(ierr);
	ierr = PetibmVorticityZComputeGrid(gridux, griduy, gridwz); CHKERRQ(ierr);
	if (rank == 0)
	{
#ifdef PETIBM_0_3
		ierr = PetibmGridHDF5Write(outdir+"/grid.h5", "wz", gridwz); CHKERRQ(ierr);
#elif PETIBM_0_2
		mkdir((outdir+"/grids").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		ierr = PetibmGridHDF5Write(outdir+"/grids/wz.h5", gridwz); CHKERRQ(ierr);
#endif
	}
	// create base DMDA object
	ierr = PetibmFieldDMDACreate2d(gridCtx, fieldCtx, da); CHKERRQ(ierr);
	// initialize field for ux velocity
	ierr = PetibmFieldDMDACreate2d("ux", da, ux.da); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(ux.da, &ux.global); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(ux.da, &ux.local); CHKERRQ(ierr);
	// initialize field for uy velocity
	ierr = PetibmFieldDMDACreate2d("uy", da, uy.da); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(uy.da, &uy.global); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(uy.da, &uy.local); CHKERRQ(ierr);
	// initialize field for wz vorticity
	ierr = PetibmFieldDMDACreate3d("wz", da, wz.da); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(wz.da, &wz.global); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(wz.da, &wz.local); CHKERRQ(ierr);

	// loop over the time steps to compute the z-vorticity
	for (ite=stepCtx.start; ite<=stepCtx.end; ite+=stepCtx.step)
	{
		ierr = PetscPrintf(
			PETSC_COMM_WORLD, "[time-step %d]\n", ite); CHKERRQ(ierr);
#ifdef PETIBM_0_3
		// get name of time-step file
		std::stringstream ss;
		ss << std::setfill('0') << std::setw(7) << ite << ".h5";
		std::string filename(ss.str());
		// read velocity field
		ierr = PetibmFieldHDF5Read(datadir+"/"+filename, "v", uy); CHKERRQ(ierr);
		ierr = PetibmFieldHDF5Read(datadir+"/"+filename, "u", ux); CHKERRQ(ierr);
		ierr = PetibmVorticityZComputeField(
			gridux, griduy, ux, uy, wz); CHKERRQ(ierr);
		ierr = PetibmFieldHDF5Write(outdir+"/"+filename, "wz", wz); CHKERRQ(ierr);
#elif PETIBM_0_2
		// get name of time-step directory
		std::stringstream ss;
		ss << directory << "/" << std::setfill('0') << std::setw(7) << ite;
		std::string folder(ss.str());
		// define if binary format or HDF5
		std::string extension = ".h5";
		PetscViewerType viewerType = PETSCVIEWERHDF5;
		if (binary_format == PETSC_TRUE)
		{
			extension = ".dat";
			viewerType = PETSCVIEWERBINARY;
		}
		// get time-step directory to save
		std::stringstream ssout;
		ssout << outdir << "/" << std::setfill('0') << std::setw(7) << ite;
		std::string outfolder(ssout.str());
		mkdir(outfolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		// read velocity field
		ierr = PetibmFieldRead(
			folder+"/uy"+extension, "uy", viewerType, uy); CHKERRQ(ierr);
		ierr = PetibmFieldRead(
			folder+"/ux"+extension, "ux", viewerType, ux); CHKERRQ(ierr);
		ierr = PetibmVorticityZComputeField(
			gridux, griduy, ux, uy, wz); CHKERRQ(ierr);
		ierr = PetibmFieldHDF5Write(outfolder+"/wz.h5", "wz", wz); CHKERRQ(ierr);
#endif
	}

	ierr = PetibmGridDestroy(gridux); CHKERRQ(ierr);
	ierr = PetibmGridDestroy(griduy); CHKERRQ(ierr);
	ierr = PetibmGridDestroy(gridwz); CHKERRQ(ierr);
	ierr = DMDestroy(&da); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(ux); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(uy); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(wz); CHKERRQ(ierr);
	
	ierr = PetscFinalize(); CHKERRQ(ierr);

	return 0;
} // main
