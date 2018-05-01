/*! Computes the vorticity field from the 3D velocity vector field.
 * \file vorticity3d.cpp
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


static char help[] = "petibm-vorticity3d (0.1.0)\n\n" \
"Computes the x- and/or z-components " \
"of the 3D vorticity field.\n\n" \
"Usage: petibm-vorticity3d [arg]\n" \
"Options and arguments:\n" \
"  -config_file <path>\tInsert options and arguments from a given file\n" \
"  -data_directory <path>\tData directory [default='.']\n" \
"  -output_directory <path>\tOutput directory [default='output']\n" \
"  -grid_directory <path>\tDirectory with HDF5 grids (PetIBM-0.2) [default='grids']\n" \
"  -grid_path <path>\tPath of the HDF5 grid file (PetIBM-0.3) [default='grid.h5']\n" \
"  -compute_wx\tCompute the x-component of the vorticity\n" \
"  -compute_wz\tCompute the z-component of the vorticity\n" \
"  -input_binary\tRead the velocity written in PETSc binary format (PetIBM-0.2)\n" \
"  -output_binary\tWrite the vorticity in PETSc binary format (PetIBM-0.2)\n" \
"  -nstart <int>\tStarting time-step\n" \
"  -nend <int>\tEnding time-step\n" \
"  -nstep <int>\tTime-step increment\n" \
"  -nx <int>\tNumber of cells in the x-direction\n" \
"  -ny <int>\tNumber of cells in the y-direction\n" \
"  -nz <int>\tNumber of cells in the z-direction\n" \
"  -periodic_x\tUse periodic conditions in the x-direction\n" \
"  -periodic_y\tUse periodic conditions in the y-direction\n" \
"  -periodic_z\tUse periodic conditions in the z-direction\n" \
"\n"
;


int main(int argc, char **argv)
{
	PetscErrorCode ierr;
	std::string datadir, outdir, griddir, gridpath;
	PetibmGrid gridux, griduy, griduz, gridwx, gridwz;
	PetibmGridCtx gridCtx, griduxCtx, griduyCtx, griduzCtx, gridwxCtx, gridwzCtx;
	PetibmField ux, uy, uz, wx, wz;
	PetibmFieldCtx fieldCtx;
	PetibmTimeStepCtx stepCtx;
	DM da;
	PetscInt ite;
	PetscMPIInt rank;
	PetscBool found = PETSC_FALSE,
	          compute_wx = PETSC_FALSE,
	          compute_wz = PETSC_FALSE;

	ierr = PetscInitialize(&argc, &argv, nullptr, help); CHKERRQ(ierr);

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

	// parse command-line options
	ierr = PetibmOptionsInsertFile(nullptr); CHKERRQ(ierr);
	ierr = PetibmGetDirectory(
		&datadir, "-data_directory", "."); CHKERRQ(ierr);
	ierr = PetibmGetDirectory(
		&outdir, "-output_directory", "output", PETSC_TRUE); CHKERRQ(ierr);
	ierr = PetibmGridGetOptions(nullptr, &gridCtx); CHKERRQ(ierr);
	ierr = PetibmFieldGetOptions(nullptr, &fieldCtx); CHKERRQ(ierr);
	ierr = PetibmTimeStepGetOptions(nullptr, &stepCtx); CHKERRQ(ierr);
#ifdef PETIBM_0_2
	ierr = PetibmGetDirectory(
		&griddir, "-grid_directory", datadir+"/grids"); CHKERRQ(ierr);
#else
	ierr = PetibmGetFilePath(
		&gridpath, "-grid_path", datadir+"/grid.h5"); CHKERRQ(ierr);
#endif
	ierr = PetscOptionsGetBool(
		nullptr, nullptr, "-compute_wx", &compute_wx, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(
		nullptr, nullptr, "-compute_wz", &compute_wz, &found); CHKERRQ(ierr);
#ifdef PETIBM_0_2
	// get the input and output formats
	PetscBool input_binary = PETSC_FALSE,
	          output_binary = PETSC_FALSE;
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

	// read staggered gridline stations for x-velocity
	griduxCtx.nx = (fieldCtx.periodic_x) ? gridCtx.nx : gridCtx.nx-1;
	griduxCtx.ny = gridCtx.ny;
	griduxCtx.nz = gridCtx.nz;
	ierr = PetibmGridCreateSeq(griduxCtx, gridux); CHKERRQ(ierr);
#ifdef PETIBM_0_2
	gridpath = griddir+"/staggered-x.h5";
#endif
	ierr = PetibmGridHDF5Read(
		PETSC_COMM_SELF, gridpath, "u", gridux); CHKERRQ(ierr);
	// read staggered gridline stations for y-velocity
	griduyCtx.nx = gridCtx.nx;
	griduyCtx.ny = (fieldCtx.periodic_y) ? gridCtx.ny : gridCtx.ny-1;
	griduyCtx.nz = gridCtx.nz;
	ierr = PetibmGridCreateSeq(griduyCtx, griduy); CHKERRQ(ierr);
#ifdef PETIBM_0_2
	gridpath = griddir+"/staggered-y.h5";
#endif
	ierr = PetibmGridHDF5Read(
		PETSC_COMM_SELF, gridpath, "v", griduy); CHKERRQ(ierr);
	// read staggered gridline stations for z-velocity
	griduzCtx.nx = gridCtx.nx;
	griduzCtx.ny = gridCtx.ny;
	griduzCtx.nz = (fieldCtx.periodic_z) ? gridCtx.nz : gridCtx.nz-1;
	ierr = PetibmGridCreateSeq(griduzCtx, griduz); CHKERRQ(ierr);
#ifdef PETIBM_0_2
	gridpath = griddir+"/staggered-z.h5";
#endif
	ierr = PetibmGridHDF5Read(
		PETSC_COMM_SELF, gridpath, "w", griduz); CHKERRQ(ierr);
	// create grid for x-vorticity
	if (compute_wx)
	{
		gridwxCtx.nx = gridCtx.nx;
		gridwxCtx.ny = gridCtx.ny - 1;
		gridwxCtx.nz = gridCtx.nz - 1;
		ierr = PetibmGridCreateSeq(gridwxCtx, gridwx); CHKERRQ(ierr);
		ierr = PetibmVorticityXComputeGrid(griduy, griduz, gridwx); CHKERRQ(ierr);
		if (rank == 0)
		{
			std::string filepath = outdir+"/grid.h5";
#ifdef PETIBM_0_2
			filepath = outdir+"/grids/wx.h5";
			ierr = PetibmCreateDirectory(outdir+"/grids"); CHKERRQ(ierr);
#endif
			ierr = PetibmGridHDF5Write(
				PETSC_COMM_SELF, filepath, "wx", gridwx); CHKERRQ(ierr);
		}
	}
	// create grid for z-vorticity
	if (compute_wz)
	{
		gridwzCtx.nx = gridCtx.nx - 1;
		gridwzCtx.ny = gridCtx.ny - 1;
		gridwzCtx.nz = gridCtx.nz;
		ierr = PetibmGridCreateSeq(gridwzCtx, gridwz); CHKERRQ(ierr);
		ierr = PetibmVorticityZComputeGrid(gridux, griduy, gridwz); CHKERRQ(ierr);
		if (rank == 0)
		{
			std::string filepath = outdir+"/grid.h5";
#ifdef PETIBM_0_2
			filepath = outdir+"/grids/wz.h5";
			ierr = PetibmCreateDirectory(outdir+"/grids"); CHKERRQ(ierr);
#endif
			ierr = PetibmGridHDF5Write(
				PETSC_COMM_SELF, filepath, "wz", gridwz); CHKERRQ(ierr);
		}
	}
	// create base DMDA object
	ierr = PetibmFieldDMDACreate3d(gridCtx, fieldCtx, da); CHKERRQ(ierr);
	// initialize field for ux velocity
	ierr = PetibmFieldDMDACreate3d("ux", da, ux.da); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(ux.da, &ux.global); CHKERRQ(ierr);
	// initialize field for uy velocity
	ierr = PetibmFieldDMDACreate3d("uy", da, uy.da); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(uy.da, &uy.global); CHKERRQ(ierr);
	// initialize field for uz velocity
	ierr = PetibmFieldDMDACreate3d("uz", da, uz.da); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(uz.da, &uz.global); CHKERRQ(ierr);
	// initialize field for wx vorticity
	if (compute_wx)
	{
		ierr = PetibmFieldDMDACreate3d("wx", da, wx.da); CHKERRQ(ierr);
		ierr = DMCreateGlobalVector(wx.da, &wx.global); CHKERRQ(ierr);
	}
	// initialize field for wz vorticity
	if (compute_wz)
	{
		ierr = PetibmFieldDMDACreate3d("wz", da, wz.da); CHKERRQ(ierr);
		ierr = DMCreateGlobalVector(wz.da, &wz.global); CHKERRQ(ierr);
	}
	ierr = DMDestroy(&da); CHKERRQ(ierr);

	// loop over the time steps to compute the z-vorticity
	for (ite=stepCtx.start; ite<=stepCtx.end; ite+=stepCtx.step)
	{
		ierr = PetscPrintf(
			PETSC_COMM_WORLD, "[time-step %d]\n", ite); CHKERRQ(ierr);
#ifndef PETIBM_0_2
		// get name of time-step file
		std::stringstream ss;
		ss << std::setfill('0') << std::setw(7) << ite << ".h5";
		std::string filename(ss.str());
		// read velocity field
		ierr = PetibmFieldHDF5Read(
			PETSC_COMM_WORLD, datadir+"/"+filename, "v", uy); CHKERRQ(ierr);
		if (compute_wx)
		{
			ierr = PetibmFieldHDF5Read(
				PETSC_COMM_WORLD, datadir+"/"+filename, "w", uz); CHKERRQ(ierr);
			ierr = PetibmVorticityXComputeField(
				griduy, griduz, uy, uz, wx); CHKERRQ(ierr);
			ierr = PetibmFieldHDF5Write(
				PETSC_COMM_WORLD, outdir+"/"+filename, "wx", wx); CHKERRQ(ierr);
		}
		if (compute_wz)
		{
			ierr = PetibmFieldHDF5Read(
				PETSC_COMM_WORLD, datadir+"/"+filename, "u", ux); CHKERRQ(ierr);
			ierr = PetibmVorticityZComputeField(
				gridux, griduy, ux, uy, wz); CHKERRQ(ierr);
			ierr = PetibmFieldHDF5Write(
				PETSC_COMM_WORLD, outdir+"/"+filename, "wz", wz); CHKERRQ(ierr);
		}
#else
		// get name of time-step directory
		std::stringstream ss;
		ss << datadir << "/" << std::setfill('0') << std::setw(7) << ite;
		std::string folder(ss.str());
		// get time-step directory to save
		std::stringstream ssout;
		ssout << outdir << "/" << std::setfill('0') << std::setw(7) << ite;
		std::string outfolder(ssout.str());
		ierr = PetibmCreateDirectory(outfolder); CHKERRQ(ierr);
		// read velocity field
		ierr = PetibmFieldRead(PETSC_COMM_WORLD, folder+"/uy"+inExt,
		                       "uy", inViewerType, uy); CHKERRQ(ierr);
		if (compute_wx)
		{
			ierr = PetibmFieldRead(PETSC_COMM_WORLD, folder+"/uz"+inExt,
			                       "uz", inViewerType, uz); CHKERRQ(ierr);
			ierr = PetibmVorticityXComputeField(
				griduy, griduz, uy, uz, wx); CHKERRQ(ierr);
			ierr = PetibmFieldWrite(PETSC_COMM_WORLD, outfolder+"/wx"+outExt,
			                        "wx", outViewerType, wx); CHKERRQ(ierr);
		}
		if (compute_wz)
		{
			ierr = PetibmFieldRead(PETSC_COMM_WORLD, folder+"/ux"+inExt,
			                       "ux", inViewerType, ux); CHKERRQ(ierr);
			ierr = PetibmVorticityZComputeField(
				gridux, griduy, ux, uy, wz); CHKERRQ(ierr);
			ierr = PetibmFieldWrite(PETSC_COMM_WORLD, outfolder+"/wz"+outExt,
			                        "wz", outViewerType, wz); CHKERRQ(ierr);
		}
#endif
	}

	ierr = PetibmGridDestroy(gridux); CHKERRQ(ierr);
	ierr = PetibmGridDestroy(griduy); CHKERRQ(ierr);
	ierr = PetibmGridDestroy(griduz); CHKERRQ(ierr);
	ierr = PetibmGridDestroy(gridwx); CHKERRQ(ierr);
	ierr = PetibmGridDestroy(gridwz); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(ux); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(uy); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(uz); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(wx); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(wz); CHKERRQ(ierr);
	
	ierr = PetscFinalize(); CHKERRQ(ierr);

	return 0;
} // main
