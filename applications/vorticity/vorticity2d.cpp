/*! Computes the vorticity field from the 2D velocity vector field.
 * \file vorticity2d.cpp
 */

#include <iomanip>
#include <petscsys.h>
#include <petscdmda.h>

#include "petibm-utilities/field.hpp"
#include "petibm-utilities/vorticity.hpp"


int main(int argc, char **argv)
{
	PetscErrorCode ierr;

	ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);

	// parse command-line options
	std::string directory(".");
	char dir[PETSC_MAX_PATH_LEN];
	PetscBool found;
	ierr = PetscOptionsGetString(
		nullptr, nullptr, "-directory", dir, sizeof(dir), &found); CHKERRQ(ierr);
	if (found)
	  directory = dir;

	PetscInt nstart = 0,
	         nend = 0,
	         nstep = 1;
	ierr = PetscOptionsGetInt(
		nullptr, nullptr, "-nstart", &nstart, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(
		nullptr, nullptr, "-nend", &nend, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(
		nullptr, nullptr, "-nstep", &nstep, &found); CHKERRQ(ierr);

	PetscInt nx = 0,
	         ny = 0;
	ierr = PetscOptionsGetInt(
		nullptr, nullptr, "-nx", &nx, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(
		nullptr, nullptr, "-ny", &ny, &found); CHKERRQ(ierr);

	PetscBool isPeriodic_x, isPeriodic_y;
	ierr = PetscOptionsGetBool(
		nullptr, nullptr, "-periodic_x", &isPeriodic_x, nullptr); CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(
		nullptr, nullptr, "-periodic_y", &isPeriodic_y, nullptr); CHKERRQ(ierr);

	DMBoundaryType bType_x,
	               bType_y;
	bType_x = (isPeriodic_x) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	bType_y = (isPeriodic_y) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	
	// create DMDA for phi
	Field phi;
	ierr = DMDACreate2d(PETSC_COMM_WORLD,
	                    bType_x, bType_y,
	                    DMDA_STENCIL_STAR,
	                    nx, ny,
	                    PETSC_DECIDE, PETSC_DECIDE, 1, 1, nullptr, nullptr,
	                    &phi.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) phi.da, nullptr, "-phi_dmda_view"); CHKERRQ(ierr);

	// create DMDA objects for velocity components from DMDA object for pressure
	Field ux, uy;
	PetscInt numX, numY;
	const PetscInt *plx, *ply;
	ierr = DMDAGetOwnershipRanges(phi.da, &plx, &ply, nullptr); CHKERRQ(ierr);
	PetscInt m, n;
	ierr = DMDAGetInfo(phi.da,
	                   nullptr, nullptr, nullptr, nullptr,
	                   &m, &n,
	                   nullptr, nullptr, nullptr, nullptr,
	                   nullptr, nullptr, nullptr); CHKERRQ(ierr);
	// x-component of velocity
	PetscInt *ulx, *uly;
	ierr = PetscMalloc(m*sizeof(*ulx), &ulx); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*uly), &uly); CHKERRQ(ierr);
	ierr = PetscMemcpy(ulx, plx, m*sizeof(*ulx)); CHKERRQ(ierr);
	ierr = PetscMemcpy(uly, ply, n*sizeof(*uly)); CHKERRQ(ierr);
	numX = nx;
	numY = ny;
	if (!isPeriodic_x)
	{
	  ulx[m-1]--;
	  numX--;
	}
	ierr = DMDACreate2d(PETSC_COMM_WORLD, 
	                    bType_x, bType_y,
	                    DMDA_STENCIL_BOX, 
	                    numX, numY, m, n, 1, 1, ulx, uly,
	                    &ux.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) ux.da, nullptr, "-ux_dmda_view"); CHKERRQ(ierr);
	ierr = PetscFree(ulx); CHKERRQ(ierr);
	ierr = PetscFree(uly); CHKERRQ(ierr);
	// y-component of velocity
	PetscInt *vlx, *vly;
	ierr = PetscMalloc(m*sizeof(*vlx), &vlx); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*vly), &vly); CHKERRQ(ierr);
	ierr = PetscMemcpy(vlx, plx, m*sizeof(*vlx)); CHKERRQ(ierr);
	ierr = PetscMemcpy(vly, ply, n*sizeof(*vly)); CHKERRQ(ierr);
	numX = nx;
	numY = ny;
	if (!isPeriodic_y)
	{
	  vly[n-1]--;
	  numY--;
	}
	ierr = DMDACreate2d(PETSC_COMM_WORLD, 
	                    bType_x, bType_y,
	                    DMDA_STENCIL_BOX, 
	                    numX, numY, m, n, 1, 1, vlx, vly,
	                    &uy.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) uy.da, nullptr, "-uy_dmda_view"); CHKERRQ(ierr);
	ierr = PetscFree(vlx); CHKERRQ(ierr);
	ierr = PetscFree(vly); CHKERRQ(ierr);

	ierr = FieldInitialize(ux); CHKERRQ(ierr);
	ierr = FieldInitialize(uy); CHKERRQ(ierr);
	ierr = FieldInitialize(phi); CHKERRQ(ierr);

	// read grids
	std::string gridsDirectory(directory + "/grids");
	ierr = FieldReadGrid(gridsDirectory + "/staggered-x.h5", ux); CHKERRQ(ierr);
	ierr = FieldReadGrid(gridsDirectory + "/staggered-y.h5", uy); CHKERRQ(ierr);
	ierr = FieldReadGrid(gridsDirectory + "/cell-centered.h5", phi); CHKERRQ(ierr);

	// create z-vorticity field
	Field wz;
	ierr = DMDACreate2d(PETSC_COMM_WORLD,
	                    bType_x, bType_y,
	                    DMDA_STENCIL_STAR,
	                    nx-1, ny-1,
	                    PETSC_DECIDE, PETSC_DECIDE, 1, 1, nullptr, nullptr,
	                    &wz.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) wz.da, nullptr, "-wz_dmda_view"); CHKERRQ(ierr);
	ierr = FieldInitialize(wz); CHKERRQ(ierr);
	ierr = ComputeGridVorticityZ(ux, uy, wz); CHKERRQ(ierr);

	PetscMPIInt rank;
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
	if (rank == 0)
	{
	  ierr = FieldWriteGrid(gridsDirectory + "/wz.h5", wz); CHKERRQ(ierr);
	}

	for (PetscInt ite=nstart; ite<=nend; ite+=nstep)
	{
	  ierr = PetscPrintf(PETSC_COMM_WORLD, "[time-step %D]\n", ite); CHKERRQ(ierr);

	  std::stringstream ss;
	  ss << directory << "/" << std::setfill('0') << std::setw(7) << ite;
	  std::string folder(ss.str());

	  // read values
	  ierr = FieldReadValues(folder + "/ux.h5", "ux", ux); CHKERRQ(ierr);
	  ierr = FieldReadValues(folder + "/uy.h5", "uy", uy); CHKERRQ(ierr);

	  // compute the z-vorticity
	  ierr = ComputeVorticityZ(ux, uy, wz); CHKERRQ(ierr);
	  ierr = FieldWriteValues(folder + "/wz.h5", "wz", wz); CHKERRQ(ierr);
	}

	ierr = FieldDestroy(ux); CHKERRQ(ierr);
	ierr = FieldDestroy(uy); CHKERRQ(ierr);
	ierr = FieldDestroy(phi); CHKERRQ(ierr);
	ierr = FieldDestroy(wz); CHKERRQ(ierr);
	ierr = PetscFinalize(); CHKERRQ(ierr);
	return 0;
} // main
