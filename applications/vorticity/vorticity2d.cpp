/*! Computes the vorticity field from the 2D velocity vector field.
 * \file vorticity2d.cpp
 */

#include <iomanip>
#include <iostream>

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

	ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);

	// parse command-line options
	std::string directory;
	ierr = PetibmGetDirectory(&directory); CHKERRQ(ierr);
	PetibmTimeStepCtx stepCtx;
	ierr = PetibmTimeStepsGetOptions(nullptr, &stepCtx); CHKERRQ(ierr);
	PetibmGridCtx gridCtx;
	ierr = PetibmGridGetOptions(nullptr, &gridCtx); CHKERRQ(ierr);

	DMBoundaryType bType_x, bType_y;
	bType_x = (gridCtx.periodic_x) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	bType_y = (gridCtx.periodic_y) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;

	// create DMDA for phi
	PetibmField phi;
	ierr = DMDACreate2d(PETSC_COMM_WORLD,
	                    bType_x, bType_y,
	                    DMDA_STENCIL_STAR,
	                    gridCtx.nx, gridCtx.ny,
	                    PETSC_DECIDE, PETSC_DECIDE, 1, 1, nullptr, nullptr,
	                    &phi.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) phi.da, nullptr, "-phi_dmda_view"); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(phi); CHKERRQ(ierr);
	ierr = PetibmGridInitialize(gridCtx, phi.grid); CHKERRQ(ierr);
	ierr = PetibmGridHDF5Read(
		directory + "/grids/cell-centered.h5", phi.grid); CHKERRQ(ierr);

	// create DMDA objects for velocity components from DMDA object for pressure
	const PetscInt *plx, *ply;
	PetscInt *lx, *ly;
	ierr = DMDAGetOwnershipRanges(phi.da, &plx, &ply, nullptr); CHKERRQ(ierr);
	PetscInt m, n;
	ierr = DMDAGetInfo(phi.da,
	                   nullptr, nullptr, nullptr, nullptr,
	                   &m, &n,
	                   nullptr, nullptr, nullptr, nullptr,
	                   nullptr, nullptr, nullptr); CHKERRQ(ierr);
	// x-component of velocity
	PetibmField ux;
	ierr = PetscMalloc(m*sizeof(*lx), &lx); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*ly), &ly); CHKERRQ(ierr);
	ierr = PetscMemcpy(lx, plx, m*sizeof(*lx)); CHKERRQ(ierr);
	ierr = PetscMemcpy(ly, ply, n*sizeof(*ly)); CHKERRQ(ierr);
	if (!gridCtx.periodic_x)
	{
	  lx[m-1]--;
	  gridCtx.nx--;
	}
	ierr = DMDACreate2d(PETSC_COMM_WORLD, 
	                    bType_x, bType_y,
	                    DMDA_STENCIL_BOX, 
	                    gridCtx.nx, gridCtx.ny, m, n, 1, 1, lx, ly,
	                    &ux.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) ux.da, nullptr, "-ux_dmda_view"); CHKERRQ(ierr);
	ierr = PetscFree(lx); CHKERRQ(ierr);
	ierr = PetscFree(ly); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(ux); CHKERRQ(ierr);
	ierr = PetibmGridInitialize(gridCtx, ux.grid); CHKERRQ(ierr);
	if (!gridCtx.periodic_x)
		gridCtx.nx++;
	ierr = PetibmGridHDF5Read(
		directory + "/grids/staggered-x.h5", ux.grid); CHKERRQ(ierr);
	
	// y-component of velocity
	PetibmField uy;
	ierr = PetscMalloc(m*sizeof(*lx), &lx); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*ly), &ly); CHKERRQ(ierr);
	ierr = PetscMemcpy(lx, plx, m*sizeof(*lx)); CHKERRQ(ierr);
	ierr = PetscMemcpy(ly, ply, n*sizeof(*ly)); CHKERRQ(ierr);
	if (!gridCtx.periodic_y)
	{
	  ly[n-1]--;
	  gridCtx.ny--;
	}
	ierr = DMDACreate2d(PETSC_COMM_WORLD, 
	                    bType_x, bType_y,
	                    DMDA_STENCIL_BOX, 
	                    gridCtx.nx, gridCtx.ny, m, n, 1, 1, lx, ly,
	                    &uy.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) uy.da, nullptr, "-uy_dmda_view"); CHKERRQ(ierr);
	ierr = PetscFree(lx); CHKERRQ(ierr);
	ierr = PetscFree(ly); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(ux); CHKERRQ(ierr);
	ierr = PetibmGridInitialize(gridCtx, uy.grid); CHKERRQ(ierr);
	if (!gridCtx.periodic_y)
		gridCtx.ny++;
	ierr = PetibmGridHDF5Read(
		directory + "/grids/staggered-y.h5", uy.grid); CHKERRQ(ierr);

	// create z-vorticity field
	PetibmField wz;
	ierr = PetscMalloc(m*sizeof(*lx), &lx); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*ly), &ly); CHKERRQ(ierr);
	ierr = PetscMemcpy(lx, plx, m*sizeof(*lx)); CHKERRQ(ierr);
	ierr = PetscMemcpy(ly, ply, n*sizeof(*ly)); CHKERRQ(ierr);
	lx[m-1]--;
	ly[n-1]--;
	gridCtx.nx--;
	gridCtx.ny--;
	ierr = DMDACreate2d(PETSC_COMM_WORLD,
	                    bType_x, bType_y,
	                    DMDA_STENCIL_STAR,
	                    gridCtx.nx, gridCtx.ny, m, n, 1, 1, lx, ly,
	                    &wz.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) wz.da, nullptr, "-wz_dmda_view"); CHKERRQ(ierr);
	ierr = PetscFree(lx); CHKERRQ(ierr);
	ierr = PetscFree(ly); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(wz); CHKERRQ(ierr);
	ierr = PetibmGridInitialize(gridCtx, wz.grid); CHKERRQ(ierr);
	gridCtx.nx++;
	gridCtx.ny++;
	ierr = PetibmVorticityZComputeGrid(ux.grid, uy.grid, wz.grid); CHKERRQ(ierr);
	ierr = PetibmGridHDF5Write(
		directory + "/grids/wz.h5", wz.grid); CHKERRQ(ierr);

	for (PetscInt ite=stepCtx.start; ite<=stepCtx.end; ite+=stepCtx.step)
	{
	  ierr = PetscPrintf(
			PETSC_COMM_WORLD, "[time-step %D]\n", ite); CHKERRQ(ierr);
	  // get time-step directory
	  std::stringstream ss;
	  ss << directory << "/" << std::setfill('0') << std::setw(7) << ite;
	  std::string folder(ss.str());
	  // read values
	  ierr = PetibmFieldHDF5Read(folder + "/ux.h5", "ux", ux); CHKERRQ(ierr);
	  ierr = PetibmFieldHDF5Read(folder + "/uy.h5", "uy", uy); CHKERRQ(ierr);
	  // compute the z-vorticity
	  ierr = PetibmVorticityZComputeField(ux, uy, wz); CHKERRQ(ierr);
	  ierr = PetibmFieldHDF5Write(folder + "/wz.h5", "wz", wz); CHKERRQ(ierr);
	}

	ierr = PetibmFieldDestroy(ux); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(uy); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(phi); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(wz); CHKERRQ(ierr);
	ierr = PetscFinalize(); CHKERRQ(ierr);
	return 0;
} // main
