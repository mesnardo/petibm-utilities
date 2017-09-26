/*! Computes the vorticity field from the 3D velocity vector field.
 * \file vorticity3d.cpp
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

	DMBoundaryType bType_x, bType_y, bType_z;
	bType_x = (gridCtx.periodic_x) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	bType_y = (gridCtx.periodic_y) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	bType_z = (gridCtx.periodic_z) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;

	// create DMDA for phi
	PetibmField phi;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
	                    bType_x, bType_y, bType_z,
	                    DMDA_STENCIL_STAR,
	                    gridCtx.nx, gridCtx.ny, gridCtx.nz,
	                    PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 1, 1,
	                    nullptr, nullptr, nullptr,
	                    &phi.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) phi.da, nullptr, "-phi_dmda_view"); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(phi); CHKERRQ(ierr);
	ierr = PetibmGridInitialize(gridCtx, phi.grid); CHKERRQ(ierr);
	ierr = PetibmGridHDF5Read(
		directory + "/grids/cell-centered.h5", phi.grid); CHKERRQ(ierr);

	// create DMDA objects for velocity components from DMDA object for pressure
	const PetscInt *plx, *ply, *plz;
	PetscInt *lx, *ly, *lz;
	ierr = DMDAGetOwnershipRanges(phi.da, &plx, &ply, &plz); CHKERRQ(ierr);
	PetscInt m, n, p;
	ierr = DMDAGetInfo(phi.da,
	                   nullptr, nullptr, nullptr, nullptr,
	                   &m, &n, &p,
	                   nullptr, nullptr, nullptr,
	                   nullptr, nullptr, nullptr); CHKERRQ(ierr);
	
	// x-component of velocity
	PetibmField ux;
	ierr = PetscMalloc(m*sizeof(*lx), &lx); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*ly), &ly); CHKERRQ(ierr);
	ierr = PetscMalloc(p*sizeof(*lz), &lz); CHKERRQ(ierr);
	ierr = PetscMemcpy(lx, plx, m*sizeof(*lx)); CHKERRQ(ierr);
	ierr = PetscMemcpy(ly, ply, n*sizeof(*ly)); CHKERRQ(ierr);
	ierr = PetscMemcpy(lz, plz, p*sizeof(*lz)); CHKERRQ(ierr);
	if (!gridCtx.periodic_x)
	{
	  lx[m-1]--;
	  gridCtx.nx--;
	}
	ierr = DMDACreate3d(PETSC_COMM_WORLD, 
	                    bType_x, bType_y, bType_z, 
	                    DMDA_STENCIL_BOX, 
	                    gridCtx.nx, gridCtx.ny, gridCtx.nz,
	                    m, n, p, 1, 1, lx, ly, lz, 
	                    &ux.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) ux.da, nullptr, "-ux_dmda_view"); CHKERRQ(ierr);
	ierr = PetscFree(lx); CHKERRQ(ierr);
	ierr = PetscFree(ly); CHKERRQ(ierr);
	ierr = PetscFree(lz); CHKERRQ(ierr);
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
	ierr = PetscMalloc(p*sizeof(*lz), &lz); CHKERRQ(ierr);
	ierr = PetscMemcpy(lx, plx, m*sizeof(*lx)); CHKERRQ(ierr);
	ierr = PetscMemcpy(ly, ply, n*sizeof(*ly)); CHKERRQ(ierr);
	ierr = PetscMemcpy(lz, plz, p*sizeof(*lz)); CHKERRQ(ierr);
	if (!gridCtx.periodic_y)
	{
	  ly[n-1]--;
	  gridCtx.ny--;
	}
	ierr = DMDACreate3d(PETSC_COMM_WORLD, 
	                    bType_x, bType_y, bType_z, 
	                    DMDA_STENCIL_BOX, 
	                    gridCtx.nx, gridCtx.ny, gridCtx.nz,
	                    m, n, p, 1, 1, lx, ly, lz, 
	                    &uy.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) uy.da, nullptr, "-uy_dmda_view"); CHKERRQ(ierr);
	ierr = PetscFree(lx); CHKERRQ(ierr);
	ierr = PetscFree(ly); CHKERRQ(ierr);
	ierr = PetscFree(lz); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(uy); CHKERRQ(ierr);
	ierr = PetibmGridInitialize(gridCtx, uy.grid); CHKERRQ(ierr);
	if (!gridCtx.periodic_y)
		gridCtx.ny++;
	ierr = PetibmGridHDF5Read(
		directory + "/grids/staggered-y.h5", uy.grid); CHKERRQ(ierr);
	
	// z-component of velocity
	PetibmField uz;
	ierr = PetscMalloc(m*sizeof(*lx), &lx); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*ly), &ly); CHKERRQ(ierr);
	ierr = PetscMalloc(p*sizeof(*lz), &lz); CHKERRQ(ierr);
	ierr = PetscMemcpy(lx, plx, m*sizeof(*lx)); CHKERRQ(ierr);
	ierr = PetscMemcpy(ly, ply, n*sizeof(*ly)); CHKERRQ(ierr);
	ierr = PetscMemcpy(lz, plz, p*sizeof(*lz)); CHKERRQ(ierr);
	if (!gridCtx.periodic_z)
	{
	  lz[p-1]--;
	  gridCtx.nz--;
	}
	ierr = DMDACreate3d(PETSC_COMM_WORLD, 
	                    bType_x, bType_y, bType_z, 
	                    DMDA_STENCIL_BOX, 
	                    gridCtx.nx, gridCtx.ny, gridCtx.nz,
	                    m, n, p, 1, 1, lx, ly, lz, 
	                    &uz.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) uz.da, nullptr, "-uz_dmda_view"); CHKERRQ(ierr);
	ierr = PetscFree(lx); CHKERRQ(ierr);
	ierr = PetscFree(ly); CHKERRQ(ierr);
	ierr = PetscFree(lz); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(uz); CHKERRQ(ierr);
	ierr = PetibmGridInitialize(gridCtx, uz.grid); CHKERRQ(ierr);
	if (!gridCtx.periodic_z)
		gridCtx.nz++;
	ierr = PetibmGridHDF5Read(
		directory + "/grids/staggered-z.h5", uz.grid); CHKERRQ(ierr);

	// create z-vorticity field
	PetibmField wz;
	ierr = PetscMalloc(m*sizeof(*lx), &lx); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*ly), &ly); CHKERRQ(ierr);
	ierr = PetscMalloc(p*sizeof(*lz), &lz); CHKERRQ(ierr);
	ierr = PetscMemcpy(lx, plx, m*sizeof(*lx)); CHKERRQ(ierr);
	ierr = PetscMemcpy(ly, ply, n*sizeof(*ly)); CHKERRQ(ierr);
	ierr = PetscMemcpy(lz, plz, p*sizeof(*lz)); CHKERRQ(ierr);
	lx[m-1]--;
	ly[n-1]--;
	gridCtx.nx--;
	gridCtx.ny--;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
	                    bType_x, bType_y, bType_z,
	                    DMDA_STENCIL_STAR,
	                    gridCtx.nx, gridCtx.ny, gridCtx.nz,
	                    m, n, p, 1, 1, lx, ly, lz, 
	                    &wz.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) wz.da, nullptr, "-wz_dmda_view"); CHKERRQ(ierr);
	ierr = PetscFree(lx); CHKERRQ(ierr);
	ierr = PetscFree(ly); CHKERRQ(ierr);
	ierr = PetscFree(lz); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(wz); CHKERRQ(ierr);
	ierr = PetibmGridInitialize(gridCtx, uz.grid); CHKERRQ(ierr);
	gridCtx.nx++;
	gridCtx.ny++;
	ierr = PetibmVorticityZComputeGrid(ux.grid, uy.grid, wz.grid); CHKERRQ(ierr);
	ierr = PetibmGridHDF5Write(
		directory + "/grids/wz.h5", wz.grid); CHKERRQ(ierr);

	// create x-vorticity field
	PetibmField wx;
	ierr = PetscMalloc(m*sizeof(*lx), &lx); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*ly), &ly); CHKERRQ(ierr);
	ierr = PetscMalloc(p*sizeof(*lz), &lz); CHKERRQ(ierr);
	ierr = PetscMemcpy(lx, plx, m*sizeof(*lx)); CHKERRQ(ierr);
	ierr = PetscMemcpy(ly, ply, n*sizeof(*ly)); CHKERRQ(ierr);
	ierr = PetscMemcpy(lz, plz, p*sizeof(*lz)); CHKERRQ(ierr);
	ly[n-1]--;
	lz[p-1]--;
	gridCtx.ny--;
	gridCtx.nz--;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
	                    bType_x, bType_y, bType_z,
	                    DMDA_STENCIL_STAR,
	                    gridCtx.nx, gridCtx.ny, gridCtx.nz,
	                    m, n, p, 1, 1, lx, ly, lz, 
	                    &wz.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) wx.da, nullptr, "-wx_dmda_view"); CHKERRQ(ierr);
	ierr = PetscFree(lx); CHKERRQ(ierr);
	ierr = PetscFree(ly); CHKERRQ(ierr);
	ierr = PetscFree(lz); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(wx); CHKERRQ(ierr);
	ierr = PetibmGridInitialize(gridCtx, ux.grid); CHKERRQ(ierr);
	gridCtx.ny++;
	gridCtx.nz++;
	ierr = PetibmVorticityXComputeGrid(uy.grid, uz.grid, wx.grid); CHKERRQ(ierr);
	ierr = PetibmGridHDF5Write(
		directory + "/grids/wx.h5", wx.grid); CHKERRQ(ierr);

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
	  ierr = PetibmFieldHDF5Read(folder + "/uz.h5", "uz", uz); CHKERRQ(ierr);
	  // compute the z-vorticity
	  ierr = PetibmVorticityZComputeField(ux, uy, wz); CHKERRQ(ierr);
	  ierr = PetibmFieldHDF5Write(folder + "/wz.h5", "wz", wz); CHKERRQ(ierr);
	  // compute the x-vorticity
	  ierr = PetibmVorticityXComputeField(uy, uz, wx); CHKERRQ(ierr);
	  ierr = PetibmFieldHDF5Write(folder + "/wx.h5", "wx", wx); CHKERRQ(ierr);
	}

	ierr = PetibmFieldDestroy(ux); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(uy); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(uz); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(phi); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(wz); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(wx); CHKERRQ(ierr);
	ierr = PetscFinalize(); CHKERRQ(ierr);
	return 0;
} // main
