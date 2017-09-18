/*! Computes the vorticity field from the 3D velocity vector field.
 * \file vorticity3d.cpp
 */

#include <iomanip>
#include <iostream>

#include <petscsys.h>
#include <petscdmda.h>

#include "petibm-utilities/misc.hpp"
#include "petibm-utilities/field.hpp"
#include "petibm-utilities/vorticity.hpp"


int main(int argc, char **argv)
{
	PetscErrorCode ierr;

	ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);

	// parse command-line options
	std::string directory;
	ierr = PetibmGetDirectory(&directory); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,
	                   "[INFO] directory: %s\n",
	                   directory.c_str()); CHKERRQ(ierr);

	PetscInt nstart = 0,
	         nend = 0,
	         nstep = 1;
	PetscBool found;
	ierr = PetscOptionsGetInt(
		nullptr, nullptr, "-nstart", &nstart, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(
		nullptr, nullptr, "-nend", &nend, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(
		nullptr, nullptr, "-nstep", &nstep, &found); CHKERRQ(ierr);

	PetscInt nx = 0,
	         ny = 0,
	         nz = 0;
	ierr = PetscOptionsGetInt(NULL, NULL, "-nx", &nx, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-ny", &ny, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-nz", &nz, &found); CHKERRQ(ierr);  

	PetscBool isPeriodic_x = PETSC_FALSE,
	          isPeriodic_y = PETSC_FALSE,
	          isPeriodic_z = PETSC_FALSE;
	ierr = PetscOptionsGetBool(NULL, NULL, "-periodic_x", &isPeriodic_x, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL, NULL, "-periodic_y", &isPeriodic_y, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL, NULL, "-periodic_z", &isPeriodic_z, NULL); CHKERRQ(ierr);

	DMBoundaryType bType_x = (isPeriodic_x) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED,
	               bType_y = (isPeriodic_y) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED,
	               bType_z = (isPeriodic_z) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;

	// create DMDA for phi
	PetibmField phi;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
	                    bType_x, bType_y, bType_z,
	                    DMDA_STENCIL_STAR,
	                    nx, ny, nz,
	                    PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, NULL,
	                    &phi.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions((PetscObject) phi.da, NULL, "-phi_dmda_view"); CHKERRQ(ierr);

	// create DMDA objects for velocity components from DMDA object for pressure
	PetibmField ux, uy, uz;
	PetscInt numX, numY, numZ;
	const PetscInt *plx, *ply, *plz;
	ierr = DMDAGetOwnershipRanges(phi.da, &plx, &ply, &plz); CHKERRQ(ierr);
	PetscInt m, n, p;
	ierr = DMDAGetInfo(phi.da, NULL, NULL, NULL, NULL, &m, &n, &p, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
	// x-component of velocity
	PetscInt *ulx, *uly, *ulz;
	ierr = PetscMalloc(m*sizeof(*ulx), &ulx); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*uly), &uly); CHKERRQ(ierr);
	ierr = PetscMalloc(p*sizeof(*ulz), &ulz); CHKERRQ(ierr);
	ierr = PetscMemcpy(ulx, plx, m*sizeof(*ulx)); CHKERRQ(ierr);
	ierr = PetscMemcpy(uly, ply, n*sizeof(*uly)); CHKERRQ(ierr);
	ierr = PetscMemcpy(ulz, plz, p*sizeof(*ulz)); CHKERRQ(ierr);
	numX = nx;
	numY = ny;
	numZ = nz;
	if (!isPeriodic_x)
	{
	  ulx[m-1]--;
	  numX--;
	}
	ierr = DMDACreate3d(PETSC_COMM_WORLD, 
	                    bType_x, bType_y, bType_z, 
	                    DMDA_STENCIL_BOX, 
	                    numX, numY, numZ, m, n, p, 1, 1, ulx, uly, ulz, 
	                    &ux.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions((PetscObject) ux.da, NULL, "-ux_dmda_view"); CHKERRQ(ierr);
	ierr = PetscFree(ulx); CHKERRQ(ierr);
	ierr = PetscFree(uly); CHKERRQ(ierr);
	ierr = PetscFree(ulz); CHKERRQ(ierr);
	// y-component of velocity
	PetscInt *vlx, *vly, *vlz;
	ierr = PetscMalloc(m*sizeof(*vlx), &vlx); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*vly), &vly); CHKERRQ(ierr);
	ierr = PetscMalloc(p*sizeof(*vlz), &vlz); CHKERRQ(ierr);
	ierr = PetscMemcpy(vlx, plx, m*sizeof(*vlx)); CHKERRQ(ierr);
	ierr = PetscMemcpy(vly, ply, n*sizeof(*vly)); CHKERRQ(ierr);
	ierr = PetscMemcpy(vlz, plz, p*sizeof(*vlz)); CHKERRQ(ierr);
	numX = nx;
	numY = ny;
	numZ = nz;
	if (!isPeriodic_y)
	{
	  vly[n-1]--;
	  numY--;
	}
	ierr = DMDACreate3d(PETSC_COMM_WORLD, 
	                    bType_x, bType_y, bType_z, 
	                    DMDA_STENCIL_BOX, 
	                    numX, numY, numZ, m, n, p, 1, 1, vlx, vly, vlz, 
	                    &uy.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions((PetscObject) uy.da, NULL, "-uy_dmda_view"); CHKERRQ(ierr);
	ierr = PetscFree(vlx); CHKERRQ(ierr);
	ierr = PetscFree(vly); CHKERRQ(ierr);
	ierr = PetscFree(vlz); CHKERRQ(ierr);
	// z-component of velocity
	PetscInt *wlx, *wly, *wlz;
	ierr = PetscMalloc(m*sizeof(*wlx), &wlx); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*wly), &wly); CHKERRQ(ierr);
	ierr = PetscMalloc(p*sizeof(*wlz), &wlz); CHKERRQ(ierr);
	ierr = PetscMemcpy(wlx, plx, m*sizeof(*wlx)); CHKERRQ(ierr);
	ierr = PetscMemcpy(wly, ply, n*sizeof(*wly)); CHKERRQ(ierr);
	ierr = PetscMemcpy(wlz, plz, p*sizeof(*wlz)); CHKERRQ(ierr);
	numX = nx;
	numY = ny;
	numZ = nz;
	if (!isPeriodic_z)
	{
	  wlz[p-1]--;
	  numZ--;
	}
	ierr = DMDACreate3d(PETSC_COMM_WORLD, 
	                    bType_x, bType_y, bType_z, 
	                    DMDA_STENCIL_BOX, 
	                    numX, numY, numZ, m, n, p, 1, 1, wlx, wly, wlz, 
	                    &uz.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions((PetscObject) uz.da, NULL, "-uz_dmda_view"); CHKERRQ(ierr);
	ierr = PetscFree(wlx); CHKERRQ(ierr);
	ierr = PetscFree(wly); CHKERRQ(ierr);
	ierr = PetscFree(wlz); CHKERRQ(ierr);

	ierr = PetibmFieldInitialize(ux); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(uy); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(uz); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(phi); CHKERRQ(ierr);

	// read grids
	std::string gridsDirectory(directory + "/grids");
	ierr = PetibmFieldReadGrid(gridsDirectory + "/staggered-x.h5", ux); CHKERRQ(ierr);
	ierr = PetibmFieldReadGrid(gridsDirectory + "/staggered-y.h5", uy); CHKERRQ(ierr);
	ierr = PetibmFieldReadGrid(gridsDirectory + "/staggered-z.h5", uz); CHKERRQ(ierr);
	ierr = PetibmFieldReadGrid(gridsDirectory + "/cell-centered.h5", phi); CHKERRQ(ierr);

	// create z-vorticity field
	PetibmField wz;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
	                    bType_x, bType_y, bType_z,
	                    DMDA_STENCIL_STAR,
	                    nx-1, ny-1, nz,
	                    PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, NULL,
	                    &wz.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions((PetscObject) wz.da, NULL, "-wz_dmda_view"); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(wz); CHKERRQ(ierr);
	ierr = PetibmComputeGridVorticityZ(ux, uy, wz); CHKERRQ(ierr);
	// create x-vorticity field
	PetibmField wx;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
	                    bType_x, bType_y, bType_z,
	                    DMDA_STENCIL_STAR,
	                    nx, ny-1, nz-1,
	                    PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, NULL,
	                    &wx.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions((PetscObject) wx.da, NULL, "-wx_dmda_view"); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(wx); CHKERRQ(ierr);
	ierr = PetibmComputeGridVorticityX(uy, uz, wx); CHKERRQ(ierr);

	PetscMPIInt rank;
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
	if (rank == 0)
	{
	  ierr = PetibmFieldWriteGrid(gridsDirectory + "/wz.h5", wz); CHKERRQ(ierr);
	  ierr = PetibmFieldWriteGrid(gridsDirectory + "/wx.h5", wx); CHKERRQ(ierr);
	}

	for (PetscInt ite=nstart; ite<=nend; ite+=nstep)
	{
	  ierr = PetscPrintf(PETSC_COMM_WORLD, "[time-step %D]\n", ite); CHKERRQ(ierr);

	  std::stringstream ss;
	  ss << directory << "/" << std::setfill('0') << std::setw(7) << ite;
	  std::string folder(ss.str());

	  // read values
	  ierr = PetibmFieldReadValues(folder + "/ux.h5", "ux", ux); CHKERRQ(ierr);
	  ierr = PetibmFieldReadValues(folder + "/uy.h5", "uy", uy); CHKERRQ(ierr);
	  ierr = PetibmFieldReadValues(folder + "/uz.h5", "uz", uz); CHKERRQ(ierr);

	  // compute the z-vorticity
	  ierr = PetibmComputeFieldVorticityZ(ux, uy, wz); CHKERRQ(ierr);
	  ierr = PetibmFieldWriteValues(folder + "/wz.h5", "wz", wz); CHKERRQ(ierr);
	  // compute the x-vorticity
	  ierr = PetibmComputeFieldVorticityX(uy, uz, wx); CHKERRQ(ierr);
	  ierr = PetibmFieldWriteValues(folder + "/wx.h5", "wx", wx); CHKERRQ(ierr);
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
