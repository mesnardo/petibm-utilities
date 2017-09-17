/*! Implementation of the functions to compute the 2D vorticity field.
 * \file vorticity.cpp
 */

#include "petibm-utilities/vorticity.hpp"


/*! Computes the gridlines for the vorticity in the z-direction.
 *
 * \param ux Field structure; velocity in the x-direction.
 * \param uy Field structure; velocity in the y-direction.
 * \param wz Field of the z-vorticity (passed by reference).
 */
PetscErrorCode ComputeGridVorticityZ(Field ux, Field uy, Field &wz)
{
  PetscErrorCode ierr;
  PetscReal *x, *y;
  PetscReal *x_m, *y_m;
  PetscInt nx, ny;

  PetscFunctionBeginUser;

  ierr = VecGetSize(wz.x, &nx); CHKERRQ(ierr);
  ierr = VecGetArray(wz.x, &x); CHKERRQ(ierr);
  ierr = VecGetArray(ux.x, &x_m); CHKERRQ(ierr);
  for (PetscInt i=0; i<nx; i++)
    x[i] = x_m[i];
  ierr = VecRestoreArray(ux.x, &x_m); CHKERRQ(ierr);
  ierr = VecRestoreArray(wz.x, &x); CHKERRQ(ierr);
  
  ierr = VecGetSize(wz.y, &ny); CHKERRQ(ierr);
  ierr = VecGetArray(wz.y, &y); CHKERRQ(ierr);
  ierr = VecGetArray(uy.y, &y_m); CHKERRQ(ierr);
  for (PetscInt j=0; j<ny; j++)
    y[j] = y_m[j];
  ierr = VecRestoreArray(wz.y, &y); CHKERRQ(ierr);
  ierr = VecRestoreArray(uy.y, &y_m); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // ComputeGridVorticityZ


/*! Computes the vorticity in the z-direction.
 *
 * First-order.
 *
 * \param ux Field structure; velocity in the x-direction.
 * \param uy Field structure; velocity in the y-direction.
 * \param wz Field of the z-vorticity (passed by reference).
 */
PetscErrorCode ComputeVorticityZ(Field ux, Field uy, Field &wz)
{
  PetscErrorCode ierr;
  DMDALocalInfo info;
  PetscInt i, j;
  PetscReal **wz_a, **ux_a, **uy_a;
  PetscReal *x_a, *y_a;
  PetscReal dx, dy, du, dv;

  PetscFunctionBeginUser;

  // create local velocity vectors from global ones
  ierr = DMGlobalToLocalBegin(ux.da, ux.global, INSERT_VALUES, ux.local); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ux.da, ux.global, INSERT_VALUES, ux.local); CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(uy.da, uy.global, INSERT_VALUES, uy.local); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(uy.da, uy.global, INSERT_VALUES, uy.local); CHKERRQ(ierr);

  // get access to vector values in multi-dimensional fashion
  ierr = DMDAVecGetArray(wz.da, wz.global, &wz_a); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ux.da, ux.local, &ux_a); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(uy.da, uy.local, &uy_a); CHKERRQ(ierr);
  ierr = VecGetArray(uy.x, &x_a); CHKERRQ(ierr);
  ierr = VecGetArray(ux.y, &y_a); CHKERRQ(ierr);
  
  ierr = DMDAGetLocalInfo(wz.da, &info); CHKERRQ(ierr);
  for (j=info.ys; j<info.ys+info.ym; j++)
  {
    for (i=info.xs; i<info.xs+info.xm; i++)
    {
      dx = x_a[i] - x_a[i-1];
      dy = y_a[j] - y_a[j-1];
      du = ux_a[j][i] - ux_a[j-1][i];
      dv = uy_a[j][i] - uy_a[j][i-1];
      wz_a[j][i] = dv/dx - du/dy;
    }
  }

  ierr = VecRestoreArray(uy.x, &x_a); CHKERRQ(ierr);
  ierr = VecRestoreArray(ux.y, &y_a); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ux.da, ux.local, &ux_a); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(uy.da, uy.local, &uy_a); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(wz.da, wz.global, &wz_a); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // FieldComputeVorticityZ
