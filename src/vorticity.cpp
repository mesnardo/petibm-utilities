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
	PetscReal *x, *y, *z;
	PetscReal *x_m, *y_m, *z_m;
	PetscInt nx, ny, nz;
	DMDALocalInfo info;

	PetscFunctionBeginUser;

	ierr = DMDAGetLocalInfo(wz.da, &info); CHKERRQ(ierr);

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

	if (info.dim == 3)
	{
		ierr = VecGetSize(wz.z, &nz); CHKERRQ(ierr);
		ierr = VecGetArray(wz.z, &z); CHKERRQ(ierr);
		ierr = VecGetArray(ux.z, &z_m); CHKERRQ(ierr);
		for (PetscInt k=0; k<nz; k++)
			z[k] = z_m[k];
		ierr = VecRestoreArray(wz.z, &z); CHKERRQ(ierr);
		ierr = VecRestoreArray(ux.z, &z_m); CHKERRQ(ierr);
	}

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
	PetscInt i, j, k;
	PetscReal *x_a, *y_a;
	PetscReal dx, dy, du, dv;

	PetscFunctionBeginUser;

	// create local velocity vectors from global ones
	ierr = DMGlobalToLocalBegin(ux.da, ux.global, INSERT_VALUES, ux.local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ux.da, ux.global, INSERT_VALUES, ux.local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(uy.da, uy.global, INSERT_VALUES, uy.local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(uy.da, uy.global, INSERT_VALUES, uy.local); CHKERRQ(ierr);

	ierr = DMDAGetLocalInfo(wz.da, &info); CHKERRQ(ierr);

	if (info.dim == 2)
	{
		PetscReal **wz_a, **ux_a, **uy_a;
		// get access to vector values in multi-dimensional fashion
		ierr = DMDAVecGetArray(wz.da, wz.global, &wz_a); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(ux.da, ux.local, &ux_a); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(uy.da, uy.local, &uy_a); CHKERRQ(ierr);
		ierr = VecGetArray(uy.x, &x_a); CHKERRQ(ierr);
		ierr = VecGetArray(ux.y, &y_a); CHKERRQ(ierr);

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
	}
	else if (info.dim == 3)
	{
		PetscReal ***wz_a, ***ux_a, ***uy_a;
		// get access to vector values in multi-dimensional fashion
		ierr = DMDAVecGetArray(wz.da, wz.global, &wz_a); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(ux.da, ux.local, &ux_a); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(uy.da, uy.local, &uy_a); CHKERRQ(ierr);
		ierr = VecGetArray(uy.x, &x_a); CHKERRQ(ierr);
		ierr = VecGetArray(ux.y, &y_a); CHKERRQ(ierr);

		for (k=info.zs; k<info.zs+info.zm; k++)
		{
			for (j=info.ys; j<info.ys+info.ym; j++)
			{
				for (i=info.xs; i<info.xs+info.xm; i++)
				{
					dx = x_a[i] - x_a[i-1];
					dy = y_a[j] - y_a[j-1];
					du = ux_a[k][j][i] - ux_a[k][j-1][i];
					dv = uy_a[k][j][i] - uy_a[k][j][i-1];
					wz_a[k][j][i] = dv/dx - du/dy;
				}
			}
		}
		ierr = VecRestoreArray(uy.x, &x_a); CHKERRQ(ierr);
		ierr = VecRestoreArray(ux.y, &y_a); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(ux.da, ux.local, &ux_a); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(uy.da, uy.local, &uy_a); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(wz.da, wz.global, &wz_a); CHKERRQ(ierr);
	}
	else
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,
		        "Function only supports 2D or 3D fields");

	PetscFunctionReturn(0);
} // FieldComputeVorticityZ


/*! Computes the gridlines for the vorticity in the x-direction.
 *
 * \param uy Field structure; velocity in the y-direction.
 * \param uz Field structure; velocity in the z-direction.
 * \param wx Field of the x-vorticity (passed by reference).
 */
PetscErrorCode ComputeGridVorticityX(Field uy, Field uz, Field &wx)
{
	PetscErrorCode ierr;
	PetscReal *x, *y, *z;
	PetscReal *x_m, *y_m, *z_m;
	PetscInt nx, ny, nz;
	DMDALocalInfo info;

	PetscFunctionBeginUser;

	ierr = DMDAGetLocalInfo(wx.da, &info); CHKERRQ(ierr);
	if (info.dim != 3)
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,
		        "Function only supports 3D fields");

	ierr = VecGetSize(wx.x, &nx); CHKERRQ(ierr);
	ierr = VecGetArray(wx.x, &x); CHKERRQ(ierr);
	ierr = VecGetArray(uy.x, &x_m); CHKERRQ(ierr);
	for (PetscInt i=0; i<nx; i++)
		x[i] = x_m[i];
	ierr = VecRestoreArray(uy.x, &x_m); CHKERRQ(ierr);
	ierr = VecRestoreArray(wx.x, &x); CHKERRQ(ierr);

	ierr = VecGetSize(wx.y, &ny); CHKERRQ(ierr);
	ierr = VecGetArray(wx.y, &y); CHKERRQ(ierr);
	ierr = VecGetArray(uy.y, &y_m); CHKERRQ(ierr);
	for (PetscInt j=0; j<ny; j++)
		y[j] = y_m[j];
	ierr = VecRestoreArray(wx.y, &y); CHKERRQ(ierr);
	ierr = VecRestoreArray(uy.y, &y_m); CHKERRQ(ierr);

	ierr = VecGetSize(wx.z, &nz); CHKERRQ(ierr);
	ierr = VecGetArray(wx.z, &z); CHKERRQ(ierr);
	ierr = VecGetArray(uz.z, &z_m); CHKERRQ(ierr);
	for (PetscInt k=0; k<nz; k++)
		z[k] = z_m[k];
	ierr = VecRestoreArray(wx.z, &z); CHKERRQ(ierr);
	ierr = VecRestoreArray(uz.z, &z_m); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // ComputeGridVorticityX


/*! Computes the vorticity in the x-direction.
 *
 * First-order.
 *
 * \param uy Field structure; velocity in the y-direction.
 * \param uz Field structure; velocity in the z-direction.
 * \param wx Field of the x-vorticity (passed by reference).
 */
PetscErrorCode ComputeVorticityX(Field uy, Field uz, Field &wx)
{
	PetscErrorCode ierr;
	DMDALocalInfo info;
	PetscInt i, j, k;
	PetscReal ***wx_a, ***uy_a, ***uz_a;
	PetscReal *y_a, *z_a;
	PetscReal dy, dz, dv, dw;

	PetscFunctionBeginUser;

	ierr = DMDAGetLocalInfo(wx.da, &info); CHKERRQ(ierr);
	if (info.dim != 3)
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,
		        "Function only supports 3D fields");

	// create local velocity vectors from global ones
	ierr = DMGlobalToLocalBegin(uy.da, uy.global, INSERT_VALUES, uy.local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(uy.da, uy.global, INSERT_VALUES, uy.local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(uz.da, uz.global, INSERT_VALUES, uz.local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(uz.da, uz.global, INSERT_VALUES, uz.local); CHKERRQ(ierr);

	// get access to vector values in multi-dimensional fashion
	ierr = DMDAVecGetArray(wx.da, wx.global, &wx_a); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(uy.da, uy.local, &uy_a); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(uz.da, uz.local, &uz_a); CHKERRQ(ierr);
	ierr = VecGetArray(uz.y, &y_a); CHKERRQ(ierr);
	ierr = VecGetArray(uy.z, &z_a); CHKERRQ(ierr);

	for (k=info.zs; k<info.zs+info.zm; k++)
	{
		for (j=info.ys; j<info.ys+info.ym; j++)
		{
			for (i=info.xs; i<info.xs+info.xm; i++)
			{
				dy = y_a[j] - y_a[j-1];
				dz = z_a[k] - z_a[k-1];
				dv = uy_a[k][j][i] - uy_a[k][j][i-1];
				dw = uz_a[k][j][i] - uz_a[k-1][j][i];
				wx_a[k][j][i] = dw/dy - dv/dz;
			}
		}
	}

	ierr = VecRestoreArray(uz.y, &y_a); CHKERRQ(ierr);
	ierr = VecRestoreArray(uy.z, &z_a); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(uy.da, uy.local, &uy_a); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(uz.da, uz.local, &uz_a); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(wx.da, wx.global, &wx_a); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // FieldComputeVorticityX
