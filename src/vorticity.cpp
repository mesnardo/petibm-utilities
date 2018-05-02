/*! Implementation of the functions to compute the vorticity field.
 * \file vorticity.cpp
 */

#include "petibm-utilities/vorticity.h"


/*! Computes the gridlines for the vorticity in the z-direction.
 *
 * \param ux [in] The grid for the x-velocity.
 * \param uy [in] The grid for the y-velocity.
 * \param wz [out] The grid for the z-vorticity (passed by reference).
 */
PetscErrorCode PetibmVorticityZComputeGrid(
	const PetibmGrid ux, const PetibmGrid uy, PetibmGrid &wz)
{
	PetscErrorCode ierr;
	PetscReal *x, *y, *z;
	PetscReal *x_m, *y_m, *z_m;
	PetscInt nx, ny, nz;
	PetscInt i, j, k;

	PetscFunctionBeginUser;

	ierr = VecGetSize(wz.x.coords, &nx); CHKERRQ(ierr);
	ierr = VecGetArray(wz.x.coords, &x); CHKERRQ(ierr);
	ierr = VecGetArray(ux.x.coords, &x_m); CHKERRQ(ierr);
	for (i=0; i<nx; i++)
		x[i] = x_m[i];
	ierr = VecRestoreArray(ux.x.coords, &x_m); CHKERRQ(ierr);
	ierr = VecRestoreArray(wz.x.coords, &x); CHKERRQ(ierr);

	ierr = VecGetSize(wz.y.coords, &ny); CHKERRQ(ierr);
	ierr = VecGetArray(wz.y.coords, &y); CHKERRQ(ierr);
	ierr = VecGetArray(uy.y.coords, &y_m); CHKERRQ(ierr);
	for (j=0; j<ny; j++)
		y[j] = y_m[j];
	ierr = VecRestoreArray(wz.y.coords, &y); CHKERRQ(ierr);
	ierr = VecRestoreArray(uy.y.coords, &y_m); CHKERRQ(ierr);

	if (wz.dim == 3)
	{
		ierr = VecGetSize(wz.z.coords, &nz); CHKERRQ(ierr);
		ierr = VecGetArray(wz.z.coords, &z); CHKERRQ(ierr);
		ierr = VecGetArray(ux.z.coords, &z_m); CHKERRQ(ierr);
		for (k=0; k<nz; k++)
			z[k] = z_m[k];
		ierr = VecRestoreArray(wz.z.coords, &z); CHKERRQ(ierr);
		ierr = VecRestoreArray(ux.z.coords, &z_m); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
} // PetibmVorticityZComputeGrid


/*! Computes the vorticity in the z-direction.
 *
 * First-order.
 *
 * \param gridux [in] The grid for the velocity in the x-direction.
 * \param griduy [in] The grid for the velocity in the y-direction.
 * \param ux [in] The velocity field in the x-direction.
 * \param uy [in] The velocity field in the y-direction.
 * \param wz [out] The vorticity field in the z-direction (passed by reference).
 */
PetscErrorCode PetibmVorticityZComputeField(
	const PetibmGrid gridux, const PetibmGrid griduy,
	PetibmField ux, PetibmField uy, PetibmField &wz)
{
	PetscErrorCode ierr;
	DMDALocalInfo info;
	PetscInt i, j, k;
	PetscReal *x_a, *y_a;
	PetscReal dx, dy, du, dv;
	Vec ux_local, uy_local;

	PetscFunctionBeginUser;

	ierr = DMDAGetLocalInfo(wz.da, &info); CHKERRQ(ierr);

	if (info.dim == 2)
	{
		PetscReal **wz_a, **ux_a, **uy_a;
		ierr = DMDAVecGetArray(wz.da, wz.global, &wz_a); CHKERRQ(ierr);
		ierr = DMGetLocalVector(ux.da, &ux_local); CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin(
			ux.da, ux.global, INSERT_VALUES, ux_local); CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd(
			ux.da, ux.global, INSERT_VALUES, ux_local); CHKERRQ(ierr);
		ierr = DMGetLocalVector(uy.da, &uy_local); CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin(
			uy.da, uy.global, INSERT_VALUES, uy_local); CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd(
			uy.da, uy.global, INSERT_VALUES, uy_local); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(ux.da, ux_local, &ux_a); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(uy.da, uy_local, &uy_a); CHKERRQ(ierr);
		ierr = VecGetArray(griduy.x.coords, &x_a); CHKERRQ(ierr);
		ierr = VecGetArray(gridux.y.coords, &y_a); CHKERRQ(ierr);
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
		ierr = VecRestoreArray(griduy.x.coords, &x_a); CHKERRQ(ierr);
		ierr = VecRestoreArray(gridux.y.coords, &y_a); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(ux.da, ux_local, &ux_a); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(uy.da, uy_local, &uy_a); CHKERRQ(ierr);
		ierr = DMRestoreLocalVector(ux.da, &ux_local); CHKERRQ(ierr);
		ierr = DMRestoreLocalVector(uy.da, &uy_local); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(wz.da, wz.global, &wz_a); CHKERRQ(ierr);
	}
	else if (info.dim == 3)
	{
		PetscReal ***wz_a, ***ux_a, ***uy_a;
		ierr = DMDAVecGetArray(wz.da, wz.global, &wz_a); CHKERRQ(ierr);
		ierr = DMGetLocalVector(ux.da, &ux_local); CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin(
			ux.da, ux.global, INSERT_VALUES, ux_local); CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd(
			ux.da, ux.global, INSERT_VALUES, ux_local); CHKERRQ(ierr);
		ierr = DMGetLocalVector(uy.da, &uy_local); CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin(
			uy.da, uy.global, INSERT_VALUES, uy_local); CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd(
			uy.da, uy.global, INSERT_VALUES, uy_local); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(ux.da, ux_local, &ux_a); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(uy.da, uy_local, &uy_a); CHKERRQ(ierr);
		ierr = VecGetArray(griduy.x.coords, &x_a); CHKERRQ(ierr);
		ierr = VecGetArray(gridux.y.coords, &y_a); CHKERRQ(ierr);
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
		ierr = VecRestoreArray(griduy.x.coords, &x_a); CHKERRQ(ierr);
		ierr = VecRestoreArray(gridux.y.coords, &y_a); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(ux.da, ux_local, &ux_a); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(uy.da, uy_local, &uy_a); CHKERRQ(ierr);
		ierr = DMRestoreLocalVector(ux.da, &ux_local); CHKERRQ(ierr);
		ierr = DMRestoreLocalVector(uy.da, &uy_local); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(wz.da, wz.global, &wz_a); CHKERRQ(ierr);
	}
	else
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,
		        "Function only supports 2D or 3D fields");

	PetscFunctionReturn(0);
} // PetibmVorticityZComputeField


/*! Computes the gridlines for the vorticity in the x-direction.
 *
 * \param uy [in] The grid for the y-velocity.
 * \param uz [in] The grid for the z-velocity.
 * \param wx [out] The grid for the x-vorticity (passed by reference).
 */
PetscErrorCode PetibmVorticityXComputeGrid(
	const PetibmGrid uy, const PetibmGrid uz, PetibmGrid &wx)
{
	PetscErrorCode ierr;
	PetscReal *x, *y, *z;
	PetscReal *x_m, *y_m, *z_m;
	PetscInt i, j, k;
	PetscInt nx, ny, nz;

	PetscFunctionBeginUser;

	if (wx.dim != 3)
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,
		        "Function only supports 3D grids");

	ierr = VecGetSize(wx.x.coords, &nx); CHKERRQ(ierr);
	ierr = VecGetArray(wx.x.coords, &x); CHKERRQ(ierr);
	ierr = VecGetArray(uy.x.coords, &x_m); CHKERRQ(ierr);
	for (i=0; i<nx; i++)
		x[i] = x_m[i];
	ierr = VecRestoreArray(uy.x.coords, &x_m); CHKERRQ(ierr);
	ierr = VecRestoreArray(wx.x.coords, &x); CHKERRQ(ierr);

	ierr = VecGetSize(wx.y.coords, &ny); CHKERRQ(ierr);
	ierr = VecGetArray(wx.y.coords, &y); CHKERRQ(ierr);
	ierr = VecGetArray(uy.y.coords, &y_m); CHKERRQ(ierr);
	for (j=0; j<ny; j++)
		y[j] = y_m[j];
	ierr = VecRestoreArray(wx.y.coords, &y); CHKERRQ(ierr);
	ierr = VecRestoreArray(uy.y.coords, &y_m); CHKERRQ(ierr);

	ierr = VecGetSize(wx.z.coords, &nz); CHKERRQ(ierr);
	ierr = VecGetArray(wx.z.coords, &z); CHKERRQ(ierr);
	ierr = VecGetArray(uz.z.coords, &z_m); CHKERRQ(ierr);
	for (k=0; k<nz; k++)
		z[k] = z_m[k];
	ierr = VecRestoreArray(wx.z.coords, &z); CHKERRQ(ierr);
	ierr = VecRestoreArray(uz.z.coords, &z_m); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmVorticityXComputeGrid


/*! Computes the vorticity in the x-direction.
 *
 * First-order.
 *
 * \param griduy [in] The grid for the velocity in the y-direction.
 * \param griduz [in] The grid for the velocity in the z-direction.
 * \param uy [in] The velocity field in the y-direction.
 * \param uz [in] The velocity field in the z-direction.
 * \param wx [out] The vorticity field in the x-direction (passed by reference).
 */
PetscErrorCode PetibmVorticityXComputeField(
	const PetibmGrid griduy, const PetibmGrid griduz,
	PetibmField uy, PetibmField uz, PetibmField &wx)
{
	PetscErrorCode ierr;
	DMDALocalInfo info;
	PetscInt i, j, k;
	PetscReal ***wx_a, ***uy_a, ***uz_a;
	PetscReal *y_a, *z_a;
	PetscReal dy, dz, dv, dw;
	Vec uy_local, uz_local;

	PetscFunctionBeginUser;

	ierr = DMDAGetLocalInfo(wx.da, &info); CHKERRQ(ierr);
	if (info.dim != 3)
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,
		        "Function only supports 3D fields");

	ierr = DMDAVecGetArray(wx.da, wx.global, &wx_a); CHKERRQ(ierr);
	ierr = DMGetLocalVector(uy.da, &uy_local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(
		uy.da, uy.global, INSERT_VALUES, uy_local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(
		uy.da, uy.global, INSERT_VALUES, uy_local); CHKERRQ(ierr);
	ierr = DMGetLocalVector(uz.da, &uz_local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(
		uz.da, uz.global, INSERT_VALUES, uz_local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(
		uz.da, uz.global, INSERT_VALUES, uz_local); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(uy.da, uy_local, &uy_a); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(uz.da, uz_local, &uz_a); CHKERRQ(ierr);
	ierr = VecGetArray(griduz.y.coords, &y_a); CHKERRQ(ierr);
	ierr = VecGetArray(griduy.z.coords, &z_a); CHKERRQ(ierr);
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
	ierr = VecRestoreArray(griduz.y.coords, &y_a); CHKERRQ(ierr);
	ierr = VecRestoreArray(griduy.z.coords, &z_a); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(uy.da, uy_local, &uy_a); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(uz.da, uz_local, &uz_a); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(uy.da, &uy_local); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(uz.da, &uz_local); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(wx.da, wx.global, &wx_a); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmVorticityXComputeField
