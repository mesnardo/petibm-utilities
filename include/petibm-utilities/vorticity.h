/* Definition of the functions to compute the 2D vorticity field.
 * \file vorticity.h
 */

#pragma once

#include <petscsys.h>
#include <petscdmda.h>

#include "petibm-utilities/field.h"
#include "petibm-utilities/grid.h"


/*! Computes the gridlines for the vorticity in the z-direction.
 *
 * \param ux [in] The grid for the x-velocity.
 * \param uy [in] The grid for the y-velocity.
 * \param wz [out] The grid for the z-vorticity (passed by reference).
 */
PetscErrorCode PetibmVorticityZComputeGrid(
	const PetibmGrid ux, const PetibmGrid uy, PetibmGrid &wz);


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
	PetibmField ux, PetibmField uy, PetibmField &wz);


/*! Computes the gridlines for the vorticity in the x-direction.
 *
 * \param uy [in] The grid for the y-velocity.
 * \param uz [in] The grid for the z-velocity.
 * \param wx [out] The grid for the x-vorticity (passed by reference).
 */
PetscErrorCode PetibmVorticityXComputeGrid(
	const PetibmGrid uy, const PetibmGrid uz, PetibmGrid &wx);


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
	PetibmField uy, PetibmField uz, PetibmField &wx);
