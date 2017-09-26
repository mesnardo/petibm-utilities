/* Definition of the functions to compute the 2D vorticity field.
 * \file vorticity.hpp
 */

#pragma once

#include <petscsys.h>
#include <petscdmda.h>

#include "petibm-utilities/field.h"
#include "petibm-utilities/grid.h"


/*! Computes the gridlines for the vorticity in the z-direction.
 *
 * \param ux PetibmGrid structure for the x-velocity.
 * \param uy PetibmGrid structure for the y-velocity.
 * \param wz PetibmGrid structure for the z-vorticity (passed by reference).
 */
PetscErrorCode PetibmVorticityZComputeGrid(
	const PetibmGrid ux, const PetibmGrid uy, PetibmGrid &wz);

/*! Computes the vorticity in the z-direction.
 *
 * First-order.
 *
 * \param ux PetibmField structure for the x-velocity.
 * \param uy PetibmField structure for the y-velocity.
 * \param wz PetibmField structure for the z-vorticity (passed by reference).
 */
PetscErrorCode PetibmVorticityZComputeField(
	const PetibmField ux, const PetibmField uy, PetibmField &wz);

/*! Computes the gridlines for the vorticity in the x-direction.
 *
 * \param uy PetibmGrid structure for the y-velocity.
 * \param uz PetibmGrid structure for the z-velocity.
 * \param wx PetibmGrid structure for the x-vorticity (passed by reference).
 */
PetscErrorCode PetibmVorticityXComputeGrid(
	const PetibmGrid uy, const PetibmGrid uz, PetibmGrid &wx);

/*! Computes the vorticity in the x-direction.
 *
 * First-order.
 *
 * \param uy PetibmField structure for the y-velocity.
 * \param uz PetibmField structure for the z-velocity.
 * \param wx PetibmField structure for the x-vorticity (passed by reference).
 */
PetscErrorCode PetibmVorticityXComputeField(
	const PetibmField uy, const PetibmField uz, PetibmField &wx);
