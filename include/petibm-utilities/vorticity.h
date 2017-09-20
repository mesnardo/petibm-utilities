/* Definition of the functions to compute the 2D vorticity field.
 * \file vorticity.hpp
 */

#pragma once

#include <petscsys.h>
#include <petscdmda.h>

#include "petibm-utilities/field.h"


/*! Computes the gridlines for the vorticity in the z-direction.
 *
 * \param ux PetibmField structure; velocity in the x-direction.
 * \param uy PetibmField structure; velocity in the y-direction.
 * \param wz PetibmField of the z-vorticity (passed by reference).
 */
PetscErrorCode PetibmComputeGridVorticityZ(
	PetibmField ux, PetibmField uy, PetibmField &wz);

/*! Computes the vorticity in the z-direction.
 *
 * First-order.
 *
 * \param ux PetibmField structure; velocity in the x-direction.
 * \param uy PetibmField structure; velocity in the y-direction.
 * \param wz PetibmField of the z-vorticity (passed by reference).
 */
PetscErrorCode PetibmComputeFieldVorticityZ(
	PetibmField ux, PetibmField uy, PetibmField &wz);

/*! Computes the gridlines for the vorticity in the x-direction.
 *
 * \param uy PetibmField structure; velocity in the y-direction.
 * \param uz PetibmField structure; velocity in the z-direction.
 * \param wx PetibmField of the x-vorticity (passed by reference).
 */
PetscErrorCode PetibmComputeGridVorticityX(
	PetibmField uy, PetibmField uz, PetibmField &wx);

/*! Computes the vorticity in the x-direction.
 *
 * First-order.
 *
 * \param uy PetibmField structure; velocity in the y-direction.
 * \param uz PetibmField structure; velocity in the z-direction.
 * \param wx PetibmField of the x-vorticity (passed by reference).
 */
PetscErrorCode PetibmComputeFieldVorticityX(
	PetibmField uy, PetibmField uz, PetibmField &wx);
