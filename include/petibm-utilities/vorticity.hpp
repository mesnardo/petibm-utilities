/* Definition of the functions to compute the 2D vorticity field.
 * \file vorticity.hpp
 */

#pragma once

#include <petscsys.h>
#include <petscdmda.h>

#include "petibm-utilities/field.hpp"


/*! Computes the gridlines for the vorticity in the z-direction.
 *
 * \param ux Field structure; velocity in the x-direction.
 * \param uy Field structure; velocity in the y-direction.
 * \param wz Field of the z-vorticity (passed by reference).
 */
PetscErrorCode ComputeGridVorticityZ(Field ux, Field uy, Field &wz);

/*! Computes the vorticity in the z-direction.
 *
 * First-order.
 *
 * \param ux Field structure; velocity in the x-direction.
 * \param uy Field structure; velocity in the y-direction.
 * \param wz Field of the z-vorticity (passed by reference).
 */
PetscErrorCode ComputeVorticityZ(Field ux, Field uy, Field &wz);
