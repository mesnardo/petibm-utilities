/*! Definition of the structure and functions related to the time-steps.
 * \file timestep.h
 */

#pragma once

#include <petscsys.h>


/*! Structure holding information regarding the time steps.
 */
struct PetibmTimeStepCtx
{
	PetscInt start = 0,  /// starting time step
	         end = 0,  /// ending time step
	         step = 1;  /// time-step increment
}; // PetibmTimeStepCtx


/*! Gets options from command-line or config file.
 *
 * \param prefix String to prepend the name of the options.
 * \param ctx The PetibmTimeStepCtx object to fill (passed by pointer).
 */
PetscErrorCode PetibmTimeStepGetOptions(
	const char prefix[], PetibmTimeStepCtx *ctx);
