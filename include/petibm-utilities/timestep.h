/*! Definition of the structure and functions related to the time-steps.
 * \file timestep.h
 */

#pragma once

#include <petscsys.h>


typedef struct
{
	PetscInt start = 0,
	         end = 0,
	         step = 1;
} PetibmTimeStepCtx;


/*! Gets options from command-line or config file.
 *
 * \param prefix String to prepend the name of the options.
 * \param timesteps The PetibmTimeSteps object to fill.
 */
PetscErrorCode PetibmTimeStepsGetOptions(
	const char prefix[], PetibmTimeStepCtx *ctx);
