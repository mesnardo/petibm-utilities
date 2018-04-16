/*! Implementation of the functions related to the time-steps.
 * \file timestep.cpp
 */

#include "petibm-utilities/timestep.h"
#include "petibm-utilities/misc.h"


/*! Gets options from command-line or configuration file.
 *
 * \param prefix [in] String to prepend the name of the options.
 * \param ctx [out] The PetibmTimeStepCtx object to fill (passed by pointer).
 */
PetscErrorCode PetibmTimeStepGetOptions(
	const char prefix[], PetibmTimeStepCtx *ctx)
{
	PetscErrorCode ierr;
	PetscBool found;

	PetscFunctionBeginUser;

	// get starting time step
	ierr = PetscOptionsGetInt(
		nullptr, prefix, "-nstart", &ctx->start, &found); CHKERRQ(ierr);
	// get ending time step
	ierr = PetscOptionsGetInt(
		nullptr, prefix, "-nend", &ctx->end, &found); CHKERRQ(ierr);
	// get time-step increment
	ierr = PetscOptionsGetInt(
		nullptr, prefix, "-nstep", &ctx->step, &found); CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
} // PetibmTimeStepGetOptions
