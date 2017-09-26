/*! Implementation of the functions related to the time-steps.
 * \file timestep.cpp
 */

#include "petibm-utilities/timestep.h"


/*! Gets options from command-line or config file.
 *
 * \param prefix String to prepend the name of the options.
 * \param ctx The PetibmTimeStepCtx object to fill (passed by pointer).
 */
PetscErrorCode PetibmTimeStepGetOptions(
	const char prefix[], PetibmTimeStepCtx *ctx)
{
	PetscErrorCode ierr;
	char path[PETSC_MAX_PATH_LEN];
	PetscBool found;

	PetscFunctionBeginUser;

	// get path of the configuration file
	ierr = PetscOptionsGetString(nullptr, prefix, "-config_file",
	                             path, sizeof(path), &found); CHKERRQ(ierr);
	if (found)
	{
		ierr = PetscOptionsInsertFile(
			PETSC_COMM_WORLD, nullptr, path, PETSC_FALSE); CHKERRQ(ierr);
	}
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
