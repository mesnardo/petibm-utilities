/*! Implementation of the functions related to the time-steps.
 * \file timestep.cpp
 */

#include "petibm-utilities/timestep.h"


/*! Gets options from command-line or config file.
 *
 * \param prefix String to prepend the name of the options.
 * \param ctx The PetibmTimeStepCtx object to fill.
 */
PetscErrorCode PetibmTimeStepsGetOptions(
	const char prefix[], PetibmTimeStepCtx *ctx)
{
	PetscErrorCode ierr;
	char path[PETSC_MAX_PATH_LEN];
	PetscBool found;

	PetscFunctionBeginUser;

	ierr = PetscOptionsGetString(nullptr, prefix, "-config_file",
	                             path, sizeof(path), &found); CHKERRQ(ierr);
	if (found)
	{
		ierr = PetscOptionsInsertFile(
			PETSC_COMM_WORLD, nullptr, path, PETSC_FALSE); CHKERRQ(ierr);
	}
	ierr = PetscOptionsGetInt(
		nullptr, prefix, "-start", &ctx->start, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(
		nullptr, prefix, "-end", &ctx->end, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(
		nullptr, prefix, "-step", &ctx->step, &found); CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
} // PetibmTimeStepsGetOptions
