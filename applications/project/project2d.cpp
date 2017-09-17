/*! Projects the 2D PetIBM solution from one grid to another.
 * \file project2d.cpp
 */

#include <petscsys.h>


int main(int argc, char **argv)
{
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);

  // parse command-line to get simulation directory
  char dir[PETSC_MAX_PATH_LEN];
  PetscBool found;
  ierr = PetscOptionsGetString(NULL, NULL, "-directory", dir, sizeof(dir), &found); CHKERRQ(ierr);
  std::string directory(".");
  if (found)
    directory = dir;
  ierr = PetscPrintf(PETSC_COMM_WORLD, "directory: %s\n", directory.c_str()); CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;
} // main
