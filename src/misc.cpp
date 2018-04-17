/*! Implementation of miscellaneous functions.
 * \file misc.cpp
 */

#include <sys/stat.h>
#include <experimental/filesystem>

#include "petibm-utilities/misc.h"


namespace fs = std::experimental::filesystem;


/*! Gets a directory from the command-line arguments.
 *
 * \param directory [out] Directory (passed by pointer).
 * \param key [in] The key option one is seeking.
 * \param defval [in] Default value to use if not found.
 * \param create [in] Create the directory if PETSC_TRUE (default if PETSC_FALSE).
 */
PetscErrorCode PetibmGetDirectory(
	std::string *directory, const char key[],
	const std::string defval, const PetscBool create)
{
	PetscErrorCode ierr;
	char dir[PETSC_MAX_PATH_LEN];
	PetscBool found;

	PetscFunctionBeginUser;

	ierr = PetscOptionsGetString(
		nullptr, nullptr, key, dir, sizeof(dir), &found); CHKERRQ(ierr);
	*directory = (found)? dir : defval;
	if (create)
	{
		ierr = PetibmCreateDirectory(*directory); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
} // PetibmGetDirectory


/*! Gets a file path from the command-line arguments.
 *
 * \param fielpath [out] File path (passed by pointer).
 * \param key [in] The key option one is seeking.
 * \param defval [in] Default value to use if not found.
 */
PetscErrorCode PetibmGetFilePath(
	std::string *filepath, const char key[], const std::string defval)
{
	PetscErrorCode ierr;
	char path[PETSC_MAX_PATH_LEN];
	PetscBool found;

	PetscFunctionBeginUser;

	ierr = PetscOptionsGetString(
		nullptr, nullptr, key, path, sizeof(path), &found); CHKERRQ(ierr);
	*filepath = (found)? path : defval;

	PetscFunctionReturn(0);
} // PetibmGetDirectory


/*! Creates a directory.
 *
 * \param directory [in] Path of the directory to create.
 */
PetscErrorCode PetibmCreateDirectory(const std::string directory)
{
	PetscErrorCode ierr;
	PetscMPIInt rank;

	PetscFunctionBeginUser;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

	if (rank == 0)
	{
		if (!fs::is_directory(directory.c_str()) or !fs::exists(directory.c_str()))
		{
			fs::create_directories(directory.c_str());
		}
	}

	PetscFunctionReturn(0);
} // PetibmCreateDirectory


PetscErrorCode PetibmGetParentDirectory(
	const std::string filepath, std::string &directory)
{
	PetscFunctionBeginUser;

	fs::path p = filepath;
	directory = p.parent_path();

	PetscFunctionReturn(0);
} // PetibmGetParentDirectory


/*! Loads options from configuration file.
 *
 * \param prefix [in] Prefix of "-config_file".
 */
PetscErrorCode PetibmOptionsInsertFile(const char prefix[])
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

	PetscFunctionReturn(0);
} // PetibmOptionsInsertFile


/*! Counts the number of points within provided boundaries.
 *
 * \param x [in] The sequential vector containing the points.
 * \param start [in] The starting boundary.
 * \param end [in] The ending boundary.
 * \param num [out] The number of points (passed by pointer).
 */
PetscErrorCode PetibmGetNumPoints1D(
	const Vec x, const PetscReal start, const PetscReal end, PetscInt *num)
{
	PetscErrorCode ierr;
	PetscInt n;
	PetscReal *x_arr;
	PetscInt count;

	PetscFunctionBeginUser;

	ierr = VecGetSize(x, &n); CHKERRQ(ierr);
	ierr = VecGetArray(x, &x_arr); CHKERRQ(ierr);
	count = 0;
	for (PetscInt i=0; i<n; i++)
		if (start <= x_arr[i] and x_arr[i] < end)
			count++;
	ierr = VecRestoreArray(x, &x_arr); CHKERRQ(ierr);
	*num = count;

	PetscFunctionReturn(0);
} // PetibmGetNumPoints1D


/*! Get the global index of the inferior neighbor from a gridline for each
 * station along another gridline.
 *
 * \param lineB [in] The gridline where each station will be associated with a neighbor.
 * \param lineA [in] The reference gridline where to find neighbors.
 * \param Iv [out] Standard vector that will contains the indices (passed by reference).
 */
PetscErrorCode PetibmGetNeighbors1D(
	PetibmGridline lineB, PetibmGridline lineA, std::vector<PetscInt> &Iv)
{
	PetscErrorCode ierr;
	DMDALocalInfo infoA, infoB;
	PetscInt i, I, Istart;
	PetscReal *xA, *xB;

	PetscFunctionBeginUser;

	ierr = DMDAGetLocalInfo(lineA.da, &infoA); CHKERRQ(ierr);
	ierr = DMDAGetLocalInfo(lineB.da, &infoB); CHKERRQ(ierr);

	Iv.reserve(infoB.xm);
	ierr = DMDAVecGetArray(lineA.da, lineA.local, &xA); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lineB.da, lineB.local, &xB); CHKERRQ(ierr);
	Istart = infoA.gxs;
	for (i=infoB.xs; i<infoB.xs+infoB.xm; i++)
	{
		for (I=Istart; I<infoA.xs+infoA.xm; I++)
		{
			if (xA[I] <= xB[i] and xB[i] < xA[I+1])
			{
				Iv.push_back(I);
				break;
			}
		}
		Istart = I;
	}
	ierr = DMDAVecRestoreArray(lineA.da, lineA.local, &xA); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lineB.da, lineB.local, &xB); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // PetibmGetNeighbors1D


/*! Helper function to print on-process vector in a sequential manner.
 *
 * \param v [in] The vector to print.
 */
PetscErrorCode PetibmVecView(const Vec v)
{
	PetscErrorCode ierr;
	PetscMPIInt r, rank, size;

	PetscFunctionBeginUser;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
	ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

	for (r=0; r<size; r++)
	{
		if (r == rank)
		{
			ierr = PetscPrintf(
				PETSC_COMM_SELF, "[%d/%d]\n", rank+1, size); CHKERRQ(ierr);
			ierr = VecView(v, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
		}
		ierr = MPI_Barrier(PETSC_COMM_WORLD);
	}

	PetscFunctionReturn(0);
} // PetibmVecView


/*! Helper function to print on-process array of integers in a sequential manner.
 *
 * \param N [in] The number of elements.
 * \param v [in] The array of integers to print.
 */
PetscErrorCode PetibmIntView(const PetscInt N, const PetscInt v[])
{
	PetscErrorCode ierr;
	PetscMPIInt r, rank, size;

	PetscFunctionBeginUser;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
	ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

	for (r=0; r<size; r++)
	{
		if (r == rank)
		{
			ierr = PetscPrintf(
				PETSC_COMM_SELF, "[%d/%d]\n", rank+1, size); CHKERRQ(ierr);
			ierr = PetscIntView(N, v, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
		}
		ierr = MPI_Barrier(PETSC_COMM_WORLD);
	}

	PetscFunctionReturn(0);
} // PetibmIntStdVecView
