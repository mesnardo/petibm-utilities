/*! Definition of the structure and functions related to the grid.
 * \file grid.h
 */

#pragma once

#include <string>

#include <petscsys.h>
#include <petscvec.h>
#include <petscdmda.h>


/*! Structure holding information about the grid.
 */
struct PetibmGridCtx
{
	char path[PETSC_MAX_PATH_LEN];  /// path of the file containing the grid
	PetscInt nx = 0,  /// number of points in the x-direction
	         ny = 0,  /// number of points in the y-direction
	         nz = 0;  /// number of points in the z-direction
	PetscReal starts[3],  /// stating point in each direction
	          ends[3];  /// end point in each direction
}; // PetibmGridCtx


/*! Structure holding the decomposition and vectors associated with a gridline.
 */
struct PetibmGridline
{
	DM da;  /// 1D DMDA object
	Vec coords,  /// vector containing the gridline stations
	    local;  /// local ghosted vector for stations on process
}; // PetibmGridline


/*! Structure holding the grid.
 */
struct PetibmGrid
{
	PetscInt dim = 2;  /// dimension of the grid
	PetibmGridline x,  /// gridline in the x-direction
	               y,  /// gridline in the y-direction
	               z;  /// gridline in the z-direction
}; // PetibmGrid


/*! Gets options from command-line or config file.
 *
 * \param prefix String to prepend the name of the options.
 * \param ctx The PetibmGridCtx structure to fill.
 */
PetscErrorCode PetibmGridGetOptions(
	const char prefix[], PetibmGridCtx *ctx);


/*! Initializes the grid based on the context.
 *
 * Creates a 1D DMDA object for each direction and creates the local global
 * vectors associated with each DMDA.
 * The domain is decomposed along the y-direction only.
 *
 * \param ctx The grid context.
 * \param grid The grid to initialize (passed by reference).
 */
PetscErrorCode PetibmGridInitialize(
	const PetibmGridCtx ctx, PetibmGrid &grid);


/*! Initializes a grid, providing the coordinates, based on a reference grid.
 *
 * Creates a 1D DMDA objects for each direction and creates the local and global
 * vectors associated with each DMDA.
 * The coordinates are provided as an array of sequential vectors, each vectors
 * containing the stations along a gridline.
 * The decomposition followed the decomposition of a reference grid so that
 * physically closed points are located on the same process.
 *
 * \param other The grid used as a reference.
 * \param coords The stations along each direction.
 * \param grid The grid to initialize (passed by reference).
 */
PetscErrorCode PetibmGridInitialize(
	const PetibmGrid other, const Vec coords[], PetibmGrid &grid);


/*! Initializes a gridline, providing the coordinates, based on a reference
 * gridline.
 *
 * Creates a 1D DMDA object for the direction and creates the local and global
 * vectors associated with it.
 * The coordinates are provided as a sequential vector containing the stations
 * along the gridline.
 * The decomposition followed the decomposition of a reference gridline so that
 * physically closed points are located on the same process.
 *
 * \param other The gridline used as a reference.
 * \param coords The stations along the direction.
 * \param line The gridline to initialize (passed by reference).
 */
PetscErrorCode PetibmGridlineInitialize(
	PetibmGridline other, const Vec coords, PetibmGridline &line);


/*! Sets the stations at external boundary points for all directions.
 *
 * \param starts List of starting points.
 * \param ends List of ending points.
 * \param grid The grid to modify (passed by reference).
 */
PetscErrorCode PetibmGridSetBoundaryPoints(
	const PetscReal starts[], const PetscReal ends[],PetibmGrid &grid);


/*! Sets the starting end ending points for a gridline.
 *
 * \param start The starting point.
 * \param end The ending point.
 * \param line The gridline to modify (passed by reference).
 */
PetscErrorCode PetibmGridlineSetBoundaryPoints(
	const PetscReal start, const PetscReal end, PetibmGridline &line);


/*! Inserts global values into local vectors for the grid.
 *
 * \param grid The grid to work on (passed by reference).
 */
PetscErrorCode PetibmGridGlobalToLocal(PetibmGrid &grid);


/*! Inserts global values into local vector for the gridline.
 *
 * \param line The gridline to work on (passed by reference).
 */
PetscErrorCode PetibmGridlineGlobalToLocal(PetibmGridline &line);


/*! Destroys a PetibmGrid structure.
 *
 * \param grid The PetibmGrid structure to destroy (passed by reference).
 */
PetscErrorCode PetibmGridDestroy(PetibmGrid &grid);


/*! Destroys a PetibmGridline structure.
 *
 * \param line The PetibmGridline structure to destroy (passed by reference).
 */
PetscErrorCode PetibmGridlineDestroy(PetibmGridline &line);


/*! Reads the gridlines from file.
 *
 * The gridlines should be stored in HDF5 format in the same file.
 *
 * \param filepath Path of the input file.
 * \param grid The PetibmGrid object to fill (passed by reference).
 */
PetscErrorCode PetibmGridHDF5Read(
	const std::string filepath, PetibmGrid &grid);


/*! Reads the gridline stations from a file.
 *
 * The stations along the gridline should be stored in HDF5 format.
 *
 * \param filepath The path of the input file.
 * \param name The name of the gridline (the direction).
 * \param line The sequential vector to fill (passed by reference).
 */
PetscErrorCode PetibmGridlineHDF5Read(
	std::string filepath, std::string name, Vec &line);

/*! Writes the gridlines into file in HDF5 format.
 *
 * \param filepath Path of the output file.
 * \param grid The PetibmGrid structure containing the gridlines.
 */
PetscErrorCode PetibmGridHDF5Write(
	const std::string filepath, const PetibmGrid grid);
