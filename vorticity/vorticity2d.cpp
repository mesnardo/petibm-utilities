/*! Computes the vorticity field from the 2D velocity vector field.
 * \file vorticity2d.cpp
 */

#include <iomanip>
#include <petscsys.h>
#include <petscdmda.h>
#include <petscviewerhdf5.h>


typedef struct
{
  DM da;
  Vec x, y;
  Vec global, local;
} Field;


/*! Initializes a Field structure.
 *
 * Creates the vectors based on the DMDA object.
 *
 * \param field The Field structure to initialize (passed by reference).
 */
PetscErrorCode FieldInitialize(Field &field)
{
  PetscErrorCode ierr;
  DMDALocalInfo info;

  PetscFunctionBeginUser;

  ierr = DMCreateGlobalVector(field.da, &field.global); CHKERRQ(ierr);
  ierr = DMCreateLocalVector(field.da, &field.local);
  ierr = DMDAGetLocalInfo(field.da, &info);
  ierr = VecCreateSeq(PETSC_COMM_SELF, info.mx, &field.x); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, info.my, &field.y); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // FieldInitialize


/*! Destroys the PETSc objects of a Field structure.
 *
 * \param field The Field structure (passed by reference).
 */
PetscErrorCode FieldDestroy(Field &field)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = VecDestroy(&field.x); CHKERRQ(ierr);
  ierr = VecDestroy(&field.y); CHKERRQ(ierr);
  ierr = VecDestroy(&field.global); CHKERRQ(ierr);
  ierr = VecDestroy(&field.local); CHKERRQ(ierr);
  ierr = DMDestroy(&field.da); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // FieldDestroy


/*! Reads the field values from a HDF5 file.
 *
 * \param filePath Path of the file to read.
 * \param name The name of the field.
 * \param field The field structure (passed by reference).
 */
PetscErrorCode FieldReadValues(std::string filePath, std::string name, Field &field)
{
  PetscErrorCode ierr;
  PetscViewer viewer;

  PetscFunctionBeginUser;

  ierr = PetscObjectSetName((PetscObject) field.global, name.c_str()); CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
  ierr = PetscViewerSetType(viewer, PETSCVIEWERHDF5); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
  ierr = VecLoad(field.global, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // FieldReadValues


/*! Reads a HDF5 grid file for a given field.
 *
 * \param filePath Path of the grid file to read.
 * \param field The field structure (passed by reference).
 */
PetscErrorCode FieldReadGrid(std::string filePath, Field &field)
{
  PetscErrorCode ierr;
  PetscViewer viewer;
  DMDALocalInfo info;

  PetscFunctionBeginUser;

  ierr = DMDAGetLocalInfo(field.da, &info); CHKERRQ(ierr);

  ierr = PetscViewerHDF5Open(PETSC_COMM_SELF, filePath.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
  // read x-stations
  ierr = PetscObjectSetName((PetscObject) field.x, "x"); CHKERRQ(ierr);
  ierr = VecLoad(field.x, viewer); CHKERRQ(ierr);
  // read y-stations
  ierr = PetscObjectSetName((PetscObject) field.y, "y"); CHKERRQ(ierr);
  ierr = VecLoad(field.y, viewer); CHKERRQ(ierr);

  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // FieldReadGrid


/*! Writes a grid into a HDF5 file.
 *
 * \param filePath Path of the file to write into.
 * \param field Field structure that contains the grid.
 */
PetscErrorCode FieldWriteGrid(std::string filePath, Field field)
{
  PetscErrorCode ierr;
  PetscViewer viewer;
  DMDALocalInfo info;

  PetscFunctionBeginUser;

  ierr = DMDAGetLocalInfo(field.da, &info); CHKERRQ(ierr);

  ierr = PetscViewerHDF5Open(PETSC_COMM_SELF, filePath.c_str(), FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
  // read x-stations
  ierr = PetscObjectSetName((PetscObject) field.x, "x"); CHKERRQ(ierr);
  ierr = VecView(field.x, viewer); CHKERRQ(ierr);
  // read y-stations
  ierr = PetscObjectSetName((PetscObject) field.y, "y"); CHKERRQ(ierr);
  ierr = VecView(field.y, viewer); CHKERRQ(ierr);

  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // FieldWriteGrid


/*! Writes the field values into a HDF5 file.
 *
 * \param filePath Path of the file to write into.
 * \param name Name of the field.
 * \param field Field structure.
 */
PetscErrorCode FieldWriteValues(std::string filePath, std::string name, Field field)
{
  PetscErrorCode ierr;
  PetscViewer viewer;

  PetscFunctionBeginUser;
  
  ierr = PetscObjectSetName((PetscObject) field.global, name.c_str()); CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
  ierr = PetscViewerSetType(viewer, PETSCVIEWERHDF5); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
  ierr = VecView(field.global, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // FieldWriteValues


/*! Computes the gridlines for the vorticity in the z-direction.
 *
 * \param ux Field structure; velocity in the x-direction.
 * \param uy Field structure; velocity in the y-direction.
 * \param wz Field of the z-vorticity (passed by reference).
 */
PetscErrorCode ComputeGridVorticityZ(Field ux, Field uy, Field &wz)
{
  PetscErrorCode ierr;
  PetscReal *x, *y;
  PetscReal *x_m, *y_m;
  PetscInt nx, ny;

  PetscFunctionBeginUser;

  ierr = VecGetSize(wz.x, &nx); CHKERRQ(ierr);
  ierr = VecGetArray(wz.x, &x); CHKERRQ(ierr);
  ierr = VecGetArray(ux.x, &x_m); CHKERRQ(ierr);
  for (PetscInt i=0; i<nx; i++)
    x[i] = x_m[i];
  ierr = VecRestoreArray(ux.x, &x_m); CHKERRQ(ierr);
  ierr = VecRestoreArray(wz.x, &x); CHKERRQ(ierr);
  
  ierr = VecGetSize(wz.y, &ny); CHKERRQ(ierr);
  ierr = VecGetArray(wz.y, &y); CHKERRQ(ierr);
  ierr = VecGetArray(uy.y, &y_m); CHKERRQ(ierr);
  for (PetscInt j=0; j<ny; j++)
    y[j] = y_m[j];
  ierr = VecRestoreArray(wz.y, &y); CHKERRQ(ierr);
  ierr = VecRestoreArray(uy.y, &y_m); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // ComputeGridVorticityZ


/*! Computes the vorticity in the z-direction.
 *
 * First-order.
 *
 * \param ux Field structure; velocity in the x-direction.
 * \param uy Field structure; velocity in the y-direction.
 * \param wz Field of the z-vorticity (passed by reference).
 */
PetscErrorCode ComputeVorticityZ(Field ux, Field uy, Field &wz)
{
  PetscErrorCode ierr;
  DMDALocalInfo info;
  PetscInt i, j;
  PetscReal **wz_a, **ux_a, **uy_a;
  PetscReal *x_a, *y_a;
  PetscReal dx, dy, du, dv;

  PetscFunctionBeginUser;

  // create local velocity vectors from global ones
  ierr = DMGlobalToLocalBegin(ux.da, ux.global, INSERT_VALUES, ux.local); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ux.da, ux.global, INSERT_VALUES, ux.local); CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(uy.da, uy.global, INSERT_VALUES, uy.local); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(uy.da, uy.global, INSERT_VALUES, uy.local); CHKERRQ(ierr);

  // get access to vector values in multi-dimensional fashion
  ierr = DMDAVecGetArray(wz.da, wz.global, &wz_a); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ux.da, ux.local, &ux_a); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(uy.da, uy.local, &uy_a); CHKERRQ(ierr);
  ierr = VecGetArray(uy.x, &x_a); CHKERRQ(ierr);
  ierr = VecGetArray(ux.y, &y_a); CHKERRQ(ierr);
  
  ierr = DMDAGetLocalInfo(wz.da, &info); CHKERRQ(ierr);
  for (j=info.ys; j<info.ys+info.ym; j++)
  {
    for (i=info.xs; i<info.xs+info.xm; i++)
    {
      dx = x_a[i] - x_a[i-1];
      dy = y_a[j] - y_a[j-1];
      du = ux_a[j][i] - ux_a[j-1][i];
      dv = uy_a[j][i] - uy_a[j][i-1];
      wz_a[j][i] = dv/dx - du/dy;
    }
  }

  ierr = VecRestoreArray(uy.x, &x_a); CHKERRQ(ierr);
  ierr = VecRestoreArray(ux.y, &y_a); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ux.da, ux.local, &ux_a); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(uy.da, uy.local, &uy_a); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(wz.da, wz.global, &wz_a); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // FieldComputeVorticityZ


int main(int argc, char **argv)
{
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);

  // parse command-line options
  std::string directory(".");
  char dir[PETSC_MAX_PATH_LEN];
  PetscBool found;
  ierr = PetscOptionsGetString(NULL, NULL, "-directory", dir, sizeof(dir), &found); CHKERRQ(ierr);
  if (found)
    directory = dir;

  PetscInt nstart = 0,
           nend = 0,
           nstep = 1;
  ierr = PetscOptionsGetInt(NULL, NULL, "-nstart", &nstart, &found); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL, NULL, "-nend", &nend, &found); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL, NULL, "-nstep", &nstep, &found); CHKERRQ(ierr);

  PetscInt nx = 0,
           ny = 0;
  ierr = PetscOptionsGetInt(NULL, NULL, "-nx", &nx, &found); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL, NULL, "-ny", &ny, &found); CHKERRQ(ierr);
  
  PetscBool isPeriodic_x, isPeriodic_y;
  ierr = PetscOptionsGetBool(NULL, NULL, "-periodic_x", &isPeriodic_x, NULL); CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL, NULL, "-periodic_y", &isPeriodic_y, NULL); CHKERRQ(ierr);

  DMBoundaryType bType_x = (isPeriodic_x) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED,
                 bType_y = (isPeriodic_y) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;

  // create DMDA for phi
  Field phi;
  ierr = DMDACreate2d(PETSC_COMM_WORLD,
                      bType_x, bType_y,
                      DMDA_STENCIL_STAR,
                      nx, ny,
                      PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL,
                      &phi.da); CHKERRQ(ierr);
  ierr = PetscObjectViewFromOptions((PetscObject) phi.da, NULL, "-phi_dmda_view"); CHKERRQ(ierr);

  // create DMDA objects for velocity components from DMDA object for pressure
  Field ux, uy;
  PetscInt numX, numY;
  const PetscInt *plx, *ply;
  ierr = DMDAGetOwnershipRanges(phi.da, &plx, &ply, NULL); CHKERRQ(ierr);
  PetscInt m, n;
  ierr = DMDAGetInfo(phi.da, NULL, NULL, NULL, NULL, &m, &n, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  // x-component of velocity
  PetscInt *ulx, *uly;
  ierr = PetscMalloc(m*sizeof(*ulx), &ulx); CHKERRQ(ierr);
  ierr = PetscMalloc(n*sizeof(*uly), &uly); CHKERRQ(ierr);
  ierr = PetscMemcpy(ulx, plx, m*sizeof(*ulx)); CHKERRQ(ierr);
  ierr = PetscMemcpy(uly, ply, n*sizeof(*uly)); CHKERRQ(ierr);
  numX = nx;
  numY = ny;
  if (!isPeriodic_x)
  {
    ulx[m-1]--;
    numX--;
  }
  ierr = DMDACreate2d(PETSC_COMM_WORLD, 
                      bType_x, bType_y,
                      DMDA_STENCIL_BOX, 
                      numX, numY, m, n, 1, 1, ulx, uly,
                      &ux.da); CHKERRQ(ierr);
  ierr = PetscObjectViewFromOptions((PetscObject) ux.da, NULL, "-ux_dmda_view"); CHKERRQ(ierr);
  ierr = PetscFree(ulx); CHKERRQ(ierr);
  ierr = PetscFree(uly); CHKERRQ(ierr);
  // y-component of velocity
  PetscInt *vlx, *vly;
  ierr = PetscMalloc(m*sizeof(*vlx), &vlx); CHKERRQ(ierr);
  ierr = PetscMalloc(n*sizeof(*vly), &vly); CHKERRQ(ierr);
  ierr = PetscMemcpy(vlx, plx, m*sizeof(*vlx)); CHKERRQ(ierr);
  ierr = PetscMemcpy(vly, ply, n*sizeof(*vly)); CHKERRQ(ierr);
  numX = nx;
  numY = ny;
  if (!isPeriodic_y)
  {
    vly[n-1]--;
    numY--;
  }
  ierr = DMDACreate2d(PETSC_COMM_WORLD, 
                      bType_x, bType_y,
                      DMDA_STENCIL_BOX, 
                      numX, numY, m, n, 1, 1, vlx, vly,
                      &uy.da); CHKERRQ(ierr);
  ierr = PetscObjectViewFromOptions((PetscObject) uy.da, NULL, "-uy_dmda_view"); CHKERRQ(ierr);
  ierr = PetscFree(vlx); CHKERRQ(ierr);
  ierr = PetscFree(vly); CHKERRQ(ierr);

  ierr = FieldInitialize(ux); CHKERRQ(ierr);
  ierr = FieldInitialize(uy); CHKERRQ(ierr);
  ierr = FieldInitialize(phi); CHKERRQ(ierr);

  // read grids
  std::string gridsDirectory(directory + "/grids");
  ierr = FieldReadGrid(gridsDirectory + "/staggered-x.h5", ux); CHKERRQ(ierr);
  ierr = FieldReadGrid(gridsDirectory + "/staggered-y.h5", uy); CHKERRQ(ierr);
  ierr = FieldReadGrid(gridsDirectory + "/cell-centered.h5", phi); CHKERRQ(ierr);

  // create z-vorticity field
  Field wz;
  ierr = DMDACreate2d(PETSC_COMM_WORLD,
                      bType_x, bType_y,
                      DMDA_STENCIL_STAR,
                      nx-1, ny-1,
                      PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL,
                      &wz.da); CHKERRQ(ierr);
  ierr = PetscObjectViewFromOptions((PetscObject) wz.da, NULL, "-wz_dmda_view"); CHKERRQ(ierr);
  ierr = FieldInitialize(wz); CHKERRQ(ierr);
  ierr = ComputeGridVorticityZ(ux, uy, wz); CHKERRQ(ierr);

  PetscMPIInt rank;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  if (rank == 0)
  {
    ierr = FieldWriteGrid(gridsDirectory + "/wz.h5", wz); CHKERRQ(ierr);
  }

  for (PetscInt ite=nstart; ite<=nend; ite+=nstep)
  {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "[time-step %D]\n", ite); CHKERRQ(ierr);

    std::stringstream ss;
    ss << directory << "/" << std::setfill('0') << std::setw(7) << ite;
    std::string folder(ss.str());
  
    // read values
    ierr = FieldReadValues(folder + "/ux.h5", "ux", ux); CHKERRQ(ierr);
    ierr = FieldReadValues(folder + "/uy.h5", "uy", uy); CHKERRQ(ierr);
  
    // compute the z-vorticity
    ierr = ComputeVorticityZ(ux, uy, wz); CHKERRQ(ierr);
    ierr = FieldWriteValues(folder + "/wz.h5", "wz", wz); CHKERRQ(ierr);
  }

  ierr = FieldDestroy(ux); CHKERRQ(ierr);
  ierr = FieldDestroy(uy); CHKERRQ(ierr);
  ierr = FieldDestroy(phi); CHKERRQ(ierr);
  ierr = FieldDestroy(wz); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
} // main
