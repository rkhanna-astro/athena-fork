//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file linear_wave.cpp
//! \brief Linear wave problem generator for 1D/2D/3D problems.
//!
//! In 1D, the problem is setup along one of the three coordinate axes (specified by
//! setting [ang_2,ang_3] = 0.0 or PI/2 in the input file).  In 2D/3D this routine
//! automatically sets the wavevector along the domain diagonal.
//========================================================================================

// C headers

// C++ headers
#include <algorithm>  // min, max
#include <cmath>      // sqrt()
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <random>


// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../reconstruct/reconstruction.hpp"


#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// void DMRInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
//   Real time, Real dt,
//   int il, int iu, int jl, int ju, int kl, int ku, int ngh);

namespace {
// Parameters which define initial solution -- made global so that they can be shared
// with functions A1,2,3 which compute vector potentials
Real d0, p0, u0, bx0, by0, bz0, dby, dbz;
int wave_flag;
Real ang_2, ang_3; // Rotation angles about the y and z' axis
bool ang_2_vert, ang_3_vert; // Switches to set ang_2 and/or ang_3 to pi/2
Real sin_a2, cos_a2, sin_a3, cos_a3;
Real amp, lambda, k_par; // amplitude, Wavelength, 2*PI/wavelength
Real gam,gm1,iso_cs,vflow;
Real ev[NWAVE], rem[NWAVE][NWAVE], lem[NWAVE][NWAVE];

// functions to compute vector potential to initialize the solution
Real A1(const Real x1, const Real x2, const Real x3);
Real A2(const Real x1, const Real x2, const Real x3);
Real A3(const Real x1, const Real x2, const Real x3);

// function to compute eigenvectors of linear waves
void Eigensystem(const Real d, const Real v1, const Real v2, const Real v3,
                 const Real h, const Real b1, const Real b2, const Real b3,
                 const Real x, const Real y, Real eigenvalues[(NWAVE)],
                 Real right_eigenmatrix[(NWAVE)][(NWAVE)],
                 Real left_eigenmatrix[(NWAVE)][(NWAVE)]);

Real MaxV2(MeshBlock *pmb, int iout);
} // namespace

// AMR refinement condition
int RefinementCondition(MeshBlock *pmb);

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // read global parameters
  wave_flag = pin->GetInteger("problem", "wave_flag");
  amp = pin->GetReal("problem", "amp");
  vflow = pin->GetOrAddReal("problem", "vflow", 0.0);
  ang_2 = pin->GetOrAddReal("problem", "ang_2", -999.9);
  ang_3 = pin->GetOrAddReal("problem", "ang_3", -999.9);

  ang_2_vert = pin->GetOrAddBoolean("problem", "ang_2_vert", false);
  ang_3_vert = pin->GetOrAddBoolean("problem", "ang_3_vert", false);

  // initialize global variables  
  if (NON_BAROTROPIC_EOS) {
    gam   = pin->GetReal("hydro", "gamma");
    gm1 = (gam - 1.0);
  } else {
    iso_cs = pin->GetReal("hydro", "iso_sound_speed");
  }

  // For wavevector along coordinate axes, set desired values of ang_2/ang_3.
  //    For example, for 1D problem use ang_2 = ang_3 = 0.0
  //    For wavevector along grid diagonal, do not input values for ang_2/ang_3.
  // Code below will automatically calculate these imposing periodicity and exactly one
  // wavelength along each grid direction
  Real x1size = mesh_size.x1max - mesh_size.x1min;
  Real x2size = mesh_size.x2max - mesh_size.x2min;
  Real x3size = mesh_size.x3max - mesh_size.x3min;

  // User should never input -999.9 in angles
  if (ang_3 == -999.9) ang_3 = std::atan(x1size/x2size);
  sin_a3 = std::sin(ang_3);
  cos_a3 = std::cos(ang_3);

  // Override ang_3 input and hardcode vertical (along x2 axis) wavevector
  if (ang_3_vert) {
    sin_a3 = 1.0;
    cos_a3 = 0.0;
    ang_3 = 0.5*M_PI;
  }

  if (ang_2 == -999.9) ang_2 = std::atan(0.5*(x1size*cos_a3 + x2size*sin_a3)/x3size);
  sin_a2 = std::sin(ang_2);
  cos_a2 = std::cos(ang_2);

  // Override ang_2 input and hardcode vertical (along x3 axis) wavevector
  if (ang_2_vert) {
    sin_a2 = 1.0;
    cos_a2 = 0.0;
    ang_2 = 0.5*M_PI;
  }

  Real x1 = x1size*cos_a2*cos_a3;
  Real x2 = x2size*cos_a2*sin_a3;
  Real x3 = x3size*sin_a2;

  // For lambda choose the smaller of the 3
  lambda = 2 * x1;
  if (f2 && ang_3 != 0.0) lambda = std::min(lambda,x2);
  if (f3 && ang_2 != 0.0) lambda = std::min(lambda,x3);

  // If cos_a2 or cos_a3 = 0, need to override lambda
  if (ang_3_vert)
    lambda = x2;
  if (ang_2_vert)
    lambda = x3;

  // Initialize k_parallel
  k_par = 2.0*(PI)/lambda;

  // Compute eigenvectors, where the quantities u0 and bx0 are parallel to the
  // wavevector, and v0,w0,by0,bz0 are perpendicular.
  d0 = 1.0;
  p0 = 0.0;
  u0 = vflow;
  Real v0 = 0.0;
  Real w0 = 0.0;
  bx0 = 0.0;
  by0 = 1.58;
  bz0 = 0.0;
  Real xfact = 0.0;
  Real yfact = 1.0;
  Real h0 = 0.0;

  if (NON_BAROTROPIC_EOS) {
    p0 = 1.0/gam;
    h0 = ((p0/gm1 + 0.5*d0*(u0*u0 + v0*v0+w0*w0)) + p0)/d0;
    if (MAGNETIC_FIELDS_ENABLED) h0 += (bx0*bx0 + by0*by0 + bz0*bz0)/d0;
  }

  Eigensystem(d0, u0, v0, w0, h0, bx0, by0, bz0, xfact, yfact, ev, rem, lem);

  if (pin->GetOrAddBoolean("problem", "test", false) && ncycle==0) {
    // reinterpret tlim as the number of orbital periods
    Real ntlim = lambda/std::abs(ev[wave_flag])*tlim;
    tlim = ntlim;
    pin->SetReal("time", "tlim", ntlim);
  }

  if (adaptive)
    EnrollUserRefinementCondition(RefinementCondition);

  // primarily used for tests of decaying linear waves (might conditionally enroll):
  AllocateUserHistoryOutput(1);
  EnrollUserHistoryOutput(0, MaxV2, "max-v2", UserHistoryOperation::max);
  // EnrollUserBoundaryFunction(BoundaryFace::inner_x2, DMRInnerX2);
  return;
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief Compute L1 error in linear waves and output to file
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  if (!pin->GetOrAddBoolean("problem", "compute_error", false)) return;

  // Initialize errors to zero
  Real l1_err[NHYDRO+NFIELD]{}, max_err[NHYDRO+NFIELD]{};

  for (int b=0; b<nblocal; ++b) {
    MeshBlock *pmb = my_blocks(b);
    BoundaryValues *pbval = pmb->pbval;
    int il = pmb->is, iu = pmb->ie, jl = pmb->js, ju = pmb->je,
        kl = pmb->ks, ku = pmb->ke;
    // adjust loop limits for fourth order error calculation
    //------------------------------------------------
    if (pmb->precon->correct_err) {
      // Expand loop limits on all sides by one
      if (pbval->nblevel[1][1][0] != -1) il -= 1;
      if (pbval->nblevel[1][1][2] != -1) iu += 1;
      if (pbval->nblevel[1][0][1] != -1) jl -= 1;
      if (pbval->nblevel[1][2][1] != -1) ju += 1;
      if (pbval->nblevel[0][1][1] != -1) kl -= 1;
      if (pbval->nblevel[2][1][1] != -1) ku += 1;
    }
    // Save analytic solution of conserved variables in 4D scratch array
    AthenaArray<Real> cons_;
    // Even for MHD, there are only cell-centered mesh variables
    int ncells4 = NHYDRO + NFIELD;
    int nl = 0;
    int nu = ncells4 - 1;
    cons_.NewAthenaArray(ncells4, pmb->ncells3, pmb->ncells2, pmb->ncells1);

    //  Compute errors at cell centers
    for (int k=kl; k<=ku; k++) {
      for (int j=jl; j<=ju; j++) {
        for (int i=il; i<=iu; i++) {
          Real x = cos_a2*(pmb->pcoord->x1v(i)*cos_a3 + pmb->pcoord->x2v(j)*sin_a3)
                   + pmb->pcoord->x3v(k)*sin_a2;
          Real sn = std::sin(k_par*x);

          Real d1 = d0 + amp*sn*rem[0][wave_flag];
          Real mx = d0*vflow;
          Real my = 0.0;
          Real mz = 0.0;
          Real m1 = mx*cos_a2*cos_a3 - my*sin_a3 - mz*sin_a2*cos_a3;
          Real m2 = mx*cos_a2*sin_a3 + my*cos_a3 - mz*sin_a2*sin_a3;
          Real m3 = mx*sin_a2                    + mz*cos_a2;

          // Store analytic solution at cell-centers
          cons_(IDN,k,j,i) = d1;
          cons_(IM1,k,j,i) = m1;
          cons_(IM2,k,j,i) = m2;
          cons_(IM3,k,j,i) = m3;

          if (NON_BAROTROPIC_EOS) {
            Real e0 = p0/gm1 + 0.5*d0*u0*u0 + amp*sn*rem[4][wave_flag];
            if (MAGNETIC_FIELDS_ENABLED) {
              e0 += 0.5*(bx0*bx0 + by0*by0 + bz0*bz0);
              Real bx = bx0;
              Real by = by0 + amp*sn*rem[5][wave_flag];
              Real bz = bz0 + amp*sn*rem[6][wave_flag];
              Real b1 = bx*cos_a2*cos_a3 - by*sin_a3 - bz*sin_a2*cos_a3;
              Real b2 = bx*cos_a2*sin_a3 + by*cos_a3 - bz*sin_a2*sin_a3;
              Real b3 = bx*sin_a2                    + bz*cos_a2;
              cons_(NHYDRO+IB1,k,j,i) = b1;
              cons_(NHYDRO+IB2,k,j,i) = b2;
              cons_(NHYDRO+IB3,k,j,i) = b3;
            }
            cons_(IEN,k,j,i) = e0;
          }
        }
      }
    }
    // begin fourth-order error correction
    // -------------------------------
    if (pmb->precon->correct_err) {
      // Restore loop limits to real cells only
      il = pmb->is, iu = pmb->ie, jl = pmb->js, ju = pmb->je, kl = pmb->ks, ku = pmb->ke;

      // Compute and store Laplacian of cell-centered conserved variables, Hydro and Bcc
      AthenaArray<Real> delta_cons_;
      delta_cons_.NewAthenaArray(ncells4, pmb->ncells3, pmb->ncells2, pmb->ncells1);
      pmb->pcoord->Laplacian(cons_, delta_cons_, il, iu, jl, ju, kl, ku, nl, nu);

      // TODO(felker): assuming uniform mesh with dx1f=dx2f=dx3f, so this factors out
      // TODO(felker): also, this may need to be dx1v, since Laplacian is cell-centered
      Real h = pmb->pcoord->dx1f(il);  // pco->dx1f(i); inside loop
      Real C = (h*h)/24.0;

      // Compute fourth-order approximation to cell-averaged conserved variables
      for (int n=nl; n<=nu; ++n) {
        for (int k=kl; k<=ku; ++k) {
          for (int j=jl; j<=ju; ++j) {
            for (int i=il; i<=iu; ++i) {
              cons_(n,k,j,i) = cons_(n,k,j,i) + C*delta_cons_(n,k,j,i);
            }
          }
        }
      }
    } // end if (pmb->precon->correct_err)
    // ------- end fourth-order error calculation

    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=il; i<=iu; ++i) {
          // Load cell-averaged <U>, either midpoint approx. or fourth-order approx
          Real d1 = cons_(IDN,k,j,i);
          Real m1 = cons_(IM1,k,j,i);
          Real m2 = cons_(IM2,k,j,i);
          Real m3 = cons_(IM3,k,j,i);
          // Weight l1 error by cell volume
          Real vol = pmb->pcoord->GetCellVolume(k, j, i);

          l1_err[IDN] += std::abs(d1 - pmb->phydro->u(IDN,k,j,i))*vol;
          max_err[IDN] = std::max(
              static_cast<Real>(std::abs(d1 - pmb->phydro->u(IDN,k,j,i))),
              max_err[IDN]);
          l1_err[IM1] += std::abs(m1 - pmb->phydro->u(IM1,k,j,i))*vol;
          l1_err[IM2] += std::abs(m2 - pmb->phydro->u(IM2,k,j,i))*vol;
          l1_err[IM3] += std::abs(m3 - pmb->phydro->u(IM3,k,j,i))*vol;
          max_err[IM1] = std::max(
              static_cast<Real>(std::abs(m1 - pmb->phydro->u(IM1,k,j,i))),
              max_err[IM1]);
          max_err[IM2] = std::max(
              static_cast<Real>(std::abs(m2 - pmb->phydro->u(IM2,k,j,i))),
              max_err[IM2]);
          max_err[IM3] = std::max(
              static_cast<Real>(std::abs(m3 - pmb->phydro->u(IM3,k,j,i))),
              max_err[IM3]);

          if (NON_BAROTROPIC_EOS) {
            Real e0 = cons_(IEN,k,j,i);
            l1_err[IEN] += std::abs(e0 - pmb->phydro->u(IEN,k,j,i))*vol;
            max_err[IEN] = std::max(
                static_cast<Real>(std::abs(e0-pmb->phydro->u(IEN,k,j,i))),
                max_err[IEN]);
          }

          if (MAGNETIC_FIELDS_ENABLED) {
            Real b1 = cons_(NHYDRO+IB1,k,j,i);
            Real b2 = cons_(NHYDRO+IB2,k,j,i);
            Real b3 = cons_(NHYDRO+IB3,k,j,i);
            Real db1 = std::abs(b1 - pmb->pfield->bcc(IB1,k,j,i));
            Real db2 = std::abs(b2 - pmb->pfield->bcc(IB2,k,j,i));
            Real db3 = std::abs(b3 - pmb->pfield->bcc(IB3,k,j,i));

            l1_err[NHYDRO + IB1] += db1*vol;
            l1_err[NHYDRO + IB2] += db2*vol;
            l1_err[NHYDRO + IB3] += db3*vol;
            max_err[NHYDRO + IB1] = std::max(db1, max_err[NHYDRO+IB1]);
            max_err[NHYDRO + IB2] = std::max(db2, max_err[NHYDRO+IB2]);
            max_err[NHYDRO + IB3] = std::max(db3, max_err[NHYDRO+IB3]);
          }
        }
      }
    }
  }
  Real rms_err = 0.0, max_max_over_l1 = 0.0;

#ifdef MPI_PARALLEL
  if (Globals::my_rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, &l1_err, (NHYDRO+NFIELD), MPI_ATHENA_REAL, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &max_err, (NHYDRO+NFIELD), MPI_ATHENA_REAL, MPI_MAX, 0,
               MPI_COMM_WORLD);
  } else {
    MPI_Reduce(&l1_err, &l1_err, (NHYDRO+NFIELD), MPI_ATHENA_REAL, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&max_err, &max_err, (NHYDRO+NFIELD), MPI_ATHENA_REAL, MPI_MAX, 0,
               MPI_COMM_WORLD);
  }
#endif

  // only the root process outputs the data
  if (Globals::my_rank == 0) {
    // normalize errors by number of cells
    Real vol= (mesh_size.x1max - mesh_size.x1min)*(mesh_size.x2max - mesh_size.x2min)
              *(mesh_size.x3max - mesh_size.x3min);
    for (int i=0; i<(NHYDRO+NFIELD); ++i) l1_err[i] = l1_err[i]/vol;
    // compute rms error
    for (int i=0; i<(NHYDRO+NFIELD); ++i) {
      rms_err += SQR(l1_err[i]);
      max_max_over_l1 = std::max(max_max_over_l1, (max_err[i]/l1_err[i]));
    }
    rms_err = std::sqrt(rms_err);

    // open output file and write out errors
    std::string fname;
    fname.assign("linearwave-errors.dat");
    std::stringstream msg;
    FILE *pfile;

    // The file exists -- reopen the file in append mode
    if ((pfile = std::fopen(fname.c_str(), "r")) != nullptr) {
      if ((pfile = std::freopen(fname.c_str(), "a", pfile)) == nullptr) {
        msg << "### FATAL ERROR in function Mesh::UserWorkAfterLoop"
            << std::endl << "Error output file could not be opened" <<std::endl;
        ATHENA_ERROR(msg);
      }

      // The file does not exist -- open the file in write mode and add headers
    } else {
      if ((pfile = std::fopen(fname.c_str(), "w")) == nullptr) {
        msg << "### FATAL ERROR in function Mesh::UserWorkAfterLoop"
            << std::endl << "Error output file could not be opened" <<std::endl;
        ATHENA_ERROR(msg);
      }
      std::fprintf(pfile, "# Nx1  Nx2  Nx3  Ncycle  ");
      std::fprintf(pfile, "RMS-L1-Error  d_L1  M1_L1  M2_L1  M3_L1  E_L1 ");
      if (MAGNETIC_FIELDS_ENABLED) std::fprintf(pfile, "  B1c_L1  B2c_L1  B3c_L1");
      std::fprintf(pfile, "  Largest-Max/L1  d_max  M1_max  M2_max  M3_max  E_max ");
      if (MAGNETIC_FIELDS_ENABLED) std::fprintf(pfile, "  B1c_max  B2c_max  B3c_max");
      std::fprintf(pfile, "\n");
    }

    // write errors
    std::fprintf(pfile, "%d  %d", mesh_size.nx1, mesh_size.nx2);
    std::fprintf(pfile, "  %d  %d", mesh_size.nx3, ncycle);
    std::fprintf(pfile, "  %e  %e", rms_err, l1_err[IDN]);
    std::fprintf(pfile, "  %e  %e  %e", l1_err[IM1], l1_err[IM2], l1_err[IM3]);
    if (NON_BAROTROPIC_EOS)
      std::fprintf(pfile, "  %e", l1_err[IEN]);
    if (MAGNETIC_FIELDS_ENABLED) {
      std::fprintf(pfile, "  %e", l1_err[NHYDRO+IB1]);
      std::fprintf(pfile, "  %e", l1_err[NHYDRO+IB2]);
      std::fprintf(pfile, "  %e", l1_err[NHYDRO+IB3]);
    }
    std::fprintf(pfile, "  %e  %e  ", max_max_over_l1, max_err[IDN]);
    std::fprintf(pfile, "%e  %e  %e", max_err[IM1], max_err[IM2], max_err[IM3]);
    if (NON_BAROTROPIC_EOS)
      std::fprintf(pfile, "  %e", max_err[IEN]);
    if (MAGNETIC_FIELDS_ENABLED) {
      std::fprintf(pfile, "  %e", max_err[NHYDRO+IB1]);
      std::fprintf(pfile, "  %e", max_err[NHYDRO+IB2]);
      std::fprintf(pfile, "  %e", max_err[NHYDRO+IB3]);
    }
    std::fprintf(pfile, "\n");
    std::fclose(pfile);
  }

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Linear wave problem generator for 1D/2D/3D problems.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Initialize the magnetic fields.  Note wavevector, eigenvectors, and other variables
  // are set in InitUserMeshData

  if (MAGNETIC_FIELDS_ENABLED) {
    AthenaArray<Real> a1, a2, a3;
    // nxN != ncellsN, in general. Allocate to extend through 2*ghost, regardless # dim
    int nx1 = block_size.nx1 + 2*NGHOST;
    int nx2 = block_size.nx2 + 2*NGHOST;
    int nx3 = block_size.nx3 + 2*NGHOST;
    a1.NewAthenaArray(nx3, nx2, nx1);
    a2.NewAthenaArray(nx3, nx2, nx1);
    a3.NewAthenaArray(nx3, nx2, nx1);

    // wave amplitudes
    dby = amp*rem[NWAVE-2][wave_flag];
    dbz = amp*rem[NWAVE-1][wave_flag];

    int level = loc.level;
    // Initialize components of the vector potential
    if (block_size.nx3 > 1) {
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie+1; i++) {
            if ((pbval->nblevel[1][0][1]>level && j==js)
                || (pbval->nblevel[1][2][1]>level && j==je+1)
                || (pbval->nblevel[0][1][1]>level && k==ks)
                || (pbval->nblevel[2][1][1]>level && k==ke+1)
                || (pbval->nblevel[0][0][1]>level && j==js   && k==ks)
                || (pbval->nblevel[0][2][1]>level && j==je+1 && k==ks)
                || (pbval->nblevel[2][0][1]>level && j==js   && k==ke+1)
                || (pbval->nblevel[2][2][1]>level && j==je+1 && k==ke+1)) {
              Real x1l = pcoord->x1f(i)+0.25*pcoord->dx1f(i);
              Real x1r = pcoord->x1f(i)+0.75*pcoord->dx1f(i);
              a1(k,j,i) = 0.5*(A1(x1l, pcoord->x2f(j), pcoord->x3f(k)) +
                               A1(x1r, pcoord->x2f(j), pcoord->x3f(k)));
            } else {
              a1(k,j,i) = A1(pcoord->x1v(i), pcoord->x2f(j), pcoord->x3f(k));
            }

            if ((pbval->nblevel[1][1][0]>level && i==is)
                || (pbval->nblevel[1][1][2]>level && i==ie+1)
                || (pbval->nblevel[0][1][1]>level && k==ks)
                || (pbval->nblevel[2][1][1]>level && k==ke+1)
                || (pbval->nblevel[0][1][0]>level && i==is   && k==ks)
                || (pbval->nblevel[0][1][2]>level && i==ie+1 && k==ks)
                || (pbval->nblevel[2][1][0]>level && i==is   && k==ke+1)
                || (pbval->nblevel[2][1][2]>level && i==ie+1 && k==ke+1)) {
              Real x2l = pcoord->x2f(j)+0.25*pcoord->dx2f(j);
              Real x2r = pcoord->x2f(j)+0.75*pcoord->dx2f(j);
              a2(k,j,i) = 0.5*(A2(pcoord->x1f(i), x2l, pcoord->x3f(k)) +
                               A2(pcoord->x1f(i), x2r, pcoord->x3f(k)));
            } else {
              a2(k,j,i) = A2(pcoord->x1f(i), pcoord->x2v(j), pcoord->x3f(k));
            }

            if ((pbval->nblevel[1][1][0]>level && i==is)
                || (pbval->nblevel[1][1][2]>level && i==ie+1)
                || (pbval->nblevel[1][0][1]>level && j==js)
                || (pbval->nblevel[1][2][1]>level && j==je+1)
                || (pbval->nblevel[1][0][0]>level && i==is   && j==js)
                || (pbval->nblevel[1][0][2]>level && i==ie+1 && j==js)
                || (pbval->nblevel[1][2][0]>level && i==is   && j==je+1)
                || (pbval->nblevel[1][2][2]>level && i==ie+1 && j==je+1)) {
              Real x3l = pcoord->x3f(k)+0.25*pcoord->dx3f(k);
              Real x3r = pcoord->x3f(k)+0.75*pcoord->dx3f(k);
              a3(k,j,i) = 0.5*(A3(pcoord->x1f(i), pcoord->x2f(j), x3l) +
                               A3(pcoord->x1f(i), pcoord->x2f(j), x3r));
            } else {
              a3(k,j,i) = A3(pcoord->x1f(i), pcoord->x2f(j), pcoord->x3v(k));
            }
          }
        }
      }
    } else {
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie+1; i++) {
            if (i != ie+1)
              a1(k,j,i) = A1(pcoord->x1v(i), pcoord->x2f(j), pcoord->x3f(k));
            if (j != je+1)
              a2(k,j,i) = A2(pcoord->x1f(i), pcoord->x2v(j), pcoord->x3f(k));
            if (k != ke+1)
              a3(k,j,i) = A3(pcoord->x1f(i), pcoord->x2f(j), pcoord->x3v(k));
          }
        }
      }
    }

    // Initialize interface fields
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie+1; i++) {
          Real y = pcoord->x2v(j);

          Real sn = std::sin(k_par*y);
          // Real sn2 = std::sin(2.0*k_par*x);
          // Real sn3 = std::sin(4.0*k_par*x);

          // Real amp3 = amp / 2.0;
          // Real amp2 = amp / 3.0;
          // Real amp1 = amp / 6.0;

          // pfield->b.x1f(k,j,i) =  amp1 * sn1 + amp2 * sn2 + amp3 * sn3;

          pfield->b.x1f(k,j,i) =  0;
        }
      }
    }

    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie; i++) {
          pfield->b.x2f(k,j,i) = by0;
        }
      }
    }

    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          Real x = pcoord->x2v(j);
          Real sn = std::sin(k_par*x);
          pfield->b.x3f(k,j,i) = 0;
        }
      }
    }
  }

  // initialize conserved variables
  for (int k=ks; k<=ke; k++) {
    // Real linear_density = 1.0;
    // Real y_grid = 0.0;
    // Real step = 1.0 / 256;

    // Real scale_height = (iso_cs * iso_cs) / 1.0;

    for (int j=js; j<=je; j++) {
      // linear_density = 1.0;

      for (int i=is; i<=ie; i++) {
        Real x = cos_a2*(pcoord->x1v(i)*cos_a3 + pcoord->x2v(j)*sin_a3) +
                 pcoord->x3v(k)*sin_a2;
        Real sn = std::sin(k_par*x);
        
        // Uinform Density
        std::random_device rd;  // Non-deterministic random seed (hardware-based if possible)
        std::mt19937 gen(rd()); // Mersenne Twister random number engine (fast and high-quality)
        std::uniform_real_distribution<> dis(-1.0, 1.0);

        phydro->u(IDN, k, j, i) = d0 + amp * dis(gen);
        // std::cout << "Random number between -1 and 1: " << phydro->u(IDN, k, j, i) << std::endl;

        // Linear Density:
        // linear_density -= step;
        // phydro->u(IDN, k, j, i) = linear_density;

        // Exponential Density:
        // phydro->u(IDN, k, j, i) = std::exp( - y_grid / scale_height);

        // Hyperbolic Secant Density Variation:
        // phydro->u(IDN, k, j, i) = 1.0 / std::cosh(M_PI * exp_value/ 2.0);

        // std::cout << "Density: " << phydro->u(IDN, k, j, i) << std::endl;

        Real mx = d0*vflow;
        Real my = 0.0;
        Real mz = 0.0;

        phydro->u(IM1,k,j,i) = mx*cos_a2*cos_a3 - my*sin_a3 - mz*sin_a2*cos_a3;
        phydro->u(IM2,k,j,i) = mx*cos_a2*sin_a3 + my*cos_a3 - mz*sin_a2*sin_a3;
        phydro->u(IM3,k,j,i) = mx*sin_a2                    + mz*cos_a2;

        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = p0/gm1 + 0.5*d0*u0*u0;
          if (MAGNETIC_FIELDS_ENABLED) {
            phydro->u(IEN,k,j,i) += 0.5*(bx0*bx0+by0*by0+bz0*bz0);
          }
        }
      }
      // linear_density -= step;
      // y_grid += step;
    }
  }
  return;
}

// void DMRInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
//   Real time, Real dt,
//   int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
//     Real amplitude = 0.01;             // Max perturbation strength
//     Real frequency = 2.0 * M_PI / 1.0;  // One cycle per time unit

//     Real perturb = amplitude * std::sin(frequency * time / 10.0);

//     // copy face-centered magnetic fields into ghost zones
//     for (int k=kl; k<=ku; ++k) {
//       for (int j=1; j<=ngh; ++j) {
// #pragma omp simd
//         for (int i=il; i<=iu; ++i) {
//           b.x1f(k, jl-j, i) = perturb;
//         }
//       }
//     }
//   return;
// }

namespace {
//----------------------------------------------------------------------------------------
//! \fn Real A1(const Real x1,const Real x2,const Real x3)
//  \brief A1: 1-component of vector potential, using a gauge such that Ax = 0, and Ay,
//  Az are functions of x and y alone.

Real A1(const Real x1, const Real x2, const Real x3) {
  Real x =  x1*cos_a2*cos_a3 + x2*cos_a2*sin_a3 + x3*sin_a2;
  Real y = -x1*sin_a3        + x2*cos_a3;
  Real Ay =  bz0*x - (dbz/k_par)*std::cos(k_par*(x));
  Real Az = -by0*x + (dby/k_par)*std::cos(k_par*(x)) + bx0*y;

  return -Ay*sin_a3 - Az*sin_a2*cos_a3;
}

//----------------------------------------------------------------------------------------
//! \fn Real A2(const Real x1,const Real x2,const Real x3)
//  \brief A2: 2-component of vector potential

Real A2(const Real x1, const Real x2, const Real x3) {
  Real x =  x1*cos_a2*cos_a3 + x2*cos_a2*sin_a3 + x3*sin_a2;
  Real y = -x1*sin_a3        + x2*cos_a3;
  Real Ay =  bz0*x - (dbz/k_par)*std::cos(k_par*(x));
  Real Az = -by0*x + (dby/k_par)*std::cos(k_par*(x)) + bx0*y;

  return Ay*cos_a3 - Az*sin_a2*sin_a3;
}

//----------------------------------------------------------------------------------------
//! \fn Real A3(const Real x1,const Real x2,const Real x3)
//  \brief A3: 3-component of vector potential

Real A3(const Real x1, const Real x2, const Real x3) {
  Real x =  x1*cos_a2*cos_a3 + x2*cos_a2*sin_a3 + x3*sin_a2;
  Real y = -x1*sin_a3        + x2*cos_a3;
  Real Az = -by0*x + (dby/k_par)*std::cos(k_par*(x)) + bx0*y;

  return Az*cos_a2;
}

//----------------------------------------------------------------------------------------
//! \fn void Eigensystem()
//  \brief computes eigenvectors of linear waves

void Eigensystem(const Real d, const Real v1, const Real v2, const Real v3,
                 const Real h, const Real b1, const Real b2, const Real b3,
                 const Real x, const Real y, Real eigenvalues[(NWAVE)],
                 Real right_eigenmatrix[(NWAVE)][(NWAVE)],
                 Real left_eigenmatrix[(NWAVE)][(NWAVE)]) {
  if (MAGNETIC_FIELDS_ENABLED) {
    //--- Adiabatic MHD ---
    if (NON_BAROTROPIC_EOS) {
      Real vsq,btsq,bt_starsq,vaxsq,hp,twid_asq,cfsq,cf,cssq,cs;
      Real bt,bt_star,bet2,bet3,bet2_star,bet3_star,bet_starsq,vbet,alpha_f,alpha_s;
      Real isqrtd,sqrtd,s,twid_a,qf,qs,af_prime,as_prime,afpbb,aspbb,vax;
      Real norm,cff,css,af,as,afpb,aspb,q2_star,q3_star,vqstr;
      Real ct2,tsum,tdif,cf2_cs2;
      Real qa,qb,qc,qd;
      vsq = v1*v1 + v2*v2 + v3*v3;
      btsq = b2*b2 + b3*b3;
      bt_starsq = (gm1 - (gm1 - 1.0)*y)*btsq;
      vaxsq = b1*b1/d;
      hp = h - (vaxsq + btsq/d);
      twid_asq = std::max((gm1*(hp-0.5*vsq)-(gm1-1.0)*x), TINY_NUMBER);

      // Compute fast- and slow-magnetosonic speeds (eq. B18)
      ct2 = bt_starsq/d;
      tsum = vaxsq + ct2 + twid_asq;
      tdif = vaxsq + ct2 - twid_asq;
      cf2_cs2 = std::sqrt(tdif*tdif + 4.0*twid_asq*ct2);

      cfsq = 0.5*(tsum + cf2_cs2);
      cf = std::sqrt(cfsq);

      cssq = twid_asq*vaxsq/cfsq;
      cs = std::sqrt(cssq);

      // Compute beta(s) (eqs. A17, B20, B28)
      bt = std::sqrt(btsq);
      bt_star = std::sqrt(bt_starsq);
      if (bt == 0.0) {
        bet2 = 1.0;
        bet3 = 0.0;
      } else {
        bet2 = b2/bt;
        bet3 = b3/bt;
      }
      bet2_star = bet2/std::sqrt(gm1 - (gm1-1.0)*y);
      bet3_star = bet3/std::sqrt(gm1 - (gm1-1.0)*y);
      bet_starsq = bet2_star*bet2_star + bet3_star*bet3_star;
      vbet = v2*bet2_star + v3*bet3_star;

      // Compute alpha(s) (eq. A16)
      if ((cfsq - cssq) == 0.0) {
        alpha_f = 1.0;
        alpha_s = 0.0;
      } else if ( (twid_asq - cssq) <= 0.0) {
        alpha_f = 0.0;
        alpha_s = 1.0;
      } else if ( (cfsq - twid_asq) <= 0.0) {
        alpha_f = 1.0;
        alpha_s = 0.0;
      } else {
        alpha_f = std::sqrt((twid_asq - cssq)/(cfsq - cssq));
        alpha_s = std::sqrt((cfsq - twid_asq)/(cfsq - cssq));
      }

      // Compute Q(s) and A(s) (eq. A14-15), etc.
      sqrtd = std::sqrt(d);
      isqrtd = 1.0/sqrtd;
      s = SIGN(b1);
      twid_a = std::sqrt(twid_asq);
      qf = cf*alpha_f*s;
      qs = cs*alpha_s*s;
      af_prime = twid_a*alpha_f*isqrtd;
      as_prime = twid_a*alpha_s*isqrtd;
      afpbb = af_prime*bt_star*bet_starsq;
      aspbb = as_prime*bt_star*bet_starsq;

      // Compute eigenvalues (eq. B17)
      vax = std::sqrt(vaxsq);
      eigenvalues[0] = v1 - cf;
      eigenvalues[1] = v1 - vax;
      eigenvalues[2] = v1 - cs;
      eigenvalues[3] = v1;
      eigenvalues[4] = v1 + cs;
      eigenvalues[5] = v1 + vax;
      eigenvalues[6] = v1 + cf;

      // Right-eigenvectors, stored as COLUMNS (eq. B21) */
      right_eigenmatrix[0][0] = alpha_f;
      right_eigenmatrix[0][1] = 0.0;
      right_eigenmatrix[0][2] = alpha_s;
      right_eigenmatrix[0][3] = 1.0;
      right_eigenmatrix[0][4] = alpha_s;
      right_eigenmatrix[0][5] = 0.0;
      right_eigenmatrix[0][6] = alpha_f;

      right_eigenmatrix[1][0] = alpha_f*eigenvalues[0];
      right_eigenmatrix[1][1] = 0.0;
      right_eigenmatrix[1][2] = alpha_s*eigenvalues[2];
      right_eigenmatrix[1][3] = v1;
      right_eigenmatrix[1][4] = alpha_s*eigenvalues[4];
      right_eigenmatrix[1][5] = 0.0;
      right_eigenmatrix[1][6] = alpha_f*eigenvalues[6];

      qa = alpha_f*v2;
      qb = alpha_s*v2;
      qc = qs*bet2_star;
      qd = qf*bet2_star;
      right_eigenmatrix[2][0] = qa + qc;
      right_eigenmatrix[2][1] = -bet3;
      right_eigenmatrix[2][2] = qb - qd;
      right_eigenmatrix[2][3] = v2;
      right_eigenmatrix[2][4] = qb + qd;
      right_eigenmatrix[2][5] = bet3;
      right_eigenmatrix[2][6] = qa - qc;

      qa = alpha_f*v3;
      qb = alpha_s*v3;
      qc = qs*bet3_star;
      qd = qf*bet3_star;
      right_eigenmatrix[3][0] = qa + qc;
      right_eigenmatrix[3][1] = bet2;
      right_eigenmatrix[3][2] = qb - qd;
      right_eigenmatrix[3][3] = v3;
      right_eigenmatrix[3][4] = qb + qd;
      right_eigenmatrix[3][5] = -bet2;
      right_eigenmatrix[3][6] = qa - qc;

      right_eigenmatrix[4][0] = alpha_f*(hp - v1*cf) + qs*vbet + aspbb;
      right_eigenmatrix[4][1] = -(v2*bet3 - v3*bet2);
      right_eigenmatrix[4][2] = alpha_s*(hp - v1*cs) - qf*vbet - afpbb;
      right_eigenmatrix[4][3] = 0.5*vsq + (gm1-1.0)*x/gm1;
      right_eigenmatrix[4][4] = alpha_s*(hp + v1*cs) + qf*vbet - afpbb;
      right_eigenmatrix[4][5] = -right_eigenmatrix[4][1];
      right_eigenmatrix[4][6] = alpha_f*(hp + v1*cf) - qs*vbet + aspbb;

      right_eigenmatrix[5][0] = as_prime*bet2_star;
      right_eigenmatrix[5][1] = -bet3*s*isqrtd;
      right_eigenmatrix[5][2] = -af_prime*bet2_star;
      right_eigenmatrix[5][3] = 0.0;
      right_eigenmatrix[5][4] = right_eigenmatrix[5][2];
      right_eigenmatrix[5][5] = right_eigenmatrix[5][1];
      right_eigenmatrix[5][6] = right_eigenmatrix[5][0];

      right_eigenmatrix[6][0] = as_prime*bet3_star;
      right_eigenmatrix[6][1] = bet2*s*isqrtd;
      right_eigenmatrix[6][2] = -af_prime*bet3_star;
      right_eigenmatrix[6][3] = 0.0;
      right_eigenmatrix[6][4] = right_eigenmatrix[6][2];
      right_eigenmatrix[6][5] = right_eigenmatrix[6][1];
      right_eigenmatrix[6][6] = right_eigenmatrix[6][0];

      // Left-eigenvectors, stored as ROWS (eq. B29)
      // Normalize by 1/2a^{2}: quantities denoted by \hat{f}
      norm = 0.5/twid_asq;
      cff = norm*alpha_f*cf;
      css = norm*alpha_s*cs;
      qf *= norm;
      qs *= norm;
      af = norm*af_prime*d;
      as = norm*as_prime*d;
      afpb = norm*af_prime*bt_star;
      aspb = norm*as_prime*bt_star;

      // Normalize by (gamma-1)/2a^{2}: quantities denoted by \bar{f}
      norm *= gm1;
      alpha_f *= norm;
      alpha_s *= norm;
      q2_star = bet2_star/bet_starsq;
      q3_star = bet3_star/bet_starsq;
      vqstr = (v2*q2_star + v3*q3_star);
      norm *= 2.0;

      left_eigenmatrix[0][0] = alpha_f*(vsq-hp) + cff*(cf+v1) - qs*vqstr - aspb;
      left_eigenmatrix[0][1] = -alpha_f*v1 - cff;
      left_eigenmatrix[0][2] = -alpha_f*v2 + qs*q2_star;
      left_eigenmatrix[0][3] = -alpha_f*v3 + qs*q3_star;
      left_eigenmatrix[0][4] = alpha_f;
      left_eigenmatrix[0][5] = as*q2_star - alpha_f*b2;
      left_eigenmatrix[0][6] = as*q3_star - alpha_f*b3;

      left_eigenmatrix[1][0] = 0.5*(v2*bet3 - v3*bet2);
      left_eigenmatrix[1][1] = 0.0;
      left_eigenmatrix[1][2] = -0.5*bet3;
      left_eigenmatrix[1][3] = 0.5*bet2;
      left_eigenmatrix[1][4] = 0.0;
      left_eigenmatrix[1][5] = -0.5*sqrtd*bet3*s;
      left_eigenmatrix[1][6] = 0.5*sqrtd*bet2*s;

      left_eigenmatrix[2][0] = alpha_s*(vsq-hp) + css*(cs+v1) + qf*vqstr + afpb;
      left_eigenmatrix[2][1] = -alpha_s*v1 - css;
      left_eigenmatrix[2][2] = -alpha_s*v2 - qf*q2_star;
      left_eigenmatrix[2][3] = -alpha_s*v3 - qf*q3_star;
      left_eigenmatrix[2][4] = alpha_s;
      left_eigenmatrix[2][5] = -af*q2_star - alpha_s*b2;
      left_eigenmatrix[2][6] = -af*q3_star - alpha_s*b3;

      left_eigenmatrix[3][0] = 1.0 - norm*(0.5*vsq - (gm1-1.0)*x/gm1);
      left_eigenmatrix[3][1] = norm*v1;
      left_eigenmatrix[3][2] = norm*v2;
      left_eigenmatrix[3][3] = norm*v3;
      left_eigenmatrix[3][4] = -norm;
      left_eigenmatrix[3][5] = norm*b2;
      left_eigenmatrix[3][6] = norm*b3;

      left_eigenmatrix[4][0] = alpha_s*(vsq-hp) + css*(cs-v1) - qf*vqstr + afpb;
      left_eigenmatrix[4][1] = -alpha_s*v1 + css;
      left_eigenmatrix[4][2] = -alpha_s*v2 + qf*q2_star;
      left_eigenmatrix[4][3] = -alpha_s*v3 + qf*q3_star;
      left_eigenmatrix[4][4] = alpha_s;
      left_eigenmatrix[4][5] = left_eigenmatrix[2][5];
      left_eigenmatrix[4][6] = left_eigenmatrix[2][6];

      left_eigenmatrix[5][0] = -left_eigenmatrix[1][0];
      left_eigenmatrix[5][1] = 0.0;
      left_eigenmatrix[5][2] = -left_eigenmatrix[1][2];
      left_eigenmatrix[5][3] = -left_eigenmatrix[1][3];
      left_eigenmatrix[5][4] = 0.0;
      left_eigenmatrix[5][5] = left_eigenmatrix[1][5];
      left_eigenmatrix[5][6] = left_eigenmatrix[1][6];

      left_eigenmatrix[6][0] = alpha_f*(vsq-hp) + cff*(cf-v1) + qs*vqstr - aspb;
      left_eigenmatrix[6][1] = -alpha_f*v1 + cff;
      left_eigenmatrix[6][2] = -alpha_f*v2 - qs*q2_star;
      left_eigenmatrix[6][3] = -alpha_f*v3 - qs*q3_star;
      left_eigenmatrix[6][4] = alpha_f;
      left_eigenmatrix[6][5] = left_eigenmatrix[0][5];
      left_eigenmatrix[6][6] = left_eigenmatrix[0][6];

      //--- Isothermal MHD ---

    } else {
      Real btsq,bt_starsq,vaxsq,twid_csq,cfsq,cf,cssq,cs;
      Real bt,bt_star,bet2,bet3,bet2_star,bet3_star,bet_starsq,alpha_f,alpha_s;
      Real sqrtd,s,twid_c,qf,qs,af_prime,as_prime,vax;
      Real norm,cff,css,af,as,afpb,aspb,q2_star,q3_star,vqstr;
      Real ct2,tsum,tdif,cf2_cs2;
      Real di = 1.0/d;
      btsq = b2*b2 + b3*b3;
      bt_starsq = btsq*y;
      vaxsq = b1*b1*di;
      twid_csq = (iso_cs*iso_cs) + x;

      // Compute fast- and slow-magnetosonic speeds (eq. B39)
      ct2 = bt_starsq*di;
      tsum = vaxsq + ct2 + twid_csq;
      tdif = vaxsq + ct2 - twid_csq;
      cf2_cs2 = std::sqrt(tdif*tdif + 4.0*twid_csq*ct2);

      cfsq = 0.5*(tsum + cf2_cs2);
      cf = std::sqrt(cfsq);

      cssq = twid_csq*vaxsq/cfsq;
      cs = std::sqrt(cssq);

      // Compute beta's (eqs. A17, B28, B40)
      bt = std::sqrt(btsq);
      bt_star = std::sqrt(bt_starsq);
      if (bt == 0.0) {
        bet2 = 1.0;
        bet3 = 0.0;
      } else {
        bet2 = b2/bt;
        bet3 = b3/bt;
      }
      bet2_star = bet2/std::sqrt(y);
      bet3_star = bet3/std::sqrt(y);
      bet_starsq = bet2_star*bet2_star + bet3_star*bet3_star;

      // Compute alpha's (eq. A16)
      if ((cfsq-cssq) == 0.0) {
        alpha_f = 1.0;
        alpha_s = 0.0;
      } else if ((twid_csq - cssq) <= 0.0) {
        alpha_f = 0.0;
        alpha_s = 1.0;
      } else if ((cfsq - twid_csq) <= 0.0) {
        alpha_f = 1.0;
        alpha_s = 0.0;
      } else {
        alpha_f = std::sqrt((twid_csq - cssq)/(cfsq - cssq));
        alpha_s = std::sqrt((cfsq - twid_csq)/(cfsq - cssq));
      }

      // Compute Q's (eq. A14-15), etc.
      sqrtd = std::sqrt(d);
      s = SIGN(b1);
      twid_c = std::sqrt(twid_csq);
      qf = cf*alpha_f*s;
      qs = cs*alpha_s*s;
      af_prime = twid_c*alpha_f/sqrtd;
      as_prime = twid_c*alpha_s/sqrtd;

      // Compute eigenvalues (eq. B38)
      vax  = std::sqrt(vaxsq);
      eigenvalues[0] = v1 - cf;
      eigenvalues[1] = v1 - vax;
      eigenvalues[2] = v1 - cs;
      eigenvalues[3] = v1 + cs;
      eigenvalues[4] = v1 + vax;
      eigenvalues[5] = v1 + cf;

      // Right-eigenvectors, stored as COLUMNS (eq. B21)
      right_eigenmatrix[0][0] = alpha_f;
      right_eigenmatrix[1][0] = alpha_f*(v1 - cf);
      right_eigenmatrix[2][0] = alpha_f*v2 + qs*bet2_star;
      right_eigenmatrix[3][0] = alpha_f*v3 + qs*bet3_star;
      right_eigenmatrix[4][0] = as_prime*bet2_star;
      right_eigenmatrix[5][0] = as_prime*bet3_star;

      right_eigenmatrix[0][1] = 0.0;
      right_eigenmatrix[1][1] = 0.0;
      right_eigenmatrix[2][1] = -bet3;
      right_eigenmatrix[3][1] = bet2;
      right_eigenmatrix[4][1] = -bet3*s/sqrtd;
      right_eigenmatrix[5][1] = bet2*s/sqrtd;

      right_eigenmatrix[0][2] = alpha_s;
      right_eigenmatrix[1][2] = alpha_s*(v1 - cs);
      right_eigenmatrix[2][2] = alpha_s*v2 - qf*bet2_star;
      right_eigenmatrix[3][2] = alpha_s*v3 - qf*bet3_star;
      right_eigenmatrix[4][2] = -af_prime*bet2_star;
      right_eigenmatrix[5][2] = -af_prime*bet3_star;

      right_eigenmatrix[0][3] = alpha_s;
      right_eigenmatrix[1][3] = alpha_s*(v1 + cs);
      right_eigenmatrix[2][3] = alpha_s*v2 + qf*bet2_star;
      right_eigenmatrix[3][3] = alpha_s*v3 + qf*bet3_star;
      right_eigenmatrix[4][3] = right_eigenmatrix[4][2];
      right_eigenmatrix[5][3] = right_eigenmatrix[5][2];

      right_eigenmatrix[0][4] = 0.0;
      right_eigenmatrix[1][4] = 0.0;
      right_eigenmatrix[2][4] = bet3;
      right_eigenmatrix[3][4] = -bet2;
      right_eigenmatrix[4][4] = right_eigenmatrix[4][1];
      right_eigenmatrix[5][4] = right_eigenmatrix[5][1];

      right_eigenmatrix[0][5] = alpha_f;
      right_eigenmatrix[1][5] = alpha_f*(v1 + cf);
      right_eigenmatrix[2][5] = alpha_f*v2 - qs*bet2_star;
      right_eigenmatrix[3][5] = alpha_f*v3 - qs*bet3_star;
      right_eigenmatrix[4][5] = right_eigenmatrix[4][0];
      right_eigenmatrix[5][5] = right_eigenmatrix[5][0];

      // Left-eigenvectors, stored as ROWS (eq. B41)
      // Normalize by 1/2a^{2}: quantities denoted by \hat{f}
      norm = 0.5/twid_csq;
      cff = norm*alpha_f*cf;
      css = norm*alpha_s*cs;
      qf *= norm;
      qs *= norm;
      af = norm*af_prime*d;
      as = norm*as_prime*d;
      afpb = norm*af_prime*bt_star;
      aspb = norm*as_prime*bt_star;

      q2_star = bet2_star/bet_starsq;
      q3_star = bet3_star/bet_starsq;
      vqstr = (v2*q2_star + v3*q3_star);

      left_eigenmatrix[0][0] = cff*(cf+v1) - qs*vqstr - aspb;
      left_eigenmatrix[0][1] = -cff;
      left_eigenmatrix[0][2] = qs*q2_star;
      left_eigenmatrix[0][3] = qs*q3_star;
      left_eigenmatrix[0][4] = as*q2_star;
      left_eigenmatrix[0][5] = as*q3_star;

      left_eigenmatrix[1][0] = 0.5*(v2*bet3 - v3*bet2);
      left_eigenmatrix[1][1] = 0.0;
      left_eigenmatrix[1][2] = -0.5*bet3;
      left_eigenmatrix[1][3] = 0.5*bet2;
      left_eigenmatrix[1][4] = -0.5*sqrtd*bet3*s;
      left_eigenmatrix[1][5] = 0.5*sqrtd*bet2*s;

      left_eigenmatrix[2][0] = css*(cs+v1) + qf*vqstr + afpb;
      left_eigenmatrix[2][1] = -css;
      left_eigenmatrix[2][2] = -qf*q2_star;
      left_eigenmatrix[2][3] = -qf*q3_star;
      left_eigenmatrix[2][4] = -af*q2_star;
      left_eigenmatrix[2][5] = -af*q3_star;

      left_eigenmatrix[3][0] = css*(cs-v1) - qf*vqstr + afpb;
      left_eigenmatrix[3][1] = css;
      left_eigenmatrix[3][2] = -left_eigenmatrix[2][2];
      left_eigenmatrix[3][3] = -left_eigenmatrix[2][3];
      left_eigenmatrix[3][4] = left_eigenmatrix[2][4];
      left_eigenmatrix[3][5] = left_eigenmatrix[2][5];

      left_eigenmatrix[4][0] = -left_eigenmatrix[1][0];
      left_eigenmatrix[4][1] = 0.0;
      left_eigenmatrix[4][2] = -left_eigenmatrix[1][2];
      left_eigenmatrix[4][3] = -left_eigenmatrix[1][3];
      left_eigenmatrix[4][4] = left_eigenmatrix[1][4];
      left_eigenmatrix[4][5] = left_eigenmatrix[1][5];

      left_eigenmatrix[5][0] = cff*(cf-v1) + qs*vqstr - aspb;
      left_eigenmatrix[5][1] = cff;
      left_eigenmatrix[5][2] = -left_eigenmatrix[0][2];
      left_eigenmatrix[5][3] = -left_eigenmatrix[0][3];
      left_eigenmatrix[5][4] = left_eigenmatrix[0][4];
      left_eigenmatrix[5][5] = left_eigenmatrix[0][5];
    }
  } else {
    //--- Adiabatic Hydrodynamics ---
    if (NON_BAROTROPIC_EOS) {
      Real vsq = v1*v1 + v2*v2 + v3*v3;
      Real asq = gm1*std::max((h-0.5*vsq), TINY_NUMBER);
      Real a = std::sqrt(asq);

      // Compute eigenvalues (eq. B2)
      eigenvalues[0] = v1 - a;
      eigenvalues[1] = v1;
      eigenvalues[2] = v1;
      eigenvalues[3] = v1;
      eigenvalues[4] = v1 + a;

      // Right-eigenvectors, stored as COLUMNS (eq. B3)
      right_eigenmatrix[0][0] = 1.0;
      right_eigenmatrix[1][0] = v1 - a;
      right_eigenmatrix[2][0] = v2;
      right_eigenmatrix[3][0] = v3;
      right_eigenmatrix[4][0] = h - v1*a;

      right_eigenmatrix[0][1] = 0.0;
      right_eigenmatrix[1][1] = 0.0;
      right_eigenmatrix[2][1] = 1.0;
      right_eigenmatrix[3][1] = 0.0;
      right_eigenmatrix[4][1] = v2;

      right_eigenmatrix[0][2] = 0.0;
      right_eigenmatrix[1][2] = 0.0;
      right_eigenmatrix[2][2] = 0.0;
      right_eigenmatrix[3][2] = 1.0;
      right_eigenmatrix[4][2] = v3;

      right_eigenmatrix[0][3] = 1.0;
      right_eigenmatrix[1][3] = v1;
      right_eigenmatrix[2][3] = v2;
      right_eigenmatrix[3][3] = v3;
      right_eigenmatrix[4][3] = 0.5*vsq;

      right_eigenmatrix[0][4] = 1.0;
      right_eigenmatrix[1][4] = v1 + a;
      right_eigenmatrix[2][4] = v2;
      right_eigenmatrix[3][4] = v3;
      right_eigenmatrix[4][4] = h + v1*a;

      // Left-eigenvectors, stored as ROWS (eq. B4)
      Real na = 0.5/asq;
      left_eigenmatrix[0][0] = na*(0.5*gm1*vsq + v1*a);
      left_eigenmatrix[0][1] = -na*(gm1*v1 + a);
      left_eigenmatrix[0][2] = -na*gm1*v2;
      left_eigenmatrix[0][3] = -na*gm1*v3;
      left_eigenmatrix[0][4] = na*gm1;

      left_eigenmatrix[1][0] = -v2;
      left_eigenmatrix[1][1] = 0.0;
      left_eigenmatrix[1][2] = 1.0;
      left_eigenmatrix[1][3] = 0.0;
      left_eigenmatrix[1][4] = 0.0;

      left_eigenmatrix[2][0] = -v3;
      left_eigenmatrix[2][1] = 0.0;
      left_eigenmatrix[2][2] = 0.0;
      left_eigenmatrix[2][3] = 1.0;
      left_eigenmatrix[2][4] = 0.0;

      Real qa = gm1/asq;
      left_eigenmatrix[3][0] = 1.0 - na*gm1*vsq;
      left_eigenmatrix[3][1] = qa*v1;
      left_eigenmatrix[3][2] = qa*v2;
      left_eigenmatrix[3][3] = qa*v3;
      left_eigenmatrix[3][4] = -qa;

      left_eigenmatrix[4][0] = na*(0.5*gm1*vsq - v1*a);
      left_eigenmatrix[4][1] = -na*(gm1*v1 - a);
      left_eigenmatrix[4][2] = left_eigenmatrix[0][2];
      left_eigenmatrix[4][3] = left_eigenmatrix[0][3];
      left_eigenmatrix[4][4] = left_eigenmatrix[0][4];

      //--- Isothermal Hydrodynamics ---

    } else {
      // Compute eigenvalues (eq. B6)
      eigenvalues[0] = v1 - iso_cs;
      eigenvalues[1] = v1;
      eigenvalues[2] = v1;
      eigenvalues[3] = v1 + iso_cs;

      // Right-eigenvectors, stored as COLUMNS (eq. B3)
      right_eigenmatrix[0][0] = 1.0;
      right_eigenmatrix[1][0] = v1 - iso_cs;
      right_eigenmatrix[2][0] = v2;
      right_eigenmatrix[3][0] = v3;

      right_eigenmatrix[0][1] = 0.0;
      right_eigenmatrix[1][1] = 0.0;
      right_eigenmatrix[2][1] = 1.0;
      right_eigenmatrix[3][1] = 0.0;

      right_eigenmatrix[0][2] = 0.0;
      right_eigenmatrix[1][2] = 0.0;
      right_eigenmatrix[2][2] = 0.0;
      right_eigenmatrix[3][2] = 1.0;

      right_eigenmatrix[0][3] = 1.0;
      right_eigenmatrix[1][3] = v1 + iso_cs;
      right_eigenmatrix[2][3] = v2;
      right_eigenmatrix[3][3] = v3;

      // Left-eigenvectors, stored as ROWS (eq. B7)

      left_eigenmatrix[0][0] = 0.5*(1.0 + v1/iso_cs);
      left_eigenmatrix[0][1] = -0.5/iso_cs;
      left_eigenmatrix[0][2] = 0.0;
      left_eigenmatrix[0][3] = 0.0;

      left_eigenmatrix[1][0] = -v2;
      left_eigenmatrix[1][1] = 0.0;
      left_eigenmatrix[1][2] = 1.0;
      left_eigenmatrix[1][3] = 0.0;

      left_eigenmatrix[2][0] = -v3;
      left_eigenmatrix[2][1] = 0.0;
      left_eigenmatrix[2][2] = 0.0;
      left_eigenmatrix[2][3] = 1.0;

      left_eigenmatrix[3][0] = 0.5*(1.0 - v1/iso_cs);
      left_eigenmatrix[3][1] = 0.5/iso_cs;
      left_eigenmatrix[3][2] = 0.0;
      left_eigenmatrix[3][3] = 0.0;
    }
  }
}

Real MaxV2(MeshBlock *pmb, int iout) {
  Real max_v2 = 0.0;
  int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je, ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> &w = pmb->phydro->w;
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        max_v2 = std::max(std::abs(w(IVY,k,j,i)), max_v2);
      }
    }
  }
  return max_v2;
}
} // namespace

// refinement condition: density curvature
int RefinementCondition(MeshBlock *pmb) {
  AthenaArray<Real> &w = pmb->phydro->w;
  Real dmax = 0.0, dmin = 2.0*d0;  // max and min densities
  for (int k=pmb->ks; k<=pmb->ke; k++) {
    for (int j=pmb->js; j<=pmb->je; j++) {
      for (int i=pmb->is; i<=pmb->ie; i++) {
        if (w(IDN,k,j,i) > dmax) dmax = w(IDN,k,j,i);
        if (w(IDN,k,j,i) < dmin) dmin = w(IDN,k,j,i);
      }
    }
  }
  // refine : delta rho > 0.9*amp
  if (dmax-d0 > 0.9*amp*rem[0][wave_flag]) return 1;
  //  Real a=std::max(dmax-d0,d0-dmin);
  //  if (a > 0.9*amp*rem[0][wave_flag]) return 1;
  // derefinement: else
  return -1;
}


void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(1);
  return;
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++)
        user_out_var(0,k,j,i) = phydro->w(IDN,k,j,i)-d0;
    }
  }
  return;
}
