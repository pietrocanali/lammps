// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)
                         Aidan Thompson (SNL)

   Please cite the related publication:
   Tranchida, J., Plimpton, S. J., Thibaudeau, P., & Thompson, A. P. (2018).
   Massively parallel symplectic algorithm for coupled magnetic spin dynamics
   and molecular dynamics. Journal of Computational Physics.
------------------------------------------------------------------------- */

#include "pair_spin_magelec.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "info.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>
#include <utility>
#include <vector>

using namespace LAMMPS_NS;

namespace {
constexpr double EANGSTROM_TO_C_PER_M2 = 16.02176634;

struct PolCandidate {
  int index;
  tagint tag;
  tagint molecule;
  double rsq;
};
}

/* ---------------------------------------------------------------------- */

PairSpinMagelec::~PairSpinMagelec()
{
  delete[] pvector;
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cut_spin_magelec);
    memory->destroy(ME);
    memory->destroy(ME_mech);
    memory->destroy(v_mex);
    memory->destroy(v_mey);
    memory->destroy(v_mez);
    memory->destroy(dynamic_pair);
    memory->destroy(gamma_me);
    memory->destroy(cutsq); // to be deteled
    memory->destroy(emag);
  }
}

/* ----------------------------------------------------------------------
   request the regular spin list and, for dynamic polarization, a second
   full all-type list with the independently specified polarization cutoff
------------------------------------------------------------------------- */

void PairSpinMagelec::init_style()
{
  PairSpin::init_style();

  if (!dynamic_configured) return;
  if (!atom->q_flag)
    error->all(FLERR,"Pair spin/magelec polarization requires atom attribute q");
  if (!atom->molecule_flag)
    error->all(FLERR,"Pair spin/magelec polarization requires molecule IDs to pair cores and shells");
  if (!atom->tag_enable)
    error->all(FLERR,"Pair spin/magelec polarization requires atom IDs");

  auto *request = neighbor->add_request(this,NeighConst::REQ_FULL);
  request->set_id(1);
  request->set_cutoff_fixed(pol_cutoff);
}

/* ---------------------------------------------------------------------- */

void PairSpinMagelec::init_list(int id, NeighList *ptr)
{
  if (id == 0)
    list = ptr;
  else if (id == 1)
    pol_list = ptr;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSpinMagelec::settings(int narg, char **arg)
{

  PairSpin::settings(narg,arg);

  cut_spin_magelec_global = utils::numeric(FLERR,arg[0],false,lmp);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut_spin_magelec[i][j] = cut_spin_magelec_global;
        }
  }

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type spin pairs (only one for now)
------------------------------------------------------------------------- */

void PairSpinMagelec::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  int mode;
  double rij, magelec = 0.0, gamma = 0.0;
  double mex = 0.0, mey = 0.0, mez = 0.0;

  if (strcmp(arg[2],"magelec") == 0) {
    if (narg != 8) error->all(FLERR,"Incorrect args in pair_coeff spin/magelec command");
    mode = FIXED_COUPLING;
    rij = utils::numeric(FLERR,arg[3],false,lmp);
    magelec = utils::numeric(FLERR,arg[4],false,lmp);
    mex = utils::numeric(FLERR,arg[5],false,lmp);
    mey = utils::numeric(FLERR,arg[6],false,lmp);
    mez = utils::numeric(FLERR,arg[7],false,lmp);

    // Preserve the historical normalization exactly for legacy inputs.
    const double inorm = 1.0/(mex*mex+mey*mey+mez*mez);
    mex *= inorm;
    mey *= inorm;
    mez *= inorm;

  } else if (strcmp(arg[2],"polarization") == 0) {
    mode = POLARIZATION_COUPLING;
    if (narg < 12 || strcmp(arg[5],"center") != 0 || strcmp(arg[8],"neighbors") != 0)
      error->all(FLERR,"Incorrect args in pair_coeff spin/magelec polarization command");

    rij = utils::numeric(FLERR,arg[3],false,lmp);
    gamma = utils::numeric(FLERR,arg[4],false,lmp);
    const int new_center_core = utils::inumeric(FLERR,arg[6],false,lmp);
    const int new_center_shell = utils::inumeric(FLERR,arg[7],false,lmp);
    const int nspecies = utils::inumeric(FLERR,arg[9],false,lmp);
    if (nspecies <= 0 || narg != 12+4*nspecies)
      error->all(FLERR,"Incorrect args in pair_coeff spin/magelec polarization command");
    if (new_center_core < 1 || new_center_core > atom->ntypes ||
        new_center_shell < 1 || new_center_shell > atom->ntypes ||
        new_center_core == new_center_shell)
      error->all(FLERR,"Invalid center types in pair_coeff spin/magelec polarization command");

    std::vector<PolarizationSpecies> new_species;
    int iarg = 10;
    for (int m = 0; m < nspecies; m++) {
      PolarizationSpecies species;
      species.core_type = utils::inumeric(FLERR,arg[iarg++],false,lmp);
      species.shell_type = utils::inumeric(FLERR,arg[iarg++],false,lmp);
      species.count = utils::inumeric(FLERR,arg[iarg++],false,lmp);
      species.weight = utils::numeric(FLERR,arg[iarg++],false,lmp);
      if (species.core_type < 1 || species.core_type > atom->ntypes ||
          species.shell_type < 1 || species.shell_type > atom->ntypes ||
          species.core_type == species.shell_type || species.count <= 0)
        error->all(FLERR,"Invalid neighbor specification in pair_coeff spin/magelec polarization command");
      new_species.push_back(species);
    }
    if (strcmp(arg[iarg],"polcut") != 0)
      error->all(FLERR,"Expected polcut in pair_coeff spin/magelec polarization command");
    const double new_pol_cutoff = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (rij <= 0.0 || new_pol_cutoff <= 0.0)
      error->all(FLERR,"Spin/magelec cutoffs must be positive");

    if (dynamic_configured) {
      int same = new_center_core == center_core_type && new_center_shell == center_shell_type &&
        new_pol_cutoff == pol_cutoff && new_species.size() == pol_species.size();
      for (std::size_t m = 0; same && m < new_species.size(); m++)
        same = new_species[m].core_type == pol_species[m].core_type &&
          new_species[m].shell_type == pol_species[m].shell_type &&
          new_species[m].count == pol_species[m].count &&
          new_species[m].weight == pol_species[m].weight;
      if (!same)
        error->all(FLERR,"All spin/magelec polarization coefficients must use the same cell definition");
    } else {
      center_core_type = new_center_core;
      center_shell_type = new_center_shell;
      pol_cutoff = new_pol_cutoff;
      pol_species = new_species;
      dynamic_configured = 1;
    }
  } else
    error->all(FLERR,"Expected magelec or polarization in pair_coeff spin/magelec command");

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut_spin_magelec[i][j] = rij;
      ME[i][j] = magelec/hbar;
      ME_mech[i][j] = magelec;
      v_mex[i][j] = mex;
      v_mey[i][j] = mey;
      v_mez[i][j] = mez;
      dynamic_pair[i][j] = mode;
      gamma_me[i][j] = gamma;
      setflag[i][j] = 1;
      count++;
    }
  }
  if (count == 0)
    error->all(FLERR,"Incorrect args in pair_style command");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSpinMagelec::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    error->all(FLERR, Error::NOLASTLINE,
               "All pair coeffs are not set. Status\n" + Info::get_pair_coeff_status(lmp));

  ME[j][i] = ME[i][j];
  ME_mech[j][i] = ME_mech[i][j];
  v_mex[j][i] = v_mex[i][j];
  v_mey[j][i] = v_mey[i][j];
  v_mez[j][i] = v_mez[i][j];
  dynamic_pair[j][i] = dynamic_pair[i][j];
  gamma_me[j][i] = gamma_me[i][j];
  cut_spin_magelec[j][i] = cut_spin_magelec[i][j];

  return cut_spin_magelec_global;
}

/* ----------------------------------------------------------------------
   extract the larger cutoff
------------------------------------------------------------------------- */

void *PairSpinMagelec::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"cut") == 0) return (void *) &cut_spin_magelec_global;
  if (strcmp(str,"polarization") == 0) {
    dim = 1;
    return (void *) polarization;
  }
  return nullptr;
}

/* ----------------------------------------------------------------------
   Build Bi-centered local dipoles.  Cores and shells are paired by their
   molecule ID; every periodic atom ID is retained only at its nearest image.
------------------------------------------------------------------------- */

void PairSpinMagelec::compute_polarization()
{
  if (!pol_list)
    error->all(FLERR,"Spin/magelec polarization neighbor list is not initialized");

  // pair_style hybrid may turn a sub-style request into a type-skipped child.
  // Walk to its unskipped parent; distances are filtered to pol_cutoff below.
  NeighList *all_type_list = pol_list;
  while (all_type_list->listcopy) all_type_list = all_type_list->listcopy;
  if (all_type_list->listskip) all_type_list = all_type_list->listskip;

  double **x = atom->x;
  double *q = atom->q;
  int *type = atom->type;
  tagint *tag = atom->tag;
  tagint *molecule = atom->molecule;
  int *numneigh = all_type_list->numneigh;
  int **firstneigh = all_type_list->firstneigh;
  const int nlocal = atom->nlocal;
  const int nall = atom->nlocal+atom->nghost;
  const double cutsq = pol_cutoff*pol_cutoff;

  double local_dipole[3] = {0.0,0.0,0.0};
  int local_centers = 0;
  pol_contributions.clear();

  for (int center = 0; center < nlocal; center++) {
    if (type[center] != center_core_type) continue;
    if (molecule[center] <= 0)
      error->one(FLERR,"Spin/magelec polarization center atom {} has invalid molecule ID {}",
                 tag[center],molecule[center]);
    local_centers++;

    std::map<tagint,PolCandidate> nearest_by_tag;
    auto add_candidate = [&](int j) {
      // A periodic ghost of the center retains the center's atom ID.  Exclude
      // it before distance sorting so it cannot replace a genuine neighbor.
      if (tag[j] == tag[center]) return;

      // Local and ghost coordinates already contain LAMMPS's periodic image
      // translation.  Reapplying minimum_image() here can collapse a periodic
      // self-image onto the center and can select an image inconsistent with
      // the neighbor list.
      const double delx = x[j][0]-x[center][0];
      const double dely = x[j][1]-x[center][1];
      const double delz = x[j][2]-x[center][2];
      const double rsq = delx*delx+dely*dely+delz*delz;
      if (rsq > cutsq) return;

      auto iter = nearest_by_tag.find(tag[j]);
      if (iter == nearest_by_tag.end() || rsq < iter->second.rsq)
        nearest_by_tag[tag[j]] = {j,tag[j],molecule[j],rsq};
    };

    int *jlist = firstneigh[center];
    for (int jj = 0; jj < numneigh[center]; jj++) {
      const int j = jlist[jj] & NEIGHMASK;
      add_candidate(j);
    }

    // Core-shell bonds can be excluded from ordinary pair lists by
    // special_bonds.  Their atoms are nevertheless communicated as ghosts;
    // add shell types from that local+ghost set while retaining Rpol filtering.
    for (int j = 0; j < nall; j++) {
      int is_shell = type[j] == center_shell_type;
      for (const auto &species : pol_species)
        if (type[j] == species.shell_type) is_shell = 1;
      if (!is_shell) continue;
      add_candidate(j);
    }

    std::vector<PolCandidate> candidates;
    candidates.reserve(nearest_by_tag.size());
    for (const auto &entry : nearest_by_tag) candidates.push_back(entry.second);

    auto by_distance = [](const PolCandidate &a, const PolCandidate &b) {
      if (a.rsq != b.rsq) return a.rsq < b.rsq;
      return a.tag < b.tag;
    };

    auto add_contribution = [&](int index, double weight) {
      PolarizationContribution contribution;
      contribution.center = center;
      contribution.atom = index;
      contribution.weighted_charge = weight*q[index];
      contribution.displacement[0] = x[index][0]-x[center][0];
      contribution.displacement[1] = x[index][1]-x[center][1];
      contribution.displacement[2] = x[index][2]-x[center][2];
      for (int k = 0; k < 3; k++)
        local_dipole[k] += contribution.weighted_charge*contribution.displacement[k];
      pol_contributions.push_back(contribution);
    };

    std::vector<PolCandidate> center_shells;
    for (const auto &candidate : candidates)
      if (type[candidate.index] == center_shell_type && candidate.molecule == molecule[center])
        center_shells.push_back(candidate);
    std::sort(center_shells.begin(),center_shells.end(),by_distance);
    if (center_shells.empty())
      error->one(FLERR,"Spin/magelec polarization center atom {} expected 1 shell of type {} "
                 "but observed 0 within cutoff {}",tag[center],center_shell_type,pol_cutoff);
    add_contribution(center_shells.front().index,1.0);

    for (const auto &species : pol_species) {
      std::vector<PolCandidate> cores;
      for (const auto &candidate : candidates)
        if (type[candidate.index] == species.core_type) cores.push_back(candidate);
      std::sort(cores.begin(),cores.end(),by_distance);
      if (static_cast<int>(cores.size()) < species.count)
        error->one(FLERR,"Spin/magelec polarization center atom {} expected {} cores of type {} "
                   "but observed {} within cutoff {}",tag[center],species.count,species.core_type,
                   cores.size(),pol_cutoff);

      for (int n = 0; n < species.count; n++) {
        const auto &core = cores[n];
        if (core.molecule <= 0)
          error->one(FLERR,"Spin/magelec polarization core atom {} has invalid molecule ID {}",
                     core.tag,core.molecule);
        std::vector<PolCandidate> shells;
        for (const auto &candidate : candidates)
          if (type[candidate.index] == species.shell_type &&
              candidate.molecule == core.molecule)
            shells.push_back(candidate);
        std::sort(shells.begin(),shells.end(),by_distance);
        if (shells.empty())
          error->one(FLERR,"Spin/magelec polarization center atom {} found core atom {} but no "
                     "paired shell of type {} within cutoff {}",tag[center],core.tag,
                     species.shell_type,pol_cutoff);
        add_contribution(core.index,species.weight);
        add_contribution(shells.front().index,species.weight);
      }
    }
  }

  int global_centers;
  MPI_Allreduce(&local_centers,&global_centers,1,MPI_INT,MPI_SUM,world);
  if (global_centers == 0)
    error->all(FLERR,"Spin/magelec polarization found no center atoms of type {}",
               center_core_type);

  double global_dipole[3];
  MPI_Allreduce(local_dipole,global_dipole,3,MPI_DOUBLE,MPI_SUM,world);
  const double volume = domain->xprd*domain->yprd*domain->zprd;
  for (int k = 0; k < 3; k++) {
    pvector[k] = EANGSTROM_TO_C_PER_M2*local_dipole[k]/volume;
    polarization[k] = EANGSTROM_TO_C_PER_M2*global_dipole[k]/volume;
  }
}


/* ---------------------------------------------------------------------- */

void PairSpinMagelec::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double evdwl, ecoul;
  double xi[3], eij[3];
  double delx,dely,delz;
  double spi[3], spj[3];
  double fi[3], fmi[3];
  double local_cut2;
  double rsq, inorm;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;
  ev_init(eflag,vflag);

  pvector[0] = pvector[1] = pvector[2] = 0.0;
  if (dynamic_configured) compute_polarization();
  else polarization[0] = polarization[1] = polarization[2] = 0.0;

  double local_pol_gradient[3] = {0.0,0.0,0.0};
  double local_dynamic_energy = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  double **fm = atom->fm;
  double **sp = atom->sp;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // checking size of emag

  if (nlocal_max < nlocal) {                    // grow emag lists if necessary
    nlocal_max = nlocal;
    memory->grow(emag,nlocal_max,"pair/spin:emag");
  }

  // magneto-electric computation
  // loop over atoms and their neighbors

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    xi[0] = x[i][0];
    xi[1] = x[i][1];
    xi[2] = x[i][2];
    spi[0] = sp[i][0];
    spi[1] = sp[i][1];
    spi[2] = sp[i][2];
    emag[i] = 0.0;

    // loop on neighbors

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];

      spj[0] = sp[j][0];
      spj[1] = sp[j][1];
      spj[2] = sp[j][2];

      evdwl = 0.0;
      fi[0] = fi[1] = fi[2] = 0.0;
      fmi[0] = fmi[1] = fmi[2] = 0.0;

      delx = xi[0] - x[j][0];
      dely = xi[1] - x[j][1];
      delz = xi[2] - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq == 0.0) continue;
      inorm = 1.0/sqrt(rsq);
      eij[0] = -inorm*delx;
      eij[1] = -inorm*dely;
      eij[2] = -inorm*delz;

      local_cut2 = cut_spin_magelec[itype][jtype]*cut_spin_magelec[itype][jtype];

      // compute me interaction

      if (rsq <= local_cut2) {
        const int is_dynamic = dynamic_pair[itype][jtype] == POLARIZATION_COUPLING;
        if (is_dynamic) {
          double pair_gradient[3];
          compute_dynamic_spin(i,j,eij,fmi,spj);
          compute_dynamic_mechanical(i,j,eij,1.0/inorm,spi,spj,fi,pair_gradient,evdwl);
          for (int k = 0; k < 3; k++) local_pol_gradient[k] += pair_gradient[k];
          local_dynamic_energy += evdwl;
          if (!lattice_flag) fi[0] = fi[1] = fi[2] = 0.0;
          if (eflag) emag[i] += evdwl;
          else if (!vflag) evdwl = 0.0;
        } else {
          compute_magelec(i,j,eij,fmi,spj);

          if (lattice_flag)
            compute_magelec_mech(i,j,fi,spi,spj);

          if (eflag) {
            evdwl -= (spi[0]*fmi[0] + spi[1]*fmi[1] + spi[2]*fmi[2]);
            evdwl *= 0.5*hbar;
            emag[i] += evdwl;
          } else evdwl = 0.0;
        }

        f[i][0] += fi[0];
        f[i][1] += fi[1];
        f[i][2] += fi[2];
        if (newton_pair || j < nlocal) {
          f[j][0] -= fi[0];
          f[j][1] -= fi[1];
          f[j][2] -= fi[2];
        }
        fm[i][0] += fmi[0];
        fm[i][1] += fmi[1];
        fm[i][2] += fmi[2];

        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
            evdwl,ecoul,fi[0],fi[1],fi[2],delx,dely,delz);

        // Dynamic spin-lattice forces are generally non-central.  LAMMPS
        // stores six symmetric virial components, so replace each off-diagonal
        // x_a F_b entry tallied above by 1/2 (x_a F_b + x_b F_a).
        if (is_dynamic && vflag_either) {
          const double correction[3] = {
            0.5*(dely*fi[0]-delx*fi[1]),
            0.5*(delz*fi[0]-delx*fi[2]),
            0.5*(delz*fi[1]-dely*fi[2])};
          if (vflag_global) {
            virial[3] += correction[0];
            virial[4] += correction[1];
            virial[5] += correction[2];
          }
          if (vflag_atom) {
            for (int k = 0; k < 3; k++) {
              vatom[i][k+3] += 0.5*correction[k];
              vatom[j][k+3] += 0.5*correction[k];
            }
          }
        }
      }
    }
  }

  if (dynamic_configured && lattice_flag) {
    double global_pol_gradient[3];
    MPI_Allreduce(local_pol_gradient,global_pol_gradient,3,MPI_DOUBLE,MPI_SUM,world);
    apply_polarization_forces(global_pol_gradient,evflag);

    // Atomic virials contain the derivative of the cell dipole.  The explicit
    // 1/V normalization contributes +E to each diagonal virial component.
    if (vflag_global) {
      virial[0] += local_dynamic_energy;
      virial[1] += local_dynamic_energy;
      virial[2] += local_dynamic_energy;
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   update the pair interactions fmi acting on the spin ii
------------------------------------------------------------------------- */

void PairSpinMagelec::compute_single_pair(int ii, double fmi[3])
{
  int *type = atom->type;
  double **x = atom->x;
  double **sp = atom->sp;
  double local_cut2;
  double xi[3], eij[3];
  double delx,dely,delz;
  double spj[3];

  int j,jnum,itype,jtype,ntypes;
  int k,locflag;
  int *jlist,*numneigh,**firstneigh;

  double rsq, inorm;

  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // check if interaction applies to type of ii

  itype = type[ii];
  ntypes = atom->ntypes;
  locflag = 0;
  k = 1;
  while (k <= ntypes) {
    if (k <= itype) {
      if (setflag[k][itype] == 1) {
        locflag =1;
        break;
      }
      k++;
    } else if (k > itype) {
      if (setflag[itype][k] == 1) {
        locflag =1;
        break;
      }
      k++;
    } else error->all(FLERR,"Wrong type number");
  }

  // if interaction applies to type ii,
  // locflag = 1 and compute pair interaction

  if (locflag == 1) {

    xi[0] = x[ii][0];
    xi[1] = x[ii][1];
    xi[2] = x[ii][2];

    jlist = firstneigh[ii];
    jnum = numneigh[ii];

    for (int jj = 0; jj < jnum; jj++) {

      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];
      local_cut2 = cut_spin_magelec[itype][jtype]*cut_spin_magelec[itype][jtype];

      spj[0] = sp[j][0];
      spj[1] = sp[j][1];
      spj[2] = sp[j][2];

      delx = xi[0] - x[j][0];
      dely = xi[1] - x[j][1];
      delz = xi[2] - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq == 0.0) continue;
      inorm = 1.0/sqrt(rsq);
      eij[0] = -inorm*delx;
      eij[1] = -inorm*dely;
      eij[2] = -inorm*delz;

      if (rsq <= local_cut2) {
        if (dynamic_pair[itype][jtype] == POLARIZATION_COUPLING)
          compute_dynamic_spin(ii,j,eij,fmi,spj);
        else
          compute_magelec(ii,j,eij,fmi,spj);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairSpinMagelec::compute_dynamic_spin(int i, int j, const double eij[3], double fmi[3],
                                           const double spj[3])
{
  int *type = atom->type;
  const double gamma_over_hbar = gamma_me[type[i]][type[j]]/hbar;

  const double d[3] = {
    gamma_over_hbar*(polarization[1]*eij[2]-polarization[2]*eij[1]),
    gamma_over_hbar*(polarization[2]*eij[0]-polarization[0]*eij[2]),
    gamma_over_hbar*(polarization[0]*eij[1]-polarization[1]*eij[0])};

  // fm = -(1/hbar) dH/ds_i for H = gamma (P x e_ij).(s_i x s_j)
  fmi[0] -= spj[1]*d[2]-spj[2]*d[1];
  fmi[1] -= spj[2]*d[0]-spj[0]*d[2];
  fmi[2] -= spj[0]*d[1]-spj[1]*d[0];
}

/* ----------------------------------------------------------------------
   Half-weighted directed contribution for the full spin neighbor list.
   The reverse directed entry supplies the other half of force and energy.
------------------------------------------------------------------------- */

void PairSpinMagelec::compute_dynamic_mechanical(int i, int j, const double eij[3], double r,
                                                 const double spi[3], const double spj[3],
                                                 double fi[3], double gradient[3], double &energy)
{
  int *type = atom->type;
  const double gamma = gamma_me[type[i]][type[j]];
  const double cross_spin[3] = {
    spi[1]*spj[2]-spi[2]*spj[1],
    spi[2]*spj[0]-spi[0]*spj[2],
    spi[0]*spj[1]-spi[1]*spj[0]};
  const double p_cross_e[3] = {
    polarization[1]*eij[2]-polarization[2]*eij[1],
    polarization[2]*eij[0]-polarization[0]*eij[2],
    polarization[0]*eij[1]-polarization[1]*eij[0]};

  energy = 0.5*gamma*(p_cross_e[0]*cross_spin[0]+p_cross_e[1]*cross_spin[1]+
                      p_cross_e[2]*cross_spin[2]);

  // dH/dP = gamma (e x (s_i x s_j)); half for this directed entry.
  gradient[0] = 0.5*gamma*(eij[1]*cross_spin[2]-eij[2]*cross_spin[1]);
  gradient[1] = 0.5*gamma*(eij[2]*cross_spin[0]-eij[0]*cross_spin[2]);
  gradient[2] = 0.5*gamma*(eij[0]*cross_spin[1]-eij[1]*cross_spin[0]);

  // dH/de = gamma ((s_i x s_j) x P), projected perpendicular to e.
  double dedirection[3] = {
    gamma*(cross_spin[1]*polarization[2]-cross_spin[2]*polarization[1]),
    gamma*(cross_spin[2]*polarization[0]-cross_spin[0]*polarization[2]),
    gamma*(cross_spin[0]*polarization[1]-cross_spin[1]*polarization[0])};
  const double parallel = dedirection[0]*eij[0]+dedirection[1]*eij[1]+
    dedirection[2]*eij[2];
  for (int k = 0; k < 3; k++) fi[k] += 0.5*(dedirection[k]-parallel*eij[k])/r;
}

/* ----------------------------------------------------------------------
   Apply -dH/dP dP/dr to every selected contributor and the opposite force
   to its cell center.  The same center-relative force pairs tally the dipole
   part of the virial without depending on wrapped coordinates.
------------------------------------------------------------------------- */

void PairSpinMagelec::apply_polarization_forces(const double gradient[3], int tally)
{
  double **f = atom->f;
  const int nlocal = atom->nlocal;
  const int newton_pair = force->newton_pair;
  const double volume = domain->xprd*domain->yprd*domain->zprd;
  const double conversion = EANGSTROM_TO_C_PER_M2/volume;

  for (const auto &contribution : pol_contributions) {
    const double scale = conversion*contribution.weighted_charge;
    const double center_force[3] =
      {scale*gradient[0],scale*gradient[1],scale*gradient[2]};
    for (int k = 0; k < 3; k++) {
      f[contribution.center][k] += center_force[k];
      f[contribution.atom][k] -= center_force[k];
    }
    if (tally) {
      ev_tally_xyz(contribution.center,contribution.atom,nlocal,newton_pair,0.0,0.0,
                   center_force[0],center_force[1],center_force[2],
                   -contribution.displacement[0],-contribution.displacement[1],
                   -contribution.displacement[2]);
      if (vflag_either) {
        const double dx = -contribution.displacement[0];
        const double dy = -contribution.displacement[1];
        const double dz = -contribution.displacement[2];
        const double correction[3] = {
          0.5*(dy*center_force[0]-dx*center_force[1]),
          0.5*(dz*center_force[0]-dx*center_force[2]),
          0.5*(dz*center_force[1]-dy*center_force[2])};
        if (vflag_global) {
          virial[3] += correction[0];
          virial[4] += correction[1];
          virial[5] += correction[2];
        }
        if (vflag_atom) {
          for (int k = 0; k < 3; k++) {
            vatom[contribution.center][k+3] += 0.5*correction[k];
            vatom[contribution.atom][k+3] += 0.5*correction[k];
          }
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairSpinMagelec::compute_magelec(int i, int j, double eij[3], double fmi[3], double spj[3])
{
  int *type = atom->type;
  int itype, jtype;
  double meix,meiy,meiz;
  double vx,vy,vz;
  itype = type[i];
  jtype = type[j];

  vx = v_mex[itype][jtype];
  vy = v_mey[itype][jtype];
  vz = v_mez[itype][jtype];

  meix = (vy*eij[2] - vz*eij[1]);
  meiy = (vz*eij[0] - vx*eij[2]);
  meiz = (vx*eij[1] - vy*eij[0]);

  meix *= ME[itype][jtype];
  meiy *= ME[itype][jtype];
  meiz *= ME[itype][jtype];

  fmi[0] += (spj[1]*meiz - spj[2]*meiy);
  fmi[1] += (spj[2]*meix - spj[0]*meiz);
  fmi[2] += (spj[0]*meiy - spj[1]*meix);
}

/* ---------------------------------------------------------------------- */

void PairSpinMagelec::compute_magelec_mech(int i, int j, double fi[3], double spi[3], double spj[3])
{
  int *type = atom->type;
  int itype, jtype;
  itype = type[i];
  jtype = type[j];

  double meix,meiy,meiz;
  double vx, vy, vz;

  vx = v_mex[itype][jtype];
  vy = v_mey[itype][jtype];
  vz = v_mez[itype][jtype];

  meix = (spi[1]*spi[2] - spi[2]*spj[1]);
  meiy = (spi[2]*spi[0] - spi[0]*spj[2]);
  meiz = (spi[0]*spi[1] - spi[1]*spj[0]);

  meix *= ME_mech[itype][jtype];
  meiy *= ME_mech[itype][jtype];
  meiz *= ME_mech[itype][jtype];

  fi[0] += 0.5*(meiy*vz - meiz*vy);
  fi[1] += 0.5*(meiz*vx - meix*vz);
  fi[2] += 0.5*(meix*vy - meiy*vx);

}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairSpinMagelec::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cut_spin_magelec,n+1,n+1,"pair/spin/me:cut_spin_magelec");
  memory->create(ME,n+1,n+1,"pair/spin/me:ME");
  memory->create(ME_mech,n+1,n+1,"pair/spin/me:ME_mech");
  memory->create(v_mex,n+1,n+1,"pair/spin/me:ME_vector_x");
  memory->create(v_mey,n+1,n+1,"pair/spin/me:ME_vector_y");
  memory->create(v_mez,n+1,n+1,"pair/spin/me:ME_vector_z");
  memory->create(dynamic_pair,n+1,n+1,"pair/spin/me:dynamic_pair");
  memory->create(gamma_me,n+1,n+1,"pair/spin/me:gamma");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  for (int i = 0; i <= n; i++)
    for (int j = 0; j <= n; j++) {
      dynamic_pair[i][j] = FIXED_COUPLING;
      gamma_me[i][j] = 0.0;
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSpinMagelec::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&dynamic_pair[i][j],sizeof(int),1,fp);
        fwrite(&ME[i][j],sizeof(double),1,fp);
        fwrite(&ME_mech[i][j],sizeof(double),1,fp);
        fwrite(&gamma_me[i][j],sizeof(double),1,fp);
        fwrite(&v_mex[i][j],sizeof(double),1,fp);
        fwrite(&v_mey[i][j],sizeof(double),1,fp);
        fwrite(&v_mez[i][j],sizeof(double),1,fp);
        fwrite(&cut_spin_magelec[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSpinMagelec::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&dynamic_pair[i][j],sizeof(int),1,fp,nullptr,error);
          utils::sfread(FLERR,&ME[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&ME_mech[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&gamma_me[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&v_mex[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&v_mey[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&v_mez[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_spin_magelec[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&dynamic_pair[i][j],1,MPI_INT,0,world);
        MPI_Bcast(&ME[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&ME_mech[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&gamma_me[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&v_mex[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&v_mey[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&v_mez[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_spin_magelec[i][j],1,MPI_DOUBLE,0,world);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSpinMagelec::write_restart_settings(FILE *fp)
{
  fwrite(&cut_spin_magelec_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&dynamic_configured,sizeof(int),1,fp);
  if (dynamic_configured) {
    fwrite(&center_core_type,sizeof(int),1,fp);
    fwrite(&center_shell_type,sizeof(int),1,fp);
    fwrite(&pol_cutoff,sizeof(double),1,fp);
    const int nspecies = static_cast<int>(pol_species.size());
    fwrite(&nspecies,sizeof(int),1,fp);
    for (const auto &species : pol_species) {
      fwrite(&species.core_type,sizeof(int),1,fp);
      fwrite(&species.shell_type,sizeof(int),1,fp);
      fwrite(&species.count,sizeof(int),1,fp);
      fwrite(&species.weight,sizeof(double),1,fp);
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSpinMagelec::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&cut_spin_magelec_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&dynamic_configured,sizeof(int),1,fp,nullptr,error);
    if (dynamic_configured) {
      utils::sfread(FLERR,&center_core_type,sizeof(int),1,fp,nullptr,error);
      utils::sfread(FLERR,&center_shell_type,sizeof(int),1,fp,nullptr,error);
      utils::sfread(FLERR,&pol_cutoff,sizeof(double),1,fp,nullptr,error);
      int nspecies;
      utils::sfread(FLERR,&nspecies,sizeof(int),1,fp,nullptr,error);
      pol_species.resize(nspecies);
      for (auto &species : pol_species) {
        utils::sfread(FLERR,&species.core_type,sizeof(int),1,fp,nullptr,error);
        utils::sfread(FLERR,&species.shell_type,sizeof(int),1,fp,nullptr,error);
        utils::sfread(FLERR,&species.count,sizeof(int),1,fp,nullptr,error);
        utils::sfread(FLERR,&species.weight,sizeof(double),1,fp,nullptr,error);
      }
    }
  }
  MPI_Bcast(&cut_spin_magelec_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&dynamic_configured,1,MPI_INT,0,world);
  if (dynamic_configured) {
    MPI_Bcast(&center_core_type,1,MPI_INT,0,world);
    MPI_Bcast(&center_shell_type,1,MPI_INT,0,world);
    MPI_Bcast(&pol_cutoff,1,MPI_DOUBLE,0,world);
    int nspecies = static_cast<int>(pol_species.size());
    MPI_Bcast(&nspecies,1,MPI_INT,0,world);
    if (comm->me != 0) pol_species.resize(nspecies);
    for (auto &species : pol_species) {
      MPI_Bcast(&species.core_type,1,MPI_INT,0,world);
      MPI_Bcast(&species.shell_type,1,MPI_INT,0,world);
      MPI_Bcast(&species.count,1,MPI_INT,0,world);
      MPI_Bcast(&species.weight,1,MPI_DOUBLE,0,world);
    }
  }
  polarization[0] = polarization[1] = polarization[2] = 0.0;
}
