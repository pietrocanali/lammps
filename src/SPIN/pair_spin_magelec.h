/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(spin/magelec,PairSpinMagelec);
// clang-format on
#else

#ifndef LMP_PAIR_SPIN_MAGELEC_H
#define LMP_PAIR_SPIN_MAGELEC_H

#include "pair_spin.h"

#include <vector>

namespace LAMMPS_NS {

class PairSpinMagelec : public PairSpin {
 public:
  PairSpinMagelec(LAMMPS *lmp) :
      PairSpin(lmp), ME(nullptr), ME_mech(nullptr), v_mex(nullptr), v_mey(nullptr), v_mez(nullptr),
      cut_spin_magelec(nullptr), dynamic_pair(nullptr), gamma_me(nullptr), pol_list(nullptr),
      dynamic_configured(0), center_core_type(0), center_shell_type(0), pol_cutoff(0.0)
  {
    polarization[0] = polarization[1] = polarization[2] = 0.0;
    nextra = 3;
    pvector = new double[nextra]();
  }
  ~PairSpinMagelec() override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void *extract(const char *, int &) override;
  void init_style() override;
  void init_list(int, class NeighList *) override;

  void compute(int, int) override;
  void compute_single_pair(int, double *) override;

  void compute_magelec(int, int, double *, double *, double *);
  void compute_magelec_mech(int, int, double *, double *, double *);

  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;

  double cut_spin_magelec_global;    // global me cutoff

 protected:
  enum { FIXED_COUPLING, POLARIZATION_COUPLING };

  struct PolarizationSpecies {
    int core_type;
    int shell_type;
    int count;
    double weight;
  };

  struct PolarizationContribution {
    int center;
    int atom;
    double weighted_charge;
    double displacement[3];
  };

  double **ME, **ME_mech;              // magelec coeff in eV
  double **v_mex, **v_mey, **v_mez;    // magelec direction
  double **cut_spin_magelec;           // magelec cutoff distance
  int **dynamic_pair;                   // fixed or polarization coupling
  double **gamma_me;                    // dynamic gamma in eV/(C/m^2)

  class NeighList *pol_list;            // all-type polarization neighbor list
  int dynamic_configured;
  int center_core_type, center_shell_type;
  double pol_cutoff;
  double polarization[3];               // globally reduced, cached C/m^2
  std::vector<PolarizationSpecies> pol_species;
  std::vector<PolarizationContribution> pol_contributions;

  void allocate() override;
  void compute_polarization();
  void compute_dynamic_spin(int, int, const double *, double *, const double *);
  void compute_dynamic_mechanical(int, int, const double *, double, const double *,
                                  const double *, double *, double *, double &);
  void apply_polarization_forces(const double *, int);
};

}    // namespace LAMMPS_NS

#endif
#endif
