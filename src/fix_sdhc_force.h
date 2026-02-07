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

#ifdef FIX_CLASS
// clang-format off
FixStyle(sdhc/force,FixSDHCForce);
// clang-format on
#else

#ifndef LMP_FIX_SDHC_FORCE_H
#define LMP_FIX_SDHC_FORCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSDHCForce : public Fix {
 public:
  FixSDHCForce(class LAMMPS *, int, char **);
  ~FixSDHCForce() override;
  int setmask() override;
  void init() override;
  void pre_force(int) override;
  void pre_force_respa(int, int, int) override;
  void min_pre_force(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double compute_scalar() override;
  double memory_usage() override;
  void tally_pair(int, int, const double *, const double *);
  void tally_triplet(int, int, int, const double *, const double *, const double *);

 private:
  enum { MODE_NONE = 0, MODE_GROUP2 = 1, MODE_REGION = 2 };
  enum { SCALAR_ALL = 0, SCALAR_REGION2 = 1, SCALAR_HALF = 2 };
  int nlevels_respa;
  int nmax;
  double **f12;
  double *h;
  int mode;
  int scalar_mode;
  int groupbit2;
  char *idregion;
  class Region *region;
};

}    // namespace LAMMPS_NS

#endif
#endif
