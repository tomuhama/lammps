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

#include "fix_sdhc_force.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "region.h"
#include "respa.h"
#include "update.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSDHCForce::FixSDHCForce(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), f12(nullptr), h(nullptr), mode(MODE_NONE),
    scalar_mode(SCALAR_REGION2), groupbit2(0), idregion(nullptr), region(nullptr)
{
  if (narg < 3) error->all(FLERR, "Illegal fix sdhc/force command");

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "group2") == 0) {
      if (mode != MODE_NONE) error->all(FLERR, "Illegal fix sdhc/force command");
      if (iarg + 1 >= narg) error->all(FLERR, "Illegal fix sdhc/force command");
      int igroup2 = group->find(arg[iarg+1]);
      if (igroup2 == -1) error->all(FLERR, "Group ID {} does not exist", arg[iarg+1]);
      groupbit2 = group->bitmask[igroup2];
      mode = MODE_GROUP2;
      iarg += 2;
      continue;
    }
    if (strcmp(arg[iarg], "region") == 0) {
      if (mode != MODE_NONE) error->all(FLERR, "Illegal fix sdhc/force command");
      if (iarg + 1 >= narg) error->all(FLERR, "Illegal fix sdhc/force command");
      idregion = utils::strdup(arg[iarg+1]);
      region = domain->get_region_by_id(idregion);
      if (!region)
        error->all(FLERR, "Region {} for fix sdhc/force does not exist", idregion);
      mode = MODE_REGION;
      iarg += 2;
      continue;
    }
    if (strcmp(arg[iarg], "scalar") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR, "Illegal fix sdhc/force command");
      if (strcmp(arg[iarg+1], "all") == 0) {
        scalar_mode = SCALAR_ALL;
      } else if (strcmp(arg[iarg+1], "region2") == 0) {
        scalar_mode = SCALAR_REGION2;
      } else if (strcmp(arg[iarg+1], "half") == 0) {
        scalar_mode = SCALAR_HALF;
      } else {
        error->all(FLERR, "Illegal fix sdhc/force command");
      }
      iarg += 2;
      continue;
    }
    error->all(FLERR, "Illegal fix sdhc/force command");
  }

  peratom_flag = 1;
  size_peratom_cols = 3;
  peratom_freq = 1;
  scalar_flag = 1;
  extscalar = 1;
  global_freq = 1;
  dynamic_group_allow = 1;
  comm_reverse = 3;

  nmax = atom->nmax;
  memory->create(f12, nmax, 3, "sdhc/force:f12");
  memory->create(h, nmax, "sdhc/force:h");
  array_atom = f12;

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    f12[i][0] = f12[i][1] = f12[i][2] = 0.0;
    h[i] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

FixSDHCForce::~FixSDHCForce()
{
  memory->destroy(f12);
  memory->destroy(h);
  delete[] idregion;
}

/* ---------------------------------------------------------------------- */

int FixSDHCForce::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= PRE_FORCE_RESPA;
  mask |= MIN_PRE_FORCE;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSDHCForce::init()
{
  if (utils::strmatch(update->integrate_style, "^respa"))
    nlevels_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels;

  if (mode == MODE_REGION) {
    region = domain->get_region_by_id(idregion);
    if (!region)
      error->all(FLERR, "Region {} for fix sdhc/force does not exist", idregion);
    region->init();
  }
}

/* ---------------------------------------------------------------------- */

void FixSDHCForce::pre_force(int /*vflag*/)
{
  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->destroy(f12);
    memory->destroy(h);
    memory->create(f12, nmax, 3, "sdhc/force:f12");
    memory->create(h, nmax, "sdhc/force:h");
    array_atom = f12;
  }

  if (mode == MODE_REGION && region) region->prematch();

  double **x = atom->x;
  int *mask = atom->mask;
  int nall = atom->nlocal + atom->nghost;
  for (int i = 0; i < nall; i++) {
    f12[i][0] = f12[i][1] = f12[i][2] = 0.0;
    if (mode == MODE_GROUP2)
      h[i] = (mask[i] & groupbit2) ? 1.0 : 0.0;
    else if (mode == MODE_REGION)
      h[i] = region->match(x[i][0], x[i][1], x[i][2]) ? 1.0 : 0.0;
    else
      h[i] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

void FixSDHCForce::pre_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa - 1) pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSDHCForce::min_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSDHCForce::post_force(int /*vflag*/)
{
  if (force->newton_pair) comm->reverse_comm(this);
}

/* ---------------------------------------------------------------------- */

void FixSDHCForce::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa - 1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSDHCForce::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

int FixSDHCForce::pack_reverse_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++) {
    buf[m++] = f12[i][0];
    buf[m++] = f12[i][1];
    buf[m++] = f12[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixSDHCForce::unpack_reverse_comm(int n, int *list, double *buf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    const int j = list[i];
    f12[j][0] += buf[m++];
    f12[j][1] += buf[m++];
    f12[j][2] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

double FixSDHCForce::compute_scalar()
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double qlocal = 0.0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (scalar_mode == SCALAR_REGION2 && h[i] < 0.5) continue;
      qlocal += f12[i][0] * v[i][0] + f12[i][1] * v[i][1] + f12[i][2] * v[i][2];
    }
  }

  double qglobal = 0.0;
  MPI_Allreduce(&qlocal, &qglobal, 1, MPI_DOUBLE, MPI_SUM, world);
  if (scalar_mode == SCALAR_HALF) qglobal *= 0.5;
  return qglobal;
}

/* ---------------------------------------------------------------------- */

void FixSDHCForce::tally_pair(int i, int j, const double *dFi, const double *dFj)
{
  const double hi = h[i];
  const double hj = h[j];
  const double h0 = 0.5 * (hi + hj);
  const double ci = hi - h0;
  const double cj = hj - h0;

  f12[i][0] += ci * dFi[0];
  f12[i][1] += ci * dFi[1];
  f12[i][2] += ci * dFi[2];
  f12[j][0] += cj * dFj[0];
  f12[j][1] += cj * dFj[1];
  f12[j][2] += cj * dFj[2];
}

/* ---------------------------------------------------------------------- */

void FixSDHCForce::tally_triplet(int i, int j, int k, const double *dFi, const double *dFj,
                                 const double *dFk)
{
  const double hi = h[i];
  const double hj = h[j];
  const double hk = h[k];
  const double h0 = (hi + hj + hk) / 3.0;
  const double ci = hi - h0;
  const double cj = hj - h0;
  const double ck = hk - h0;

  f12[i][0] += ci * dFi[0];
  f12[i][1] += ci * dFi[1];
  f12[i][2] += ci * dFi[2];
  f12[j][0] += cj * dFj[0];
  f12[j][1] += cj * dFj[1];
  f12[j][2] += cj * dFj[2];
  f12[k][0] += ck * dFk[0];
  f12[k][1] += ck * dFk[1];
  f12[k][2] += ck * dFk[2];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixSDHCForce::memory_usage()
{
  double bytes = (double) atom->nmax * 4 * sizeof(double);
  return bytes;
}
