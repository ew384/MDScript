/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef DIHEDRAL_CLASS

DihedralStyle(trappea,DihedralTRAPPEA)
//DihedralStyle(opls,DihedralOPLS)

#else

#ifndef LMP_DIHEDRAL_TRAPPEA_H
#define LMP_DIHEDRAL_TRAPPEA_H

#include "stdio.h"
#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralTRAPPEA : public Dihedral {
 public:
  DihedralTRAPPEA(class LAMMPS *);
  virtual ~DihedralTRAPPEA();
  virtual void compute(int, int);
  void coeff(int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);

 protected:
  double *k0,*k1,*k2,*k3,*k4,*k5,*k6,*k7;

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

W: Dihedral problem: %d %ld %d %d %d %d

Conformation of the 4 listed dihedral atoms is extreme; you may want
to check your simulation geometry.

E: Incorrect args for dihedral coefficients

Self-explanatory.  Check the input script or data file.

*/
