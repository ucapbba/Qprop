#ifndef gfunc
#define gfunc gfunc

#include<wavefunction.h>
#include<grid.h>
#include<fluid.h>

class grid;
class hamop;
class wavefunction;

//
// Applies a Green's function operator 1/(H-Ek) to the wavefunction
//
void gfunc(wavefunction &psi_at_RI, wavefunction &d_psi_dr_at_RI,
  double energy, wavefunction staticpot, fluid hartree_potential_zero, double nuclear_charge, grid g, wavefunction wf, double R_surff);

#endif
