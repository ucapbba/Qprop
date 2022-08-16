#ifndef Header_HH
#define Header_HH

#include <iostream>
#include <complex>
#include <functional>
#include <string>
#include <time.h>

#include <bar.h>
#include <fluid.h>
#include <grid.h>
#include <hamop.h>
#include <wavefunction.h>
#include <powers.hh>
#include <parameter.hh>
#include <smallHelpers.hh>
#include <tsurffSpectrum.hh>
#include <gfunc.h>

// Functions determining the potentials are defined in potentials.hh
#include <potentials.hh>


void print_banner(long tsurff_version)
{
    fprintf(stdout, " ---------------------------------------------------\n");
    if (tsurff_version == 1)
        fprintf(stdout, "                 t-SURFF initialized              \n");
    if (tsurff_version == 2)
        fprintf(stdout, "                 i-SURFV initialized              \n");
    fprintf(stdout, " ---------------------------------------------------\n");
};

void ell_m_consistency(long ell, long m, grid g) {
    if (ell < labs(m)) {
        cerr << "|m| is greater than ell" << endl;
        exit(-1);
    };
    if (g.ngps_y() < ell + 1) {
        cerr << "grid too small for ell" << endl;
        exit(-1);
    };
};

#endif