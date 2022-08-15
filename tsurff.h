#include "Header.h"

void print_banner_t()
{
    fprintf(stdout, " ---------------------------------------------------\n");
    fprintf(stdout, "                 t-SURFF initialized              \n");
    fprintf(stdout, " ---------------------------------------------------\n");
};


void tsurff()
{
    clock_t t = clock();
    parameterListe para_ini("initial.param");
    parameterListe para_prop("propagate.param");
    parameterListe para_tsurff("tsurff.param");

    const long tsurff_version = para_tsurff.getLong("tsurff-version");
    const long     qprop_dim = para_ini.getLong("qprop-dim");


    print_banner_t();

    double omega = para_prop.getDouble("omega");
    double E_0 = para_prop.getDouble("max-electric-field");
    double phase = para_prop.getDouble("phase");
    double n_c = para_prop.getDouble("num-cycles");


    if ((omega == 0.0) && (E_0 != 0.0))
    {
        cout << " -----------------------------------------------------------\n";
        cout << " err: Zero frequency for a non-zero electric field detected.\n Exiting program...\n";
        cout << " -----------------------------------------------------------\n";
        exit(1);
    };

    double Ex, Ey, Ez, nx, ny, nz;
    // Only A_z in linear polarization case, renaming Ex1,2 -> Ez1,2
    if (qprop_dim == 34)
    {
        Ex = 0.0;
        Ey = 0.0;
        Ez = E_0;
        nx = 0.0;
        ny = 0.0;
        nz = n_c;

    }
    // Only A_x and A_y in XY-plane polarization case
    else if (qprop_dim == 44)
    {
        Ex = E_0;
        Ey = E_0;
        Ez = 0.0;
        nx = n_c;
        ny = n_c;
        nz = 0.0;
    };
    vecpot vecpot_x(omega, nx, Ex, phase);
    vecpot vecpot_y(omega, ny, Ey, phase + 0.5 * M_PI);
    vecpot vecpot_z(omega, nz, Ez, phase);

    tsurffSpectrum<vecpot, vecpot, vecpot> tsurff(para_ini, para_prop, para_tsurff, vecpot_x, vecpot_y, vecpot_z);
    tsurff.time_integration();
    // tsurff.print_int_dt_psi();
    // tsurff.print_wigner(0);
    // tsurff.print_Plms();
    // tsurff.print_bessel();
    tsurff.polar_spectrum();
    tsurff.print_partial_amplitudes();

    float sec = ((float)(clock() - t)) / CLOCKS_PER_SEC;
    cout << "Timer: " << long(sec) / 3600 << "h " << long(sec) % 3600 / 60 << "min " << (long(sec * 100) % 6000) / 100.0 << "sec\n";
    fprintf(stdout, " ---------------------------------------------------\n");
    if (tsurff_version == 1)
        fprintf(stdout, "                 t-SURFF finished              \n");
    if (tsurff_version == 2)
        fprintf(stdout, "                 i-SURFV finished              \n");
    fprintf(stdout, " ---------------------------------------------------\n");


}
