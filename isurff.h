#pragma once
void print_banner_i() {
    fprintf(stdout, " -----------------------------------------------------------\n");
    fprintf(stdout, "  Applying the Green's function 1/(H-Ek) to Psi and dPsi/dr\n");
    fprintf(stdout, "  (C) Copyright by Tulsky V A  and Bauer D, Rostock  (2019)\n");
    fprintf(stdout, " -----------------------------------------------------------\n");
};

void isurff()
{
    clock_t t = clock();
    print_banner_i();
    // Set verbosity
    const int iv = 1; // 1 leads to common output, 0 is silent regime

    // Get parameters
    parameterListe para_ini("initial.param");
    parameterListe para_prop("propagate.param");
    parameterListe para_tsurff("tsurff.param");

    const long   delta_k_scheme = para_tsurff.getLong("delta-k-scheme");
    const long   num_energy = para_tsurff.getLong("num-k-surff");
    const double k_max_surff = para_tsurff.getDouble("k-max-surff");
    const double energy_max = 0.5 * k_max_surff * k_max_surff;
    const double delta_r = para_ini.getDouble("delta-r");
    const long   qprop_dim = para_ini.getLong("qprop-dim");
    const long   ell_grid_size = para_prop.getLong("ell-grid-size");


    // This is the grid for the wavefunction that is to be analyzed
    const double R_surff = para_prop.getDouble("R-tsurff");
    const double omega = para_prop.getDouble("omega");
    const double E_0 = para_prop.getDouble("max-electric-field");

    if ((omega == 0.0) && (E_0 != 0.0))
    {
        cout << " -----------------------------------------------------------\n";
        cout << " err: Zero frequency for a non-zero electric field detected.\n Exiting program...\n";
        cout << " -----------------------------------------------------------\n";
        exit(1);
    };
    const double quiver_amplitude = E_0 / (omega * omega);
    const double w_imag_fact = para_tsurff.getDouble("isurfv-imag-width-factor");
    const double imag_potential_width = w_imag_fact * para_prop.getDouble("imag-width");
    const double grid_size = imag_potential_width + R_surff + quiver_amplitude;
    grid g;
    g.set_dim(qprop_dim);
    g.set_ngps(long(grid_size / delta_r), ell_grid_size, 1);
    g.set_delt(delta_r);
    g.set_offs(0, 0, 0);

    const double imag_potential_width_small = para_prop.getDouble("imag-width");
    const double grid_size_small = imag_potential_width_small + R_surff + quiver_amplitude;
    grid g_small;
    g_small.set_dim(qprop_dim);
    g_small.set_ngps(long(grid_size_small / delta_r), ell_grid_size, 1);
    g_small.set_delt(delta_r);
    g_small.set_offs(0, 0, 0);

    wavefunction psi_at_RI, d_psi_dr_at_RI;
    if (qprop_dim == 34)
    {
        psi_at_RI.init(g.ngps_y());
        d_psi_dr_at_RI.init(g.ngps_y());
    }
    else if (qprop_dim == 44)
    {
        psi_at_RI.init(g.ngps_y() * g.ngps_y());
        d_psi_dr_at_RI.init(g.ngps_y() * g.ngps_y());
    }
    else
    {
        fprintf(stderr, "err: wrong qprop-dim!=34 and qprop-dim!=44!\n"); exit(12);
    };

    fluid  V_ee_0;
    V_ee_0.init(g.ngps_x());

    //
    // The Hamiltonian
    //
    const double nuclear_charge = para_ini.getDouble("nuclear-charge");
    const double Rco = para_ini.getDouble("pot-cutoff");
    scalarpot scalarpotx(nuclear_charge, Rco);
    hamop hamilton;
    const double im_ampl = para_prop.getDouble("imag-ampl");
    const long imag_width_ngps = long(imag_potential_width / delta_r);
    imagpot imaginarypot(imag_width_ngps, im_ampl);
    hamilton.init(g, always_zero2, always_zero2, always_zero2, scalarpotx, always_zero5, always_zero5, imaginarypot, always_zero2);


    //
    // Specify where the wavefunction is located
    //
    string str_fname_wf = string("./dat/real_prop_wf.dat");
    FILE* file_wf = fopen_with_check(str_fname_wf, "r");
    // The wavefunction arrays
    wavefunction  wf, wf_load;
    wf.init(g.size());
    wf_load.init(g_small.size());

    wf_load.init(g_small, file_wf, 0, 0);
    fclose(file_wf);


    wf.regrid(g, g_small, wf_load);
    // Display the grid parameters
    if (iv == 1)
    {
        fprintf(stdout, "Grid: \n");
        fprintf(stdout, "g.ngps_x()=%ld\n", g.ngps_x());
        fprintf(stdout, "g.ngps_y()=%ld\n", g.ngps_y());
        fprintf(stdout, "g.ngps_z()=%ld\n", g.ngps_z());
        fprintf(stdout, "g.dimens()=%d\n", g.dimens());
        fprintf(stdout, "g.delt_x()=%20.15le\n", g.delt_x());
        fprintf(stdout, "nuclear_charge   =%20.15le\n", nuclear_charge);
        fprintf(stdout, "str_fname_wf=%s\n", str_fname_wf.c_str());
        fflush(stdout);
    };

    wavefunction staticpot;
    staticpot.init(g.size());
    staticpot.calculate_staticpot(g, hamilton);

    // Output files
    ofstream file_res(string("./dat/isurfv.dat"));
    file_res.precision(17);


    // Energy loop
    for (long i_energy = 0; i_energy < num_energy; i_energy++)
    {
        if (i_energy % 100 == 0)
        {
            cout << i_energy << " energies of " << num_energy << " complete." << endl;
        };
        double energy, delta_energy;
        if (delta_k_scheme == 1)
        {
            delta_energy = energy_max / (double(num_energy) * double(num_energy));
            energy = double(i_energy) * double(i_energy) * delta_energy;
        }
        else if (delta_k_scheme == 2)
        {
            delta_energy = energy_max / double(num_energy);
            energy = i_energy * delta_energy;
        };

        gfunc(psi_at_RI, d_psi_dr_at_RI, energy, staticpot, V_ee_0, nuclear_charge, g, wf, R_surff);

        // write the partial results (i.e., for individual l (and m)) to file
        file_res << energy << " " << sqrt(2.0 * energy) << " ";
        for (long l_index = 0; l_index < g.ngps_y(); l_index++)
        {
            if (qprop_dim == 34) // Linear polarization and fixed value of m
            {
                file_res << real(psi_at_RI[l_index]) << " " << imag(psi_at_RI[l_index]) << " " << real(d_psi_dr_at_RI[l_index]) << " " << imag(d_psi_dr_at_RI[l_index]) << " ";
            }
            else if (qprop_dim == 44) // Polarization in the XY plane
            {
                const long m_limit = l_index;
                for (long m_index = -m_limit; m_index <= m_limit; m_index++)
                {
                    long lmindex = l_index * (l_index + 1) + m_index;
                    file_res << real(psi_at_RI[lmindex]) << " " << imag(psi_at_RI[lmindex]) << " " << real(d_psi_dr_at_RI[lmindex]) << " " << imag(d_psi_dr_at_RI[lmindex]) << " ";
                };
            }
            else
            {
                cerr << "Unknown propagation mode " << g.dimens() << endl;
            };
        };
        file_res << endl; // Break line for every new energy
    }; // End of energy loop
    float sec = ((float)(clock() - t)) / CLOCKS_PER_SEC;
    cout << "Timer: " << long(sec) / 3600 << "h " << long(sec) % 3600 / 60 << "min " << (long(sec * 100) % 6000) / 100.0 << "sec\n";
    fprintf(stdout, " -----------------------------------------------------------\n");
    fprintf(stdout, "            The Green's function trick performed            \n");
    fprintf(stdout, " -----------------------------------------------------------\n");
}
