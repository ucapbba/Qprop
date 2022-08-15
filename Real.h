#include "Header.h"

void print_banner()
{
    fprintf(stdout, " --------------------------------------------------------\n");
    fprintf(stdout, "                  Real time propagation                  \n");
    fprintf(stdout, "  (C) Copyright by Bauer D and Koval P, Heidelberg (2005)\n");
    fprintf(stdout, " --------------------------------------------------------\n");
};

void RealProp()
{
    clock_t t = clock();
    grid g_prop, g_load;
    wavefunction staticpot, wf, wf_load;

    print_banner();

    int iv = 1; // verbosity of stdout

    // get parameters
    parameterListe para_ini("initial.param");
    parameterListe para_prop("propagate.param");
    parameterListe para_tsurff("tsurff.param");

    const long gen_tsurff_data = para_prop.getLong("generate-tsurff-data");
    const long gen_hhg_data = para_prop.getLong("generate-hhg-data");

    // Define if wavefunction is to be stored at intermediate steps
    long wf_saving_interval = para_prop.getLong("wf-saving-interval");
    // Define if several observables are to be stored at intermediate steps
    long obser_saving_interval = para_prop.getLong("obser-saving-interval");

    //  Input (result from imaginary propagation generated file or other)
    string str_fname_wf_ini = para_prop.getString("init-wf-filename");
    FILE* file_wf_ini = fopen_with_check(str_fname_wf_ini, "r");

    // *** Declare the grid for load ***
    const long qprop_dim = para_ini.getLong("qprop-dim");
    g_load.set_dim(qprop_dim);
    const double delta_r = para_ini.getDouble("delta-r");
    g_load.set_ngps(long(para_ini.getDouble("radial-grid-size") / delta_r), para_ini.getLong("ell-grid-size"), 1);
    g_load.set_delt(delta_r, 0.0, 0.0);
    g_load.set_offs(0, 0, 0);

    const int my_m_quantum_num = para_ini.getLong("initial-m");
    const double nuclear_charge = para_ini.getDouble("nuclear-charge");
    const double Rco = para_ini.getDouble("pot-cutoff");

    // Everything related to the laser
    const double omega = para_prop.getDouble("omega");
    const double n_c = para_prop.getDouble("num-cycles");
    const double E_0 = para_prop.getDouble("max-electric-field");
    const double phase = para_prop.getDouble("phase");

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
        Ex = 0.0;

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
    double U_p, pulse_duration;
    U_p = (qprop_dim == 34) ? vecpot_z.get_Up() : vecpot_x.get_Up();
    pulse_duration = (qprop_dim == 34) ? vecpot_z.get_duration() : vecpot_x.get_duration();
    double I_p = nuclear_charge * nuclear_charge / 2.0; // Warn: hydrogenic ionization energy!
    const double gamma = sqrt(I_p / 2.0 / U_p);

    // How long do the slowest electrons have time to reach the t-SURFF boundary

    const double real_timestep = para_prop.getDouble("delta-t");
    int tsurff_method;
    const double additional_time = para_prop.getDouble("additional-time");
    double duration;
    if (gen_tsurff_data == 0)
    {
        duration = pulse_duration + additional_time;
        cout << "  No t-SURFF/i-SURFF data will be generated\n";
        cout << " --------------------------------------------------------\n";
    }
    else
    {
        tsurff_method = para_tsurff.getLong("tsurff-version");
        if (tsurff_method == 1)
        {
            const double time_surff = para_prop.getDouble("R-tsurff") / para_tsurff.getDouble("k-min-tsurff");
            duration = pulse_duration + additional_time + time_surff;
            cout << "  Output data suitable for t-SURFF will be generated\n";
            cout << " --------------------------------------------------------\n";
        }
        else if (tsurff_method == 2)
        {
            duration = pulse_duration + additional_time;
            cout << "  Output data suitable for i-SURFV will be generated\n";
            cout << " --------------------------------------------------------\n";
        };
    };
    duration = double(long(duration / real_timestep) + 1) * real_timestep;


    // *** Declare the grid for propagation ***
    g_prop.set_dim(qprop_dim);
    const double quiver_amplitude = E_0 / (omega * omega);
    const double grid_size = para_prop.getDouble("imag-width") + para_prop.getDouble("R-tsurff") + quiver_amplitude;
    g_prop.set_ngps(long(grid_size / delta_r), para_prop.getLong("ell-grid-size"), 1);
    g_prop.set_delt(delta_r);
    g_prop.set_offs(0, 0, 0);


    // Output that will be created by this program
    string common_prefix("./dat/real_prop");
    string str_fname_logfi = common_prefix + string(".log");
    FILE* file_logfi = fopen_with_check(str_fname_logfi, "w");
    string str_fname_yield = common_prefix + string("_yield.dat");
    FILE* file_yield = fopen_with_check(str_fname_yield, "w");
    string str_fname_obser = common_prefix + string("_obser.dat");
    FILE* file_obser = fopen_with_check(str_fname_obser, "w");

    if (iv != 0)
    {
        fprintf(stdout, "%s will be (re)written.\n", str_fname_logfi.c_str());
        fprintf(stdout, "%s will be (re)written.\n", str_fname_yield.c_str());
        fprintf(stdout, "%s will be (re)written.\n", str_fname_obser.c_str());
    };

    // Create an instance of the class for doing the tsurff related work
    tsurffSaveWF       tsurff_save_wf(para_ini, para_prop, para_tsurff, g_prop);

    // The absorbing imaginary potential
    const double im_ampl = para_prop.getDouble("imag-ampl");
    const long imag_width_ngps = long(para_prop.getDouble("imag-width") / delta_r);
    imagpot imaginarypot(imag_width_ngps, im_ampl);
    // Set the binding potential and the hamiltonian
    scalarpot scalarpotx(nuclear_charge, Rco);
    hamop hamilton;
    hamilton.init(g_prop, vecpot_x, vecpot_y, vecpot_z, scalarpotx, always_zero5, always_zero5, imaginarypot, always_zero2);

    // This is the linear and constant part of the Hamiltonian
    staticpot.init(g_prop.size());
    staticpot.calculate_staticpot(g_prop, hamilton);

    // *** Wavefunction array 
    wf.init(g_prop.size());
    wf_load.init(g_load.size());

    cout << "Imag time prop grid with size: g_load.size() = " << g_load.size() << endl;
    cout << "  is fit into" << endl;
    cout << "Real time prop grid with size: g_prop.size() = " << g_prop.size() << endl;


    // *** Wavefunction initialization ***
    wf_load.init(g_load, file_wf_ini, 0, iv);
    wf.regrid(g_prop, g_load, wf_load);
    fclose(file_wf_ini);

    string str_fname_wf = common_prefix + string("_wf.dat");

    fprintf(stdout, "Norm of orbital: %le\n", wf.norm(g_prop));

    long lno_of_ts = long(duration * 1.0 / real_timestep);// + 1;

    //
    // Write to log file
    //
    fprintf(file_logfi, "Real-time propagation\n");
    fprintf(file_logfi, "Grid: \n");
    fprintf(file_logfi, "  g_prop.dimens() = %d\n\n", g_prop.dimens());
    fprintf(file_logfi, "  g_prop.ngps_x() = %ld\n", g_prop.ngps_x());
    fprintf(file_logfi, "  g_prop.ngps_y() = %ld\n", g_prop.ngps_y());
    fprintf(file_logfi, "  g_prop.ngps_z() = %ld\n", g_prop.ngps_z());
    fprintf(file_logfi, "  g_prop.delt_x() = %15.10le\n", g_prop.delt_x());

    fprintf(file_logfi, "  real_timestep     = %15.10le\n", real_timestep);
    fprintf(file_logfi, "  lno_of_ts         = %ld\n", lno_of_ts);
    fprintf(file_logfi, "  nuclear_charge    = %15.10le\n", nuclear_charge);
    fprintf(file_logfi, "  str_fname_wf_ini = %s\n", str_fname_wf_ini.c_str());
    fprintf(file_logfi, "  str_fname_obser  = %s\n", str_fname_obser.c_str());
    if (gen_tsurff_data == 1)
        fprintf(file_logfi, "  str_fname_wf = %s\n", str_fname_wf.c_str());

    fprintf(file_logfi, "Laser: \n");
    fprintf(file_logfi, "  E_0         = %15.10le\n", E_0);
    fprintf(file_logfi, "  omega      = %15.10le\n", omega);
    fprintf(file_logfi, "  n_cycles      = %15.10le\n", n_c);
    fprintf(file_logfi, "Total propagation time duration = %15.10le a.u.\n", duration);

    fprintf(file_logfi, "Keldysh gamma    = %15.10le ", gamma);
    if (gamma > 1.0)
    {
        fprintf(file_logfi, "(multi-photon regime)\n");
    }
    else
    {
        fprintf(file_logfi, "(tunneling regime)\n");
    };
    fflush(file_logfi);
    fclose(file_logfi);

    long ldumpwidth = obser_saving_interval;  // How often to save observables &/or HHG data
    int me = 0; // dummy here
    // Write vector potential to file
    string str_fname_vpot;
    if (qprop_dim == 34)
    {
        str_fname_vpot = common_prefix + string("_vpot_z.dat");
        FILE* file_vpot = fopen_with_check(str_fname_vpot, "w");
        for (long ts = 0; ts < lno_of_ts; ts++)
        {
            const double time = real_timestep * double(ts);
            if (ts % ldumpwidth == 0)
            {
                fprintf(file_vpot, "%15.10le %15.10le\n", time, vecpot_z(time, me));
            };
        };
        fclose(file_vpot);
    }
    else if (qprop_dim == 44)
    {
        str_fname_vpot = common_prefix + string("_vpot_xy.dat");
        FILE* file_vpot = fopen_with_check(str_fname_vpot, "w");
        for (long ts = 0; ts < lno_of_ts; ts++)
        {
            const double time = real_timestep * double(ts);
            if (ts % ldumpwidth == 0)
            {
                fprintf(file_vpot, "%15.10le %15.10le %15.10le\n", time, vecpot_x(time, me), vecpot_y(time, me));
            };
        };
        fclose(file_vpot);
    };

    // *********************************** //
    // ****** Real time propagation ****** //
    // *********************************** //
    cplxd timestep = cplxd(real_timestep, 0.0);
    cplxd P;
    double N;
    for (long ts = 0; ts < lno_of_ts; ts++)
    {
        const double time = real_timestep * double(ts);

        if (gen_tsurff_data == 1)
        {
            // Save the orbitals \varphi_{\ell}(\RI) and the derivative \partial_r\varphi_{\ell}(r)|_{r=\RI}
            tsurff_save_wf(wf);
        };
        if (ts % ldumpwidth == 0)
        {
            // Calculate total energy, projection onto initial state, norm, and <z>. Also <d2z/dt2> or <d2(x+iy)/dt2> if HHG spectrum required
            double E_tot = real(wf.energy(0.0, g_prop, hamilton, me, staticpot, nuclear_charge));
            P = wf.project(g_prop, g_load, wf_load, 0);
            N = wf.norm(g_prop);
            double z_expect = real(wf.expect_z(g_prop, my_m_quantum_num));
            if (gen_hhg_data == 1)
            {
                cplxd accel;
                if (qprop_dim == 34)
                {
                    accel = real(wf.accel_z(g_prop, my_m_quantum_num, hamilton, time, me, real_timestep)); //d^2<z>/dt^2
                }
                else if (qprop_dim == 44)
                {
                    accel = wf.accel_cycl_pol_plus(g_prop, hamilton, time, me, real_timestep); // d^2<x+iy>/dt^2
                };
                cplxd imagi(0.0, 1.0);
                fprintf(file_obser, "%15.17le %15.10le %15.10le %15.10le %15.10le %15.10le %15.10le\n", time, E_tot, real(conj(P) * P), N, z_expect, real(accel), imag(accel));
            }
            else
            {
                fprintf(file_obser, "%15.10le %15.10le %15.10le %15.10le %15.10le\n", time, E_tot, real(conj(P) * P), N, z_expect);
            };
        };
        //
        // Propagate one step forward in (real) time.
        //
        wf.propagate(timestep, time, g_prop, hamilton, me, staticpot, my_m_quantum_num, nuclear_charge);
        //if (ts % (ldumpwidth * 10) == 0)
        //{
        cout << "\r" << "timestep " << ts << " of " << lno_of_ts << ", Norm of wave function: " << N << "                 ";
        //};

        // Save the full wavefunction after each %(...) timesteps. 
        // Don't forget to create a folder in which to save this data.
        if (wf_saving_interval != 0)
        {
            if (ts % wf_saving_interval == 0)
            {
                string fname_wf = string("./wf/real_prop_wf_") + to_string(ts) + string(".dat");
                FILE* file_wf_t = fopen_with_check(fname_wf, "w");
                wf.dump_to_file_sh(g_prop, file_wf_t, 1); // wf at timestep ts is saved
                fclose(file_wf_t);
            };
        };

    }; // End of real-time-propagation loop

    fclose(file_obser);

    double yield_N = (1.0 - N);
    double yield_P = (1.0 - real(conj(P) * P));
    fprintf(file_yield, "%15.10le %15.10le\n", yield_N, yield_P);
    fclose(file_yield);

    if (gen_tsurff_data == 1)
    {
        FILE* file_wf;
        file_wf = fopen_with_check(str_fname_wf, "w");
        wf.dump_to_file_sh(g_prop, file_wf, 1, iv); // Final wf is saved
        fclose(file_wf);
    };

    if (iv != 0)
    {
        fprintf(stdout, "%s was read.\n", str_fname_wf_ini.c_str());
        if (gen_tsurff_data == 1)
            fprintf(stdout, "%s is written.\n", str_fname_wf.c_str());
        fprintf(stdout, "%s is written.\n", str_fname_obser.c_str());
        fprintf(stdout, "%s is written.\n", str_fname_vpot.c_str());
        fprintf(stdout, "%s is written.\n", str_fname_logfi.c_str());
        fprintf(stdout, "%s is written.\n", str_fname_yield.c_str());
        //         fprintf(stdout, "Hasta la vista...\n");
    };
    float sec = ((float)(clock() - t)) / CLOCKS_PER_SEC;
    cout << "Timer: " << long(sec) / 3600 << "h " << long(sec) % 3600 / 60 << "min " << (long(sec * 100) % 6000) / 100.0 << "sec\n";
    cout << " -------------------------------------------------------\n";
    cout << "             Real time propagation finished             \n";
}

