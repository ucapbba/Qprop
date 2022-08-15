#ifndef _potentials_H
#define _potentials_H
// Definitions of potentials
#define M_PI       3.14159265358979323846   // pi
// Definitions of potentials

double always_zero2(double t, int me) {
    return 0;
};

double always_zero5(double x, double y, double z, double t, int me) {
    return 0;
};

// The vector potential with sine squared envelope
class vecpot
{
    double omega;
    double  E_0;
    double   n_c;
    double    phi_cep;
    double       ww;
    double  duration;
public:
    vecpot(double om, double n_cyc, double E_max, double cep) : omega(om), n_c(n_cyc), E_0(E_max), phi_cep(cep)
    {
        duration = n_c * 2 * M_PI / omega; // Total duration of the (overlaping) pulses
        // Angular frequency of the envelopes
        ww = omega / (2.0 * n_c);
    };
    double operator()(double time, int me) const
    {
        double result = 0.0;
        if (time < duration)
        {
            // Here's the shape of the laser pulse
            result += E_0 / omega * pow2(sin(ww * time)) * sin(omega * time + phi_cep);
        };
        return result;
    };
    double get_duration()
    {
        return duration;
    };
    double get_Up()
    {
        return 0.25 * (E_0 * E_0) / (omega * omega);
    };
    double integral(double Time)
    {
        const double pulse_dur = n_c * 2.0 * M_PI / omega;
        // const double pulse_dur = (pulse_dur1<pulse_dur2)?pulse_dur2:pulse_dur1; // Total duration of the (overlaping) pulses
        const double excursion_ampl = E_0 / (omega * omega);
        // This is the time integral of A(t), i.e., the free-electron excursion alpha(t)
        double result = 0.0;
        if (Time < pulse_dur)
        {
            result += 0.5 * excursion_ampl * (
                cos(phi_cep) * (1.0 - 0.5 * n_c / (n_c + 1.0) - 0.5 * n_c / (n_c - 1.0))
                - cos(omega * Time + phi_cep)
                + cos(omega * Time + omega * Time / n_c + phi_cep) * 0.5 * n_c / (n_c + 1.0)
                + cos(omega * Time - omega * Time / n_c + phi_cep) * 0.5 * n_c / (n_c - 1.0));
        };
        return result;
    };
};

class scalarpot
{
    double nuclear_charge;
    double R_co;
public:
    scalarpot(double charge, double co) : nuclear_charge(charge), R_co(co) {
        //
    };
    double operator()(double x, double y, double z, double time, int me) const {
        // Scrinzi potential
        // double result=(x<R_co)?nuclear_charge*(-1.0/x-pow2(x)/(2*pow3(R_co))+3.0/(2.0*R_co)):0.0;
        // Simple Volker potential; first -1/r then linear then zero
        const double m = 1.0 / pow2(R_co);
        double result = (x < R_co) ? -1.0 / x : ((x < 2 * R_co) ? -1.0 / R_co + m * (x - R_co) : 0.0);
        return nuclear_charge * result;
    };
    double get_nuclear_charge() { return nuclear_charge; };
};

class imagpot {
    long imag_potential_width;
    double ampl_im;  // Amplitude of imaginary absorbing potential  <--------------  100.0 imag pot on,  0.0 off
public:
    imagpot(long width, double ampl = 100.0) : ampl_im(ampl), imag_potential_width(width) {
        //
    };
    double operator()(long xindex, long yindex, long zindex, double time, grid g) {
        if (ampl_im > 1.0) {
            const long imag_start = g.ngps_x() - imag_potential_width;
            if (xindex < imag_start)
                return 0;
            else {
                const double r = double(xindex - imag_start) / double(imag_potential_width);
                return ampl_im * pow2(pow2(pow2(pow2(r))));
            };
        }
        else {
            return 0.0;
        };
    };
};
#endif