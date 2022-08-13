#ifndef POWERS_HH
#define POWERS_HH

#include <complex>

typedef std::complex<double> cplxd;

/// calculate \f$x^2\f$
inline double pow2(double x) {return x*x;};

/// calculate \f$x^3\f$
inline double pow3(double x) {return x*x*x;};

/// calculate \f$i^l\f$ with integer \f$l\f$
inline cplxd pow_i(long l) {
	cplxd p_i;
  if (l%4==0) p_i = cplxd(1.0, 0.0);
  if (l%4==1) p_i = cplxd(0.0, 1.0);
  if (l%4==2) p_i = cplxd(-1.0, 0.0);
  if (l%4==3) p_i = cplxd(0.0, -1.0);
  if (l%4==-0) p_i = cplxd(1.0, 0.0);
  if (l%4==-1) p_i = cplxd(0.0, -1.0);
  if (l%4==-2) p_i = cplxd(-1.0, 0.0);
  if (l%4==-3) p_i = cplxd(0.0, 1.0);
  return p_i;
};

/// calculate \f${(-i)}^l\f$ with integer \f$l\f$
inline cplxd pow_neg_i(long l) {
  cplxd p_i;
  if (l%4==0) p_i = cplxd(1.0, 0.0); //1
  if (l%4==1) p_i = cplxd(0.0, -1.0); //-i
  if (l%4==2) p_i = cplxd(-1.0, 0.0); //-1
  if (l%4==3) p_i = cplxd(0.0, 1.0); //i
  if (l%4==-0) p_i = cplxd(1.0, 0.0);
  if (l%4==-1) p_i = cplxd(0.0, 1.0);
  if (l%4==-2) p_i = cplxd(-1.0, 0.0);
  if (l%4==-3) p_i = cplxd(0.0, -1.0);
  return p_i;
};

/// calculate \f${(-1)}^p\f$ with integer \f$p\f$
inline double pow_neg_1(long p) {
  return (labs(p%2)==1)?-1.0:1.0;
};

#endif
