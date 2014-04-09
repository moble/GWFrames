// Copyright (c) 2014, Michael Boyle
// See LICENSE file for details

#include <iostream>
#include <cmath>
#include "Utilities.hpp"
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include "Quaternions.hpp"
#include "Errors.hpp"
using GWFrames::Matrix;
using GWFrames::MatrixC;
using Quaternions::Quaternion;
using std::vector;
using std::complex;
using std::cerr;
using std::endl;


#define Utilities_Epsilon 1.0e-14

const std::complex<double> ComplexI(0.0, 1.0);


double GWFrames::abs(const std::vector<double>& v) {
  double n=0.0;
  for(unsigned int i=0; i<v.size(); ++i) {
    n += v[i]*v[i];
  }
  return std::sqrt(n);
}

std::vector<double> GWFrames::abs(const std::vector<std::vector<double> >& v) {
  const unsigned int vsize = v.size();
  vector<double> n(vsize);
  for(unsigned int i=0; i<vsize; ++i) {
    n[i] = GWFrames::abs(v[i]);
  }
  return n;
}

std::vector<double> GWFrames::operator+(const std::vector<double>& a, const std::vector<double>& b) {
  const unsigned int size = a.size();
  // if(b.size() != size) {
  //   cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": a.size()=" << a.size() << "; b.size()=" << b.size() << endl;
  //   throw(GWFrames_VectorSizeMismatch);
  // }
  vector<double> c(a);
  for(unsigned int i=0; i<size; ++i) {
    c[i] += b[i];
  }
  return c;
}

std::vector<double> GWFrames::operator-(const std::vector<double>& a, const std::vector<double>& b) {
  const unsigned int size = a.size();
  // if(b.size() != size) {
  //   cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": a.size()=" << a.size() << "; b.size()=" << b.size() << endl;
  //   throw(GWFrames_VectorSizeMismatch);
  // }
  vector<double> c(a);
  for(unsigned int i=0; i<size; ++i) {
    c[i] -= b[i];
  }
  return c;
}

std::vector<double> GWFrames::operator+(const std::vector<double>& a, const double b) {
  const unsigned int size = a.size();
  vector<double> c(a);
  for(unsigned int i=0; i<size; ++i) {
    c[i] += b;
  }
  return c;
}

std::vector<double> GWFrames::operator-(const std::vector<double>& a, const double b) {
  const unsigned int size = a.size();
  vector<double> c(a);
  for(unsigned int i=0; i<size; ++i) {
    c[i] -= b;
  }
  return c;
}

std::vector<double> GWFrames::operator/(const std::vector<double>& a, const double b) {
  const unsigned int size = a.size();
  vector<double> c(a);
  for(unsigned int i=0; i<size; ++i) {
    c[i] /= b;
  }
  return c;
}

std::vector<double> GWFrames::operator-(const std::vector<double>& a) {
  const unsigned int size = a.size();
  vector<double> c(a);
  for(unsigned int i=0; i<size; ++i) {
    c[i] *= -1;
  }
  return c;
}

std::vector<std::vector<double> > GWFrames::operator/(const std::vector<std::vector<double> >& a, const std::vector<double>& b) {
  const unsigned int size1 = a.size();
  if(size1<1) { return vector<vector<double> >(0); }
  const unsigned int size2 = a[0].size();
  vector<vector<double> > c(a);
  for(unsigned int i=0; i<size1; ++i) {
    for(unsigned int j=0; j<size2; ++j) {
      c[i][j] /= b[i];
    }
  }
  return c;
}

std::vector<std::vector<double> > GWFrames::operator*(const std::vector<double>& a, const std::vector<std::vector<double> >& b) {
  const unsigned int size1 = b.size();
  if(size1<1) { return vector<vector<double> >(0); }
  const unsigned int size2 = b[0].size();
  vector<vector<double> > c(b);
  for(unsigned int i=0; i<size1; ++i) {
    for(unsigned int j=0; j<size2; ++j) {
      c[i][j] *= a[i];
    }
  }
  return c;
}

std::vector<std::vector<double> > GWFrames::operator-(const std::vector<std::vector<double> >& a, const std::vector<std::vector<double> >& b) {
  const unsigned int size1 = a.size();
  if(size1==0) { return std::vector<std::vector<double> >(0); }
  const unsigned int size2 = a[0].size();
  // if(b.size() != size) {
  //   cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": a.size()=" << a.size() << "; b.size()=" << b.size() << endl;
  //   throw(GWFrames_VectorSizeMismatch);
  // }
  vector<vector<double> > c(size1, vector<double>(size2));
  for(unsigned int i=0; i<size1; ++i) {
    c[i] = a[i]-b[i];
  }
  return c;
}

std::vector<double> GWFrames::pow(const std::vector<double>& a, const double b) {
  const unsigned int size = a.size();
  vector<double> c(size);
  for(unsigned int i=0; i<size; ++i) {
    c[i] = std::pow(a[i], b);
  }
  return c;
}

/// Unwrap phase so that it is (roughly) continuous.
std::vector<double> GWFrames::Unwrap(const std::vector<double>& Arg) {
  // Compare Matlab's unwrap.m file
  vector<double> ArgUnwrapped = Arg;
  double Dp = 0.0;
  double Dps = 0.0;
  double CumCorr = 0.0;

  // Dp will contain the incremental phase variations;
  // Dps will contain the equivalents, confined to [-pi,pi)
  // CumCorr will contain the incremental phase corrections
  for(unsigned int i=1; i<Arg.size(); ++i) {
    Dp = Arg[i]-Arg[i-1];
    // C++'s fmod is unlike Matlab's for negative values, so:
    if(Dp+M_PI<0) {
      Dps = M_PI - fmod(-Dp-M_PI, 2.0*M_PI);
    } else {
      Dps = fmod(Dp+M_PI, 2.0*M_PI) - M_PI;
    }
    if(Dps==-M_PI && Dp>0) { Dps = M_PI; }
    CumCorr += Dps - Dp;
    ArgUnwrapped[i] += CumCorr;
  }

  return ArgUnwrapped;
}


/// Integrate scalar function by simple trapezoidal rule.
std::vector<double> GWFrames::ScalarIntegral(const std::vector<double>& fdot, const std::vector<double>& t) {
  ///
  /// \param fdot Vector of scalars.
  /// \param t Vector of corresponding time steps.
  if(fdot.size() != t.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": fdot.size()=" << fdot.size() << " != t.size()=" << t.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  const unsigned int Size=fdot.size();
  vector<double> f(Size);
  f[0] = 0.0;
  for(unsigned int i=1; i<Size; ++i) {
    f[i] = f[i-1] + (t[i]-t[i-1])*(fdot[i]+fdot[i-1])/2.0;
  }
  return f;
}

/// Integrate scalar function by simple trapezoidal rule.
double GWFrames::CumulativeScalarIntegral(const std::vector<double>& fdot, const std::vector<double>& t) {
  ///
  /// \param fdot Vector of scalars.
  /// \param t Vector of corresponding time steps.
  if(fdot.size() != t.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": fdot.size()=" << fdot.size() << " != t.size()=" << t.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  const unsigned int Size=fdot.size();
  double f = 0.0;
  for(unsigned int i=1; i<Size; ++i) {
    f += (t[i]-t[i-1])*(fdot[i]+fdot[i-1])/2.0;
  }
  return f;
}

inline double SQR(const double a) { return a*a; }

/// Five-point finite-differencing of vector of doubles.
std::vector<double> GWFrames::ScalarDerivative(const std::vector<double>& f, const std::vector<double>& t) {
  ///
  /// \param f Vector of doubles.
  /// \param t Vector of corresponding time steps.
  ///
  /// The formula for this finite difference comes from Eq. (A 5b) of
  /// "Derivative formulas and errors for non-uniformly spaced points"
  /// by M. K. Bowen and Ronald Smith.  As explained in their Eqs. (B
  /// 9b) and (B 10b), this is a fourth-order formula.
  ///
  /// If there are fewer than five points, the function reverts to
  /// simpler formulas.  If there are fewer than two points, or there
  /// are different numbers of points in the two input vectors, an
  /// exception is thrown.

  if(f.size() != t.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": f.size()=" << f.size() << " != t.size()=" << t.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  if(f.size()<2) { cerr << "\n" << __FILE__ << ":" << __LINE__ << ": size=" << f.size() << endl; throw(GWFrames_NotEnoughPointsForDerivative); }

  vector<double> D(f.size());
  const unsigned int i_f = f.size()-1;

  if(f.size()==2) {
    D[0] = (f[1]-f[0])/(t[1]-t[0]);
    D[1] = D[0];
    return D;
  }

  if(f.size()==3 || f.size()==4) {
    double hprev = t[1]-t[0];
    { // Compute first point
      const double hnext = t[2]-t[1];
      D[0] = -((2*hprev+hnext)/(hprev*(hprev+hnext)))*f[0] + ((hnext+hprev)/(hnext*hprev))*f[1] - (hprev/(hnext*(hnext+hprev)))*f[2];
    }
    for(unsigned int i=1; i<i_f; ++i) { // Compute intermediate points
      const double hnext = t[i+1]-t[i];
      /// Sundqvist and Veronis, Tellus XXII (1970), 1
      D[i] = (f[i+1] - f[i-1]*SQR(hnext/hprev) - f[i]*(1-SQR(hnext/hprev))) / (hnext*(1+hnext/hprev));
      hprev = hnext;
    }
    { // Compute final point
      const double hnext = t[i_f]  -t[i_f-1];
      const double hprev = t[i_f-1]-t[i_f-2];
      D[i_f] = (hnext/(hprev*(hprev+hnext)))*f[i_f-2] - ((hnext+hprev)/(hnext*hprev))*f[i_f-1] + ((hprev+2*hnext)/(hnext*(hnext+hprev)))*f[i_f];
    }
    return D;
  }

  for(unsigned int i=0; i<2; ++i) {
    const double x = t[i];
    const double& f1 = f[0];
    const double& f2 = f[1];
    const double& f3 = f[2];
    const double& f4 = f[3];
    const double& f5 = f[4];
    const double& x1 = t[0];
    const double& x2 = t[1];
    const double& x3 = t[2];
    const double& x4 = t[3];
    const double& x5 = t[4];
    const double h1 = x1 - x;
    const double h2 = x2 - x;
    const double h3 = x3 - x;
    const double h4 = x4 - x;
    const double h5 = x5 - x;
    const double h12 = x1 - x2;
    const double h13 = x1 - x3;
    const double h14 = x1 - x4;
    const double h15 = x1 - x5;
    const double h23 = x2 - x3;
    const double h24 = x2 - x4;
    const double h25 = x2 - x5;
    const double h34 = x3 - x4;
    const double h35 = x3 - x5;
    const double h45 = x4 - x5;
    D[i] =
      (-(h2*h3*h4 +h2*h3*h5 +h2*h4*h5 +h3*h4*h5)*f1/((h12)*(h13)*(h14)*(h15))
       +(h1*h3*h4 + h1*h3*h5 + h1*h4*h5 + h3*h4*h5)*f2/((h12)*(h23)*(h24)*(h25))
       -(h1*h2*h4 + h1*h2*h5 + h1*h4*h5 + h2*h4*h5)*f3/((h13)*(h23)*(h34)*(h35))
       +(h1*h2*h3 + h1*h2*h5 + h1*h3*h5 + h2*h3*h5)*f4/((h14)*(h24)*(h34)*(h45))
       -(h1*h2*h3 + h1*h2*h4 + h1*h3*h4 + h2*h3*h4)*f5/((h15)*(h25)*(h35)*(h45)));
  }
  for(unsigned int i=2; i<i_f-1; ++i) {
    const double& f1 = f[i-2];
    const double& f2 = f[i-1];
    const double& f3 = f[i];
    const double& f4 = f[i+1];
    const double& f5 = f[i+2];
    const double& x1 = t[i-2];
    const double& x2 = t[i-1];
    const double& x3 = t[i];
    const double& x4 = t[i+1];
    const double& x5 = t[i+2];
    const double h1 = x1 - x3;
    const double h2 = x2 - x3;
    const double h4 = x4 - x3;
    const double h5 = x5 - x3;
    const double h12 = x1 - x2;
    const double h13 = x1 - x3;
    const double h14 = x1 - x4;
    const double h15 = x1 - x5;
    const double h23 = x2 - x3;
    const double h24 = x2 - x4;
    const double h25 = x2 - x5;
    const double h34 = x3 - x4;
    const double h35 = x3 - x5;
    const double h45 = x4 - x5;
    D[i] =
      (-(h2*h4*h5)*f1/(h12*h13*h14*h15)
       +(h1*h4*h5)*f2/(h12*h23*h24*h25)
       -(h1*h2*h4 + h1*h2*h5 + h1*h4*h5 + h2*h4*h5)*f3/((h13)*(h23)*(h34)*(h35))
       +(h1*h2*h5)*f4/(h14*h24*h34*h45)
       -(h1*h2*h4)*f5/(h15*h25*h35*h45));
  }
  for(unsigned int i=i_f-1; i<i_f+1; ++i) {
    const double x = t[i];
    const double& f1 = f[i_f-4];
    const double& f2 = f[i_f-3];
    const double& f3 = f[i_f-2];
    const double& f4 = f[i_f-1];
    const double& f5 = f[i_f];
    const double& x1 = t[i_f-4];
    const double& x2 = t[i_f-3];
    const double& x3 = t[i_f-2];
    const double& x4 = t[i_f-1];
    const double& x5 = t[i_f];
    const double h1 = x1 - x;
    const double h2 = x2 - x;
    const double h3 = x3 - x;
    const double h4 = x4 - x;
    const double h5 = x5 - x;
    const double h12 = x1 - x2;
    const double h13 = x1 - x3;
    const double h14 = x1 - x4;
    const double h15 = x1 - x5;
    const double h23 = x2 - x3;
    const double h24 = x2 - x4;
    const double h25 = x2 - x5;
    const double h34 = x3 - x4;
    const double h35 = x3 - x5;
    const double h45 = x4 - x5;
    D[i] =
      (-(h2*h3*h4 +h2*h3*h5 +h2*h4*h5 +h3*h4*h5)*f1/((h12)*(h13)*(h14)*(h15))
       +(h1*h3*h4 + h1*h3*h5 + h1*h4*h5 + h3*h4*h5)*f2/((h12)*(h23)*(h24)*(h25))
       -(h1*h2*h4 + h1*h2*h5 + h1*h4*h5 + h2*h4*h5)*f3/((h13)*(h23)*(h34)*(h35))
       +(h1*h2*h3 + h1*h2*h5 + h1*h3*h5 + h2*h3*h5)*f4/((h14)*(h24)*(h34)*(h45))
       -(h1*h2*h3 + h1*h2*h4 + h1*h3*h4 + h2*h3*h4)*f5/((h15)*(h25)*(h35)*(h45)));
  }
  return D;

}
/// Five-point finite-differencing of vector of complex numbers.
std::vector<std::complex<double> > GWFrames::ComplexDerivative(const std::vector<std::complex<double> >& f, const std::vector<double>& t) {
  ///
  /// \param f Vector of std::complex<double>s.
  /// \param t Vector of corresponding time steps.
  ///
  /// The formula for this finite difference comes from Eq. (A 5b) of
  /// "Derivative formulas and errors for non-uniformly spaced points"
  /// by M. K. Bowen and Ronald Smith.  As explained in their Eqs. (B
  /// 9b) and (B 10b), this is a fourth-order formula.
  ///
  /// If there are fewer than five points, the function reverts to
  /// simpler formulas.  If there are fewer than two points, or there
  /// are different numbers of points in the two input vectors, an
  /// exception is thrown.

  if(f.size() != t.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": f.size()=" << f.size() << " != t.size()=" << t.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  if(f.size()<2) { cerr << "\n" << __FILE__ << ":" << __LINE__ << ": size=" << f.size() << endl; throw(GWFrames_NotEnoughPointsForDerivative); }

  vector<std::complex<double> > D(f.size());
  const unsigned int i_f = f.size()-1;

  if(f.size()==2) {
    D[0] = (f[1]-f[0])/(t[1]-t[0]);
    D[1] = D[0];
    return D;
  }

  if(f.size()==3 || f.size()==4) {
    double hprev = t[1]-t[0];
    { // Compute first point
      const double hnext = t[2]-t[1];
      D[0] = -((2*hprev+hnext)/(hprev*(hprev+hnext)))*f[0] + ((hnext+hprev)/(hnext*hprev))*f[1] - (hprev/(hnext*(hnext+hprev)))*f[2];
    }
    for(unsigned int i=1; i<i_f; ++i) { // Compute intermediate points
      const double hnext = t[i+1]-t[i];
      /// Sundqvist and Veronis, Tellus XXII (1970), 1
      D[i] = (f[i+1] - f[i-1]*SQR(hnext/hprev) - f[i]*(1-SQR(hnext/hprev))) / (hnext*(1+hnext/hprev));
      hprev = hnext;
    }
    { // Compute final point
      const double hnext = t[i_f]  -t[i_f-1];
      const double hprev = t[i_f-1]-t[i_f-2];
      D[i_f] = (hnext/(hprev*(hprev+hnext)))*f[i_f-2] - ((hnext+hprev)/(hnext*hprev))*f[i_f-1] + ((hprev+2*hnext)/(hnext*(hnext+hprev)))*f[i_f];
    }
    return D;
  }

  for(unsigned int i=0; i<2; ++i) {
    const double x = t[i];
    const std::complex<double>& f1 = f[0];
    const std::complex<double>& f2 = f[1];
    const std::complex<double>& f3 = f[2];
    const std::complex<double>& f4 = f[3];
    const std::complex<double>& f5 = f[4];
    const double& x1 = t[0];
    const double& x2 = t[1];
    const double& x3 = t[2];
    const double& x4 = t[3];
    const double& x5 = t[4];
    const double h1 = x1 - x;
    const double h2 = x2 - x;
    const double h3 = x3 - x;
    const double h4 = x4 - x;
    const double h5 = x5 - x;
    const double h12 = x1 - x2;
    const double h13 = x1 - x3;
    const double h14 = x1 - x4;
    const double h15 = x1 - x5;
    const double h23 = x2 - x3;
    const double h24 = x2 - x4;
    const double h25 = x2 - x5;
    const double h34 = x3 - x4;
    const double h35 = x3 - x5;
    const double h45 = x4 - x5;
    D[i] =
      (-(h2*h3*h4 +h2*h3*h5 +h2*h4*h5 +h3*h4*h5)*f1/((h12)*(h13)*(h14)*(h15))
       +(h1*h3*h4 + h1*h3*h5 + h1*h4*h5 + h3*h4*h5)*f2/((h12)*(h23)*(h24)*(h25))
       -(h1*h2*h4 + h1*h2*h5 + h1*h4*h5 + h2*h4*h5)*f3/((h13)*(h23)*(h34)*(h35))
       +(h1*h2*h3 + h1*h2*h5 + h1*h3*h5 + h2*h3*h5)*f4/((h14)*(h24)*(h34)*(h45))
       -(h1*h2*h3 + h1*h2*h4 + h1*h3*h4 + h2*h3*h4)*f5/((h15)*(h25)*(h35)*(h45)));
  }
  for(unsigned int i=2; i<i_f-1; ++i) {
    const std::complex<double>& f1 = f[i-2];
    const std::complex<double>& f2 = f[i-1];
    const std::complex<double>& f3 = f[i];
    const std::complex<double>& f4 = f[i+1];
    const std::complex<double>& f5 = f[i+2];
    const double& x1 = t[i-2];
    const double& x2 = t[i-1];
    const double& x3 = t[i];
    const double& x4 = t[i+1];
    const double& x5 = t[i+2];
    const double h1 = x1 - x3;
    const double h2 = x2 - x3;
    const double h4 = x4 - x3;
    const double h5 = x5 - x3;
    const double h12 = x1 - x2;
    const double h13 = x1 - x3;
    const double h14 = x1 - x4;
    const double h15 = x1 - x5;
    const double h23 = x2 - x3;
    const double h24 = x2 - x4;
    const double h25 = x2 - x5;
    const double h34 = x3 - x4;
    const double h35 = x3 - x5;
    const double h45 = x4 - x5;
    D[i] =
      (-(h2*h4*h5)*f1/(h12*h13*h14*h15)
       +(h1*h4*h5)*f2/(h12*h23*h24*h25)
       -(h1*h2*h4 + h1*h2*h5 + h1*h4*h5 + h2*h4*h5)*f3/((h13)*(h23)*(h34)*(h35))
       +(h1*h2*h5)*f4/(h14*h24*h34*h45)
       -(h1*h2*h4)*f5/(h15*h25*h35*h45));
  }
  for(unsigned int i=i_f-1; i<i_f+1; ++i) {
    const double x = t[i];
    const std::complex<double>& f1 = f[i_f-4];
    const std::complex<double>& f2 = f[i_f-3];
    const std::complex<double>& f3 = f[i_f-2];
    const std::complex<double>& f4 = f[i_f-1];
    const std::complex<double>& f5 = f[i_f];
    const double& x1 = t[i_f-4];
    const double& x2 = t[i_f-3];
    const double& x3 = t[i_f-2];
    const double& x4 = t[i_f-1];
    const double& x5 = t[i_f];
    const double h1 = x1 - x;
    const double h2 = x2 - x;
    const double h3 = x3 - x;
    const double h4 = x4 - x;
    const double h5 = x5 - x;
    const double h12 = x1 - x2;
    const double h13 = x1 - x3;
    const double h14 = x1 - x4;
    const double h15 = x1 - x5;
    const double h23 = x2 - x3;
    const double h24 = x2 - x4;
    const double h25 = x2 - x5;
    const double h34 = x3 - x4;
    const double h35 = x3 - x5;
    const double h45 = x4 - x5;
    D[i] =
      (-(h2*h3*h4 +h2*h3*h5 +h2*h4*h5 +h3*h4*h5)*f1/((h12)*(h13)*(h14)*(h15))
       +(h1*h3*h4 + h1*h3*h5 + h1*h4*h5 + h3*h4*h5)*f2/((h12)*(h23)*(h24)*(h25))
       -(h1*h2*h4 + h1*h2*h5 + h1*h4*h5 + h2*h4*h5)*f3/((h13)*(h23)*(h34)*(h35))
       +(h1*h2*h3 + h1*h2*h5 + h1*h3*h5 + h2*h3*h5)*f4/((h14)*(h24)*(h34)*(h45))
       -(h1*h2*h3 + h1*h2*h4 + h1*h3*h4 + h2*h3*h4)*f5/((h15)*(h25)*(h35)*(h45)));
  }
  return D;

}


/// Integrate vector function by simple trapezoidal rule.
std::vector<std::vector<double> > GWFrames::VectorIntegral(const std::vector<std::vector<double> >& fdot, const std::vector<double>& t) {
  ///
  /// \param fdot Vector of vectors (first index time).
  /// \param t Vector of corresponding time steps.
  if(fdot.size() != t.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": fdot.size()=" << fdot.size() << " != t.size()=" << t.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  if(fdot.size()==0) {
    return fdot;
  }
  const unsigned int Size1=fdot.size();
  const unsigned int Size2=fdot[0].size();
  vector<vector<double> > f(Size1, vector<double>(Size2));
  f[0] = vector<double>(Size2, 0.0);
  for(unsigned int i=1; i<Size1; ++i) {
    for(unsigned int j=0; j<Size2; ++j) {
      f[i][j] = f[i-1][j] + (t[i]-t[i-1])*(fdot[i][j]+fdot[i-1][j])/2.0;
    }
  }
  return f;
}

/// Integrate vector function by simple trapezoidal rule.
std::vector<double> GWFrames::CumulativeVectorIntegral(const std::vector<std::vector<double> >& fdot, const std::vector<double>& t) {
  ///
  /// \param fdot Vector of vectors (first index time).
  /// \param t Vector of corresponding time steps.
  if(fdot.size() != t.size()) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": fdot.size()=" << fdot.size() << " != t.size()=" << t.size() << std::endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  if(fdot.size()==0) {
    return vector<double>(0);
  }
  const unsigned int Size1=fdot.size();
  const unsigned int Size2=fdot[0].size();
  vector<double> f(Size2, 0.0);
  for(unsigned int i=1; i<Size1; ++i) {
    for(unsigned int j=0; j<Size2; ++j) {
      f[j] += (t[i]-t[i-1])*(fdot[i][j]+fdot[i-1][j])/2.0;
    }
  }
  return f;
}

/// Return the intersection of two time sequences.
std::vector<double> GWFrames::Intersection(const std::vector<double>& t1, const std::vector<double>& t2, const double MinStep, const double MinTime, const double MaxTime) {
  /// The time step at each point is the minimum of the time steps in
  /// t1 and t2 at that instant, or MinStep, whichever is greater.
  /// The output starts at the earliest moment common to t1 and t2, or
  /// MinTime, whichever is greater.
  ///
  /// The input to this function is assumed to be strictly monotonic.
  if(t1.size()==0) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": t1 is empty.  Assuming this is not desired." << std::endl;
    throw(GWFrames_EmptyIntersection);
  }
  if(t2.size()==0) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": t2 is empty.  Assuming this is not desired." << std::endl;
    throw(GWFrames_EmptyIntersection);
  }
  vector<double> t(t1.size()+t2.size());
  double min1 = t1[0];
  double min2 = t2[0];
  double mint = std::max(std::max(min1, min2), MinTime);
  double max1 = t1.back();
  double max2 = t2.back();
  double maxt = std::min(std::min(max1, max2), MaxTime);
  if(mint > max1 || mint > max2) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": t1=[" << min1 << ", ..., " << max1 << "]\tt2=[" << min2 << ", ..., " << max2 << "]"
              << "\nIntersection is empty." << std::endl;
    throw(GWFrames_EmptyIntersection);
  }
  t[0] = mint;
  unsigned int I=0, I1=0, I2=0;
  const unsigned int size1=t1.size(), size2=t2.size();
  while(t[I] < maxt) {
    // adjust I1 to ensure that t[I] is in the interval ( t1[I1-1], t1[I1] ]
    if(t[I]<min1 || t[I]>max1) { // if t[I] is less than the smallest t1, or greater than the largest t1, I1=0;
      I1=0;
    } else {
      I1 = std::max(I1, (unsigned int)1);
      while(t[I]>t1[I1] && I1<size1) { I1++; }
    }
    // adjust I2 to ensure that t[I] is in the interval ( t2[I2-1], t2[I2] ]
    if(t[I]<min2 || t[I]>max2) { // if t[I] is less than the smallest t2, or greater than the largest t2, I2=0;
      I2=0;
    } else {
      I2 = std::max(I2, (unsigned int)1);
      while(t[I]>t2[I2] && I2<size2) { I2++; }
    }
    t[I+1] = t[I] + std::max(std::min(t1[I1]-t1[I1-1], t2[I2]-t2[I2-1]), MinStep);
    if(t[I+1]>maxt) { break; }
    I++;
  }
  t.erase(t.begin()+I+1, t.end()); // only take the relevant part of the reserved vector
  return t;
}

/// Return the union of two time sequences.
std::vector<double> GWFrames::Union(const std::vector<double>& t1, const std::vector<double>& t2, const double MinStep) {
  /// On the overlap between the two sequences, the time is built up
  /// by taking the smaller time step in either of the two sequences,
  /// or MinStep if that step is smaller.
  ///
  /// The input to this function is assumed to be strictly monotonic.
  if(t1.size()==0) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": t1 is empty.  Returning trivial Union." << std::endl;
    return vector<double>(t2);
  }
  if(t2.size()==0) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": t2 is empty.  Returning trivial Union." << std::endl;
    return vector<double>(t1);
  }
  vector<double> t(t1.size()+t2.size()+4);
  double min1 = t1[0];
  double min2 = t2[0];
  double mint = std::min(min1, min2);
  double max1 = t1.back();
  double max2 = t2.back();
  double maxt = std::max(max1, max2);
  if(min2 > max1 || min1 > max2) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Disjoint union.  Assuming this is not desired."
              << "\nmin(t1)=" << min1 << "\tmin(t2)=" << min2
              << "\nmax(t1)=" << max1 << "\tmax(t2)=" << max2 << endl;
    throw(GWFrames_EmptyIntersection);
  }
  t[0] = mint;
  unsigned int I=0, I1=0, I2=0;
  const unsigned int size1=t1.size(), size2=t2.size();
  while(t[I] < maxt) {
    // adjust I1 to ensure that t[I] is in the interval ( t1[I1-1], t1[I1] ]
    if(t[I]<min1 || t[I]>max1) { // if t[I] is less than the smallest t1, or greater than the largest t1, I1=0;
      I1=0;
    } else {
      I1 = std::max(I1, (unsigned int)1);
      while(t[I]>t1[I1] && I1<size1) { ++I1; }
    }
    // adjust I2 to ensure that t[I] is in the interval ( t2[I2-1], t2[I2] ]
    if(t[I]<min2 || t[I]>max2) { // if t[I] is less than the smallest t2, or greater than the largest t2, I2=0;
      I2=0;
    } else {
      I2 = std::max(I2, (unsigned int)1);
      while(t[I]>t2[I2] && I2<size2) { ++I2; }
    }
    if(I1==0) {
      t[I+1] = t[I] + std::max(t2[I2]-t2[I2-1], MinStep);
    } else if(I2==0) {
      t[I+1] = t[I] + std::max(t1[I1]-t1[I1-1], MinStep);
    } else {
      t[I+1] = t[I] + std::max(std::min(t1[I1]-t1[I1-1], t2[I2]-t2[I2-1]), MinStep);
    }
    if(t[I+1]>maxt) { break; }
    ++I;
  }
  t.erase(t.begin()+I+1, t.end()); // only take the relevant part of the reserved vector
  return t;
}


Matrix::Matrix()
  : m(NULL)
{ }

Matrix::Matrix(unsigned int rows, unsigned int cols)
  : m(gsl_matrix_calloc(rows, cols))
{ }

Matrix::Matrix(unsigned int rows, unsigned int cols, const double a)
  : m(gsl_matrix_alloc(rows, cols))
{
  gsl_matrix_set_all(m, a);
}

Matrix::Matrix(const Matrix& rhs)
  : m(gsl_matrix_alloc(rhs.nrows(), rhs.ncols()))
{
  gsl_matrix_memcpy(m, rhs.m);
}

Matrix::Matrix(const std::vector<std::vector<double> >& DataIn)
  : m(gsl_matrix_alloc(DataIn.size(), (DataIn.size()==0 ? 0 : DataIn[0].size())))
{
  for(unsigned int r=0; r<m->size1; ++r) {
    for(unsigned int c=0; c<m->size2; ++c) {
      gsl_matrix_set(m, r, c, DataIn[r][c]);
    }
  }
}

Matrix& Matrix::operator=(const Matrix& rhs) {
  if(m) { gsl_matrix_free(m); }
  m = gsl_matrix_alloc(rhs.nrows(), rhs.ncols());
  gsl_matrix_memcpy(m, rhs.m);
  return *this;
}

Matrix& Matrix::operator=(const std::vector<std::vector<double> >& DataIn) {
  if(m) { gsl_matrix_free(m); }
  m = gsl_matrix_alloc(DataIn.size(), (DataIn.size()==0 ? 0 : DataIn[0].size()));
  for(unsigned int r=0; r<m->size1; ++r) {
    for(unsigned int c=0; c<m->size2; ++c) {
      gsl_matrix_set(m, r, c, DataIn[r][c]);
    }
  }
  return *this;
}

Matrix Matrix::operator-(const Matrix& rhs) {
  if(ncols() != rhs.ncols() || nrows() != rhs.nrows()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__
         << ": ncols=" << ncols() << "; rhs.ncols()=" << rhs.ncols()
         << "\nnrows=" << nrows() << "; rhs.nrows()=" << rhs.nrows() << endl;
    throw(GWFrames_MatrixSizeMismatch);
  }
  Matrix M(*this);
  for(unsigned int r=0; r<M.nrows(); ++r) {
    for(unsigned int c=0; c<M.ncols(); ++c) {
      gsl_matrix_set(M.m, r, c, M(r,c)-rhs(r,c));
    }
  }
  return M;
}

// / \@cond
void Matrix::resize(unsigned int newNRows, unsigned int newNCols, const double a) {
  if(m) { gsl_matrix_free(m); }
  if(a!=0.0) {
    m = gsl_matrix_alloc(newNRows, newNCols);
    gsl_matrix_set_all(m, a);
  } else {
    m = gsl_matrix_calloc(newNRows, newNCols);
  }
  return;
}
// / \@endcond

void Matrix::clear() {
  if(m) { gsl_matrix_free(m); }
  return;
}

void Matrix::swap(Matrix& b) {
  gsl_matrix* tmp = m;
  m = b.m;
  b.m = tmp;
  return;
}

std::vector<double> Matrix::operator*(const std::vector<double>& b) const {
  if(ncols() != b.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": ncols=" << ncols() << "; vector.size()=" << b.size() << endl;
    throw(GWFrames_MatrixSizeMismatch);
  }
  const unsigned int C=ncols();
  const unsigned int R=nrows();
  vector<double> Result(R, 0.0);
  for(unsigned int r=0; r<R; ++r) {
    for(unsigned int c=0; c<C; ++c) {
      Result[r] += this->operator()(r,c) * b[c];
    }
  }
  return Result;
}

// Quaternion Matrix::operator*(const Quaternion& b) const {
//   if(ncols() != 3) {
//     cerr << "\n\nn" << __FILE__ << ":" << __LINE__ << ": cols=" << ncols() << "; Quaternions have 3 vector components" << endl;
//     throw(GWFrames_MatrixSizeMismatch);
//   }
//   if(nrows() != 3) {
//     cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": nrows=" << nrows() << "; Quaternions have 3 vector components" << endl;
//     throw(GWFrames_MatrixSizeMismatch);
//   }
//   const unsigned int C=3;
//   const unsigned int R=3;
//   Quaternion Result(0,0,0,0);
//   for(unsigned int r=0; r<R; ++r) {
//     for(unsigned int c=0; c<C; ++c) {
//       Result[r+1] += this->operator()(r,c) * b[c+1];
//     }
//   }
//   return Result;
// }

std::vector<double> GWFrames::operator*(const std::vector<double>& a, const Matrix& b) {
  if(b.nrows() != a.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": b.nrows=" << b.nrows() << "; vector.size()=" << a.size() << endl;
    throw(GWFrames_MatrixSizeMismatch);
  }
  const unsigned int C=b.ncols();
  const unsigned int R=b.nrows();
  vector<double> Result(C, 0.0);
  for(unsigned int c=0; c<C; ++c) {
    for(unsigned int r=0; r<R; ++r) {
      Result[c] += b(r,c) * a[r];
    }
  }
  return Result;
}

std::vector<double> GWFrames::DominantPrincipalAxis(Matrix& M) {
  if(M.nrows()!=3 || M.ncols()!=3) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": M.nrows()=" << M.nrows() << "; M.ncols()=" << M.ncols() << endl;
    throw(GWFrames_MatrixSizeAssumedToBeThree);
  }
  gsl_vector *eval = gsl_vector_alloc (3);
  gsl_matrix *evec = gsl_matrix_alloc (3, 3);
  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (3);
  gsl_eigen_symmv(M.gslobj(), eval, evec, w); // Do the work
  gsl_eigen_symmv_free(w);
  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC); // Sort by eigenvalue magnitude
  gsl_vector_const_view evec_Dominant = gsl_matrix_const_column(evec, 0);
  vector<double> Result(3);
  for(unsigned int i=0; i<3; ++i) {
    Result[i] = gsl_vector_get(&evec_Dominant.vector, i);
  }
  gsl_vector_free(eval);
  gsl_matrix_free(evec);
  return Result;
}

std::vector<double> GWFrames::Eigenvalues(Matrix& M) {
  if(M.nrows()!=3 || M.ncols()!=3) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": M.nrows()=" << M.nrows() << "; M.ncols()=" << M.ncols() << endl;
    throw(GWFrames_MatrixSizeAssumedToBeThree);
  }
  gsl_vector *eval = gsl_vector_alloc (3);
  gsl_matrix *evec = gsl_matrix_alloc (3, 3);
  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (3);
  gsl_eigen_symmv(M.gslobj(), eval, evec, w); // Do the work
  gsl_eigen_symmv_free(w);
  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC); // Sort by eigenvalue magnitude
  vector<double> Result(3);
  for(unsigned int i=0; i<3; ++i) {
    Result[i] = gsl_vector_get(eval, i);
  }
  gsl_vector_free(eval);
  gsl_matrix_free(evec);
  return Result;
}

std::vector<double> GWFrames::Eigenvectors(Matrix& M) {
  if(M.nrows()!=3 || M.ncols()!=3) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": M.nrows()=" << M.nrows() << "; M.ncols()=" << M.ncols() << endl;
    throw(GWFrames_MatrixSizeAssumedToBeThree);
  }
  gsl_vector *eval = gsl_vector_alloc (3);
  gsl_matrix *evec = gsl_matrix_alloc (3, 3);
  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (3);
  gsl_eigen_symmv(M.gslobj(), eval, evec, w); // Do the work
  gsl_eigen_symmv_free(w);
  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC); // Sort by eigenvalue magnitude
  vector<double> Result(9);
  for(unsigned int j=0; j<3; ++j) {
    gsl_vector_const_view evec_Dominant = gsl_matrix_const_column(evec, j);
    for(unsigned int i=0; i<3; ++i) {
      Result[i+3*j] = gsl_vector_get(&evec_Dominant.vector, i);
    }
  }
  gsl_vector_free(eval);
  gsl_matrix_free(evec);
  return Result;
}

std::vector<double> GWFrames::Eigensystem(Matrix& M) {
  if(M.nrows()!=3 || M.ncols()!=3) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": M.nrows()=" << M.nrows() << "; M.ncols()=" << M.ncols() << endl;
    throw(GWFrames_MatrixSizeAssumedToBeThree);
  }
  gsl_vector *eval = gsl_vector_alloc (3);
  gsl_matrix *evec = gsl_matrix_alloc (3, 3);
  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (3);
  gsl_eigen_symmv(M.gslobj(), eval, evec, w); // Do the work
  gsl_eigen_symmv_free(w);
  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC); // Sort by eigenvalue magnitude
  vector<double> Result(12);
  for(unsigned int j=0; j<3; ++j) {
    gsl_vector_const_view evec_Dominant = gsl_matrix_const_column(evec, j);
    for(unsigned int i=0; i<3; ++i) {
      Result[i+4*j] = gsl_vector_get(&evec_Dominant.vector, i);
    }
    Result[3+4*j] = gsl_vector_get(eval, j);
  }
  gsl_vector_free(eval);
  gsl_matrix_free(evec);
  return Result;
}

double GWFrames::Determinant(Matrix& M) {
  if(M.nrows()!=3 || M.ncols()!=3) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": M.nrows()=" << M.nrows() << "; M.ncols()=" << M.ncols() << endl;
    throw(GWFrames_MatrixSizeAssumedToBeThree);
  }
  return M(0,0)*M(1,1)*M(2,2) + M(0,1)*M(1,2)*M(2,0) + M(0,2)*M(1,0)*M(2,1) - M(0,2)*M(1,1)*M(2,0) - M(0,1)*M(1,0)*M(2,2) - M(0,0)*M(1,2)*M(2,1);
}






////////////////////////////////////////////////////////////////

MatrixC::MatrixC()
  : nn(0), mm(0), v(NULL)
{ }

MatrixC::MatrixC(int n, int m)
  : nn(n), mm(m), v(n>0 ? new std::complex<double>*[n] : NULL)
{
  const int nel=m*n;
  if(v) v[0] = nel>0 ? new std::complex<double>[nel] : NULL;
  for(int i=1;i<n;i++) {
    v[i] = v[i-1] + m;
  }
}

MatrixC::MatrixC(int n, int m, const std::complex<double> &a)
  : nn(n), mm(m), v(n>0 ? new std::complex<double>*[n] : NULL)
{
  const int nel=m*n;
  if (v) v[0] = nel>0 ? new std::complex<double>[nel] : NULL;
  for (int i=1; i< n; i++) v[i] = v[i-1] + m;
  for (int i=0; i< n; i++) for (int j=0; j<m; j++) v[i][j] = a;
}

MatrixC::MatrixC(int n, int m, const std::complex<double> *a)
  : nn(n), mm(m), v(n>0 ? new std::complex<double>*[n] : NULL)
{
  const int nel=m*n;
  if (v) v[0] = nel>0 ? new std::complex<double>[nel] : NULL;
  for (int i=1; i<n; i++) v[i] = v[i-1] + m;
  for (int i=0; i<n; i++) for (int j=0; j<m; j++) v[i][j] = *a++;
}

MatrixC::MatrixC(const MatrixC &rhs)
  : nn(rhs.nn), mm(rhs.mm), v(nn>0 ? new std::complex<double>*[nn] : NULL)
{
  const int nel=mm*nn;
  if (v) v[0] = nel>0 ? new std::complex<double>[nel] : NULL;
  for (int i=1; i<nn; i++) v[i] = v[i-1] + mm;
  for (int i=0; i<nn; i++) for (int j=0; j<mm; j++) v[i][j] = rhs[i][j];
}

MatrixC::MatrixC(const std::vector<std::vector<std::complex<double> > >& rhs)
  : nn(rhs.size()), mm(nn>0 ? rhs[0].size() : 0), v(nn>0 ? new std::complex<double>*[nn] : NULL)
{
  const int nel=mm*nn;
  if (v) v[0] = nel>0 ? new std::complex<double>[nel] : NULL;
  for (int i=1; i<nn; i++) v[i] = v[i-1] + mm;
  for (int i=0; i<nn; i++) for (int j=0; j<mm; j++) v[i][j] = rhs[i][j];
}

MatrixC & MatrixC::operator=(const MatrixC &rhs) {
  if (this != &rhs) {
    if (nn != rhs.nn || mm != rhs.mm) {
      if (v != NULL) {
        delete[] (v[0]);
        delete[] (v);
      }
      nn=rhs.nn;
      mm=rhs.mm;
      v = nn>0 ? new std::complex<double>*[nn] : NULL;
      const int nel = mm*nn;
      if (v) v[0] = nel>0 ? new std::complex<double>[nel] : NULL;
      for (int i=1; i< nn; i++) v[i] = v[i-1] + mm;
    }
    for (int i=0; i<nn; i++) for (int j=0; j<mm; j++) v[i][j] = rhs[i][j];
  }
  return *this;
}

void MatrixC::swap(MatrixC& b) {
  { const int n = b.nn; b.nn=nn; nn=n; }
  { const int m = b.mm; b.mm=mm; mm=m; }
  { std::complex<double>** vv=b.v; b.v=v; v=vv; }
  return;
}

// / \@cond
void MatrixC::resize(int newn, int newm) {
  if (newn != nn || newm != mm) {
    if (v != NULL) {
      delete[] (v[0]);
      delete[] (v);
    }
    nn = newn;
    mm = newm;
    v = nn>0 ? new std::complex<double>*[nn] : NULL;
    const int nel = mm*nn;
    if (v) { v[0] = nel>0 ? new std::complex<double>[nel] : NULL; }
    for(int i=1; i< nn; i++) {
      v[i] = v[i-1] + mm;
    }
  }
}
// / \@endcond

void MatrixC::assign(int newn, int newm, const std::complex<double>& a) {
  if (newn != nn || newm != mm) {
    if (v != NULL) {
      delete[] (v[0]);
      delete[] (v);
    }
    nn = newn;
    mm = newm;
    v = nn>0 ? new std::complex<double>*[nn] : NULL;
    const int nel = mm*nn;
    if(v) { v[0] = nel>0 ? new std::complex<double>[nel] : NULL; }
    for(int i=1; i< nn; i++) {
      v[i] = v[i-1] + mm;
    }
  }
  for(int i=0; i< nn; i++) {
    for(int j=0; j<mm; j++) {
      v[i][j] = a;
    }
  }
}

MatrixC::~MatrixC()
{
  if (v != NULL) {
    delete[] (v[0]);
    delete[] (v);
  }
}


///////////////////////////////////////////////////////////////////
