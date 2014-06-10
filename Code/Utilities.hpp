// Copyright (c) 2013, Michael Boyle
// See LICENSE file for details

#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <vector>
#include <complex>
#include <iostream>
#include <gsl/gsl_matrix.h>

#include "Quaternions.hpp"
#include "Errors.hpp"

namespace GWFrames {
  
  // Useful operations on vectors
  double abs(const std::vector<double>& v);
  std::vector<double> abs(const std::vector<std::vector<double> >& v);
  std::vector<double> fabs(const std::vector<double>& x);
  std::vector<double> pow(const std::vector<double>& base, const double& exponent);
  std::vector<double> log(const std::vector<double>& x);
  std::vector<double> exp(const std::vector<double>& x);
  std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b);
  std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b);
  std::vector<double> operator+(const std::vector<double>& a, const double b);
  std::vector<double> operator-(const std::vector<double>& a, const double b);
  std::vector<double> operator/(const std::vector<double>& a, const double b);
  std::vector<double> operator*(const std::vector<double>& a, const double b);
  std::vector<double> operator*(const double a, const std::vector<double>& b);
  std::vector<std::vector<double> > operator/(const std::vector<std::vector<double> > & a, const std::vector<double>& b);
  std::vector<double> Unwrap(const std::vector<double>& In);

  
  // Integrals and derivatives
  std::vector<double> ScalarIntegral(const std::vector<double>& fdot, const std::vector<double>& t);
  double CumulativeScalarIntegral(const std::vector<double>& fdot, const std::vector<double>& t);
  std::vector<double> ScalarDerivative(const std::vector<double>& f, const std::vector<double>& t);
  std::vector<std::complex<double> > ComplexDerivative(const std::vector<std::complex<double> >& f, const std::vector<double>& t);
  std::vector<std::vector<double> > VectorIntegral(const std::vector<std::vector<double> >& fdot, const std::vector<double>& t);
  std::vector<double> CumulativeVectorIntegral(const std::vector<std::vector<double> >& fdot, const std::vector<double>& t);
  
  // Common-time functions
  std::vector<double> Intersection(const std::vector<double>& t1, const std::vector<double>& t2,
				   const double MinStep=0.005, const double MinTime=-1e300, const double MaxTime=1e300);
  std::vector<double> Union(const std::vector<double>& t1, const std::vector<double>& t2, const double MinStep=0.005);
  
  /// 3x3 object wrapping GSL matrix; probably not needed directly
  class Matrix {
  private:
    gsl_matrix* m;
  public:
    Matrix();
    Matrix(unsigned int rows, unsigned int cols); // Zero-based array
    Matrix(unsigned int rows, unsigned int cols, const double a); // array of a's
    Matrix(const Matrix &rhs);
    Matrix(const std::vector<std::vector<double> >& DataIn);
    Matrix& operator=(const Matrix &rhs);
    Matrix& operator=(const std::vector<std::vector<double> >& newData);
    Matrix operator-(const Matrix &rhs);
    inline double operator()(const unsigned int row, const unsigned int col) const { return *gsl_matrix_ptr(m, row, col); }
    inline double& operator()(const unsigned int row, const unsigned int col) { return *gsl_matrix_ptr(m, row, col); }
    inline Matrix& set(const unsigned int r, const unsigned int c, const double v) { gsl_matrix_set(m, r, c, v); return *this; }
    inline unsigned int nrows() const { return m->size1; }
    inline unsigned int ncols() const { return m->size2; }
    // / \@cond
    void resize(unsigned int newNRows, unsigned int newNCols, const double=0.0);
    // / \@endcond
    void clear(); // empty contents of this Matrix
    void swap(Matrix& b); // Swap data sets
    ~Matrix() { if(m) { gsl_matrix_free(m); } }
    
    std::vector<double> operator*(const std::vector<double>& b) const;
    Quaternion operator*(const Quaternion& b) const;
    
    inline gsl_matrix* gslobj() { return m; }
    inline const gsl_matrix* gslobj() const { return m; }
  };
  std::vector<double> operator*(const std::vector<double>& a, const Matrix& b);
  Quaternion operator*(const Quaternion& a, const Matrix& b);
  std::vector<double> DominantPrincipalAxis(Matrix& M);
  std::vector<double> Eigenvalues(Matrix& M);
  std::vector<double> Eigenvectors(Matrix& M);
  std::vector<double> Eigensystem(Matrix& M);
  double Determinant(Matrix& M);
  
  /// Rectangular array of complex data; probably not needed directly
  class MatrixC {
  private:
    int nn;
    int mm;
    std::complex<double> **v;
  public:
    MatrixC();
    MatrixC(int n, int m);			// Zero-based array
    MatrixC(int n, int m, const std::complex<double> &a);	//Initialize to constant
    MatrixC(int n, int m, const std::complex<double> *a);	// Initialize to array
    MatrixC(const std::vector<std::vector<std::complex<double> > >& DataIn);
    MatrixC(const MatrixC &rhs);		// Copy constructor
    MatrixC& operator=(const MatrixC &rhs);	//assignment
    void swap(MatrixC& b);
    inline std::complex<double>* operator[](const int i) { return v[i]; }
    inline const std::complex<double>* operator[](const int i) const { return v[i]; }
    inline int nrows() const { return nn; }
    inline int ncols() const { return mm; }
    // / \@cond
    void resize(int newn, int newm); // resize (contents not preserved)
    // / \@endcond
    void assign(int newn, int newm, const std::complex<double> &a); // resize and assign a constant value
    ~MatrixC();
  };
  // std::vector<std::complex<double> > operator*(const std::vector<std::complex<double> >& b) const;
  // std::vector<std::complex<double> > operator*(const std::vector<std::complex<double> >& a, const MatrixC& b);
  
  const int ellMax_Utilities = 16;
  
  /// Object for pre-computing and retrieving factorials
  class FactorialFunctor {
  private:
    static const std::vector<double> FactorialTable;
  public:
    FactorialFunctor() { };
    inline double operator()(const unsigned int i) const {
      #ifdef DEBUG
      if(i>171) {
  	std::cerr << "\n\ni = " << i << "\tiMax = 171"
  		  << "\nFactorialFunctor is only implemented up to 171!; larger values overflow."
  		  << std::endl;
  	throw(GWFrames_IndexOutOfBounds);
      }
      #endif
      return FactorialTable[i];
    }
  };
  
  /// Object for pre-computing and retrieving binomials
  class BinomialCoefficientFunctor {
  private:
    static const std::vector<double> BinomialCoefficientTable;
  public:
    BinomialCoefficientFunctor() { };
    inline double operator()(const unsigned int n, const unsigned int k) const {
      #ifdef DEBUG
      if(n>2*ellMax_Utilities || k>n) {
  	std::cerr << "\n\n(n, k) = (" << n << ", " << k << ")\t2*ellMax_Utilities = " << 2*ellMax_Utilities
  		  << "\nBinomialCoefficientFunctor is only implemented up to n=2*ellMax_Utilities=" << 2*ellMax_Utilities
  		  << ".\nTo increase this bound, edit 'ellMax_Utilities' in " << __FILE__ << " and recompile." << std::endl;
  	throw(GWFrames_IndexOutOfBounds);
      }
      #endif
      return BinomialCoefficientTable[(n*(n+1))/2+k];
    }
  };
  
  /// Object for pre-computing and retrieving values of the ladder operators
  class LadderOperatorFactorFunctor {
  private:
    static const std::vector<double> FactorTable;
  public:
    LadderOperatorFactorFunctor() { };
    inline double operator()(const int ell, const int m) const {
      #ifdef DEBUG
      if(ell>ellMax_Utilities || std::abs(m)>ell) {
  	std::cerr << "\n\n(ell, m) = (" << ell << ", " << m << ")\tellMax_Utilities = " << ellMax_Utilities
  		  << "\nLadderOperatorFactorFunctor is only implemented up to ell=" << ellMax_Utilities
  		  << ".\nTo increase this bound, edit 'ellMax_Utilities' in " << __FILE__ << " and recompile." << std::endl;
  	throw(GWFrames_IndexOutOfBounds);
      }
      #endif
      return FactorTable[ell*ell+ell+m];
    }
  };
  
  /// Object for pre-computing and retrieving coefficients for the Wigner D matrices
  class WignerCoefficientFunctor {
  private:
    static const std::vector<double> CoefficientTable;
  public:
    WignerCoefficientFunctor() { };
    inline double operator()(const int ell, const int mp, const int m) const {
      #ifdef DEBUG
      if(ell>ellMax_Utilities || std::abs(mp)>ell || std::abs(m)>ell) {
  	std::cerr << "\n\n(ell, mp, m) = (" << ell << ", " << mp << ", " << m << ")\tellMax_Utilities = " << ellMax_Utilities
  		  << "\nWignerCoefficientFunctor is only implemented up to ell=" << ellMax_Utilities
  		  << ".\nTo increase this bound, edit 'ellMax_Utilities' in " << __FILE__ << " and recompile." << std::endl;
  	throw(GWFrames_IndexOutOfBounds);
      }
      #endif
      return CoefficientTable[int(ell*(ell*(1.3333333333333333*ell + 2) + 1.6666666666666667) + mp*(2*ell + 1) + m + 0.5)];
    }
  };
  
  /// Object for computing the Wigner D matrices as functions of quaternion rotors
  class WignerDMatrix {
  private:
    BinomialCoefficientFunctor BinomialCoefficient;
    WignerCoefficientFunctor WignerCoefficient;
    std::complex<double> Ra, Rb;
    double absRa, absRb, absRRatioSquared;
  public:
    WignerDMatrix(const Quaternion& iR=Quaternion(1,0,0,0));
    WignerDMatrix& SetRotation(const Quaternion& iR);
    std::complex<double> operator()(const int ell, const int mp, const int m) const;
  };
  
  /// Object for computing values of the spin-weighted spherical harmonics
  class SWSH {
  private:
    WignerDMatrix D;
    int spin;
    double sign;
  public:
    // / \@cond
    SWSH(const int s=-2, const Quaternion& iR=Quaternion(1,0,0,0))
      : D(iR), spin(s), sign(s%2==0 ? 1.0 : -1.0)
    { }
    // / \@endcond
    inline SWSH& SetRotation(const Quaternion& iR) { D.SetRotation(iR); return *this; }
    inline SWSH& SetAngles(const double vartheta, const double varphi) { D.SetRotation(Quaternion(vartheta, varphi)); return *this; }
    inline std::complex<double> operator()(const int ell, const int m) const {
      return sign * std::sqrt((2*ell+1)/(4*M_PI)) * D(ell, m, -spin);
    }
  };
  
  
} // namespace GWFrames

#endif // UTILITIES_HPP
