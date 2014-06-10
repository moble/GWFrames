// Copyright (c) 2014, Michael Boyle
// See LICENSE file for details

#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <vector>
#include <complex>
#include <iostream>
#include <gsl/gsl_matrix.h>

namespace GWFrames {

  // Typedefs
  typedef std::vector<double> ThreeVector; // Can be assumed to have three components
  typedef std::vector<double> FourVector; // Can be assumed to have four components

  // Useful operations on vectors
  double abs(const std::vector<double>& v);
  std::vector<double> abs(const std::vector<std::vector<double> >& v);
  std::vector<double> fabs(const std::vector<double>& x);
  std::vector<double> pow(const std::vector<double>& base, const double exponent);
  std::vector<double> log(const std::vector<double>& x);
  std::vector<double> exp(const std::vector<double>& x);
  std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b);
  std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b);
  std::vector<double> operator+(const std::vector<double>& a, const double b);
  std::vector<double> operator-(const std::vector<double>& a, const double b);
  std::vector<double> operator/(const std::vector<double>& a, const double b);
  std::vector<double> operator-(const std::vector<double>& a);
  std::vector<double> operator*(const std::vector<double>& a, const double b);
  std::vector<double> operator*(const double a, const std::vector<double>& b);
  std::vector<std::vector<double> > operator/(const std::vector<std::vector<double> > & a, const std::vector<double>& b);
  std::vector<std::vector<double> > operator*(const std::vector<double> & a, const std::vector<std::vector<double> >& b);
  std::vector<std::vector<double> > operator-(const std::vector<std::vector<double> >& a, const std::vector<std::vector<double> >& b);
  std::vector<double> Unwrap(const std::vector<double>& In);
  bool xINy(const int x, const std::vector<int>& y);

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
    // Quaternion operator*(const Quaternion& b) const;

    inline gsl_matrix* gslobj() { return m; }
    inline const gsl_matrix* gslobj() const { return m; }
  };
  std::vector<double> operator*(const std::vector<double>& a, const Matrix& b);
  // Quaternion operator*(const Quaternion& a, const Matrix& b);
  std::vector<double> DominantPrincipalValue(std::vector<Matrix>& M);
  std::vector<double> SubordinatePrincipalValue(std::vector<Matrix>& M);
  std::vector<double> DominantPrincipalAxis(Matrix& M);
  std::vector<double> SubordinatePrincipalAxis(Matrix& M);
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

  std::ostream& operator<<(std::ostream& out, const std::vector<double>& v);
  std::ostream& operator<<(std::ostream& out, const std::vector<int>& v);
  std::ostream& operator<<(std::ostream& out, const std::vector<std::vector<int> >& vv);
  std::ostream& operator<<(std::ostream& out, const std::vector<unsigned int>& v);

} // namespace GWFrames

#endif // UTILITIES_HPP
