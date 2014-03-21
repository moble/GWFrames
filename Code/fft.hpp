#ifndef FFT_HPP
#define FFT_HPP

#include <vector>
#include <complex>

namespace WaveformUtilities {
  
  /// This helper class is useful for untangling NR's fft routines.
  class WrapVecDoub {
  private:
    std::vector<double> vvec;
    std::vector<double> &v;
    int n, mask;
  public:
    WrapVecDoub(const int nn) : vvec(nn, 0.0), v(vvec), n(nn/2), mask(n-1) {validate();}
    WrapVecDoub(std::vector<double>&vec) : v(vec), n(vec.size()/2), mask(n-1) {validate();}
    void validate();
    inline std::complex<double>& operator[] (int i) {return (std::complex<double> &)v[(i&mask) << 1];}
    inline double& real(int i) {return v[(i&mask) << 1];}
    inline double& imag(int i) {return v[((i&mask) << 1)+1];}
    operator std::vector<double>&() {return v;}
  };
  
  /// This function creates a frequency vector in (0 -> positives -> negatives -> 0) order
  std::vector<double> TimeToFrequency(const std::vector<double>& Time);
  
  /// This function returns the positive half of the frequencies, so returned size is 1/2 input size + 1
  std::vector<double> TimeToPositiveFrequencies(const std::vector<double>& Time);
  
  /// The following call the Numerical Recipes fft routines (from fourier.h)
  /// Note that the returned quantities represent the bare fft sum, with no normalization constants
  void dft(std::vector<double>& data);
  void idft(std::vector<double>& data);
  void realdft(std::vector<double>& data);
  std::vector<double> convlv(const std::vector<double>& data, const std::vector<double>& respns, const int isign);
  
}

#endif // FFT_HPP
