#ifndef DOXYGEN

#ifdef __restrict
#define restrict __restrict
#endif

extern "C" {
  #include <stdlib.h>
  #include <stdio.h>
  #include <math.h>
  #include <complex.h>
  #include "fftw3.h"
  #include "alm.h"
  #include "wigner_d_halfpi.h"
  #include "spinsfast_forward.h"
  #include "spinsfast_backward.h"
}

#endif // DOXYGEN

/// Pointwise multiply this object by another Waveform object
template <typename Op>
GWFrames::Waveform GWFrames::Waveform::BinaryOp(const GWFrames::Waveform& B) const {
  const Waveform& A = *this;

  if(A.NTimes() != B.NTimes()) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nError: Asking for the product of two Waveform objects with different time data."
              << "\n       A.NTimes()=" << A.NTimes() << "\tB.NTimes()=" << B.NTimes()
              << "\n       Interpolate to a common set of times first.\n"
              << std::endl;
    throw(GWFrames_MatrixSizeMismatch);
  }

  if(A.frameType != GWFrames::Inertial || B.frameType != GWFrames::Inertial) {
    if(A.frameType != B.frameType) {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
                << "\nError: Asking for the pointwise product of Waveforms in " << GWFrames::WaveformFrameNames[A.frameType]
                << " and " << GWFrames::WaveformFrameNames[B.frameType] << " frames."
                << "\n       This should only be applied to Waveforms in the same frame.\n"
                << std::endl;
      throw(GWFrames_WrongFrameType);
    } else if(A.frame.size() != B.frame.size()) {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
                << "\nError: Asking for the pointwise product of Waveforms with " << A.frame.size() << " and " << B.frame.size() << " frame data points."
                << "\n       This should only be applied to Waveforms in the same frame.\n"
                << std::endl;
      throw(GWFrames_WrongFrameType);
    }
  }

  // This will be the new object holding the multiplied data
  GWFrames::Waveform C;

  // The new spin weight is the sum of the old ones
  C.spinweight = A.spinweight + B.spinweight;

  // Store both old histories in C's
  C.history << "### *this = A*B\n"
            << "#### A.history.str():\n" << A.history.str()
            << "#### B.history.str():\n" << B.history.str()
            << "#### End of old histories from `A*B`" << std::endl;

  // Just copy other data from A
  C.t = A.t;
  C.frame = A.frame;
  C.frameType = A.frameType;
  C.dataType = A.dataType;
  C.rIsScaledOut = A.rIsScaledOut;
  C.mIsScaledOut = A.mIsScaledOut;

  // Determine the ranges of l that the output should have
  int lMinA = std::abs(A.SpinWeight());
  int lMinB = std::abs(B.SpinWeight());
  int lMin = std::abs(C.SpinWeight());
  int lMax = lMin;
  {
    int lMaxA = lMin;
    for(unsigned int i=0; i<A.NModes(); ++i) {
      if(A.lm[i][0]>lMaxA) { lMaxA = A.lm[i][0]; }
    }
    int lMaxB = lMin;
    for(unsigned int i=0; i<B.NModes(); ++i) {
      if(B.lm[i][0]>lMaxB) { lMaxB = B.lm[i][0]; }
    }
    lMax = (lMaxA>lMaxB ? lMaxB : lMaxA); // Take the smaller
  }
  int Nlm = N_lm(lMax);

  // Set the output lm data
  C.lm = std::vector<std::vector<int> >(lMax*(2+lMax)-lMin*lMin+1, std::vector<int>(2,0));
  {
    unsigned int i=0;
    for(int l=lMin; l<=lMax; ++l) {
      for(int m=-l; m<=l; ++m) {
        C.lm[i][0] = l;
        C.lm[i][1] = m;
        ++i;
      }
    }
  }

  // These numbers determine the equi-angular grid on which we will do
  // the pointwise multiplication.  For best accuracy, have N_phi>
  // 2*lMax and N_theta > 2*lMax; but for speed, don't make them much
  // greater.
  int N_phi = 2*lMax + 1;
  int N_theta = 2*lMax + 1;

  // These will be work arrays
  const std::complex<double> I(0.0,1.0);
  const std::complex<double> zero(0.0,0.0);
  std::vector<std::complex<double> > almA(Nlm);
  std::vector<std::complex<double> > almB(Nlm);
  std::vector<std::complex<double> > almC(Nlm);

  // Now, loop through each time step doing the work
  C.data.resize(C.lm.size(), C.t.size());
  for(unsigned int i_t=0; i_t<C.t.size(); ++i_t) {
    std::vector<std::complex<double> > fA(N_phi*N_theta, zero);
    std::vector<std::complex<double> > fB(N_phi*N_theta, zero);
    std::vector<std::complex<double> > fC(N_phi*N_theta, zero);

    { // Set the a_lm coefficients of A
      unsigned int i=0;
      for(int l=0; l<lMinA; ++l) {
        for(int m=-l; m<=-l; ++m) {
          almA[i] = zero;
          ++i;
        }
      }
      for(int l=lMinA; l<=lMax; ++l) {
        for(int m=-l; m<=-l; ++m) {
          const unsigned int iA = A.FindModeIndexWithoutError(l, m);
          almA[i] = A.Data(iA, i_t);
          ++i;
        }
      }
    }

    { // Set the a_lm coefficients of B
      unsigned int i=0;
      for(int l=0; l<lMinB; ++l) {
        for(int m=-l; m<=-l; ++m) {
          almB[i] = zero;
          ++i;
        }
      }
      for(int l=lMinB; l<=lMax; ++l) {
        for(int m=-l; m<=-l; ++m) {
          const unsigned int iB = B.FindModeIndexWithoutError(l, m);
          almB[i] = B.Data(iB, i_t);
          ++i;
        }
      }
    }

    { // Transform each and multiply pointwise
      spinsfast_salm2map(reinterpret_cast<fftw_complex*>(&almA[0]),
                         reinterpret_cast<fftw_complex*>(&fA[0]),
                         A.SpinWeight(), N_theta, N_phi, lMax);
      spinsfast_salm2map(reinterpret_cast<fftw_complex*>(&almB[0]),
                         reinterpret_cast<fftw_complex*>(&fB[0]),
                         B.SpinWeight(), N_theta, N_phi, lMax);
      for(int i=0; i<N_phi*N_theta; ++i) {
        fC[i] = Op()(fA[i], fB[i]);
      }
    }

    // Transform back and record the new data in C
    spinsfast_map2salm(reinterpret_cast<fftw_complex*>(&fC[0]),
                       reinterpret_cast<fftw_complex*>(&almC[0]),
                       C.SpinWeight(), N_theta, N_phi, lMax);
    for(unsigned int i_m=0; i_m<C.NModes(); ++i_m) {
      C.SetData(i_m, i_t, almC[i_m+lMin*lMin]);
    }

  } // Finish loop over time

  return C;
}
