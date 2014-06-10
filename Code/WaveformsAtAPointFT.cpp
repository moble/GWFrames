#include "WaveformsAtAPointFT.hpp"

#include "NoiseCurves.hpp"
#include "Utilities.hpp"
#include "fft.hpp"
#include "Errors.hpp"

#include "Waveforms.hpp"
#include <complex>

namespace WU = WaveformUtilities;
using std::vector;
using std::cerr;
using std::endl;
using std::ios_base;
using std::complex;
using GWFrames::operator*;
using GWFrames::operator-;

namespace {

double sqr(const double& a) { return a*a; }
double cube(const double& a) { return a*a*a; }

double BumpFunction(const double x, const double x0, const double x1) {
  if(x<=x0) { return 0.0; }
  if(x>=x1) { return 1.0; }
  const double t = (x-x0)/(x1-x0);
  return 1.0 / (1.0 + exp(1.0/t - 1.0/(1-t)));
}

double DoubleSidedF(const unsigned int i, const unsigned int N, const double df) {
  if(i<N/2) { return i*df; }
  return i*df - N*df;
}

}


GWFrames::WaveformAtAPointFT::WaveformAtAPointFT(const GWFrames::Waveform& W,
                                       const double Dt,
                                       const double Vartheta,
                                       const double Varphi,
                                       const double TotalMass,
                                       const unsigned int WindowNCycles,
                                       const double DetectorResponseAmp,
                                       const double DetectorResponsePhase)
  : mDt(Dt), mVartheta(Vartheta), mVarphi(Varphi), mNormalized(false)
{

  // Interpolate to an even time spacing dt whose size is the next power of 2
  const unsigned int N1 = (unsigned int)(std::floor((W.T().back()-W.T(0))/Dt));
  const unsigned int N2 = (unsigned int)(std::pow(2.0,ceil(log2(N1))));
  vector<double> NewTimes(N2);
  for(unsigned int i=0; i<N2; ++i) {
    NewTimes[i] = W.T(0) + i*Dt;
  }

  // TODO: expects a waveform in inertial frame
  GWFrames::Waveform W2 = W.Interpolate(NewTimes, true);
  const vector<complex<double> > ComplexHData = W2.EvaluateAtPoint(Vartheta, Varphi);

  // Construct initial real,imag H as a function of time
  vector<double> InitRealT(ComplexHData.size());
  vector<double> InitImagT(ComplexHData.size());
  for (size_t i=0; i<ComplexHData.size(); i++) {
    InitRealT[i] = ComplexHData[i].real();
    InitImagT[i] = ComplexHData[i].imag();
  }



  // Account for detector response in H as a function of time
  vector<double> RealT;
  if(DetectorResponsePhase != 0.0) {
    RealT = InitRealT*(DetectorResponseAmp*cos(DetectorResponsePhase))
            -InitImagT*(DetectorResponseAmp*sin(DetectorResponsePhase));
  } else if(DetectorResponseAmp != 1.0) {
    RealT = InitRealT*DetectorResponseAmp;
  } else {
    RealT = InitRealT;
  }


  {
    // Zero up to the first zero crossing for continuity
    unsigned int i=0;
    const double Sign = RealT[0] / std::abs(RealT[0]);
    while(RealT[i++]*Sign>0) { }
    for(unsigned int j=0; j<i; ++j) {
      RealT[j] = 0.0;
    }
    // Now find the following 2*N zero crossings
    const unsigned int i0 = i;
    const double t0 = NewTimes[i0];
    for(unsigned int j=0; j<WindowNCycles; ++j) {
      while(RealT[i++]*Sign<0) { }
      while(RealT[i++]*Sign>0) { }
    }
    // And window the data
    const unsigned int i1 = i;
    const double t1 = NewTimes[i1];
    for(unsigned int j=i0; j<=i1; j++) {
      RealT[j] *= BumpFunction(NewTimes[j], t0, t1);
    }
  }

  // Set up the frequency domain (in Hz)
  {
    // Conversion of time to seconds from geometric units
    // Fundamental constants
    const double G    = 6.67259e-11;     // Units of m^3 kg^-1 s^-2
    const double c    = 299792458;       // Units of m / s
    const double Msol = 1.98892e30;      // Units of kg
    const double TotalMassInSeconds = TotalMass * Msol * G / cube(c);
    const vector<double> PhysicalTimes = TotalMassInSeconds * NewTimes;
    mFreqs = WU::TimeToPositiveFrequencies(PhysicalTimes);
  }

  // Construct real,imag H as a function of frequency
  // The return from realdft needs to be multiplied by dt to correspond to the continuum FT
  WU::realdft(RealT);
  if (mFreqs.size() != RealT.size()/2+1) {
    cerr << "Time and frequency data size mismatch: "
         << mFreqs.size() << "," << RealT.size()/2+1 << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  mRealF.resize(mFreqs.size());
  mImagF.resize(mFreqs.size());
  for(unsigned int i=0; i<RealT.size()/2; ++i) {
    mRealF[i] = Dt*RealT[2*i];
    mImagF[i] = Dt*RealT[2*i+1];
  }
  // Sort out some funky storage
  mRealF.back() = 0.0; // RealT[1]; // just ignore data at the Nyquist frequency
  mImagF.back() = 0.0;
  mImagF[0] = 0.0;
}

namespace GWFrames {
WaveformAtAPointFT& WaveformAtAPointFT::Normalize(const vector<double>& InversePSD)
{
  if(mNormalized) { return *this; }
  const double snr = SNR(InversePSD);
  for (size_t i=0; i<mRealF.size(); i++) {
    mRealF[i] /= snr;
    mImagF[i] /= snr;
  }
  mNormalized = true;
  return *this;
}

WaveformAtAPointFT& WaveformAtAPointFT::Normalize(const std::string& Detector)
{
  return Normalize(WU::NoiseCurve(F(), Detector, true));
}

WaveformAtAPointFT& WaveformAtAPointFT::ZeroAbove(const double Frequency)
{
  for(unsigned int f=0; f<NFreq(); ++f) {
    if(std::fabs(F(f)) > Frequency) {
      mRealF[f] = 0.0;
      mImagF[f] = 0.0;
    }
  }
  return *this;
}

// TODO: Assert same df, not just NFreq()
double WaveformAtAPointFT::InnerProduct(const WaveformAtAPointFT& B,
                                        const vector<double>& InversePSD) const
{
  if(NFreq() != B.NFreq()) {
    cerr << "\nthis->NFreq()=" << NFreq() << "\tB.NFreq()=" << B.NFreq() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  if(NFreq() != InversePSD.size()) {
    cerr << "\nthis->NFreq()=" << NFreq() << "\tInversePSD.size()=" << InversePSD.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  double InnerProduct = 0.0;
  for(unsigned int i=0; i<NFreq(); ++i) {
    InnerProduct += (Re(i)*B.Re(i)+Im(i)*B.Im(i))*InversePSD[i];
  }
  InnerProduct = 4*(F(1)-F(0))*InnerProduct; // Remember: single-sided frequency
  return InnerProduct;
}

vector<double> WaveformAtAPointFT::InversePSD(const std::string& Detector) const
{
  return WU::NoiseCurve(F(), Detector, true);
}

double WaveformAtAPointFT::SNR(const std::vector<double>& InversePSD) const
{
  /// \param[in] InversePSD Noise spectrum
  if(NFreq() != InversePSD.size()) {
    cerr << "\nthis->NFreq()=" << NFreq() << "\tInversePSD.size()=" << InversePSD.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  double SNRSquared = 0.0;
  for(unsigned int i=0; i<NFreq(); ++i) {
    SNRSquared += (Re(i)*Re(i)+Im(i)*Im(i))*InversePSD[i];
  }
  SNRSquared = 4*(F(1)-F(0))*SNRSquared; // Remember: single-sided frequency
  return std::sqrt(SNRSquared);
}

double WaveformAtAPointFT::SNR(const std::string& Detector) const
{
  /// \param[in] Detector Noise spectrum from this detector
  return SNR(WU::NoiseCurve(F(), Detector, true));
}

/// Compute the match between two WaveformAtAPointFT
void WaveformAtAPointFT::Match(const WaveformAtAPointFT& B,
                               const std::vector<double>& InversePSD,
                               double& timeOffset, double& phaseOffset,
                               double& match) const
{
  /// \param[in] B WaveformAtAPointFT to compute match with
  /// \param[in] InversePSD Spectrum used to weight contributions by frequencies to match
  /// \param[out] timeOffset Time offset used between the waveforms
  /// \param[out] phaseOffset Phase offset used between the waveforms
  /// \param[out] match Match between the two waveforms
  const unsigned int n = NFreq(); // Only positive frequencies are stored in t
  const unsigned int N = 2*(n-1);  // But this is how many there really are
  if(!IsNormalized() || !B.IsNormalized()) {
    cerr << "\n\nWARNING!!! Matching non-normalized WaveformAtAPointFT objects. WARNING!!!\n" << endl;
  }
  if(n != B.NFreq() || n != InversePSD.size()) {
    cerr << "Waveform sizes, " << n << " and " << B.NFreq()
         << ", are not compatible with InversePSD size, " << InversePSD.size() << "." << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  const double eps = 1e-8;
  const double df = F(1)-F(0);
  const double df_B = B.F(1)-B.F(0);
  const double rel_diff_df = std::fabs(1 - df/df_B);
  if(rel_diff_df > eps) {
    cerr << "Waveform frequency steps, " << df << " and " << df_B
         << ", are not compatible in Match: rel_diff="<< rel_diff_df << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  // s1 s2* = (a1 + i b1) (a2 - i b2) = (a1 a2 + b1 b2) + i(b1 a2 - a1 b2)
  WU::WrapVecDoub data(2*N);
  for(unsigned int i=0; i<n; ++i) {
    data.real(i) = (Re(i)*B.Re(i)+Im(i)*B.Im(i))*InversePSD[i];
    data.imag(i) = (Im(i)*B.Re(i)-Re(i)*B.Im(i))*InversePSD[i];
  }
  idft(data);
  unsigned int maxi=0;
  double maxmag = std::sqrt(sqr(data.real(0)) + sqr(data.imag(0)));
  for(unsigned int i=1; i<N; ++i) {
    const double mag = std::sqrt(sqr(data.real(i)) + sqr(data.imag(i)));
    if(mag>maxmag) { maxmag = mag; maxi = int(i); }
  }
  // note: assumes N is even
  timeOffset = (maxi<N/2 ? double(maxi)/(N*df) : double(maxi-int(N))/(N*df));
  phaseOffset = atan2(data.imag(maxi), data.real(maxi))/2.0;
  /// The return from ifft is just the bare FFT sum, so we multiply by df to get
  /// the continuum-analog FT.  This is correct because the input data (re,im) are
  /// the continuum-analog data, rather than just the return from the bare FFT sum.
  /// See, e.g., Eq. (A.33) [rather than Eq. (A.35)] of
  /// http://etd.caltech.edu/etd/available/etd-01122009-143851/
  match = 4.0*df*maxmag;
  return;
}

/// Compute the match between two WaveformAtAPointFT
void WaveformAtAPointFT::Match(const WaveformAtAPointFT& B, double& timeOffset,
                               double& phaseOffset, double& match,
                               const std::string& Detector) const
{
  Match(B, WU::NoiseCurve(F(), Detector, true), timeOffset, phaseOffset, match);
  return;
}

/// Compute the match between two WaveformAtAPointFT
double WaveformAtAPointFT::Match(const WaveformAtAPointFT& B,
                                 const std::vector<double>& InversePSD) const
{
  double match, timeOffset, phaseOffset;
  Match(B, InversePSD, timeOffset, phaseOffset, match);
  return match;
}

/// Compute the match between two WaveformAtAPointFT
double WaveformAtAPointFT::Match(const WaveformAtAPointFT& B,
                                 const std::string& Detector) const
{
  return Match(B, WU::NoiseCurve(F(), Detector, true));
}

}
