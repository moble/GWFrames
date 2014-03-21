#include "NumericalRecipes.hpp"

#include "WaveformAtAPointFT.hpp"

#include "NoiseCurves.hpp"
#include "VectorFunctions.hpp"
#include "fft.hpp"
#include "Fit.hpp"

namespace WU = WaveformUtilities;
using namespace WaveformUtilities;
using namespace WaveformObjects;
using std::vector;
using std::cerr;
using std::endl;
using std::ios_base;

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


WaveformAtAPointFT::WaveformAtAPointFT()
  : Normalized(false)
{ }

WaveformAtAPointFT::WaveformAtAPointFT(const WaveformAtAPoint& W, const unsigned int WindowNCycles,
				       const double DetectorResponseAmp, const double DetectorResponsePhase)
  : WaveformAtAPoint(W), Normalized(false)
{
  // Record that this is happening
  History() << "### WaveformAtAPointFT(W, " << DetectorResponseAmp << ", " << DetectorResponsePhase << ");" << endl;
  
  // Set up the frequency data
  const double dt = W.T(1)-W.T(0);
  TRef() = TimeToPositiveFrequencies(W.T());
  RRef().resize(NTimes());
  
  // Now get the actual data
  vector<double> RealT;
  if(DetectorResponsePhase!=0.0) {
    RealT = W.Re()*(DetectorResponseAmp*cos(DetectorResponsePhase)) - W.Im()*(DetectorResponseAmp*sin(DetectorResponsePhase));
  } else if(DetectorResponseAmp!=1.0) {
    RealT = DetectorResponseAmp*W.Re();
  } else {
    RealT = W.Re();
  }
  
  {
    // Zero up to the first zero crossing for continuity
    unsigned int i=0;
    const double Sign = RealT[0] / abs(RealT[0]);
    while(RealT[i++]*Sign>0) { }
    for(unsigned int j=0; j<i; ++j) {
      RealT[j] = 0.0;
    }
    // Now find the following 2*N zero crossings
    const unsigned int i0 = i;
    const double t0 = W.T(i0);
    for(unsigned int j=0; j<WindowNCycles; ++j) {
      while(RealT[i++]*Sign<0) { }
      while(RealT[i++]*Sign>0) { }
    }
    // And window the data
    const unsigned int i1 = i;
    const double t1 = W.T(i1);
    for(unsigned int j=i0; j<=i1; j++) {
      RealT[j] *= BumpFunction(W.T(j), t0, t1);
    }
  }
  
  // Do the actual work
  realdft(RealT);
  // The return from realdft needs to be multiplied by dt to correspond to the continuum FT
  ReRef().resize(NTimes());
  ImRef().resize(NTimes());
  for(unsigned int i=0; i<RealT.size()/2; ++i) {
    ReRef(i) = dt*RealT[2*i];
    ImRef(i) = dt*RealT[2*i+1];
  }
  // Sort out some funky storage
  ReRef().back() = 0.0; // RealT[1]; // just ignore data at the Nyquist frequency
  ImRef().back() = 0.0;
  ImRef(0) = 0.0; 
}

WaveformAtAPointFT& WaveformAtAPointFT::Normalize(const std::vector<double>& InversePSD) {
  if(Normalized) { return *this; }
  const double snr = SNR(InversePSD);
  ReRef() /= snr;
  ImRef() /= snr;
  Normalized = true;
  return *this;
}

WaveformAtAPointFT& WaveformAtAPointFT::ZeroAbove(const double Frequency) {
  for(unsigned int f=0; f<NTimes(); ++f) {
    if(fabs(F(f))>Frequency) { ReRef(f) = 0.0; ImRef(f) = 0.0; }
  }
  return *this;
}

double WaveformAtAPointFT::InnerProduct(const WaveformAtAPointFT& B, const std::vector<double>& InversePSD) const {
  if(NTimes() != B.NTimes()) {
    cerr << "\nthis->NTimes()=" << NTimes() << "\tB.NTimes()=" << B.NTimes() << endl;
    Throw1WithMessage("Incompatible sizes");
  }
  if(NTimes() != InversePSD.size()) {
    cerr << "\nWaveform size=" << NTimes() << "\tInversePSD.size()=" << InversePSD.size() << endl;
    Throw1WithMessage("Incompatible sizes");
  }
  double InnerProduct = 0.0;
  for(unsigned int i=0; i<NTimes(); ++i) {
    InnerProduct += (Re(i)*B.Re(i)+Im(i)*B.Im(i))*InversePSD[i];
  }
  InnerProduct = 4*(F(1)-F(0))*InnerProduct; // Remember: single-sided frequency
  return InnerProduct;
}

double WaveformAtAPointFT::SNR(const std::vector<double>& InversePSD) const {
  if(NTimes() != InversePSD.size()) {
    cerr << "\nWaveform size=" << NTimes() << "\tInversePSD.size()=" << InversePSD.size() << endl;
    Throw1WithMessage("Incompatible sizes");
  }
  double SNRSquared = 0.0;
  for(unsigned int i=0; i<NTimes(); ++i) {
//   cerr << "\nDEBUGGING!!! " << __LINE__ << " " << __FILE__ << endl;
//   for(unsigned int i=NTimes()/2; i<NTimes(); ++i) {
    SNRSquared += (Re(i)*Re(i)+Im(i)*Im(i))*InversePSD[i];
  }
  SNRSquared = 4*(F(1)-F(0))*SNRSquared; // Remember: single-sided frequency
  return sqrt(SNRSquared);
}

double WaveformAtAPointFT::Match(const WaveformAtAPointFT& B, const std::vector<double>& InversePSD) const {
  const unsigned int n = NTimes();
  const unsigned int N = 2*(n-1);
//   if(!Normalized || !B.Normalized) {
//     cerr << "\n\nWARNING!!!  Matching non-normalized WaveformAtAPointFT objects.  WARNING!!!\n" << endl;
//   }
  if(n != B.NTimes() || n != InversePSD.size()) {
    cerr << "Waveform sizes, " << n << " and " << B.NTimes()
	 << ", are not compatible with InversePSD size, " << InversePSD.size() << "." << endl;
    Throw1WithMessage("Incompatible sizes");
  }
  const double df = F(1)-F(0);
  if(df != B.F(1)-B.F(0)) {
    cerr << "Waveform frequency steps, " << df << " and " << B.F(1)-B.F(0) << ", are not compatible in Match." << endl;
    Throw1WithMessage("Incompatible resolutions");
  }
  // s1 s2* = (a1 + i b1) (a2 - i b2) = (a1 a2 + b1 b2) + i(b1 a2 - a1 b2)
  WaveformUtilities::WrapVecDoub data(2*N);
  for(unsigned int i=0; i<n; ++i) {
    data.real(i) = (Re(i)*B.Re(i)+Im(i)*B.Im(i))*InversePSD[i];
    data.imag(i) = (Im(i)*B.Re(i)-Re(i)*B.Im(i))*InversePSD[i];
  }
  idft(data);
  double maxmag=sqrt(SQR(data.real(0)) + SQR(data.imag(0)));
  for(unsigned int i=1; i<N; ++i) {
    const double mag = sqrt(SQR(data.real(i)) + SQR(data.imag(i)));
    if(mag>maxmag) { maxmag = mag; }
  }
  /// The return from ifft is just the bare FFT sum, so we multiply by df to get
  /// the continuum-analog FT.  This is correct because the input data (re,im) are
  /// the continuum-analog data, rather than just the return from the bare FFT sum.
  /// See, e.g., Eq. (A.33) [rather than Eq. (A.35)] of
  /// http://etd.caltech.edu/etd/available/etd-01122009-143851/
  return 4.0*df*maxmag;
}

double WaveformAtAPointFT::Match(const WaveformAtAPointFT& B, const std::string& Detector) const {
  return Match(B, NoiseCurve(F(), Detector, true));
}

void WaveformAtAPointFT::Match(const WaveformAtAPointFT& B, const std::vector<double>& InversePSD, double& timeOffset, double& phaseOffset, double& match) const {
  const unsigned int n = NTimes(); // Only positive frequencies are stored in t
  const unsigned int N = 2*(n-1);  // But this is how many there really are
  if(n != B.NTimes() || n != InversePSD.size()) {
    cerr << "Waveform sizes, " << n << " and " << B.NTimes()
	 << ", are not compatible with InversePSD size, " << InversePSD.size() << "." << endl;
    Throw1WithMessage("Incompatible sizes");
  }
  const double df = F(1)-F(0);
  if(df != B.F(1)-B.F(0)) {
    cerr << "Waveform frequency steps, " << df << " and " << B.F(1)-B.F(0) << ", are not compatible in Match." << endl;
    Throw1WithMessage("Incompatible resolutions");
  }
  // s1 s2* = (a1 + i b1) (a2 - i b2) = (a1 a2 + b1 b2) + i(b1 a2 - a1 b2)
  WaveformUtilities::WrapVecDoub data(2*N);
  for(unsigned int i=0; i<n; ++i) {
    data.real(i) = (Re(i)*B.Re(i)+Im(i)*B.Im(i))*InversePSD[i];
    data.imag(i) = (Im(i)*B.Re(i)-Re(i)*B.Im(i))*InversePSD[i];
  }
  idft(data);
  int maxi=0;
  double maxmag=sqrt(SQR(data.real(0)) + SQR(data.imag(0)));
  for(unsigned int i=1; i<N; ++i) {
    const double mag = sqrt(SQR(data.real(i)) + SQR(data.imag(i)));
    if(mag>maxmag) { maxmag = mag; maxi = int(i); }
  }
  timeOffset = (maxi<N/2 ? double(maxi)/(N*df) : double(maxi-int(N))/(N*df));
  // timeOffset = maxi/(N*df);
  phaseOffset = atan2(data.imag(maxi), data.real(maxi))/2.0;
  // std::cout << maxi << " of " << N << "\tdf=" << df << "\ttO=" << timeOffset << "\tpO=" << phaseOffset << "\t" << __FILE__ << ":" << __LINE__ << endl;
  /// The return from ifft is just the bare FFT sum, so we multiply by df to get
  /// the continuum-analog FT.  This is correct because the input data (re,im) are
  /// the continuum-analog data, rather than just the return from the bare FFT sum.
  /// See, e.g., Eq. (A.33) [rather than Eq. (A.35)] of
  /// http://etd.caltech.edu/etd/available/etd-01122009-143851/
  match = 4.0*df*maxmag;
  return;
}

void WaveformAtAPointFT::Match(const WaveformAtAPointFT& B, double& timeOffset, double& phaseOffset, double& match, const std::string& Detector) const {
  Match(B, NoiseCurve(F(), Detector, true), timeOffset, phaseOffset, match);
  return;
}

WaveformAtAPointFT WaveformAtAPointFT::operator-(const WaveformAtAPointFT& b) const {
  WaveformAtAPointFT c(*this);
  c.ReRef() = Re()-b.Re();
  c.ImRef() = Im()-b.Im();
  return c;
}

WaveformAtAPointFT WaveformAtAPointFT::operator*(const double b) const {
  WaveformAtAPointFT c(*this);
  c.ReRef() *= b;
  c.ImRef() *= b;
  return c;
}

std::ostream& operator<<(std::ostream& os, const WaveformAtAPointFT& a) {
  os << a.HistoryStr()
     << "# [1] = " << a.TimeScale() << endl
     << "# [2] = Re{F[" << a.Type() << "(" << a.Vartheta() << "," << a.Varphi() << ")]}" << endl
     << "# [3] = Im{F[" << a.Type() << "(" << a.Vartheta() << "," << a.Vartheta() << ")]}" << endl;
  for(unsigned int i=0; i<a.NTimes(); ++i) {
    os << a.T(i) << " " << a.Re(i) << " " << a.Im(i) << endl;
  }
  return os;
}
