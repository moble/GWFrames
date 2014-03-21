#include "NumericalRecipes.hpp"

#include "NoiseCurves.hpp"

#include "VectorFunctions.hpp"
#include "Interpolate.hpp"

namespace WU = WaveformUtilities;
using std::vector;
using std::string;
using std::max;
using std::min;
using std::cerr;
using std::endl;

vector<double> AdvLIGO_NSNSOptimal(const vector<double>& F, const bool Invert=false, const double NoiseFloor=0.0) {
  const double FMin = max(NoiseFloor, WU::AdvLIGOSeismicWall);
  const double FMax = 8192;
  vector<double> PSD(F.size(), 0.0);
  if(Invert) {
    for(unsigned int i=0; i<F.size(); ++i) {
      if(fabs(F[i])<FMin || fabs(F[i])>FMax) { PSD[i] = 0.0; }
      else {
	PSD[i] = 1.0e49 / (pow(10, -4*(fabs(F[i])-7.9)*(fabs(F[i])-7.9)+16)
			+ 2.4e-62*pow(fabs(F[i])/215.0,-50)
			+ 0.08*pow(fabs(F[i])/215.0,-4.69)
			+ 123.35*(1.0 - 0.23*(fabs(F[i])/215)*(fabs(F[i])/215) + 0.0764*(fabs(F[i])/215)*(fabs(F[i])/215)*(fabs(F[i])/215)*(fabs(F[i])/215)) / (1.0 + 0.17*(fabs(F[i])/215)*(fabs(F[i])/215)) );
      }
    }
  } else {
    for(unsigned int i=0; i<F.size(); ++i) {
      if(fabs(F[i])<FMin || fabs(F[i])>FMax) { PSD[i] = numeric_limits<double>::infinity(); }
      else {
	PSD[i] = 1.0e-49*(pow(10, -4*(fabs(F[i])-7.9)*(fabs(F[i])-7.9)+16)
			  + 2.4e-62*pow(fabs(F[i])/215.0,-50)
			  + 0.08*pow(fabs(F[i])/215.0,-4.69)
			  + 123.35*(1.0-0.23*(fabs(F[i])/215)*(fabs(F[i])/215)+ 0.0764*(fabs(F[i])/215)*(fabs(F[i])/215)*(fabs(F[i])/215)*(fabs(F[i])/215)) / (1.0 +0.17*(fabs(F[i])/215)*(fabs(F[i])/215)));
      }
    }
  }
  return PSD;
}

vector<double> AdvLIGO_ZeroDet_HighP(const vector<double>& F, const bool Invert=false, const double NoiseFloor=0.0) {
  #include "AdvLIGO_ZeroDet_HighP.ipp"
  vector<double> LogPSD;
  LogPSD = WU::Interpolate(ZERO_DET_high_PLogF, ZERO_DET_high_PLogPSD, log(fabs(F)));
  const double MinFreq(max(NoiseFloor, WU::AdvLIGOSeismicWall));
  for(unsigned int i=0; i<LogPSD.size(); ++i) if(fabs(F[i])<MinFreq) { LogPSD[i] = 500.0; }
  if(Invert) { return exp(-1.0*LogPSD); }
  return exp(LogPSD);
}
vector<double> AdvLIGO_ZeroDet_LowP(const vector<double>& F, const bool Invert=false, const double NoiseFloor=0.0) {
  #include "AdvLIGO_ZeroDet_LowP.ipp"
  vector<double> LogPSD(WU::Interpolate(ZERO_DET_low_PLogF, ZERO_DET_low_PLogPSD, log(fabs(F))));
  const double MinFreq(max(NoiseFloor, WU::AdvLIGOSeismicWall));
  for(unsigned int i=0; i<LogPSD.size(); ++i) if(fabs(F[i])<MinFreq) { LogPSD[i] = 500.0; }
  if(Invert) { return exp(-1.0*LogPSD); }
  return exp(LogPSD);
}

vector<double> IniLIGO_Approx(const vector<double>& F, const bool Invert=false, const double NoiseFloor=0.0) {
  const double FMin = max(NoiseFloor, WU::IniLIGOSeismicWall);
  const double FMax = WU::IniLIGOSamplingFreq;
  vector<double> PSD(F.size(), 0.0);
  if(Invert) {
    for(unsigned int i=0; i<F.size(); ++i) {
      const double f = fabs(F[i]);
      if(f<FMin || f>FMax) { PSD[i] = 0.0; }
      else {
 	const double x = f / 150.0;
	//PSD[i] = 1.0 / (9e-46 * (pow(4.49*x, -56) + 0.16*pow(x, -4.52) + 0.32*pow(x,2) + 0.52)); // Table IV of PRD 63 044023
	PSD[i] = 1.0 / (3.136e-46 * (pow(4.49*x, -56) + 0.16*pow(x, -4.52) + 0.32*pow(x,2) + 0.52)); // Eq. (10) of CQG 26 (2009) 114006
      }
    }
  } else {
    for(unsigned int i=0; i<F.size(); ++i) {
      const double f = fabs(F[i]);
      if(f<FMin || f>FMax) { PSD[i] = numeric_limits<double>::infinity(); }
      else {
 	const double x = f / 150.0;
	//PSD[i] = 9e-46 * (pow(4.49*x, -56) + 0.16*pow(x, -4.52) + 0.32*pow(x,2) + 0.52); // Table IV of PRD 63 044023
	PSD[i] = 3.136e-46 * (pow(4.49*x, -56) + 0.16*pow(x, -4.52) + 0.32*pow(x,2) + 0.52); // Eq. (10) of CQG 26 (2009) 114006
      }
    }
  }
  return PSD;
}

vector<double> WU::NoiseCurve(const vector<double>& F, const string& Detector, const bool Invert, const double NoiseFloor) {
  if(Detector.compare("AdvLIGO_NSNSOptimal")==0) { return AdvLIGO_NSNSOptimal(F, Invert, NoiseFloor); }
  else if(Detector.compare("AdvLIGO_ZeroDet_HighP")==0) { return AdvLIGO_ZeroDet_HighP(F, Invert, NoiseFloor); }
  else if(Detector.compare("AdvLIGO_ZeroDet_LowP")==0) { return AdvLIGO_ZeroDet_LowP(F, Invert, NoiseFloor); }
  else if(Detector.compare("IniLIGO_Approx")==0) { return IniLIGO_Approx(F, Invert, NoiseFloor); }
  else { cerr << "\nDetector type: '" << Detector << "'" << endl;  Throw1WithMessage("Unknown detector"); }
}

vector<double> WU::InverseNoiseCurve(const vector<double>& F, const string& Detector, const double NoiseFloor) {
  return NoiseCurve(F, Detector, true, NoiseFloor);
}
