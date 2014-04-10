// Copyright (c) 2014, Michael Boyle
// See LICENSE file for details

// #include <omp.h>

#include <unistd.h>
#include <sys/param.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <climits>
#include <cmath>
#include <functional>
#include <algorithm>
#include <complex>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multifit.h>
#include "Waveforms.hpp"
#include "Utilities.hpp"
#include "Quaternions.hpp"
#include "SphericalFunctions/Quaternions/Utilities.hpp"
#include "IntegrateAngularVelocity.hpp"
#include "SphericalFunctions/SWSHs.hpp"
#include "Errors.hpp"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#include "SpacetimeAlgebra.hpp"
#pragma clang diagnostic pop

using Quaternions::Quaternion;
using Quaternions::QuaternionArray;
using GWFrames::Matrix;
using SphericalFunctions::LadderOperatorFactorFunctor;
using GWFrames::abs;
using Quaternions::PrescribedRotation;
using Quaternions::FrameFromZ;
using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::flush;
using std::endl;
using std::setprecision;
using std::stringstream;
using std::istringstream;
using std::ostream;
using std::ifstream;
using std::ofstream;
using std::min;
using std::max;
using std::ios_base;
using std::complex;


const LadderOperatorFactorFunctor LadderOperatorFactor;

std::string tolower(const std::string& A) {
  string B = A;
  string::iterator it;
  for(it=B.begin(); it<B.end(); it++) {
    *it = tolower(*it);
  }
  return B;
}

std::string StringForm(const std::vector<int>& Lmodes) {
  if(Lmodes.size()==0) { return "[]"; }
  stringstream S;
  S << "[";
  for(unsigned int i=0; i<Lmodes.size()-1; ++i) {
    S << Lmodes[i] << ", ";
  }
  S << Lmodes[Lmodes.size()-1] << "]";
  return S.str();
}

/// Default constructor for an empty object
GWFrames::Waveform::Waveform() :
  spinweight(-2), boostweight(-1), history(""), t(0),frame(0), frameType(GWFrames::UnknownFrameType),
  dataType(GWFrames::UnknownDataType), rIsScaledOut(false), mIsScaledOut(false), lm(), data()
{
  {
    char path[MAXPATHLEN];
    getcwd(path, MAXPATHLEN);
    string pwd = path;
    char host[MAXHOSTNAMELEN];
    gethostname(host, MAXHOSTNAMELEN);
    string hostname = host;
    time_t rawtime;
    time ( &rawtime );
    string date = asctime ( localtime ( &rawtime ) );
    history << "# Code revision (`git rev-parse HEAD` or arXiv version) = " << CodeRevision << endl
            << "# pwd = " << pwd << endl
            << "# hostname = " << hostname << endl
            << "# date = " << date // comes with a newline
            << "Waveform(); // empty constructor" << endl;
  }
}

/// Copy constructor
GWFrames::Waveform::Waveform(const GWFrames::Waveform& a) :
  spinweight(a.spinweight), boostweight(a.boostweight), history(a.history.str()), t(a.t), frame(a.frame), frameType(a.frameType),
  dataType(a.dataType), rIsScaledOut(a.rIsScaledOut), mIsScaledOut(a.mIsScaledOut), lm(a.lm), data(a.data)
{
  /// Simply copies all fields in the input object to the constructed
  /// object, including history
  history.seekp(0, ios_base::end);
}

/// Constructor from data file
GWFrames::Waveform::Waveform(const std::string& FileName, const std::string& DataFormat) :
  spinweight(-2), boostweight(-1), history(""), t(0), frame(0), frameType(GWFrames::UnknownFrameType),
  dataType(GWFrames::UnknownDataType), rIsScaledOut(false), mIsScaledOut(false), lm(), data()
{
  ///
  /// \param FileName Relative path to data file
  /// \param DataFormat Either 'ReIm' or 'MagArg'
  ///
  /// NOTE: This function assumes that the data are stored as (ell,m)
  /// modes, starting with (2,-2), incrementing m, then incrementing
  /// ell and starting again at m=-ell.  If this is not how the modes
  /// are stored, the 'lm' data of this object needs to be reset or
  /// bad things will happen when trying to find the angular-momentum
  /// vector or rotate the waveform.
  {
    char path[MAXPATHLEN];
    getcwd(path, MAXPATHLEN);
    string pwd = path;
    char host[MAXHOSTNAMELEN];
    gethostname(host, MAXHOSTNAMELEN);
    string hostname = host;
    time_t rawtime;
    time ( &rawtime );
    string date = asctime ( localtime ( &rawtime ) );
    history << "# Code revision (`git rev-parse HEAD` or arXiv version) = " << CodeRevision << endl
            << "# pwd = " << pwd << endl
            << "# hostname = " << hostname << endl
            << "# date = " << date // comes with a newline
            << "Waveform(" << FileName << ", " << DataFormat << "); // Constructor from data file" << endl;
  }

  // Get the number of lines in the file
  char LengthChar[9];
  FILE* fp = popen(("wc -l " + FileName).c_str(), "r");
  if(!fp) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Cannot call wc to get size of '" << FileName << "'" << endl;
    throw(GWFrames_FailedSystemCall);
  }
  if(fgets(LengthChar, 9, fp) == NULL) {
    cerr << "\n\n";
    if(ferror(fp)) {
      cerr << __FILE__ << ":" << __LINE__ << ": Read error in popen while calling wc. LengthChar='" << LengthChar << "'" << endl;
    }
    cerr << __FILE__ << ":" << __LINE__ << ": Failed calling wc to get size of '" << FileName << "'" << endl;
    throw(GWFrames_FailedSystemCall);
  }
  pclose(fp);
  const int FileLength = atoi(LengthChar);

  // Open the input file stream
  ifstream ifs(FileName.c_str(), ifstream::in);
  if(!ifs.is_open()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Couldn't open '" << FileName << "'" << endl;
    throw(GWFrames_BadFileName);
  }

  // Get the header and save to 'history'
  int HeaderLines = 0;
  {
    history << "#### Begin Previous History\n";
    string Temp;
    while(ifs.peek() == '#' || ifs.peek() == '%') {
      getline(ifs, Temp);
      HeaderLines++;
      history << "#" << Temp << "\n";
    }
    history << "#### End Previous History\n";
  }

  // Read the complex data
  vector<double> Line;
  unsigned int NTimes = FileLength - HeaderLines;
  double Test = 0.0;
  string LineString;
  const complex<double> ImaginaryI(0.0,1.0);
  getline(ifs, LineString);
  istringstream LineStream(LineString);
  while(LineStream >> Test) Line.push_back(Test);
  const unsigned int NModes = (Line.size()-1)/2;
  data.resize(NModes, NTimes);
  t.resize(NTimes);
  t[0] = Line[0];
  if(tolower(DataFormat).find("reim")!=string::npos) {
    for(unsigned int i_m=0; i_m<NModes; ++i_m) {
      data[i_m][0] = Line[1+2*i_m] + ImaginaryI*Line[2+2*i_m];
    }
    double Re, Im;
    for(unsigned int i_t=1; i_t<NTimes; ++i_t) {
      ifs >> t[i_t];
      for(unsigned int i_m=0; i_m<NModes; ++i_m) {
        ifs >> Re;
        ifs >> Im;
        data[i_m][i_t] = Re + ImaginaryI*Im;
      }
    }
  } else if(tolower(DataFormat).find("magarg")!=string::npos) {
    for(unsigned int i_m=0; i_m<NModes; ++i_m) {
      data[i_m][0] = Line[1+2*i_m]*std::exp(ImaginaryI*Line[2+2*i_m]);
    }
    double Mag, Arg;
    for(unsigned int i_t=1; i_t<NTimes; ++i_t) {
      ifs >> t[i_t];
      for(unsigned int i_m=0; i_m<NModes; ++i_m) {
        ifs >> Mag;
        ifs >> Arg;
        data[i_m][i_t] = Mag*std::exp(ImaginaryI*Arg);
      }
    }
  } else {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Unknown data format '" << DataFormat << "' not yet implemented" << endl;
    throw(GWFrames_NotYetImplemented);
  }

  // Close the file stream
  ifs.close();

  // Set the (ell,m) data
  lm = vector<vector<int> >(NModes, vector<int>(2,0));
  cerr << "Warning: Waveform constructor assumes (ell,m) modes are stored as (2,-2), (2,-1), ...\n";
  {
    unsigned int i_m = 0;
    for(int ell=2; ell<10000; ++ell) { // Note ridiculous upper bound on ell to insure against infinite loops
      for(int m=-ell; m<=ell; ++m) {
        if(i_m>=NModes) { break; }
        lm[i_m][0] = ell;
        lm[i_m][1] = m;
        ++i_m;
      }
    }
  }
}

/// Assignment operator
GWFrames::Waveform& GWFrames::Waveform::operator=(const GWFrames::Waveform& a) {
  spinweight = a.spinweight;
  boostweight = a.boostweight;
  history.str(a.history.str());
  history.clear();
  history.seekp(0, ios_base::end);
  t = a.t;
  frame = a.frame;
  frameType = a.frameType;
  dataType = a.dataType;
  rIsScaledOut = a.rIsScaledOut;
  mIsScaledOut = a.mIsScaledOut;
  lm = a.lm;
  data = a.data;
  return *this;
}

/// Remove all data relating to times outside of the given range
GWFrames::Waveform& GWFrames::Waveform::DropTimesOutside(const double ta, const double tb) {
  history << "this->DropTimesOutside(" << ta << ", " << tb << ");" << std::endl;
  const unsigned int i_a = Quaternions::hunt(t, ta)+1;
  const unsigned int i_b = Quaternions::hunt(t, tb);
  vector<vector<complex<double> > > newdata(NModes(), vector<complex<double> >(t.size()));
  for(unsigned int mode=0; mode<NModes(); ++mode) {
    vector<complex<double> > ModeData = Data(mode);
    ModeData.erase(ModeData.begin()+i_b, ModeData.end());
    ModeData.erase(ModeData.begin(), ModeData.begin()+i_a);
    newdata[mode].swap(ModeData);
  }
  data = MatrixC(newdata);
  frame.erase(frame.begin()+i_b, frame.end());
  frame.erase(frame.begin(), frame.begin()+i_a);
  t.erase(t.begin()+i_b, t.end());
  t.erase(t.begin(), t.begin()+i_a);

  return *this;
}

/// Remove data relating to the given ell modes
GWFrames::Waveform& GWFrames::Waveform::DropEllModes(const std::vector<unsigned int>& EllModesToDrop) {
  history << "this->DropEllModes([" << EllModesToDrop << "]);" << std::endl;
  std::vector<unsigned int> IndicesToKeep(0);
  std::vector<std::vector<int> > newlm(0);
  for(unsigned int i_m=0; i_m<NModes(); ++i_m) {
    bool KeepThis = true;
    for(unsigned int i=0; i<EllModesToDrop.size(); ++i) {
      if(lm[i_m][0] == EllModesToDrop[i]) {
        KeepThis = false;
        break;
      }
    }
    if(KeepThis) {
      IndicesToKeep.push_back(i_m);
      newlm.push_back(lm[i_m]);
    }
  }
  lm = newlm;
  vector<vector<complex<double> > > NewData(IndicesToKeep.size(), vector<complex<double> >(NTimes()));
  for(unsigned int i_m=0; i_m<IndicesToKeep.size(); ++i_m) {
    NewData[i_m] = Data(IndicesToKeep[i_m]);
  }
  data = MatrixC(NewData);
  return *this;
}

/// Remove data relating to all but the given ell modes
GWFrames::Waveform& GWFrames::Waveform::KeepOnlyEllModes(const std::vector<unsigned int>& EllModesToKeep) {
  history << "this->KeepOnlyEllModes([" << EllModesToKeep << "]);" << std::endl;
  std::vector<unsigned int> IndicesToKeep(0);
  std::vector<std::vector<int> > newlm(0);
  for(unsigned int i_m=0; i_m<NModes(); ++i_m) {
    bool KeepThis = false;
    for(unsigned int i=0; i<EllModesToKeep.size(); ++i) {
      if(lm[i_m][0] == EllModesToKeep[i]) {
        KeepThis = true;
        break;
      }
    }
    if(KeepThis) {
      IndicesToKeep.push_back(i_m);
      newlm.push_back(lm[i_m]);
    }
  }
  lm = newlm;
  vector<vector<complex<double> > > NewData(IndicesToKeep.size(), vector<complex<double> >(NTimes()));
  for(unsigned int i_m=0; i_m<IndicesToKeep.size(); ++i_m) {
    NewData[i_m] = Data(IndicesToKeep[i_m]);
  }
  data = MatrixC(NewData);
  return *this;
}

/// Efficiently swap data between two Waveform objects.
void GWFrames::Waveform::swap(GWFrames::Waveform& b) {
  /// This function uses the std::vector method 'swap' which simply
  /// swaps pointers to data, for efficiency.

  // This call should not be recorded explicitly in the history,
  // because the histories are swapped
  { const int NewSpinWeight=b.spinweight; b.spinweight=spinweight; spinweight=NewSpinWeight; }
  { const int NewBoostWeight=b.boostweight; b.boostweight=boostweight; boostweight=NewBoostWeight; }
  { const string historyb=b.history.str(); b.history.str(history.str()); history.str(historyb); }
  history.seekp(0, ios_base::end);
  b.history.seekp(0, ios_base::end);
  t.swap(b.t);
  frame.swap(b.frame);
  { const GWFrames::WaveformFrameType bType=b.frameType; b.frameType=frameType; frameType=bType; }
  { const GWFrames::WaveformDataType bType=b.dataType; b.dataType=dataType; dataType=bType; }
  { const bool brIsScaledOut=b.rIsScaledOut; b.rIsScaledOut=rIsScaledOut; rIsScaledOut=brIsScaledOut; }
  { const bool bmIsScaledOut=b.mIsScaledOut; b.mIsScaledOut=mIsScaledOut; mIsScaledOut=bmIsScaledOut; }
  lm.swap(b.lm);
  data.swap(b.data);
  return;
}


/// Explicit constructor from data
GWFrames::Waveform::Waveform(const std::vector<double>& T, const std::vector<std::vector<int> >& LM,
                             const std::vector<std::vector<std::complex<double> > >& Data)
  : spinweight(-2), boostweight(-1), history(""), t(T), frame(), frameType(GWFrames::UnknownFrameType),
    dataType(GWFrames::UnknownDataType), rIsScaledOut(false), mIsScaledOut(false), lm(LM), data(Data)
{
  /// Arguments are T, LM, Data, which consist of the explicit data.

  // Check that dimensions match (though this is not an exhaustive check)
  if( Data.size()!=0 && ( (t.size() != Data[0].size()) || (Data.size() != lm.size()) ) ) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": t.size()=" << t.size()
         << "\nlm.shape=(" << lm.size() << "," << (lm.size()>0 ? lm[0].size() : 0) << ")"
         << "\nData.shape=(" << Data.size() << "," << (Data.size()>0 ? Data[0].size() : 0) << ")" << endl;
    throw(GWFrames_MatrixSizeMismatch);
  }
}


/// Return time derivative of data
std::vector<std::complex<double> > GWFrames::Waveform::DataDot(const unsigned int Mode) const {
  // TODO: Make this work properly for non-inertial systems
  // TODO: Why did I program it like this?
  const std::complex<double>* D = this->operator()(Mode);
  return ComplexDerivative(vector<std::complex<double> >(D,D+NTimes()), T());
}

/// Differentiate the waveform as a function of time
GWFrames::Waveform& GWFrames::Waveform::Differentiate() {
  if(frameType != GWFrames::Inertial) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nError: Asking to Differentiate a Waveform in the " << GWFrames::WaveformFrameNames[frameType] << " frame."
              << "\n       This is possible, but not yet implemented."
              << "\n       This should only be applied to Waveforms in the " << GWFrames::WaveformFrameNames[GWFrames::Inertial] << " frame.\n"
              << std::endl;
    throw(GWFrames_NotYetImplemented);
  }

  vector<vector<complex<double> > > NewData(NModes(), vector<complex<double> >(NTimes()));
  for(unsigned int i_m=0; i_m<NModes(); ++i_m) {
    NewData[i_m] = DataDot(i_m);
  }
  data = MatrixC(NewData);

  boostweight -= 1;
  if(dataType == GWFrames::h) { dataType = GWFrames::hdot; }
  if(dataType == GWFrames::hdot) { dataType = GWFrames::Psi4; }
  if(dataType == GWFrames::Psi4) { dataType = GWFrames::UnknownDataType; }
  std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
            << "\nWarning: Tracking the Waveform data type is not yet correct.  For Differentiate"
            << "\n resulting in Psi4, the real part of the data should be negated.\n\n" << std::endl;
  history << "this->Differentiate();" << std::endl;

  return *this;
}

/// Return the norm (sum of squares of modes) of the waveform
std::vector<double> GWFrames::Waveform::Norm(const bool TakeSquareRoot) const {
  ///
  /// \param TakeSquareRoot If true, the square root is taken at each instant before returning
  ///
  /// This returns the norm of the waveform, defined as the sum of the
  /// complex norms of the modes.  Note that we are calling this norm
  /// in analogy with the c++ std::complex norm, which is the square
  /// of the absolute value.  However, there is also an option to take
  /// the square root of the data at each time step, which would be
  /// the usual L2 norm of the waveform.
  ///
  /// \sa MaxNormIndex
  /// \sa MaxNormTime
  ///
  vector<double> norm(NTimes(), 0.0);
  for(unsigned int i_t=0; i_t<NTimes(); ++i_t) {
    for(unsigned int i_m=0; i_m<NModes(); ++i_m) {
      norm[i_t] += std::norm(Data(i_m, i_t));
    }
    if(TakeSquareRoot) {
      norm[i_t] = std::sqrt(norm[i_t]);
    }
  }
  return norm;
}

/// Return the data index corresponding to the time of the largest norm
unsigned int GWFrames::Waveform::MaxNormIndex(const unsigned int SkipFraction) const {
  ///
  /// \param SkipFraction Integer fraction of data to skip before looking
  ///
  /// The default value of `SkipFraction` is 4, meaning that we start
  /// looking for the maximum after 1/4th of the data, so as to cut
  /// out junk radiation.  Note that this is integer division, so an
  /// argument of `NTimes()+1` will look through all of the data.
  ///
  /// \sa Norm()
  /// \sa MaxNormTime()
  ///
  const unsigned int StartIndex = NTimes()/SkipFraction;
  const vector<double> norm = Norm(false); // don't bother taking the square root
  unsigned int index = StartIndex;
  double max = norm[index];
  for(unsigned int i_t=StartIndex+1; i_t<NTimes(); ++i_t) {
    if(norm[i_t]>max) {
      index = i_t;
      max = norm[index];
    }
  }
  return index;
}

// Return a descriptive string appropriate for a file name, like rhOverM.
std::string GWFrames::Waveform::DescriptorString() const {
  std::string Descriptor = "";
  if(RIsScaledOut()) Descriptor = "r";
  if(MIsScaledOut()) {
    if(DataType()==UnknownDataType or DataType()==h)
      Descriptor = Descriptor + DataTypeString() + "OverM";
    else if(DataType()==hdot)
      Descriptor = Descriptor + DataTypeString(); // hdot is independent of M
    else if(DataType()==Psi4)
      Descriptor = Descriptor + "M" + DataTypeString();
  } else {
    Descriptor = Descriptor + DataTypeString();
  }
  return Descriptor;
}

/// Return vector of real parts of a given mode as function of time.
std::vector<double> GWFrames::Waveform::Re(const unsigned int Mode) const {
  const unsigned int ntimes = NTimes();
  std::vector<double> re(ntimes);
  for(unsigned int i=0; i<ntimes; ++i) {
    re[i] = Re(Mode,i);
  }
  return re;
}

/// Return vector of imaginary parts of a given mode as function of time.
std::vector<double> GWFrames::Waveform::Im(const unsigned int Mode) const {
  const unsigned int ntimes = NTimes();
  std::vector<double> im(ntimes);
  for(unsigned int i=0; i<ntimes; ++i) {
    im[i] = Im(Mode,i);
  }
  return im;
}

/// Return vector of absolute value of a given mode as function of time.
std::vector<double> GWFrames::Waveform::Abs(const unsigned int Mode) const {
  const unsigned int ntimes = NTimes();
  std::vector<double> abs(ntimes);
  for(unsigned int i=0; i<ntimes; ++i) {
    abs[i] = Abs(Mode,i);
  }
  return abs;
}

/// Return vector of arg of a given mode as function of time.
std::vector<double> GWFrames::Waveform::Arg(const unsigned int Mode) const {
  /// Note that this quantity is not "unwrapped".  That is, the arg is
  /// between -pi and +pi.  To get a smooth, continuous phase in
  /// python, use numpy.unwrap.
  const unsigned int ntimes = NTimes();
  std::vector<double> arg(ntimes);
  for(unsigned int i=0; i<ntimes; ++i) {
    arg[i] = Arg(Mode,i);
  }
  return arg;
}

/// Return vector of complex data of a given mode as function of time.
std::vector<std::complex<double> > GWFrames::Waveform::Data(const unsigned int Mode) const {
  const unsigned int ntimes = NTimes();
  std::vector<std::complex<double> > dat(ntimes);
  for(unsigned int i=0; i<ntimes; ++i) {
    dat[i] = data[Mode][i];
  }
  return dat;
}


/// Return vector of vector of real parts of all modes as function of time.
std::vector<std::vector<double> > GWFrames::Waveform::Re() const {
  const unsigned int nmodes = NModes();
  const unsigned int ntimes = NTimes();
  std::vector<std::vector<double> > re(nmodes, std::vector<double>(ntimes));
  for(unsigned int i_m=0; i_m<nmodes; ++i_m) {
    for(unsigned int i_t=0; i_t<ntimes; ++i_t) {
      re[i_m][i_t] = Re(i_m,i_t);
    }
  }
  return re;
}

/// Return vector of vector of imaginary parts of all modes as function of time.
std::vector<std::vector<double> > GWFrames::Waveform::Im() const {
  const unsigned int nmodes = NModes();
  const unsigned int ntimes = NTimes();
  std::vector<std::vector<double> > im(nmodes, std::vector<double>(ntimes));
  for(unsigned int i_m=0; i_m<nmodes; ++i_m) {
    for(unsigned int i_t=0; i_t<ntimes; ++i_t) {
      im[i_m][i_t] = Im(i_m,i_t);
    }
  }
  return im;
}

/// Return vector of vector of absolute value of all modes as function of time.
std::vector<std::vector<double> > GWFrames::Waveform::Abs() const {
  const unsigned int nmodes = NModes();
  const unsigned int ntimes = NTimes();
  std::vector<std::vector<double> > abs(nmodes, std::vector<double>(ntimes));
  for(unsigned int i_m=0; i_m<nmodes; ++i_m) {
    for(unsigned int i_t=0; i_t<ntimes; ++i_t) {
      abs[i_m][i_t] = Abs(i_m,i_t);
    }
  }
  return abs;
}

/// Return vector of vector of arg of all modes as function of time.
std::vector<std::vector<double> > GWFrames::Waveform::Arg() const {
  const unsigned int nmodes = NModes();
  const unsigned int ntimes = NTimes();
  std::vector<std::vector<double> > arg(nmodes, std::vector<double>(ntimes));
  for(unsigned int i_m=0; i_m<nmodes; ++i_m) {
    for(unsigned int i_t=0; i_t<ntimes; ++i_t) {
      arg[i_m][i_t] = Arg(i_m,i_t);
    }
  }
  return arg;
}

/// Return vector of vector of arg of all modes as function of time.
std::vector<std::vector<double> > GWFrames::Waveform::ArgUnwrapped() const {
  const unsigned int nmodes = NModes();
  const unsigned int ntimes = NTimes();
  std::vector<std::vector<double> > uarg(nmodes, std::vector<double>(ntimes));
  for(unsigned int i_m=0; i_m<nmodes; ++i_m) {
    uarg[i_m] = ArgUnwrapped(i_m);
  }
  return uarg;
}

/// Return vector of vector of complex data of all modes as function of time.
std::vector<std::vector<std::complex<double> > > GWFrames::Waveform::Data() const {
  const unsigned int nmodes = NModes();
  const unsigned int ntimes = NTimes();
  std::vector<std::vector<std::complex<double> > > dat(nmodes, std::vector<std::complex<double> >(ntimes));
  for(unsigned int i_m=0; i_m<nmodes; ++i_m) {
    for(unsigned int i_t=0; i_t<ntimes; ++i_t) {
      dat[i_m][i_t] = data[i_m][i_t];
    }
  }
  return dat;
}


/// Return greatest ell value present in the data.
int GWFrames::Waveform::EllMax() const {
  int ell = lm[0][0];
  for(unsigned int i=1; i<NModes(); ++i) {
    if(lm[i][0]>ell) { ell = lm[i][0]; }
  }
  return ell;
}

/// Find index of mode with given (l,m) data.
unsigned int GWFrames::Waveform::FindModeIndex(const int l, const int m) const {
  //ORIENTATION!!! following loop
  for(unsigned int i=0; i<NModes(); ++i) {
    if(lm[i][0]==l && lm[i][1]==m) { return i; }
  }
  cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Can't find (ell,m)=(" << l << ", " << m << ")" << endl;
  throw(GWFrames_WaveformMissingLMIndex);
}

/// Find index of mode with given (l,m) data without the chance of throwing an exception.
unsigned int GWFrames::Waveform::FindModeIndexWithoutError(const int l, const int m) const {
  /// If the requested mode is not present, the returned index is 1
  /// beyond the end of the mode vector.

  // ORIENTATION!!! following loop
  unsigned int i=0;
  for(; i<NModes(); ++i) {
    if(lm[i][0]==l && lm[i][1]==m) { return i; }
  }
  ++i;
  return i;
}

/// Return the contrast in the given mode pair
std::vector<double> GWFrames::Waveform::Contrast(const int L, const int M) const {
  /// \param L \f$ell\f$ value of the mode pair
  /// \param M \f$m\f$ value of the mode pair
  ///
  /// This function just returns the value of the contrast
  /// \f$\kappa^{\ell,m}\f$ defined by Boyle et al. (2014):
  ///
  /// \f{equation*}{ kappa^{\ell,m} = 2 \frac{\lvert h^{\ell,m} \rvert
  /// - \lvert h^{\ell,-m} \rvert} {\lvert h^{\ell,m} \rvert - \lvert
  /// h^{\ell,-m} \rvert}. \f}
  ///
  /// That is, the difference between mode pairs normalized by their
  /// average.

  const unsigned int Size = NTimes();
  const unsigned int i_m1 = FindModeIndex(L,M);
  const unsigned int i_m2 = FindModeIndex(L,-M);
  std::vector<double> contrast(Size);
  for(unsigned int i_t=0; i_t<Size; ++i_t) {
    contrast[i_t] = 2 * (Abs(i_m1,i_t)-Abs(i_m2,i_t)) / (Abs(i_m1,i_t)+Abs(i_m2,i_t));
  }
  return contrast;
}


/// Rotate the physical content of the Waveform by a constant rotor.
GWFrames::Waveform& GWFrames::Waveform::RotatePhysicalSystem(const Quaternions::Quaternion& R_phys) {
  history << "this->RotatePhysicalSystem(" << R_phys << ");" << endl;

  this->TransformModesToRotatedFrame(vector<Quaternion>(1,R_phys.conjugate()));

  // Record the change of frame
  if(frame.size()==0) { // set frame data equal to input data
    frame = vector<Quaternion>(1,R_phys.conjugate());
  } else if(frame.size()==1) { // multiply frame constant by input rotation
    frame[0] = frame[0] * R_phys.conjugate();
  } else { // multiply frame data by input rotation
    frame = frame * R_phys.conjugate();
  }

  return *this;
}

/// Rotate the physical content of the Waveform.
GWFrames::Waveform& GWFrames::Waveform::RotatePhysicalSystem(std::vector<Quaternions::Quaternion> R_phys) {
  ///
  /// \param R_phys Vector of Quaternions by which to rotate
  ///
  /// This rotates the physical system, leaving the coordinates in
  /// place.
  ///
  /// The Waveform's `frame` data records the rotors needed to rotate
  /// the standard (x,y,z) basis into the (X,Y,Z) basis with respect
  /// to which the Waveform modes are decomposed.  If this is not the
  /// first rotation of the frame, we need to be careful about how we
  /// record the total rotation.  Here, we are rotating the physical
  /// system, while leaving fixed the basis with respect to which the
  /// modes are decomposed.  Therefore, the new frame must be the
  /// original `frame` data times \f$\bar{R}_{phys}\f$.
  ///
  /// Note that this function does not change the `frameType`; this is
  /// left to the calling function.
  ///

  history << "this->RotatePhysicalSystem(R_phys); // R_phys=["
          << std::setprecision(16) << R_phys[0];
  if(R_phys.size()>1) {
    history << ", " << R_phys[1] << ", ...";
  }
  history << "]" << endl;

  // We won't use R_phys from now on, so transform to R_physbar
  R_phys = Quaternions::conjugate(R_phys);
  vector<Quaternion>& R_physbar = R_phys;

  this->TransformModesToRotatedFrame(R_physbar);

  // Record the change of frame
  if(frame.size()==0) { // set frame data equal to input data
    frame = R_physbar;
  } else if(frame.size()==1) { // multiply frame constant by input rotation
    frame = frame[0] * R_physbar;
  } else { // multiply frame data by input rotation
    frame = frame * R_physbar;
  }

  return *this;
}

/// Rotate the basis in which this Waveform is measured by a constant rotor.
GWFrames::Waveform& GWFrames::Waveform::RotateDecompositionBasis(const Quaternions::Quaternion& R_frame) {
  history << "this->RotateDecompositionBasis(" << R_frame << ");" << endl;

  this->TransformModesToRotatedFrame(vector<Quaternion>(1,R_frame));

  // Record the change of frame
  if(frame.size()==0) { // set frame data equal to input data
    frame = vector<Quaternion>(1,R_frame);
  } else if(frame.size()==1) { // multiply frame constant by input rotation
    frame[0] = frame[0] * R_frame;
  } else { // multiply frame data by input rotation
    frame = frame * R_frame;
  }

  return *this;
}

/// Rotate the basis in which this Waveform is measured.
GWFrames::Waveform& GWFrames::Waveform::RotateDecompositionBasis(const std::vector<Quaternions::Quaternion>& R_frame) {
  ///
  /// \param R_frame Vector of Quaternions by which to rotate
  ///
  /// This rotates the coordinate basis, leaving the physical system
  /// in place.
  ///
  /// The Waveform's `frame` data records the rotors needed to rotate
  /// the standard (x,y,z) basis into the (X,Y,Z) basis with respect
  /// to which the Waveform modes are decomposed.  If this is not the
  /// first rotation of the frame, we need to be careful about how we
  /// record the total rotation.  Here, we are just composing
  /// rotations, so we need to store R_frame times the original frame
  /// data.
  ///
  /// Note that this function does not change the `frameType`; this is
  /// left to the calling function.
  ///

  if(R_frame.size()==0) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nWarning: You have passed an empty set of rotors to `RotateDecompositionBasis`.  I will"
              << "\n         do nothing, and simply return the Waveform as is.  But this is most unusual!"
              << std::endl;
    return *this;
  }
  if(R_frame.size()!=NTimes()) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nError: (R_frame.size()=" << R_frame.size() << ") != (NTimes()=" << NTimes() << ")" << std::endl;
    throw(GWFrames_VectorSizeMismatch);
  }

  history << "this->RotateDecompositionBasis(R_frame); // R_frame=["
          << std::setprecision(16) << R_frame[0];
  if(R_frame.size()>1) {
    history << ", " << R_frame[1] << ", ...";
  }
  history << "]" << endl;

  this->TransformModesToRotatedFrame(R_frame);

  // Record the change of frame
  if(frame.size()==0) { // set frame data equal to input data
    frame = R_frame;
  } else if(frame.size()==1) { // multiply frame constant by input rotation
    frame = frame[0] * R_frame;
  } else { // multiply frame data by input rotation
    frame = frame * R_frame;
  }

  return *this;
}

/// Rotate the basis in which this Waveform's uncertainties are measured.
GWFrames::Waveform& GWFrames::Waveform::RotateDecompositionBasisOfUncertainties(const std::vector<Quaternions::Quaternion>& R_frame) {
  ///
  /// \param R_frame Vector of Quaternions by which to rotate
  ///
  /// This rotates the coordinate basis, leaving the physical system
  /// in place.
  ///
  /// The Waveform's `frame` data records the rotors needed to rotate
  /// the standard (x,y,z) basis into the (X,Y,Z) basis with respect
  /// to which the Waveform modes are decomposed.  If this is not the
  /// first rotation of the frame, we need to be careful about how we
  /// record the total rotation.  Here, we are just composing
  /// rotations, so we need to store R_frame times the original frame
  /// data.
  ///
  /// Note that this function does not change the `frameType`; this is
  /// left to the calling function.
  ///

  history << "this->RotateDecompositionBasisOfUncertainties(R_frame); // R_frame=["
          << std::setprecision(16) << R_frame[0];
  if(R_frame.size()>1) {
    history << ", " << R_frame[1] << ", ...";
  }
  history << "]" << endl;

  this->TransformUncertaintiesToRotatedFrame(R_frame);

  // Record the change of frame
  if(frame.size()==0) { // set frame data equal to input data
    frame = R_frame;
  } else if(frame.size()==1) { // multiply frame constant by input rotation
    frame = frame[0] * R_frame;
  } else { // multiply frame data by input rotation
    frame = frame * R_frame;
  }

  return *this;
}

/// Rotate modes of the Waveform object.
GWFrames::Waveform& GWFrames::Waveform::TransformModesToRotatedFrame(const std::vector<Quaternion>& R_frame) {
  /// Given a Waveform object, alter the modes stored in this Waveform
  /// so that it measures the same physical field with respect to a
  /// basis rotated with respect to the first by \f$R_{frame}\f$.
  /// This is equivalent to rotating the physical system by
  /// \f$\bar{R}_{frame}\f$, while leaving the basis fixed.
  ///
  /// Note that this private function does not rotate the `frame`
  /// data, which should instead be done by the public functions
  /// RotatePhysicalSystem and RotateDecompositionBasis.  Also, this
  /// does not adjust the `frameType`, which is left to the calling
  /// functions.
  ///

  const int NModes = this->NModes();
  const int NTimes = this->NTimes();

  if(int(R_frame.size())!=NTimes && R_frame.size()!=1) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": R_frame.size()=" << R_frame.size() << "; NTimes()=" << NTimes << endl;
    throw(GWFrames_VectorSizeMismatch);
  }

  // Loop through each mode and do the rotation
  {
    int mode=1;
    for(int l=std::abs(SpinWeight()); l<NModes; ++l) {
      if(NModes<mode) { break; }

      // Use a vector of mode indices, in case the modes are out of
      // order.  This still assumes that we have each l from l=2 up to
      // some l_max, but it's better than assuming that, plus assuming
      // that everything is in order.
      vector<unsigned int> ModeIndices(2*l+1);
      for(int m=-l, i=0; m<=l; ++m, ++i) {
        try {
          ModeIndices[i] = FindModeIndex(l, m);
        } catch(int thrown) {
          cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Incomplete mode information in Waveform; cannot rotate." << endl;
          throw(thrown);
        }
      }

      {

        // Construct the D matrix and data storage
        SphericalFunctions::WignerDMatrix D(R_frame[0]);
        vector<vector<complex<double> > > Ds(2*l+1, vector<complex<double> >(2*l+1));
        vector<complex<double> > Data(2*l+1);

        if(R_frame.size()==1) {
          // Get the Wigner D matrix data just once
          for(int m=-l; m<=l; ++m) {
            for(int mp=-l; mp<=l; ++mp) {
              Ds[mp+l][m+l] = D(l,mp,m);
            }
          }
          // Loop through each time step
          for(int t=0; t<NTimes; ++t) {
            // Store the data for all m' modes at this time step
            for(int mp=-l, i=0; mp<=l; ++mp, ++i) {
              Data[mp+l] = this->operator()(ModeIndices[i], t);
            }
            // Compute the data at this time step for each m
            for(int m=-l, i=0; m<=l; ++m, ++i) {
              data[ModeIndices[i]][t] = 0.0;
              for(int mp=-l; mp<=l; ++mp) { // Sum over m'
                data[ModeIndices[i]][t] += Ds[mp+l][m+l]*Data[mp+l];
              }
            }
          }
        } else {
          // Loop through each time step
          for(int t=0; t<NTimes; ++t) {
            // Get the Wigner D matrix data at this time step
            D.SetRotation(R_frame[t]);
            for(int m=-l; m<=l; ++m) {
              for(int mp=-l; mp<=l; ++mp) {
                Ds[mp+l][m+l] = D(l,mp,m);
              }
            }
            // Store the data for all m' modes at this time step
            for(int mp=-l, i=0; mp<=l; ++mp, ++i) {
              Data[mp+l] = this->operator()(ModeIndices[i], t);
            }
            // Compute the data at this time step for each m
            for(int m=-l, i=0; m<=l; ++m, ++i) {
              data[ModeIndices[i]][t] = 0.0;
              for(int mp=-l; mp<=l; ++mp) { // Sum over m'
                data[ModeIndices[i]][t] += Ds[mp+l][m+l]*Data[mp+l];
              }
            }
          }
        }

      }

      mode += 2*l+1;
    }
  }

  return *this;
}

inline double SQR(const double a) { return a*a; }

/// Rotate modes of the uncertainty of a Waveform object.
GWFrames::Waveform& GWFrames::Waveform::TransformUncertaintiesToRotatedFrame(const std::vector<Quaternion>& R_frame) {
  ///
  /// Assuming the data stored in the Waveform represent uncertainties
  /// (square-root of sigma), rotate those uncertainties, using the
  /// standard method of combining by quadrature.
  ///
  /// \sa TransformModesToRotatedFrame
  ///

  const int NModes = this->NModes();
  const int NTimes = this->NTimes();

  if(int(R_frame.size())!=NTimes && R_frame.size()!=1) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": R_frame.size()=" << R_frame.size() << "; NTimes()=" << NTimes << endl;
    throw(GWFrames_VectorSizeMismatch);
  }

  // Loop through each mode and do the rotation
  {
    int mode=1;
    for(int l=2; l<NModes; ++l) {
      if(NModes<mode) { break; }

      // Use a vector of mode indices, in case the modes are out of
      // order.  This still assumes that we have each l from l=2 up to
      // some l_max, but it's better than assuming that, plus assuming
      // that everything is in order.
      vector<unsigned int> ModeIndices(2*l+1);
      for(int m=-l, i=0; m<=l; ++m, ++i) {
        try {
          ModeIndices[i] = FindModeIndex(l, m);
        } catch(int thrown) {
          cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Incomplete mode information in Waveform; cannot rotate." << endl;
          throw(thrown);
        }
      }

      {

        // Construct the D matrix and data storage
        SphericalFunctions::WignerDMatrix D(R_frame[0]);
        vector<vector<complex<double> > > Ds(2*l+1, vector<complex<double> >(2*l+1));
        vector<complex<double> > Data(2*l+1);

        if(R_frame.size()==1) {
          // Get the Wigner D matrix data just once
          for(int m=-l; m<=l; ++m) {
            for(int mp=-l; mp<=l; ++mp) {
              Ds[mp+l][m+l] = D(l,mp,m);
            }
          }
          // Loop through each time step
          for(int t=0; t<NTimes; ++t) {
            // Store the data for all m' modes at this time step
            for(int mp=-l, i=0; mp<=l; ++mp, ++i) {
              Data[mp+l] = this->operator()(ModeIndices[i], t);
            }
            // Compute the data at this time step for each m
            for(int m=-l, i=0; m<=l; ++m, ++i) {
              data[ModeIndices[i]][t] = 0.0;
              for(int mp=-l; mp<=l; ++mp) { // Sum over m'
                const double dataReSqr = SQR(std::real(Data[mp+l]));
                const double dataImSqr = SQR(std::imag(Data[mp+l]));
                const double DReSqr = SQR(std::real(Ds[mp+l][m+l]));
                const double DImSqr = SQR(std::imag(Ds[mp+l][m+l]));
                data[ModeIndices[i]][t] += complex<double>(dataReSqr*DReSqr+dataImSqr*DImSqr,
                                                           dataImSqr*DReSqr+dataReSqr*DImSqr);
              }
              data[ModeIndices[i]][t] = complex<double>(std::sqrt(std::real(data[ModeIndices[i]][t])),
                                                        std::sqrt(std::imag(data[ModeIndices[i]][t])));
            }
          }
        } else {
          // Loop through each time step
          for(int t=0; t<NTimes; ++t) {
            // Get the Wigner D matrix data at this time step
            D.SetRotation(R_frame[t]);
            for(int m=-l; m<=l; ++m) {
              for(int mp=-l; mp<=l; ++mp) {
                Ds[mp+l][m+l] = D(l,mp,m);
              }
            }
            // Store the data for all m' modes at this time step
            for(int mp=-l, i=0; mp<=l; ++mp, ++i) {
              Data[mp+l] = this->operator()(ModeIndices[i], t);
            }
            // Compute the data at this time step for each m
            for(int m=-l, i=0; m<=l; ++m, ++i) {
              data[ModeIndices[i]][t] = 0.0;
              for(int mp=-l; mp<=l; ++mp) { // Sum over m'
                const double dataReSqr = SQR(std::real(Data[mp+l]));
                const double dataImSqr = SQR(std::imag(Data[mp+l]));
                const double DReSqr = SQR(std::real(Ds[mp+l][m+l]));
                const double DImSqr = SQR(std::imag(Ds[mp+l][m+l]));
                data[ModeIndices[i]][t] += complex<double>(dataReSqr*DReSqr+dataImSqr*DImSqr,
                                                           dataImSqr*DReSqr+dataReSqr*DImSqr);
              }
              data[ModeIndices[i]][t] = complex<double>(std::sqrt(std::real(data[ModeIndices[i]][t])),
                                                        std::sqrt(std::imag(data[ModeIndices[i]][t])));
            }
          }
        }

      }

      mode += 2*l+1;
    }
  }

  return *this;
}


/// Calculate the \f$<L \partial_t>\f$ quantity defined in the paper.
vector<vector<double> > GWFrames::Waveform::LdtVector(vector<int> Lmodes) const {
  ///
  /// \param Lmodes L modes to evaluate
  ///
  /// If Lmodes is empty (default), all L modes are used.  Setting
  /// Lmodes to [2] or [2,3,4], for example, restricts the range of
  /// the sum.
  ///
  /// \f$<L \partial_t>^a = \sum_{\ell,m,m'} \Im [ \bar{f}^{\ell,m'} < \ell,m' | L_a | \ell,m > \dot{f}^{\ell,m} ]\f$

  // L+ = Lx + i Ly      Lx =    (L+ + L-) / 2     Im(Lx) =  ( Im(L+) + Im(L-) ) / 2
  // L- = Lx - i Ly      Ly = -i (L+ - L-) / 2     Im(Ly) = -( Re(L+) - Re(L-) ) / 2
  // Lz = Lz             Lz = Lz                   Im(Lz) = Im(Lz)

  LadderOperatorFactorFunctor LadderOperatorFactor;
  if(Lmodes.size()==0) {
    Lmodes.push_back(lm[0][0]);
    for(unsigned int i_m=0; i_m<NModes(); ++i_m) {
      if(std::find(Lmodes.begin(), Lmodes.end(), lm[i_m][0]) == Lmodes.end() ) {
        Lmodes.push_back(lm[i_m][0]);
      }
    }
  }
  vector<vector<double> > l(NTimes(), vector<double>(3, 0.0));
  for(unsigned int iL=0; iL<Lmodes.size(); ++iL) {
    const int L = Lmodes[iL];
    for(int M=-L; M<=L; ++M) {
      const int iMode = FindModeIndex(L,M);
      const vector<complex<double> > dDdt = DataDot(iMode);
      if(M+1<=L) { // L+
        const int iModep1 = FindModeIndex(L,M+1);
        const double c_ell_posm = LadderOperatorFactor(L, M);
        for(unsigned int iTime=0; iTime<NTimes(); ++iTime) {
          const complex<double> Lplus  = c_ell_posm * conj(data[iModep1][iTime]) * dDdt[iTime];
          l[iTime][0] += 0.5 * imag(Lplus);
          l[iTime][1] -= 0.5 * real(Lplus);
        }
      }
      { // Lz; always evaluate this one
        for(unsigned int iTime=0; iTime<NTimes(); ++iTime) {
          const complex<double> Lz = (conj(data[iMode][iTime]) * dDdt[iTime]) * double(M);
          l[iTime][2] += imag(Lz);
        }
      }
      if(M-1>=-L) { // L-
        const int iModem1 = FindModeIndex(L,M-1);
        const double c_ell_negm = LadderOperatorFactor(L, -M);
        for(unsigned int iTime=0; iTime<NTimes(); ++iTime) {
          const complex<double> Lminus  = c_ell_negm * conj(data[iModem1][iTime]) * dDdt[iTime];
          l[iTime][0] += 0.5 * imag(Lminus);
          l[iTime][1] += 0.5 * real(Lminus);
        }
      }
    }
  }
  return l;
}

/// Calculate the \f$<LL>\f$ quantity defined in the paper.
vector<Matrix> GWFrames::Waveform::LLMatrix(vector<int> Lmodes) const {
  ///
  /// \param Lmodes L modes to evaluate
  ///
  /// If Lmodes is empty (default), all L modes are used.  Setting
  /// Lmodes to [2] or [2,3,4], for example, restricts the range of
  /// the sum.
  ///
  /// \f$<LL>^{ab} = \sum_{\ell,m,m'} [\bar{f}^{\ell,m'} < \ell,m' | L_a L_b | \ell,m > f^{\ell,m} ]\f$

  // Big, bad, ugly, obvious way to do the calculation

  // L+ = Lx + i Ly      Lx =    (L+ + L-) / 2     Im(Lx) =  ( Im(L+) + Im(L-) ) / 2
  // L- = Lx - i Ly      Ly = -i (L+ - L-) / 2     Im(Ly) = -( Re(L+) - Re(L-) ) / 2
  // Lz = Lz             Lz = Lz                   Im(Lz) = Im(Lz)
  // LxLx =   (L+ + L-)(L+ + L-) / 4
  // LxLy = -i(L+ + L-)(L+ - L-) / 4
  // LxLz =   (L+ + L-)(  Lz   ) / 2
  // LyLx = -i(L+ - L-)(L+ + L-) / 4
  // LyLy =  -(L+ - L-)(L+ - L-) / 4
  // LyLz = -i(L+ - L-)(  Lz   ) / 2
  // LzLx =   (  Lz   )(L+ + L-) / 2
  // LzLy = -i(  Lz   )(L+ - L-) / 2
  // LzLz =   (  Lz   )(  Lz   )
  vector<Matrix> ll(NTimes(), Matrix(3,3));
  const complex<double> I(0.0,1.0);
  if(Lmodes.size()==0) {
    Lmodes.push_back(lm[0][0]);
    for(unsigned int i_m=0; i_m<NModes(); ++i_m) {
      if(std::find(Lmodes.begin(), Lmodes.end(), lm[i_m][0]) == Lmodes.end() ) {
        Lmodes.push_back(lm[i_m][0]);
      }
    }
  }
  for(unsigned int iL=0; iL<Lmodes.size(); ++iL) {
    const int L = Lmodes[iL];
    for(int M=-L; M<=L; ++M) {
      const int iMm2 = (M-2>=-L ? FindModeIndex(L,M-2) : 0);
      const int iMm1 = (M-1>=-L ? FindModeIndex(L,M-1) : 0);
      const int iM   = FindModeIndex(L,M);
      const int iMp1 = (M+1<=L ? FindModeIndex(L,M+1) : 0);
      const int iMp2 = (M+2<=L ? FindModeIndex(L,M+2) : 0);
      for(unsigned int iTime=0; iTime<NTimes(); ++iTime) {
        const complex<double> LpLp = (M+2<=L  ? conj(data[iMp2][iTime]) * LadderOperatorFactor(L, M+1)    * LadderOperatorFactor(L, M) * data[iM][iTime] : 0.0);
        const complex<double> LpLm = (M-1>=-L ? conj(data[iM][iTime])   * LadderOperatorFactor(L, M-1)    * LadderOperatorFactor(L, -M) * data[iM][iTime] : 0.0);
        const complex<double> LmLp = (M+1<=L  ? conj(data[iM][iTime])   * LadderOperatorFactor(L, -(M+1)) * LadderOperatorFactor(L, M) * data[iM][iTime] : 0.0);
        const complex<double> LmLm = (M-2>=-L ? conj(data[iMm2][iTime]) * LadderOperatorFactor(L, -(M-1)) * LadderOperatorFactor(L, -M) * data[iM][iTime] : 0.0);
        const complex<double> LpLz = (M+1<=L  ? conj(data[iMp1][iTime]) * LadderOperatorFactor(L, M) * double(M)      * data[iM][iTime] : 0.0);
        const complex<double> LzLp = (M+1<=L  ? conj(data[iMp1][iTime]) * double(M+1) * LadderOperatorFactor(L, M)  * data[iM][iTime] : 0.0);
        const complex<double> LmLz = (M-1>=-L ? conj(data[iMm1][iTime]) * LadderOperatorFactor(L, -M) * double(M)     * data[iM][iTime] : 0.0);
        const complex<double> LzLm = (M-1>=-L ? conj(data[iMm1][iTime]) * double(M-1) * LadderOperatorFactor(L, -M) * data[iM][iTime] : 0.0);
        const complex<double> LzLz = conj(data[iM][iTime]) * double(M) * double(M) * data[iM][iTime];
        //
        const complex<double> LxLx = 0.25 * (LpLp + LmLm + LmLp + LpLm);
        const complex<double> LxLy = -0.25 * I * (LpLp - LmLm + LmLp - LpLm);
        const complex<double> LxLz = 0.5 * (LpLz + LmLz);
        const complex<double> LyLx = -0.25 * I * (LpLp - LmLp + LpLm - LmLm);
        const complex<double> LyLy = -0.25 * (LpLp - LmLp - LpLm + LmLm);
        const complex<double> LyLz = -0.5 * I * (LpLz - LmLz);
        const complex<double> LzLx = 0.5 * (LzLp + LzLm);
        const complex<double> LzLy = -0.5 * I * (LzLp - LzLm);
        //const complex<double> LzLz = (LzLz);
        ll[iTime](0,0) += real( LxLx );
        ll[iTime](0,1) += real( LxLy + LyLx )/2.0;
        ll[iTime](0,2) += real( LxLz + LzLx )/2.0;
        ll[iTime](1,0) += real( LyLx + LxLy )/2.0;
        ll[iTime](1,1) += real( LyLy );
        ll[iTime](1,2) += real( LyLz + LzLy )/2.0;
        ll[iTime](2,0) += real( LzLx + LxLz )/2.0;
        ll[iTime](2,1) += real( LzLy + LyLz )/2.0;
        ll[iTime](2,2) += real( LzLz );
      }
    }
  }
  return ll;
}


/// Calculate the principal axis of the LL matrix, as prescribed by O'Shaughnessy et al.
std::vector<std::vector<double> > GWFrames::Waveform::OShaughnessyEtAlVector(const std::vector<int>& Lmodes) const {
  ///
  /// \param Lmodes L modes to evaluate
  ///
  /// If Lmodes is empty (default), all L modes are used.  Setting
  /// Lmodes to [2] or [2,3,4], for example, restricts the range of
  /// the sum.
  ///

  // Calculate the LL matrix at each instant
  vector<Matrix> ll = LLMatrix(Lmodes);

  // Calculate the dominant principal axis of LL at each instant
  vector<vector<double> > dpa(NTimes(), vector<double>(3));
  for(unsigned int i=0; i<ll.size(); ++i) {
    dpa[i] = GWFrames::DominantPrincipalAxis(ll[i]);
  }

  // Now, go through and make the vectors reasonably continuous.
  for(unsigned int i=1; i<dpa.size(); ++i) {
    const double x=dpa[i][0];
    const double y=dpa[i][1];
    const double z=dpa[i][2];
    const double dx=x-dpa[i-1][0];
    const double dy=y-dpa[i-1][1];
    const double dz=z-dpa[i-1][2];
    const double Norm = std::sqrt(x*x+y*y+z*z);
    const double dNorm = std::sqrt(dx*dx+dy*dy+dz*dz);
    if(dNorm>Norm) {
      dpa[i][0] = -dpa[i][0];
      dpa[i][1] = -dpa[i][1];
      dpa[i][2] = -dpa[i][2];
    }
  }

  return dpa;
}

/// Calculate the angular velocity of the Waveform.
vector<vector<double> > GWFrames::Waveform::AngularVelocityVector(const vector<int>& Lmodes) const {
  ///
  /// \param Lmodes L modes to evaluate
  ///
  /// This returns the angular velocity of the Waveform, as defined in
  /// Sec. II of <a href="http://arxiv.org/abs/1302.2919">"Angular
  /// velocity of gravitational radiation and the corotating
  /// frame"</a>.  Note that the returned vector is relative to the
  /// inertial frame.
  ///
  /// If Lmodes is empty (default), all L modes are used.  Setting
  /// Lmodes to [2] or [2,3,4], for example, restricts the range of
  /// the sum.
  ///

  // Calculate the L vector and LL matrix at each instant
  vector<vector<double> > l = LdtVector(Lmodes);
  vector<Matrix> ll = LLMatrix(Lmodes);

  // If the frame is nontrivial, include its contribution
  const bool TimeDependentFrame = (frame.size()>1);
  const vector<Quaternion> Rdot = (TimeDependentFrame
                                   ? Quaternions::QuaternionDerivative(frame, t)
                                   : vector<Quaternion>(0));
  const bool ConstantNontrivialFrame = (frame.size()==1);
  const Quaternion R0 = (ConstantNontrivialFrame
                         ? frame[0]
                         : Quaternion(1,0,0,0));

  // Construct some objects for storage
  vector<vector<double> > omega(NTimes(), vector<double>(3));
  int s;
  gsl_vector* x = gsl_vector_alloc(3);
  gsl_permutation* p = gsl_permutation_alloc(3);

  // Loop through time steps
  for(unsigned int iTime=0; iTime<omega.size(); ++iTime) {
    // Solve   -omega * LL = L   at each time step, using LU decomposition
    gsl_vector_view b = gsl_vector_view_array(&l[iTime][0], 3);
    gsl_linalg_LU_decomp(ll[iTime].gslobj(), p, &s);
    gsl_linalg_LU_solve(ll[iTime].gslobj(), p, &b.vector, x);

    // Save data for this time step
    omega[iTime][0] = -gsl_vector_get(x, 0);
    omega[iTime][1] = -gsl_vector_get(x, 1);
    omega[iTime][2] = -gsl_vector_get(x, 2);

    if(TimeDependentFrame) { // Include frame-rotation effects
      const Quaternion& R = frame[iTime];
      omega[iTime] = (R*Quaternion(omega[iTime])*R.conjugate() + 2*Rdot[iTime]*R.conjugate()).vec();
    } else if(ConstantNontrivialFrame) { // Just rotate the result
      omega[iTime] = (R0*Quaternion(omega[iTime])*R0.conjugate()).vec();
    }
  }

  // Free the memory
  gsl_permutation_free(p);
  gsl_vector_free(x);

  return omega;
}


/// Frame in which the rotation is minimal.
std::vector<Quaternions::Quaternion> GWFrames::Waveform::CorotatingFrame(const std::vector<int>& Lmodes) const {
  ///
  /// \param Lmodes L modes to evaluate
  ///
  /// This function combines the steps required to obtain the
  /// corotating frame.
  ///
  /// If Lmodes is empty (default), all L modes are used.  Setting
  /// Lmodes to [2] or [2,3,4], for example, restricts the range of
  /// the sum.
  ///

  return Quaternions::FrameFromAngularVelocity(QuaternionArray(AngularVelocityVector(Lmodes)), T());
}

/// Deduce PN-equivalent orbital angular velocity from Waveform
vector<vector<double> > GWFrames::Waveform::PNEquivalentOrbitalAV(const vector<int>& Lmodes) const {
  ///
  /// \param Lmodes L modes to evaluate
  ///
  /// This function simply takes the projection of the field's
  /// angular-velocity vector \f$\vec{\omega}\f$ along the dominant
  /// eigenvector \f$\hat{V}_f\f$ of \f$<LL>\f$.  This should be
  /// equivalent to the orbital angular velocity of the PN system.
  /// Note that the returned vector is relative to the inertial frame.
  ///
  /// If Lmodes is empty (default), all L modes are used.  Setting
  /// Lmodes to [2] or [2,3,4], for example, restricts the range of
  /// the sum.
  ///

  // If the frame is nontrivial, include its contribution
  const bool TimeDependentFrame = (frame.size()>1);
  const bool ConstantNontrivialFrame = (frame.size()==1);
  const Quaternion R0 = (ConstantNontrivialFrame
                         ? frame[0]
                         : Quaternion(1,0,0,0));

  // Get the vectors we need
  const vector<vector<double> > omega = AngularVelocityVector(Lmodes);
  vector<vector<double> > V_f = OShaughnessyEtAlVector(Lmodes);

  // Do the calculation
  for(unsigned int i=0; i<V_f.size(); ++i) {
    // Rotate V_f[i] into the inertial frame, if necessary
    if(TimeDependentFrame) {
      const Quaternion& R = frame[i];
      V_f[i] = (R*Quaternion(V_f[i])*R.conjugate()).vec();
    } else if(ConstantNontrivialFrame) {
      V_f[i] = (R0*Quaternion(V_f[i])*R0.conjugate()).vec();
    }
    // Calculate the dot product between V_f and omega
    double dotproduct = 0.0;
    for(unsigned int j=0; j<3; ++j) {
      dotproduct += V_f[i][j]*omega[i][j];
    }
    // Apply the rescaling
    for(unsigned int j=0; j<3; ++j) {
      V_f[i][j] *= dotproduct;
    }
  }
  return V_f; // [which is now V_f * (V_f \cdot omega)]
}

/// Deduce PN-equivalent precessional angular velocity from Waveform
vector<vector<double> > GWFrames::Waveform::PNEquivalentPrecessionalAV(const vector<int>& Lmodes) const {
  ///
  /// \param Lmodes L modes to evaluate
  ///
  /// This function subtracts the PN-equivalent *orbital* angular
  /// velocity (given by PNEquivalentOrbitalAV) from the field's
  /// angular velocity.  This should be equivalent to the precessional
  /// angular velocity of the PN system.  Note that the returned
  /// vector is relative to the inertial frame.
  ///
  /// This may be essentially numerical noise if there is no
  /// precession, or if precession has oscillated to zero.
  ///
  /// If Lmodes is empty (default), all L modes are used.  Setting
  /// Lmodes to [2] or [2,3,4], for example, restricts the range of
  /// the sum.
  ///
  /// \sa PNEquivalentOrbitalAV
  ///

  vector<vector<double> > omega = AngularVelocityVector(Lmodes);
  const vector<vector<double> > Omega_orb = PNEquivalentOrbitalAV(Lmodes);
  for(unsigned int i=0; i<Omega_orb.size(); ++i) {
    for(unsigned int j=0; j<3; ++j) {
      omega[i][j] -= Omega_orb[i][j];
    }
  }
  return omega;
}


/// Transform Waveform to O'Shaughnessy et al. frame.
GWFrames::Waveform& GWFrames::Waveform::TransformToOShaughnessyEtAlFrame(const std::vector<int>& Lmodes) {
  ///
  /// \param Lmodes L modes to evaluate
  ///
  /// This function combines the steps required to obtain the Waveform
  /// in the O'Shaughnessy et al. frame.
  ///
  /// If Lmodes is empty (default), all L modes are used.  Setting
  /// Lmodes to [2] or [2,3,4], for example, restricts the range of
  /// the sum.
  ///
  history << "this->TransformToOShaughnessyEtAlFrame(" << StringForm(Lmodes) << ")\n#";
  this->frameType = GWFrames::Coprecessing;
  return this->RotateDecompositionBasis(FrameFromZ(normalized(QuaternionArray(this->OShaughnessyEtAlVector(Lmodes))), T()));
}

/// Transform Waveform to frame aligned with angular-velocity vector.
GWFrames::Waveform& GWFrames::Waveform::TransformToAngularVelocityFrame(const std::vector<int>& Lmodes) {
  ///
  /// \param Lmodes L modes to evaluate
  ///
  /// This function combines the steps required to obtain the Waveform
  /// in the frame aligned with the angular-velocity vector.  Note
  /// that this frame is not the corotating frame; this frame has its
  /// z axis aligned with the angular-velocity vector.
  ///
  /// If Lmodes is empty (default), all L modes are used.  Setting
  /// Lmodes to [2] or [2,3,4], for example, restricts the range of
  /// the sum.
  ///
  history << "this->TransformToAngularVelocityFrame(" << StringForm(Lmodes) << ")\n#";
  vector<Quaternion> R_AV = normalized(QuaternionArray(this->AngularVelocityVector(Lmodes)));
  this->frameType = GWFrames::Coprecessing;
  return this->RotateDecompositionBasis(FrameFromZ(R_AV, T()));
}

/// Transform Waveform to corotating frame.
GWFrames::Waveform& GWFrames::Waveform::TransformToCorotatingFrame(const std::vector<int>& Lmodes) {
  ///
  /// \param Lmodes L modes to evaluate
  ///
  /// This function combines the steps required to obtain the Waveform
  /// in the corotating frame.  Note that this leaves an integration
  /// constant unset.  To set it, the modes should be rotated so that
  /// they are aligned with the frame using `AlignModesToFrame`.
  ///
  /// If Lmodes is empty (default), all L modes are used.  Setting
  /// Lmodes to [2] or [2,3,4], for example, restricts the range of
  /// the sum.
  ///

  if(frameType != GWFrames::Inertial) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nWarning: Asking to transform a Waveform in the " << GWFrames::WaveformFrameNames[frameType] << " frame into the corotating frame."
              << "\n         You have to think very carefully about whether or not this is what you really want.\n"
              << "\n         This should probably only be applied to Waveforms in the " << GWFrames::WaveformFrameNames[GWFrames::Inertial] << " frame.\n"
              << std::endl;
  }

  vector<Quaternion> R_corot = this->CorotatingFrame(Lmodes);
  this->frameType = GWFrames::Corotating;
  history << "this->TransformToCorotatingFrame(" << StringForm(Lmodes) << ")\n#";
  return this->RotateDecompositionBasis(R_corot);
}

/// Transform Waveform to an inertial frame.
GWFrames::Waveform& GWFrames::Waveform::TransformToInertialFrame() {
  ///
  /// This function uses the stored frame information to transform
  /// from whatever rotating frame the waveform is currently in, to a
  /// stationary, inertial frame.  This is the usual frame of scri^+,
  /// and is the frame in which GW observations should be made.
  ///

  if(frameType == GWFrames::Inertial) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nWarning: Waveform is already in the " << GWFrames::WaveformFrameNames[GWFrames::Inertial] << " frame;"
              << " this function will have no effect."
              << "\n         If it's not really in that frame, tell the Waveform first.\n"
              << std::endl;
    return *this;
  }

  if(frameType == GWFrames::UnknownFrameType) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nWarning: Waveform frame type is " << GWFrames::WaveformFrameNames[GWFrames::UnknownFrameType] << "."
              << "\n         I assume you know what you're doing...\n"
              << std::endl;
    return *this;
  }

  if(frame.size() == 0) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nWarning: frame.size()=" << frame.size() << ".  I will assume that this *is* the"
              << "           inertial frame, so this function will have no effect.\n"
              <<"            And I will continue to assume you know what you're doing...\n"
              << std::endl;
    return *this;
  }

  if(frame.size() != NTimes()) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nError: (frame.size()=" << frame.size() << ") != (NTimes()=" << NTimes() << ")."
              << "         I don't know what to do with this data, or what the inertial frame is."
              << std::endl;
    throw(GWFrames_VectorSizeMismatch);
  }

  history << "this->TransformToInertialFrame();\n#";
  this->frameType = GWFrames::Inertial; // Must come first
  this->RotateDecompositionBasis(Quaternions::conjugate(frame));
  this->SetFrame(vector<Quaternion>(0));
  return *this;
}

/// Transform Waveform uncertainties to corotating frame.
GWFrames::Waveform& GWFrames::Waveform::TransformUncertaintiesToCorotatingFrame(const std::vector<Quaternions::Quaternion>& R_frame) {
  ///
  /// \param R_frame Vector of rotors giving corotating frame of the data.
  ///

  if(frameType != GWFrames::Inertial) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nWarning: Asking to transform a Waveform in the " << GWFrames::WaveformFrameNames[frameType] << " frame into the corotating frame."
              << "\n         You have to think very carefully about whether or not this is what you really want.\n"
              << "\n         This should probably only be applied to Waveforms in the " << GWFrames::WaveformFrameNames[GWFrames::Inertial] << " frame.\n"
              << std::endl;
  }
  if(R_frame.size() != NTimes()) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nError: Asking to transform uncertainties of size " << NTimes() << " by rotors of size " << R_frame.size() << "\n"
              << std::endl;
    throw(GWFrames_VectorSizeMismatch);
  }

  this->frameType = GWFrames::Corotating;
  history << "this->TransformUncertaintiesToCorotatingFrame(R_frame); // R_frame=[" << setprecision(16) << R_frame[0] << ", " << R_frame[1] << ", ...]\n#";
  return this->RotateDecompositionBasisOfUncertainties(R_frame);
}

/// Transform Waveform to an inertial frame.
GWFrames::Waveform& GWFrames::Waveform::TransformUncertaintiesToInertialFrame() {
  ///
  /// This function uses the stored frame information to transform
  /// from whatever rotating frame the waveform is currently in, to a
  /// stationary, inertial frame.  This is the usual frame of scri^+,
  /// and is the frame in which GW observations should be made.
  ///
  history << "this->TransformUncertaintiesToInertialFrame();\n#";
  this->frameType = GWFrames::Inertial; // Must come first
  this->RotateDecompositionBasisOfUncertainties(Quaternions::conjugate(frame));
  this->SetFrame(vector<Quaternion>(0));
  return *this;
}

/// Interpolate the Waveform to a new set of time instants.
GWFrames::Waveform& GWFrames::Waveform::InterpolateInPlace(const std::vector<double>& NewTime) {
  if(NewTime.size()==0) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Asking for empty Waveform." << std::endl;
    throw(GWFrames_EmptyIntersection);
  }
  if(NewTime[0]<t[0]) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Asking for extrapolation; we only do interpolation.\n"
              << "NewTime[0]=" << NewTime[0] << "\tt[0]=" << t[0] << std::endl;
    throw(GWFrames_EmptyIntersection);
  }
  if(NewTime.back()>t.back()) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Asking for extrapolation; we only do interpolation.\n"
              << "NewTime.back()=" << NewTime.back() << "\tt.back()=" << t.back() << std::endl;
    throw(GWFrames_EmptyIntersection);
  }

  history << HistoryStr()
            << "*this = this->Interpolate(NewTime);" << std::endl;
  const vector<double> OldTime(t);
  if(frame.size()==1) { // Assume we have just a constant non-trivial frame
    frame = frame;
  } else if(frame.size()>1) { // Assume we have frame data for each time step
    frame = Squad(frame, OldTime, NewTime);
  }
  MatrixC NewData;
  NewData.resize(NModes(), NewTime.size());
  // Initialize the GSL interpolators for the data
  gsl_interp_accel* accRe = gsl_interp_accel_alloc();
  gsl_interp_accel* accIm = gsl_interp_accel_alloc();
  gsl_spline* splineRe = gsl_spline_alloc(gsl_interp_cspline, OldTime.size());
  gsl_spline* splineIm = gsl_spline_alloc(gsl_interp_cspline, OldTime.size());
  // Now loop over each mode filling in the waveform data
  for(unsigned int i_m=0; i_m<NModes(); ++i_m) {
    // Extract the real and imaginary parts of the data separately for GSL
    const vector<double> re(Re(i_m));
    const vector<double> im(Im(i_m));
    // Initialize the interpolators for this data set
    gsl_spline_init(splineRe, &(OldTime)[0], &re[0], OldTime.size());
    gsl_spline_init(splineIm, &(OldTime)[0], &im[0], OldTime.size());
    // Assign the interpolated data
    for(unsigned int i_t=0; i_t<NewTime.size(); ++i_t) {
      NewData[i_m][i_t] = complex<double>( gsl_spline_eval(splineRe, NewTime[i_t], accRe), gsl_spline_eval(splineIm, NewTime[i_t], accIm) );
    }
  }
  data.swap(NewData);
  t = NewTime;
  // Free the interpolators
  gsl_interp_accel_free(accRe);
  gsl_interp_accel_free(accIm);
  gsl_spline_free(splineRe);
  gsl_spline_free(splineIm);

  return *this;
}

/// Interpolate the Waveform to a new set of time instants.
GWFrames::Waveform GWFrames::Waveform::Interpolate(const std::vector<double>& NewTime) const {
  if(NewTime.size()==0) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Asking for empty Waveform." << std::endl;
    throw(GWFrames_EmptyIntersection);
  }
  if(NewTime[0]<t[0]) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Asking for extrapolation; we only do interpolation.\n"
              << "NewTime[0]=" << NewTime[0] << "\tt[0]=" << t[0] << std::endl;
    throw(GWFrames_EmptyIntersection);
  }
  if(NewTime.back()>t.back()) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Asking for extrapolation; we only do interpolation.\n"
              << "NewTime.back()=" << NewTime.back() << "\tt.back()=" << t.back() << std::endl;
    throw(GWFrames_EmptyIntersection);
  }

  Waveform C;
  C.history << HistoryStr()
            << "*this = this->Interpolate(NewTime);" << std::endl;
  C.t = NewTime;
  if(frame.size()==1) { // Assume we have just a constant non-trivial frame
    C.frame = frame;
  } else if(frame.size()>1) { // Assume we have frame data for each time step
    C.frame = Squad(frame, t, NewTime);
  }
  C.frameType = frameType;
  C.dataType = dataType;
  C.rIsScaledOut = rIsScaledOut;
  C.mIsScaledOut = mIsScaledOut;
  C.lm = lm;
  C.data.resize(NModes(), NewTime.size());
  // Initialize the GSL interpolators for the data
  gsl_interp_accel* accRe = gsl_interp_accel_alloc();
  gsl_interp_accel* accIm = gsl_interp_accel_alloc();
  gsl_spline* splineRe = gsl_spline_alloc(gsl_interp_cspline, NTimes());
  gsl_spline* splineIm = gsl_spline_alloc(gsl_interp_cspline, NTimes());
  // Now loop over each mode filling in the waveform data
  for(unsigned int i_m=0; i_m<C.NModes(); ++i_m) {
    // Extract the real and imaginary parts of the data separately for GSL
    const vector<double> re = Re(i_m);
    const vector<double> im = Im(i_m);
    // Initialize the interpolators for this data set
    gsl_spline_init(splineRe, &(t)[0], &re[0], NTimes());
    gsl_spline_init(splineIm, &(t)[0], &im[0], NTimes());
    // Assign the interpolated data
    for(unsigned int i_t=0; i_t<C.t.size(); ++i_t) {
      C.data[i_m][i_t] = complex<double>( gsl_spline_eval(splineRe, C.t[i_t], accRe), gsl_spline_eval(splineIm, C.t[i_t], accIm) );
    }
  }
  // Free the interpolators
  gsl_interp_accel_free(accRe);
  gsl_interp_accel_free(accIm);
  gsl_spline_free(splineRe);
  gsl_spline_free(splineIm);

  return C;
}

/// Extract a segment of a Waveform.
GWFrames::Waveform GWFrames::Waveform::Segment(const unsigned int i1, const unsigned int i2) const {
  ///
  /// \param i1 Index of initial time
  /// \param i2 Index just beyond final time
  ///
  if(i1>i2) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Requesting impossible segment"
              << "\ni1=" << i1 << "  >  i2=" << i2 << std::endl;
    throw(GWFrames_EmptyIntersection);
  }
  if(i2>NTimes()) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Requesting impossible segment"
              << "\ni2=" << i2 << "  >  NTimes()=" << NTimes() << std::endl;
    throw(GWFrames_IndexOutOfBounds);
  }
  GWFrames::Waveform C;
  C.history << HistoryStr()
            << "this->Segment(" << i1 << ", " << i2 << ");" << std::endl;
  C.t = vector<double>(&t[i1], &t[i2]);
  if(frame.size()==1) {
    C.frame = frame;
  } else if(frame.size()>1) {
    C.frame = vector<Quaternions::Quaternion>(&frame[i1], &frame[i2]);
  }
  C.frameType = frameType;
  C.dataType = dataType;
  C.rIsScaledOut = rIsScaledOut;
  C.mIsScaledOut = mIsScaledOut;
  C.lm = lm;
  C.data.resize(NModes(), i2-i1);
  for(unsigned int i_m=0; i_m<NModes(); ++i_m) {
    for(unsigned int i_t=i1; i_t<i2; ++i_t) {
      C.data[i_m][i_t-i1] = data[i_m][i_t];
    }
  }
  return C;
}

/// Find the time offset aligning this waveform to the other at the fiducial time.
void GWFrames::Waveform::GetAlignmentOfTime(const Waveform& A, const double t_fid, double& deltat) const {
  ///
  /// \param A Fixed Waveform in inertial frame to which this Waveform is aligned
  /// \param t_fid
  /// \param deltat The value to be returned
  ///
  /// This function simply finds the appropriate time offset, rather
  /// than applying it.  This is called by `AlignTime` and probably
  /// does not need to be called directly; see that function's
  /// documentation for more details.
  ///
  /// \sa AlignTime
  ///

  const Waveform& B = *this;

  // Interpolate the angular velocity magnitude of A to t_fid
  double MagAV_fid=0.0;
  {
    const vector<vector<double> > AV_A = A.AngularVelocityVector();
    std::vector<double> MagAV_A(A.NTimes());
    for(unsigned int i=0; i<A.NTimes(); ++i) {
      MagAV_A[i] = GWFrames::abs(AV_A[i]);
    }
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, A.NTimes());
    gsl_spline_init(spline, &A.t[0], &MagAV_A[0], A.NTimes());
    MagAV_fid = gsl_spline_eval(spline, t_fid, acc);
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }

  // Find that value in this waveform
  {
    const vector<vector<double> > AV_B = B.AngularVelocityVector();
    std::vector<double> MagAV_B(B.NTimes());
    for(unsigned int i=0; i<B.NTimes(); ++i) {
      MagAV_B[i] = GWFrames::abs(AV_B[i]);
    }
    unsigned int i=0;
    while(MagAV_B[i]<MagAV_fid && i<MagAV_B.size()-1) { ++i; }
    if(i==0 || i==MagAV_B.size()) {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Cannot find the appropriate angular velocity in this waveform."
                << "\nt_fid=" << t_fid << "\tAV_fid=" << MagAV_fid << "\tAV_B[0]=" << MagAV_B[0] << "\tAV_B[-1]=" << MagAV_B.back() << std::endl;
      throw(GWFrames_IndexOutOfBounds);
    }
    const unsigned int NPoints = 5;
    const unsigned int i1 = (i<NPoints ? 0 : i-NPoints);
    const unsigned int i2 = (i>MagAV_B.size()-1-NPoints ? MagAV_B.size()-1 : i+NPoints);
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, i2-i1+1);
    gsl_spline_init(spline, &MagAV_B[i1], &t[i1], i2-i1+1);
    deltat = t_fid - gsl_spline_eval(spline, MagAV_fid, acc);
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }
}

/// Change this Waveform by aligning to the other at the given time
GWFrames::Waveform& GWFrames::Waveform::AlignTime(const GWFrames::Waveform& A, const double t_fid) {
  ///
  /// \param A Fixed Waveform in inertial frame to which this Waveform is aligned
  /// \param t_fid
  ///
  /// Note that this function operates in place; the Waveform to which
  /// it is applied will change.
  ///
  /// As noted above, it is implicitly assumed that both Waveforms are
  /// in an inertial frame, so that the magnitude of the angular
  /// velocity may be properly measured.  This could be adjusted to
  /// account for the angular velocity of the frame, but hasn't been
  /// yet.
  ///
  /// To improve accuracy, the angular velocity of A is interpolated
  /// to t_fid.  The time of B is then interpolated to the
  /// interpolated angular velocity.  This assumes that B's angular
  /// velocity is strictly monotonic for roughly 5 data points to
  /// either side.
  ///

  // Call GetAlignment to find the deltat
  double deltat;
  GetAlignmentOfTime(A, t_fid, deltat);

  // Offset the time axis by that amount
  t = t + deltat;

  // Record what happened
  history << "this->AlignTime(A, " << std::setprecision(16) << t_fid << ");  # deltat=" << deltat << std::endl;

  return *this;
}

/// Find the appropriate rotation to fix the orientation of the corotating frame.
void GWFrames::Waveform::GetAlignmentOfDecompositionFrameToModes(const double t_fid, Quaternion& R_eps,
                                                                 const std::vector<int>& Lmodes) const {
  ///
  /// \param t_fid Fiducial time at which the alignment should happen
  /// \param R_eps Returned rotor
  /// \param Lmodes Lmodes to use in computing \f$<LL>\f$
  ///
  /// This function simply finds the rotation necessary to align the
  /// corotating frame to the waveform at the fiducial time, rather
  /// than applying it.  This is called by
  /// `AlignDecompositionFrameToModes` and probably does not need to
  /// be called directly; see that function's documentation for more
  /// details.
  ///
  /// \sa AlignDecompositionFrameToModes
  ///

  // We seek that R_c such that R_corot(t_fid)*R_c rotates the z axis
  // onto V_f.  V_f measured in this frame is given by
  //     V_f = R_V_f * Z * R_V_f.conjugate(),
  // (note Z rather than z) where R_V_f is found below.  But
  //     Z = R_corot * z * R_corot.conjugate(),
  // so in the (x,y,z) frame,
  //     V_f = R_V_f * R_corot * z * R_corot.conjugate() * R_V_f.conjugate().
  // Now, this is the standard composition for rotating physical
  // vectors.  However, rotation of the basis behaves oppositely, so
  // we want R_V_f as our constant rotation, applied as a rotation of
  // the decomposition basis.  We also want to rotate so that the
  // phase of the (2,2) mode is zero at t_fid.  This can be achieved
  // with an initial rotation.

  if(frameType!=GWFrames::Coprecessing && frameType!=GWFrames::Coorbital && frameType!=GWFrames::Corotating) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ":"
              << "\nError: GetAlignmentOfDecompositionFrameToModes only takes Waveforms in the "
              << GWFrames::WaveformFrameNames[GWFrames::Coprecessing] << ", "
              << GWFrames::WaveformFrameNames[GWFrames::Coorbital] << ", or "
              << GWFrames::WaveformFrameNames[GWFrames::Corotating] << " frames."
              << "\n       This Waveform is in the " << GWFrames::WaveformFrameNames[frameType] << " frame." << std::endl;
    throw(GWFrames_WrongFrameType);
  }

  if(t_fid<t[0] || t_fid>t.back()) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ":"
              << "\nError: The requested alignment time t_fid=" << t_fid << " is outside the range of times in this waveform ("
              << t[0] << ", " << t.back() << ")." << std::endl;
    throw(GWFrames_EmptyIntersection);
  }

  // Get direction of angular-velocity vector near t_fid
  Quaternion omegaHat;
  {
    int i_t_fid = 0;
    while(t[i_t_fid]<t_fid && i_t_fid<int(t.size())) { ++i_t_fid; }
    unsigned int i1 = (i_t_fid-10<0 ? 0 : i_t_fid-10);
    unsigned int i2 = (i_t_fid+11>int(t.size()) ? t.size() : i_t_fid+11);
    vector<double> tRegion(&t[i1], &t[i2]);
    const Waveform Region = (this->Interpolate(tRegion)).TransformToInertialFrame();
    omegaHat = Quaternion(Region.AngularVelocityVector()[i_t_fid-i1]).normalized();
    // omegaHat contains the components of that vector relative to the
    // inertial frame.  To get its components in this Waveform's
    // (possibly rotating) frame, we need to rotate it by the inverse
    // of this Waveform's `frame` data:
    if(Region.Frame().size()>1) {
      const Quaternion& R = Region.Frame(i_t_fid-i1);
      omegaHat = R.inverse() * omegaHat * R;
    } else if(Region.Frame().size()==1) {
      const Quaternion& R = Region.Frame(0);
      omegaHat = R.inverse() * omegaHat * R;
    }
  }

  // Interpolate the Waveform to t_fid
  Waveform Instant = this->Interpolate(vector<double>(1,t_fid));

  // V_f is the dominant eigenvector of <LL>, suggested by O'Shaughnessy et al.
  const Quaternion V_f = Quaternions::Quaternion(Instant.OShaughnessyEtAlVector(Lmodes)[0]);
  const Quaternion V_f_aligned = (omegaHat.dot(V_f.normalized()) < 0 ? -V_f : V_f);

  // R_V_f is the rotor taking the Z axis onto V_f
  const Quaternion R_V_f = Quaternions::sqrtOfRotor(-V_f_aligned*Quaternions::Quaternion(0,0,0,1));

  // Now rotate Instant so that its z axis is aligned with V_f
  Instant.RotateDecompositionBasis(R_V_f);

  // Get the phase of the (2,2) mode after rotation
  const unsigned int i_22 = Instant.FindModeIndex(2,2);
  const double phase_22 = std::atan2(Instant.Im(i_22,0),Instant.Re(i_22,0));

  // R_eps is the rotation we will be applying on the right-hand side
  R_eps = R_V_f * Quaternions::exp(Quaternions::Quaternion(0,0,0,(-phase_22/4)));
}

/// Fix the orientation of the corotating frame.
GWFrames::Waveform& GWFrames::Waveform::AlignDecompositionFrameToModes(const double t_fid, const std::vector<int>& Lmodes) {
  ///
  /// \param t_fid Fiducial time at which the alignment should happen
  /// \param Lmodes Lmodes to use in computing \f$<LL>\f$
  ///
  /// The corotating frame is only defined up to some constant rotor
  /// R_eps; if R_corot is corotating, then so is R_corot*R_eps.  This
  /// function uses that freedom to ensure that the frame is aligned
  /// with the Waveform modes at the fiducial time.  In particular, it
  /// ensures that the Z axis of the frame in which the decomposition
  /// is done is along the dominant eigenvector of \f$<LL>\f$
  /// (suggested by O'Shaughnessy et al.), and the phase of the (2,2)
  /// mode is zero.
  ///
  /// If Lmodes is empty (default), all L modes are used.  Setting
  /// Lmodes to [2] or [2,3,4], for example, restricts the range of
  /// the sum.
  ///

  // Find the appropriate rotation
  Quaternion R_eps;
  GetAlignmentOfDecompositionFrameToModes(t_fid, R_eps, Lmodes);

  // Record what happened
  history << "this->AlignDecompositionFrameToModes(" << std::setprecision(16) << t_fid << ", " << Lmodes << ");  # R_eps=" << R_eps << std::endl;

  // Now, apply the rotation
  this->RotateDecompositionBasis(R_eps);

  return *this;
}

/// Find the appropriate rotation to fix the orientation of the corotating frame over a range of time.
void GWFrames::Waveform::GetAlignmentOfDecompositionFrameToModes(const double t1, const double t2, const Quaternions::Quaternion& nHat_t1,
                                                                 Quaternions::Quaternion& R_eps, const std::vector<int>& Lmodes) const {
  ///
  /// \param t1 Beginning of time range over which the alignment should happen
  /// \param t2 End of time range over which the alignment should happen
  /// \param nHat_t1 The approximate direction of nHat at t1
  /// \param R_eps Returned rotor
  /// \param Lmodes Lmodes to use in computing \f$<LL>\f$
  ///
  /// This function simply finds the rotation necessary to align the
  /// corotating frame to the waveform at the fiducial time, rather
  /// than applying it.  This is called by
  /// `AlignDecompositionFrameToModes` and probably does not need to
  /// be called directly; see that function's documentation for more
  /// details.
  ///
  /// \sa AlignDecompositionFrameToModes
  ///

  Waveform W(*this);
  W.DropTimesOutside(t1, t2);
  std::vector<Quaternion> R_zeta(W.NTimes());
  R_zeta[0] = Quaternions::One;
  // First, do the naive alignment using a single time, to get us started
  this->GetAlignmentOfDecompositionFrameToModes(W.T(0), R_eps, Lmodes);
  W.RotateDecompositionBasis(R_eps);
  // Now, the decomposition frame will be *nearly* aligned for all times.  But we want to make it better.
  for(int i=1; i<W.NTimes(); ++i) {
    // We will use the previous time step's result as a starting point
    const Quaternion ellHat_fid = R_zeta[i-1]*Quaternions::zHat*R_zeta[i-1].inverse();
    const double t_fid = W.T(i);
    vector<double> vt_fid(1);
    vt_fid[0] = t_fid;
    Waveform W_fid = W.Interpolate(vt_fid);
    const Quaternion V_f = Quaternion(W_fid.OShaughnessyEtAlVector(Lmodes)[0]).normalized();
    const Quaternion V_f_aligned = ( ellHat_fid.dot(V_f) < 0 ? -V_f : V_f);
    // R_V_f is the rotor taking zHat onto V_f_aligned via the previous step's ellHat
    const Quaternion R_V_f = Quaternions::sqrtOfRotor(-V_f_aligned*ellHat_fid)*R_zeta[i-1];
    // Now rotate W_fid so that its Z axis is aligned with V_f
    vector<unsigned int> Two(1);
    Two[0] = 2;
    W_fid.KeepOnlyEllModes(Two);
    W_fid.RotateDecompositionBasis(R_V_f);
    // Get the phase of the (2,2) mode after rotation
    const unsigned int i_22 = W_fid.FindModeIndex(2,2);
    const double phase_22 = std::atan2(W_fid.Im(i_22,0),W_fid.Re(i_22,0));
    // R_zeta[i] is the rotation we would like to apply on the right-hand side for this time step
    R_zeta[i] = (R_V_f * Quaternions::exp(Quaternion(0,0,0,(-phase_22/2)/2))).normalized();
  }
  R_eps = R_eps * Quaternions::ApproximateMeanRotor(R_zeta, W.T());
  if(nHat_t1.dot(R_eps*Quaternions::xHat*R_eps.inverse()) < 0) {
    // std::cerr << __FILE__ << ":" << __LINE__ << ": Rotating by pi/2 about the z axis initially." << std::endl;
    R_eps = R_eps * Quaternions::exp((M_PI/2.)*Quaternions::zHat);
  }
}

/// Fix the orientation of the corotating frame by optimizing over a range of times.
GWFrames::Waveform& GWFrames::Waveform::AlignDecompositionFrameToModes(const double t1, const double t2,
                                                                       const Quaternions::Quaternion& nHat_t1,
                                                                       const std::vector<int>& Lmodes) {
  ///
  /// \param t1 Beginning of time range over which the alignment should happen
  /// \param t2 End of time range over which the alignment should happen
  /// \param nHat_t1 The approximate direction of nHat at t1
  /// \param Lmodes Lmodes to use in computing \f$<LL>\f$
  ///
  /// The corotating frame is only defined up to some constant rotor
  /// R_eps; if R_corot is corotating, then so is R_corot*R_eps.  This
  /// function uses that freedom to ensure that the frame is aligned
  /// with the Waveform modes as well as possible across the given
  /// time range.  In particular, it ensures that the Z axis of the
  /// frame in which the decomposition is done is along the dominant
  /// eigenvector of \f$<LL>\f$ (suggested by O'Shaughnessy et al.),
  /// and the phase of the (2,2) mode is zero.  These two conditions
  /// only give us axes, but we need vectors to fully specify the
  /// frame.  So we also impose the condition that the eigenvector is
  /// more parallel to the angular velocity of the waveform than
  /// anti-parallel, and the X axis of the rotated frame is more
  /// parallel to the input nHat_t1 than anti-parallel.  These
  /// conditions are imposed as accurately as possible across the
  /// range of times (t1, t2).
  ///
  /// If Lmodes is empty (default), all L modes are used.  Setting
  /// Lmodes to [2] or [2,3,4], for example, restricts the range of
  /// the sum.
  ///

  // Find the appropriate rotation
  Quaternion R_eps;
  GetAlignmentOfDecompositionFrameToModes(t1, t2, nHat_t1, R_eps, Lmodes);

  // Record what happened
  history << "this->AlignDecompositionFrameToModes(" << std::setprecision(16)
          << t1 << ", " << t2 << ", " << nHat_t1 << ", " << Lmodes << ");  # R_eps=" << R_eps << std::endl;

  // Now, apply the rotation
  this->RotateDecompositionBasis(R_eps);

  return *this;
}


/// Get the rotor needed to align this waveform's frame to the other's at the given time
void GWFrames::Waveform::GetAlignmentOfFrame(const Waveform& A, const double t_fid, Quaternion& R_delta) const {
  ///
  /// \param A Fixed Waveform in corotating frame to which this Waveform is aligned
  /// \param t_fid Fiducial time at which to equate frames
  /// \param R_delta Returned rotor
  ///
  /// This function simply finds the rotation necessary to align this
  /// waveform's frame to the other at the fiducial time, rather than
  /// applying it.  This is called by `AlignFrame` and probably does
  /// not need to be called directly; see that function's
  /// documentation for more details.
  ///
  /// \sa AlignFrame
  ///

  if(frameType != A.frameType && (frameType!=GWFrames::Corotating || frameType!=GWFrames::Coprecessing)) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nError: GetAlignmentOfFrame assumes that the two Waveforms are in the same type of frame, and that"
              << "\n       that frame is physically meaningful (either " << GWFrames::WaveformFrameNames[GWFrames::Corotating]
              << " or " << GWFrames::WaveformFrameNames[GWFrames::Coprecessing] << "),"
              << "\n       so that it makes sense to align the frames."
              << "\n       This Waveform is in the " << GWFrames::WaveformFrameNames[frameType] << " frame,"
              << "\n       The Waveform in the argument is in the " << GWFrames::WaveformFrameNames[A.frameType] << " frame.\n"
              << std::endl;
    throw(GWFrames_WrongFrameType);
  }

  const Waveform& B = *this;
  const Quaternion RA = Squad(A.Frame(), A.T(), vector<double>(1, t_fid))[0];
  const Quaternion RB = Squad(B.Frame(), B.T(), vector<double>(1, t_fid))[0];
  R_delta = RA * RB.inverse();
}

/// Change this Waveform by aligning the frame to the other's at the given time
GWFrames::Waveform& GWFrames::Waveform::AlignFrame(const GWFrames::Waveform& A, const double t_fid) {
  ///
  /// \param A Fixed Waveform in corotating frame to which this Waveform is aligned
  /// \param t_fid Fiducial time at which to equate frames
  ///
  /// Note that this function operates in place; the Waveform to which
  /// it is applied will change.  However, the modes are not altered;
  /// only the `frame` data is.
  ///
  /// As noted above, it is implicitly assumed that both Waveforms are
  /// in their corotating frames, with the modes appropriately aligned
  /// to the frames at t_fid.  The assumption is that the frames
  /// actually represent something physically meaningful, so that it
  /// is meaningful to insist that they be the same.
  ///
  /// Then, this function aligns the frames at t_fid by multiplying
  /// this->frame on the left by a constant rotor such that
  /// this->frame at t_fid is exactly A.frame at t_fid.  The resulting
  /// frame is now corotating with an angular-velocity vector that has
  /// been rotated by that constant rotor, relative to the inertial
  /// basis.
  ///
  /// \sa AlignDecompositionFrameToModes
  ///

  // Get the necessary rotation
  Quaternion R_delta;
  GetAlignmentOfFrame(A, t_fid, R_delta);

  // Apply the shift in frame
  SetFrame(R_delta * Frame());

  // Record what happened
  history << "this->AlignFrame(A, " << std::setprecision(16) << t_fid << ");  # R_delta=" << R_delta << std::endl;

  return *this;
}

/// Get time and frame offset for alignment over extended region.
void GWFrames::Waveform::GetAlignmentOfTimeAndFrame(const Waveform& A, const double t1, const double t2, double& deltat, Quaternion& R_delta) const {
  ///
  /// \param A Fixed Waveform in corotating frame to which this Waveform is aligned
  /// \param t1 Initial time of region over which differences are minimized
  /// \param t2 Final time of region over which differences are minimized
  /// \param deltat Returned time offset
  /// \param R_delta Returned rotation offset
  ///
  /// This function simply finds the time and rotation shifts
  /// necessary to align this waveform to the other at the fiducial
  /// time, rather than applying it.  This is called by
  /// `AlignTimeAndFrame` and probably does not need to be called
  /// directly; see that function's documentation for more details.
  ///
  /// In particular, note that this function is basically just a
  /// wrapper for the `Quaternions::OptimalAlignment` function.
  ///
  /// \sa AlignTimeAndFrame
  ///

  if(frameType != A.frameType || (frameType!=GWFrames::Corotating && frameType!=GWFrames::Coprecessing)) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nError: GetAlignmentOfTimeAndFrame assumes that the two Waveforms are in the type of same frame, and that"
              << "\n       that frame is physically meaningful (either " << GWFrames::WaveformFrameNames[GWFrames::Corotating]
              << " or " << GWFrames::WaveformFrameNames[GWFrames::Coprecessing] << "),"
              << "\n       so that it makes sense to align the frames."
              << "\n       This Waveform is in the " << GWFrames::WaveformFrameNames[frameType] << " frame,"
              << "\n       The Waveform in the argument is in the " << GWFrames::WaveformFrameNames[A.frameType] << " frame.\n"
              << "\n       Also note that each waveform should have its decomposition frame aligned to its modes (see `AlignDecompositionFrameToModes`).\n"
              << std::endl;
    throw(GWFrames_WrongFrameType);
  }

  const Waveform& B = *this;

  Quaternions::OptimalAlignment(t1, t2, A.Frame(), A.T(), B.Frame(), B.T(), deltat, R_delta);

  return;
}

/// Align time and frame over extended region.
GWFrames::Waveform& GWFrames::Waveform::AlignTimeAndFrame(const GWFrames::Waveform& A, const double t1, const double t2) {
  ///
  /// \param A Fixed Waveform in corotating frame to which this Waveform is aligned
  /// \param t1 Initial time of region over which differences are minimized
  /// \param t2 Final time of region over which differences are minimized
  ///
  /// Note that this function operates in place; the Waveform to which
  /// it is applied will change.  However, the modes are not altered;
  /// only the `t` and `frame` data are.
  ///
  /// The times `t1` and `t2` are measured relative to the time in
  /// Waveform `A`, and all are left fixed; only this Waveform is
  /// shifted (in time and orientation) to achieve alignment.
  ///
  /// It is implicitly assumed that both Waveforms are in their
  /// corotating frames, with the modes appropriately aligned to the
  /// frames using `AlignDecompositionFrameToModes` at some fiducial
  /// time at roughly the average of `t1` and `t2`.  The assumption is
  /// that the frames then actually represent something physically
  /// meaningful, so that it is meaningful to insist that they be the
  /// same.
  ///
  /// Also, it is assumed that the time data for the two waveforms are
  /// fairly closely aligned.  In particular, the minimization
  /// algorithm searches over time offsets of magnitude (t2-t1)/2.0 or
  /// less.  So, basically, the time data for this Waveform must be
  /// within \f$\pm (t2-t1)/2\f$ of the "correct" result.
  ///
  /// Then, this function adjust the time and orientation of this
  /// Waveform, so that the difference between the two frames is
  /// minimized.  That difference is measured by finding the rotor
  /// R_Delta required to rotate one frame into the other, taking the
  /// angle of that rotor, and integrating over the region [t1, t2].
  ///
  /// Relative to the inertial basis, the physical measurables
  /// (angular-velocity vector and dominant eigenvector of \f$<LL>\f$)
  /// of this Waveform are rotated.
  ///
  /// \sa AlignDecompositionFrameToModes
  ///

  // Apply time shift and rotation
  double deltat;
  Quaternion R_delta;
  GetAlignmentOfTimeAndFrame(A, t1, t2, deltat, R_delta);
  SetTime(T() + deltat);
  SetFrame(R_delta * Frame());

  // Record what happened
  history << "this->AlignTimeAndFrame(A, " << std::setprecision(16) << t1 << ", " << t2 << ");  # deltat=" << deltat << "; R_delta=" << R_delta << std::endl;

  return *this;
}


/// Return a Waveform with differences between the two inputs.
GWFrames::Waveform GWFrames::Waveform::Compare(const GWFrames::Waveform& A, const double MinTimeStep, const double MinTime) const {
  /// This function simply subtracts the data in this Waveform from
  /// the data in Waveform A, and finds the rotation needed to take
  /// this frame into frame A.  Note that the waveform data are stored
  /// as complex numbers, rather than as modulus and phase.

  // Make B a convenient alias for *this
  const GWFrames::Waveform& B = *this;

  // Check to see if the frameTypes are the same
  if(frameType != A.frameType) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nWarning:"
              << "\n       This Waveform is in the " << GWFrames::WaveformFrameNames[frameType] << " frame,"
              << "\n       The Waveform in the argument is in the " << GWFrames::WaveformFrameNames[A.frameType] << " frame."
              << "\n       Comparing them probably does not make sense.\n"
              << std::endl;
  }

  // Make sure we have the same number of modes in the input data
  if(NModes() != B.NModes()) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Trying to Compare Waveforms with mismatched LM data."
              << "\nA.NModes()=" << A.NModes() << "\tB.NModes()=" << B.NModes() << std::endl;
    throw(GWFrames_WaveformMissingLMIndex);
  }
  // We'll put all the data in a new Waveform C
  GWFrames::Waveform C;
  // Store both old histories in C's
  C.history << "B.Compare(A)\n"
            << "#### A.history.str():\n" << A.history.str()
            << "#### B.history.str():\n" << B.history.str()
            << "#### End of old histories from `Compare`" << std::endl;
  // The new time axis will be the intersection of the two old ones
  C.t = GWFrames::Intersection(A.t, B.t);
  // Copy the basics
  C.frameType = frameType;
  C.dataType = dataType;
  C.rIsScaledOut = rIsScaledOut;
  C.mIsScaledOut = mIsScaledOut;
  // We'll assume that {A.lm}=={B.lm} as sets, and account for
  // disordering below
  C.lm = A.lm;
  // Process the frame, depending on the sizes of the input frames
  if(A.Frame().size()>1 && B.Frame().size()>1) {
    // Find the frames interpolated to the appropriate times
    const vector<Quaternion> Aframe = Quaternions::Squad(A.frame, A.t, C.t);
    const vector<Quaternion> Bframe = Quaternions::Squad(B.frame, B.t, C.t);
    // Assign the data
    C.frame.resize(C.NTimes());
    for(unsigned int i_t=0; i_t<C.t.size(); ++i_t) {
      C.frame[i_t] = Aframe[i_t] * Quaternions::inverse(Bframe[i_t]);
    }
  } else if(A.Frame().size()==1 && B.Frame().size()>1) {
    // Find the frames interpolated to the appropriate times
    const vector<Quaternion> Bframe = Quaternions::Squad(B.frame, B.t, C.t);
    // Assign the data
    C.frame.resize(C.NTimes());
    for(unsigned int i_t=0; i_t<C.t.size(); ++i_t) {
      C.frame[i_t] = A.Frame(0) * Quaternions::inverse(Bframe[i_t]);
    }
  } else if(A.Frame().size()>1 && B.Frame().size()==1) {
    // Find the frames interpolated to the appropriate times
    const vector<Quaternion> Aframe = Quaternions::Squad(A.frame, A.t, C.t);
    // Assign the data
    C.frame.resize(C.NTimes());
    for(unsigned int i_t=0; i_t<C.t.size(); ++i_t) {
      C.frame[i_t] = Aframe[i_t] * Quaternions::inverse(B.Frame(0));
    }
  } else if(A.Frame().size()==1 && B.Frame().size()==1) {
    // Assign the data
    C.frame = vector<Quaternion>(1,A.Frame(0) * Quaternions::inverse(B.Frame(0)));
  } else if(A.Frame().size()==0 && B.Frame().size()==1) {
    // Assign the data
    C.frame = vector<Quaternion>(1,Quaternions::inverse(B.Frame(0)));
  } else if(A.Frame().size()==1 && B.Frame().size()==0) {
    // Assign the data
    C.frame = vector<Quaternion>(1,A.Frame(0));
  } // else, leave the frame data empty
  // Reserve space for the data
  C.data.resize(C.lm.size(), C.t.size());
  // Construct the GSL interpolators for the data
  gsl_interp_accel* accReA = gsl_interp_accel_alloc();
  gsl_interp_accel* accImA = gsl_interp_accel_alloc();
  gsl_interp_accel* accReB = gsl_interp_accel_alloc();
  gsl_interp_accel* accImB = gsl_interp_accel_alloc();
  gsl_spline* splineReA = gsl_spline_alloc(gsl_interp_cspline, A.NTimes());
  gsl_spline* splineImA = gsl_spline_alloc(gsl_interp_cspline, A.NTimes());
  gsl_spline* splineReB = gsl_spline_alloc(gsl_interp_cspline, B.NTimes());
  gsl_spline* splineImB = gsl_spline_alloc(gsl_interp_cspline, B.NTimes());
  // Now loop over each mode filling in the waveform data
  for(unsigned int Mode=0; Mode<A.NModes(); ++Mode) {
    // Assume that all the ell,m data are the same, but not necessarily in the same order
    const unsigned int BMode = B.FindModeIndex(A.lm[Mode][0], A.lm[Mode][1]);
    // Extract the real and imaginary parts of the data separately for GSL
    const vector<double> ReA = A.Re(Mode);
    const vector<double> ImA = A.Im(Mode);
    const vector<double> ReB = B.Re(BMode);
    const vector<double> ImB = B.Im(BMode);
    // Initialize the interpolators for this data set
    gsl_spline_init(splineReA, &(A.t)[0], &ReA[0], A.NTimes());
    gsl_spline_init(splineImA, &(A.t)[0], &ImA[0], A.NTimes());
    gsl_spline_init(splineReB, &(B.t)[0], &ReB[0], B.NTimes());
    gsl_spline_init(splineImB, &(B.t)[0], &ImB[0], B.NTimes());
    // Assign the data from the transition
    for(unsigned int i_t=0; i_t<C.t.size(); ++i_t) {
      C.data[Mode][i_t] = complex<double>( gsl_spline_eval(splineReA, C.t[i_t], accReA), gsl_spline_eval(splineImA, C.t[i_t], accImA) )
        - complex<double>( gsl_spline_eval(splineReB, C.t[i_t], accReB), gsl_spline_eval(splineImB, C.t[i_t], accImB) );
    }
  }
  gsl_interp_accel_free(accReA);
  gsl_interp_accel_free(accImA);
  gsl_interp_accel_free(accReB);
  gsl_interp_accel_free(accImB);
  gsl_spline_free(splineReA);
  gsl_spline_free(splineImA);
  gsl_spline_free(splineReB);
  gsl_spline_free(splineImB);

  return C;
}

/// Local utility function
inline double TransitionFunction_Smooth(const double x) {
  /// This smoothly transitions from 0.0 for x<=0.0 to 1.0 for x>=1.0.
  /// The function is just the usual transition function based on the
  /// familiar smooth but non-analytic function.
  return ( x<=0.0 ? 0.0 : ( x>=1.0 ? 1.0 : 1.0/(1.0+std::exp(1.0/(x-1.0) + 1.0/x)) ) );
}

/// Hybridize this Waveform with another.
GWFrames::Waveform GWFrames::Waveform::Hybridize(const GWFrames::Waveform& B, const double t1, const double t2, const double tMinStep) const {
  ///
  /// \param B Second Waveform to hybridize with
  /// \param t1 Beginning of time over which to transition
  /// \param t2 End of time over which to transition
  /// \param tMinStep Lower limit on time step appearing in the output
  ///
  /// This function simply takes two Waveforms and blends them
  /// together.  In particular, it does not align the Waveforms; that
  /// is assumed to have been done already.
  ///
  /// The transition function is a \f$C^\infty\f$ function, meaning
  /// that the output data has exactly this Waveform's data before
  /// `t1`, exactly Waveform `B`'s data after t2, and a smooth blend
  /// in between.
  ///
  /// Note that this function does NOT operate in place; a new
  /// Waveform object is constructed and returned.

  // Make A a convenient alias
  const GWFrames::Waveform& A = *this;

  // Check to see if the various type flags agree
  if(A.spinweight != B.spinweight) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nWarning:"
              << "\n       This Waveform has spin weight " << A.spinweight << "."
              << "\n       The Waveform in the argument has spin weight " << B.spinweight << "."
              << "\n       Hybridizing them probably does not make sense.\n"
              << std::endl;
  }
  if(A.boostweight != B.boostweight) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nWarning:"
              << "\n       This Waveform has boost weight " << A.boostweight << "."
              << "\n       The Waveform in the argument has boost weight " << B.boostweight << "."
              << "\n       Hybridizing them probably does not make sense.\n"
              << std::endl;
  }
  if(A.frameType != B.frameType) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nWarning:"
              << "\n       This Waveform is in the " << GWFrames::WaveformFrameNames[A.frameType] << " frame."
              << "\n       The Waveform in the argument is in the " << GWFrames::WaveformFrameNames[B.frameType] << " frame."
              << "\n       Hybridizing them probably does not make sense.\n"
              << std::endl;
  }
  if(A.dataType != B.dataType) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nWarning:"
              << "\n       This Waveform has data type " << GWFrames::WaveformDataNames[A.dataType] << "."
              << "\n       The Waveform in the argument has data type " << GWFrames::WaveformDataNames[B.dataType] << "."
              << "\n       Hybridizing them probably does not make sense.\n"
              << std::endl;
  }
  if(A.rIsScaledOut != B.rIsScaledOut) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nWarning:"
              << "\n       This Waveform claims radius is " << (A.rIsScaledOut ? "" : "not ") << "scaled out."
              << "\n       The Waveform in the argument claims radius is " << (B.rIsScaledOut ? "" : "not ") << "scaled out."
              << "\n       Hybridizing them probably does not make sense.\n"
              << std::endl;
  }
  if(A.mIsScaledOut != B.mIsScaledOut) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nWarning:"
              << "\n       This Waveform claims mass is " << (A.mIsScaledOut ? "" : "not ") << "scaled out."
              << "\n       The Waveform in the argument claims mass is " << (B.mIsScaledOut ? "" : "not ") << "scaled out."
              << "\n       Hybridizing them probably does not make sense.\n"
              << std::endl;
  }

  // Make sure we have the same number of modes in the input data
  if(NModes() != B.NModes()) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Trying to Align Waveforms with mismatched LM data."
              << "\nA.NModes()=" << A.NModes() << "\tB.NModes()=" << B.NModes() << std::endl;
    throw(GWFrames_WaveformMissingLMIndex);
  }

  // Make sure we have sufficient times for the requested hybrid
  {
    const double t1A = A.T(0);
    const double t1B = B.T(0);
    const double t2A = A.T(A.NTimes()-1);
    const double t2B = B.T(B.NTimes()-1);
    if(t1<t1A || t1<t1B || t2>t2A || t2>t2B) {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
                << ": These Waveforms do not overlap on the requested range."
                << "\nA.T(0)=" << t1A << "\tB.T(0)" << t1B << "\tt1=" << t1
                << "\nA.T(-1)=" << t2A << "\tB.T(-1)" << t2B << "\tt2=" << t2 << std::endl;
      throw(GWFrames_EmptyIntersection);
    }
  }

  // We'll put all the data in a new Waveform C
  GWFrames::Waveform C;
  C.spinweight = A.spinweight;
  C.boostweight = A.boostweight;
  C.frameType = A.frameType;
  C.dataType = A.dataType;
  C.rIsScaledOut = A.rIsScaledOut;
  C.mIsScaledOut = A.mIsScaledOut;
  C.history << "A.Hybridize(B, " << t1 << ", " << t2 << ", " << tMinStep << ")\n"
            << "#### A.history.str():\n" << A.history.str()
            << "#### B.history.str():\n" << B.history.str()
            << "#### End of old histories from `Hybridize`" << std::endl;
  C.t = GWFrames::Union(A.t, B.t, tMinStep);
  // We'll assume that A.lm==B.lm, though we'll account for disordering below
  C.lm = A.lm;
  // Make sure the time stops at the end of B's time (in case A extended further)
  int i_t=C.t.size()-1;
  while(C.T(i_t)>B.t.back() && i_t>0) { --i_t; }
  C.t.erase(C.t.begin()+i_t, C.t.end());
  // Reserve space for everything
  C.frame.resize(C.NTimes());
  C.data.resize(C.lm.size(), C.t.size());
  // Find the indices of the transition points
  unsigned int J01=0, J12=C.NTimes()-1;
  while(C.T(J01)<t1 && J01<C.NTimes()) { J01++; }
  while(C.T(J12)>t2 && J12>0) { J12--; }
  const double T01 = C.T(J01);
  const double TransitionLength = C.T(J12)-T01;
  // Find the frames interpolated to the appropriate times
  vector<Quaternion> Aframe = Quaternions::Squad(A.frame, A.t, vector<double>(&(C.t)[0], &(C.t)[J12]));
  vector<Quaternion> Bframe = Quaternions::Squad(B.frame, B.t, vector<double>(&(C.t)[J01], &(C.t)[C.NTimes()-1]));
  // Assign the frame data from earliest part
  for(unsigned int j=0; j<J01; ++j) {
    C.frame[j] = Aframe[j];
  }
  // Assign the frame data from the transition
  for(unsigned int j=J01; j<J12; ++j) {
    const double Transition = TransitionFunction_Smooth((C.T(j)-T01)/TransitionLength);
    C.frame[j] = Quaternions::Slerp(Transition, Aframe[j], Bframe[j-J01]);
  }
  // Assign the frame data from the latest part
  for(unsigned int j=J12; j<C.NTimes(); ++j) {
    C.frame[j] = Bframe[j-J01];
  }
  // Construct the GSL interpolators for the data
  gsl_interp_accel* accReA = gsl_interp_accel_alloc();
  gsl_interp_accel* accImA = gsl_interp_accel_alloc();
  gsl_interp_accel* accReB = gsl_interp_accel_alloc();
  gsl_interp_accel* accImB = gsl_interp_accel_alloc();
  gsl_spline* splineReA = gsl_spline_alloc(gsl_interp_cspline, A.NTimes());
  gsl_spline* splineImA = gsl_spline_alloc(gsl_interp_cspline, A.NTimes());
  gsl_spline* splineReB = gsl_spline_alloc(gsl_interp_cspline, B.NTimes());
  gsl_spline* splineImB = gsl_spline_alloc(gsl_interp_cspline, B.NTimes());
  // Now loop over each mode filling in the waveform data
  for(unsigned int Mode=0; Mode<A.NModes(); ++Mode) {
    // Assume that all the ell,m data are the same, but not necessarily in the same order
    const unsigned int BMode = B.FindModeIndex(A.lm[Mode][0], A.lm[Mode][1]);
    // Extract the real and imaginary parts of the data separately for GSL
    const vector<double> ReA = A.Re(Mode);
    const vector<double> ImA = A.Im(Mode);
    const vector<double> ReB = B.Re(BMode);
    const vector<double> ImB = B.Im(BMode);
    // Initialize the interpolators for this data set
    gsl_spline_init(splineReA, &(A.t)[0], &ReA[0], A.NTimes());
    gsl_spline_init(splineImA, &(A.t)[0], &ImA[0], A.NTimes());
    gsl_spline_init(splineReB, &(B.t)[0], &ReB[0], B.NTimes());
    gsl_spline_init(splineImB, &(B.t)[0], &ImB[0], B.NTimes());
    // Assign the data from earliest part
    for(unsigned int j=0; j<J01; ++j) {
      C.data[Mode][j] = complex<double>( gsl_spline_eval(splineReA, C.t[j], accReA), gsl_spline_eval(splineImA, C.t[j], accImA) );
    }
    // Assign the data from the transition
    for(unsigned int j=J01; j<J12; ++j) {
      const double Transition = TransitionFunction_Smooth((C.T(j)-T01)/TransitionLength);
      C.data[Mode][j] = complex<double>( gsl_spline_eval(splineReA, C.t[j], accReA), gsl_spline_eval(splineImA, C.t[j], accImA) ) * (1.0-Transition)
        + complex<double>( gsl_spline_eval(splineReB, C.t[j], accReB), gsl_spline_eval(splineImB, C.t[j], accImB) ) * Transition ;
    }
    // Assign the data from the latest part
    for(unsigned int j=J12; j<C.NTimes(); ++j) {
      C.data[Mode][j] = complex<double>( gsl_spline_eval(splineReB, C.t[j], accReB), gsl_spline_eval(splineImB, C.t[j], accImB) );
    }
  }
  gsl_interp_accel_free(accReA);
  gsl_interp_accel_free(accImA);
  gsl_interp_accel_free(accReB);
  gsl_interp_accel_free(accImB);
  gsl_spline_free(splineReA);
  gsl_spline_free(splineImA);
  gsl_spline_free(splineReB);
  gsl_spline_free(splineImB);

  return C;
}


/// Evaluate Waveform at a particular sky location
std::vector<std::complex<double> > GWFrames::Waveform::EvaluateAtPoint(const double vartheta, const double varphi) const {
  ///
  /// \param vartheta Polar angle of detector
  /// \param varphi Azimuthal angle of detector
  ///
  /// Note that the input angle parameters are measured relative to
  /// the binary's coordinate system.  In particular, this will make
  /// no sense if the frame type is something other than inertial, and
  /// will fail if the `FrameType` is neither `UnknownFrameType` nor
  /// `Inertial`.
  ///

  if(frameType != GWFrames::Inertial) {
    if(frameType == GWFrames::UnknownFrameType) {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
                << "\nWarning: Asking for a Waveform in the " << GWFrames::WaveformFrameNames[frameType] << " frame to be evaluated at a point."
                << "\n         This should only be applied to Waveforms in the " << GWFrames::WaveformFrameNames[GWFrames::Inertial] << " frame.\n"
                << std::endl;
    } else {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
                << "\nError: Asking for a Waveform in the " << GWFrames::WaveformFrameNames[frameType] << " frame to be evaluated at a point."
                << "\n       This should only be applied to Waveforms in the " << GWFrames::WaveformFrameNames[GWFrames::Inertial] << " frame.\n"
                << std::endl;
      throw(GWFrames_WrongFrameType);
    }
  }

  const int NT = NTimes();
  const int NM = NModes();
  vector<complex<double> > d(NT, complex<double>(0.,0.));
  SphericalFunctions::SWSH Y(SpinWeight());
  Y.SetAngles(vartheta, varphi);

  for(int i_m=0; i_m<NM; ++i_m) {
    const int ell = LM(i_m)[0];
    const int m   = LM(i_m)[1];
    const complex<double> Ylm = Y(ell,m);
    for(int i_t=0; i_t<NT; ++i_t) {
      d[i_t] += Data(i_m, i_t) * Ylm;
    }
  }

  return d;
}

/// Evaluate Waveform at a particular sky location and an instant of time
std::complex<double> GWFrames::Waveform::EvaluateAtPoint(const double vartheta, const double varphi, const unsigned int i_t) const {
  ///
  /// \param vartheta Polar angle of detector
  /// \param varphi Azimuthal angle of detector
  /// \param i_t Index of time at which to evaluate
  ///
  /// Note that the input angle parameters are measured relative to
  /// the binary's coordinate system.  In particular, this will make
  /// no sense if the frame type is something other than inertial, and
  /// will fail if the `FrameType` is neither `UnknownFrameType` nor
  /// `Inertial`.
  ///

  if(frameType != GWFrames::Inertial) {
    if(frameType == GWFrames::UnknownFrameType) {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
                << "\nWarning: Asking for a Waveform in an " << GWFrames::WaveformFrameNames[frameType] << " frame to be evaluated at a point."
                << "\n         This should probably only be applied to Waveforms in the " << GWFrames::WaveformFrameNames[GWFrames::Inertial] << " frame."
                << "\n         But I'll trust you to know what you're doing.\n"
                << std::endl;
    } else {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
                << "\nError: Asking for a Waveform in the " << GWFrames::WaveformFrameNames[frameType] << " frame to be evaluated at a point."
                << "\n       This should only be applied to Waveforms in the " << GWFrames::WaveformFrameNames[GWFrames::Inertial] << " frame.\n"
                << std::endl;
      throw(GWFrames_WrongFrameType);
    }
  }

  const int NM = NModes();
  complex<double> d(0.,0.);
  SphericalFunctions::SWSH Y(SpinWeight());
  Y.SetAngles(vartheta, varphi);

  for(int i_m=0; i_m<NM; ++i_m) {
    const int ell = LM(i_m)[0];
    const int m   = LM(i_m)[1];
    const complex<double> Ylm = Y(ell,m);
    d += Data(i_m, i_t) * Ylm;
  }

  return d;
}

/// Evaluate Waveform at a particular sky location and an instant of time
std::complex<double> GWFrames::Waveform::InterpolateToPoint(const double vartheta, const double varphi, const double t_i) const {
  ///
  /// \param vartheta Polar angle of detector
  /// \param varphi Azimuthal angle of detector
  /// \param t_i New time to interpolate to
  ///
  /// Note that the input angle parameters are measured relative to
  /// the binary's coordinate system.  In particular, this will make
  /// no sense if the frame type is something other than inertial, and
  /// will fail if the `FrameType` is neither `UnknownFrameType` nor
  /// `Inertial`.
  ///

  if(frameType != GWFrames::Inertial) {
    if(frameType == GWFrames::UnknownFrameType) {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
                << "\nWarning: Asking for a Waveform in an " << GWFrames::WaveformFrameNames[frameType] << " frame to be evaluated at a point."
                << "\n         This should probably only be applied to Waveforms in the " << GWFrames::WaveformFrameNames[GWFrames::Inertial] << " frame."
                << "\n         But I'll trust you to know what you're doing.\n"
                << std::endl;
    } else {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
                << "\nError: Asking for a Waveform in the " << GWFrames::WaveformFrameNames[frameType] << " frame to be evaluated at a point."
                << "\n       This should only be applied to Waveforms in the " << GWFrames::WaveformFrameNames[GWFrames::Inertial] << " frame.\n"
                << std::endl;
      throw(GWFrames_WrongFrameType);
    }
  }

  if(NTimes()<4) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nError: " << NTimes() << " is not enough points to interpolate.\n"
              << std::endl;
    throw(GWFrames_BadWaveformInformation);
  }

  if(t_i<T(0)) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nError: (t_i=" << t_i << ") is earlier than the earliest time in the data (T(0)=" << T(0) << ").\n"
              << std::endl;
    throw(GWFrames_ValueError);
  }

  if(t_i>T(NTimes()-1)) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nError: (t_i=" << t_i << ") is later than the latest time in the data (T(" << NTimes()-1 << ")=" << T(NTimes()-1) << ").\n"
              << std::endl;
    throw(GWFrames_ValueError);
  }

  const int NM = NModes();
  vector<double> dRe(4), dIm(dRe.size());
  SphericalFunctions::SWSH Y(SpinWeight());
  Y.SetAngles(vartheta, varphi);

  // Find a series of time steps around the requested time
  unsigned int i_t_0(std::max(int(Quaternions::hunt(t, t_i))-1, 0));
  if(i_t_0+3>=NTimes()) {
    i_t_0 = NTimes()-4;
  }

  // Evaluate at this point at a series of times around the requested time
  for(int i_m=0; i_m<NM; ++i_m) {
    const int ell = LM(i_m)[0];
    const int m   = LM(i_m)[1];
    const complex<double> Ylm = Y(ell,m);
    for(unsigned int i_t=0; i_t<dRe.size(); ++i_t) {
      const complex<double> val = Data(i_m, i_t+i_t_0) * Ylm;
      dRe[i_t] += std::real(val);
      dIm[i_t] += std::imag(val);
    }
  }

  // Now interpolate in time
  gsl_interp_accel* accRe = gsl_interp_accel_alloc();
  gsl_interp_accel* accIm = gsl_interp_accel_alloc();
  gsl_spline* splineRe = gsl_spline_alloc(gsl_interp_cspline, dRe.size());
  gsl_spline* splineIm = gsl_spline_alloc(gsl_interp_cspline, dRe.size());
  gsl_spline_init(splineRe, &(t)[i_t_0], &dRe[0], dRe.size());
  gsl_spline_init(splineIm, &(t)[i_t_0], &dIm[0], dRe.size());
  const complex<double> value( gsl_spline_eval(splineRe, t_i, accRe), gsl_spline_eval(splineIm, t_i, accIm) );
  gsl_interp_accel_free(accRe);
  gsl_interp_accel_free(accIm);
  gsl_spline_free(splineRe);
  gsl_spline_free(splineIm);

  return value;
}


/// Output Waveform object to data file.
const GWFrames::Waveform& GWFrames::Waveform::Output(const std::string& FileName, const unsigned int precision) const {
  const std::string Descriptor = DescriptorString();
  ofstream ofs(FileName.c_str(), ofstream::out);
  ofs << setprecision(precision) << flush;
  ofs << history.str() << "this->Output(" << FileName << ", " << precision << ")" << endl;
  ofs << "# [1] = Time" << endl;
  for(unsigned int i_m=0; i_m<NModes(); ++i_m) {
    ofs << "# [" << 2*i_m+2 << "] = Re{" << Descriptor << "(" << lm[i_m][0] << "," << lm[i_m][1] << ")}" << endl;
    ofs << "# [" << 2*i_m+3 << "] = Im{" << Descriptor << "(" << lm[i_m][0] << "," << lm[i_m][1] << ")}" << endl;
  }
  for(unsigned int i_t=0; i_t<NTimes(); ++i_t) {
    ofs << T(i_t) << " ";
    for(unsigned int i_m=0; i_m<NModes()-1; ++i_m) {
      ofs << data[i_m][i_t].real() << " " << data[i_m][i_t].imag() << " ";
    }
    ofs << data[NModes()-1][i_t].real() << " " << data[NModes()-1][i_t].imag() << endl;
  }
  ofs.close();
  return *this;
}


// /// Correct the error in RWZ extraction from older SpEC files
// GWFrames::Waveform& GWFrames::Waveform::HackSpECSignError() {
//   // h(ell,m) -> (-1)^m * h(ell,-m)*
//   unsigned int i_m = 0;
//   unsigned int N_m = NModes();
//   unsigned int N_t = NTimes();
//   for(int ell=2; ; ++ell) {
//     if(i_m>=N_m) {
//       break;
//     }
//     for(int m=1; m<=ell; ++m) {
//       // First, swap the m and -m modes
//       unsigned int iPositive = FindModeIndex(ell,m);
//       unsigned int iNegative = FindModeIndex(ell,-m);
//       lm[iNegative].swap(lm[iPositive]);
//       // Next, swap the signs as appropriate
//       if(m%2==0) {
//         for(unsigned int i_t=0; i_t<N_t; ++i_t) {
//           data[iPositive][i_t].imag() *= -1;
//           data[iNegative][i_t].imag() *= -1;
//         }
//       } else {
//         for(unsigned int i_t=0; i_t<N_t; ++i_t) {
//           data[iPositive][i_t].real() *= -1;
//           data[iNegative][i_t].real() *= -1;
//         }
//       }
//       ++i_m; // For positive m
//       ++i_m; // For negative m
//     }
//     // Don't forget about m=0
//     unsigned int iZero = FindModeIndex(ell,0);
//     for(unsigned int i_t=0; i_t<N_t; ++i_t) {
//       data[iZero][i_t].imag() *= -1;
//     }
//     ++i_m; // For m=0
//   }
//   return *this;
// }


// GWFrames::Waveform GWFrames::Waveform::operator+(const GWFrames::Waveform& B) const { return BinaryOp<std::plus<std::complex<double> > >(B); }
// GWFrames::Waveform GWFrames::Waveform::operator-(const GWFrames::Waveform& B) const { return BinaryOp<std::minus<std::complex<double> > >(B); }

GWFrames::Waveform GWFrames::Waveform::operator+(const GWFrames::Waveform& B) const {
  const Waveform& A = *this;

  if(A.spinweight != B.spinweight) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nError: Asking for the sum of two Waveform objects with different spin weights."
              << "\n       A.SpinWeight()=" << A.SpinWeight() << "\tB.SpinWeight()=" << B.SpinWeight()
              << std::endl;
    throw(GWFrames_BadWaveformInformation);
  }

  if(A.frameType != GWFrames::Inertial || B.frameType != GWFrames::Inertial) {
    if(A.frameType != B.frameType) {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
                << "\nError: Asking for the sum of Waveforms in " << GWFrames::WaveformFrameNames[A.frameType]
                << " and " << GWFrames::WaveformFrameNames[B.frameType] << " frames."
                << "\n       This should only be applied to Waveforms in the same frame.\n"
                << std::endl;
      throw(GWFrames_WrongFrameType);
    } else if(A.frame.size() != B.frame.size()) {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
                << "\nError: Asking for the sum of Waveforms with " << A.frame.size() << " and " << B.frame.size() << " frame data points."
                << "\n       This should only be applied to Waveforms in the same frame.\n"
                << std::endl;
      throw(GWFrames_WrongFrameType);
    }
  }

  if(A.NTimes() != B.NTimes()) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nError: Asking for the sum of two Waveform objects with different time data."
              << "\n       A.NTimes()=" << A.NTimes() << "\tB.NTimes()=" << B.NTimes()
              << "\n       Interpolate to a common set of times first.\n"
              << std::endl;
    throw(GWFrames_MatrixSizeMismatch);
  }

  // This will be the new object holding the multiplied data
  GWFrames::Waveform C(A);

  // Store both old histories in C's
  C.history << "*this = *this+B\n"
            << "#### B.history.str():\n" << B.history.str()
            << "#### End of old histories from `A+B`" << std::endl;

  // Do the work of addition
  const unsigned int ntimes = NTimes();
  const unsigned int nmodes = NModes();
  for(unsigned int i_A=0; i_A<nmodes; ++i_A) {
    const unsigned int i_B = B.FindModeIndex(lm[i_A][0], lm[i_A][1]);
    for(unsigned int i_t=0; i_t<ntimes; ++i_t) {
      C.SetData(i_A, i_t, A.Data(i_A,i_t)+B.Data(i_B,i_t));
    }
  }

  return C;
}

GWFrames::Waveform GWFrames::Waveform::operator-(const GWFrames::Waveform& B) const {
  const Waveform& A = *this;

  if(A.spinweight != B.spinweight) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nError: Asking for the difference of two Waveform objects with different spin weights."
              << "\n       A.SpinWeight()=" << A.SpinWeight() << "\tB.SpinWeight()=" << B.SpinWeight()
              << std::endl;
    throw(GWFrames_BadWaveformInformation);
  }

  if(A.frameType != GWFrames::Inertial || B.frameType != GWFrames::Inertial) {
    if(A.frameType != B.frameType) {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
                << "\nError: Asking for the difference of Waveforms in " << GWFrames::WaveformFrameNames[A.frameType]
                << " and " << GWFrames::WaveformFrameNames[B.frameType] << " frames."
                << "\n       This should only be applied to Waveforms in the same frame.\n"
                << std::endl;
      throw(GWFrames_WrongFrameType);
    } else if(A.frame.size() != B.frame.size()) {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
                << "\nError: Asking for the difference of Waveforms with " << A.frame.size() << " and " << B.frame.size() << " frame data points."
                << "\n       This should only be applied to Waveforms in the same frame.\n"
                << std::endl;
      throw(GWFrames_WrongFrameType);
    }
  }

  if(A.NTimes() != B.NTimes()) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nError: Asking for the difference of two Waveform objects with different time data."
              << "\n       A.NTimes()=" << A.NTimes() << "\tB.NTimes()=" << B.NTimes()
              << "\n       Interpolate to a common set of times first.\n"
              << std::endl;
    throw(GWFrames_MatrixSizeMismatch);
  }

  // This will be the new object holding the multiplied data
  GWFrames::Waveform C(A);

  // Store both old histories in C's
  C.history << "*this = *this-B\n"
            << "#### B.history.str():\n" << B.history.str()
            << "#### End of old histories from `A-B`" << std::endl;

  // Do the work of addition
  const unsigned int ntimes = NTimes();
  const unsigned int nmodes = NModes();
  for(unsigned int i_A=0; i_A<nmodes; ++i_A) {
    const unsigned int i_B = B.FindModeIndex(lm[i_A][0], lm[i_A][1]);
    for(unsigned int i_t=0; i_t<ntimes; ++i_t) {
      C.SetData(i_A, i_t, A.Data(i_A,i_t)-B.Data(i_B,i_t));
    }
  }

  return C;
}

GWFrames::Waveform GWFrames::Waveform::operator*(const double b) const {
  const Waveform& A = *this;

  // This will be the new object holding the multiplied data
  GWFrames::Waveform C(A);

  // Record the activity
  C.history << "*this = (*this) * " << b << std::endl;

  // Do the work of addition
  const unsigned int ntimes = NTimes();
  const unsigned int nmodes = NModes();
  for(unsigned int i_A=0; i_A<nmodes; ++i_A) {
    for(unsigned int i_t=0; i_t<ntimes; ++i_t) {
      C.SetData(i_A, i_t, A.Data(i_A,i_t)*b);
    }
  }

  return C;
}

GWFrames::Waveform GWFrames::Waveform::operator/(const double b) const {
  const Waveform& A = *this;

  // This will be the new object holding the multiplied data
  GWFrames::Waveform C(A);

  // Record the activity
  C.history << "*this = (*this) / " << b << std::endl;

  // Do the work of addition
  const unsigned int ntimes = NTimes();
  const unsigned int nmodes = NModes();
  for(unsigned int i_A=0; i_A<nmodes; ++i_A) {
    for(unsigned int i_t=0; i_t<ntimes; ++i_t) {
      C.SetData(i_A, i_t, A.Data(i_A,i_t)/b);
    }
  }

  return C;
}


GWFrames::Waveform GWFrames::Waveform::operator*(const GWFrames::Waveform& B) const { return BinaryOp<std::multiplies<std::complex<double> > >(B); }
GWFrames::Waveform GWFrames::Waveform::operator/(const GWFrames::Waveform& B) const { return BinaryOp<std::divides<std::complex<double> > >(B); }


/// Newman-Penrose edth operator
GWFrames::Waveform GWFrames::Waveform::NPEdth() const {
  /// This operator is the one defined by Newman and Penrose (1966)
  /// and further described by Goldberg et al. (1967).  It raises the
  /// spin weight of any field on the sphere by 1.  Note that this
  /// operator does not preserve boost weights in any nice way --
  /// except in special cases.  The GHP version does.  Note that, in
  /// this implementation, the only difference between the NP and GHP
  /// versions is the factor of \f$\sqrt{2}\f$.  The additional GHP
  /// term that keeps the boost weight meaningful is zero in any given
  /// frame -- though it transforms nontrivially.
  ///
  /// Note that the boost weight is set to the value of `WeightError`,
  /// which is just meant to be large enough that it will give
  /// improbable values if used.  This is not fool-proof.
  ///
  /// \sa NPEdthBar
  /// \sa GHPEdth
  /// \sa GHPEdthBar
  /// \sa IntegrateNPEdth
  /// \sa IntegrateNPEdthBar
  /// \sa IntegrateGHPEdth
  /// \sa IntegrateGHPEdthBar

  const Waveform& A = *this;
  const int NModes = A.NModes();
  const int NTimes = A.NTimes();
  const int s = A.SpinWeight();
  const int lMin = std::abs(s+1);
  Waveform EdthA(A);
  EdthA.history << "this->NPEdth();" << endl;

  for(int i_m=0; i_m<NModes; ++i_m) {
    const int l = A.LM(i_m)[0];
    if(l<lMin) {
      for(int i_t=0; i_t<NTimes; ++i_t) {
        EdthA.SetData(i_m, i_t, 0.0);
      }
    } else {
      const double factor = std::sqrt((l-s)*(l+s+1));
      for(int i_t=0; i_t<NTimes; ++i_t) {
        EdthA.SetData(i_m, i_t, EdthA.Data(i_m, i_t)*factor);
      }
    }
  }

  EdthA.SetSpinWeight(A.SpinWeight()+1);
  EdthA.SetBoostWeight(WeightError);

  return EdthA;
}

/// Newman-Penrose edth operator conjugate
GWFrames::Waveform GWFrames::Waveform::NPEdthBar() const {
  /// This operator is the one defined by Newman and Penrose (1966)
  /// and further described by Goldberg et al. (1967).  It lowers the
  /// spin weight of any field on the sphere by 1.  Note that this
  /// operator does not preserve boost weights in any nice way --
  /// except in special cases.  The GHP version does.  Note that, in
  /// this implementation, the only difference between the NP and GHP
  /// versions is the factor of \f$\sqrt{2}\f$.  The additional GHP
  /// term that keeps the boost weight meaningful is zero in any given
  /// frame -- though it transforms nontrivially.
  ///
  /// Note that the boost weight is set to the value of `WeightError`,
  /// which is just meant to be large enough that it will give
  /// improbable values if used.  This is not fool-proof.
  ///
  /// \sa NPEdth
  /// \sa GHPEdth
  /// \sa GHPEdthBar
  /// \sa IntegrateNPEdth
  /// \sa IntegrateNPEdthBar
  /// \sa IntegrateGHPEdth
  /// \sa IntegrateGHPEdthBar

  const Waveform& A = *this;
  const int NModes = A.NModes();
  const int NTimes = A.NTimes();
  const int s = A.SpinWeight();
  const int lMin = std::abs(s-1);
  Waveform EdthBarA(A);
  EdthBarA.history << "this->NPEdthBar();" << endl;

  for(int i_m=0; i_m<NModes; ++i_m) {
    const int l = A.LM(i_m)[0];
    if(l<lMin) {
      for(int i_t=0; i_t<NTimes; ++i_t) {
        EdthBarA.SetData(i_m, i_t, 0.0);
      }
    } else {
      const double factor = -std::sqrt((l+s)*(l-s+1));
      for(int i_t=0; i_t<NTimes; ++i_t) {
        EdthBarA.SetData(i_m, i_t, EdthBarA.Data(i_m, i_t)*factor);
      }
    }
  }

  EdthBarA.SetSpinWeight(A.SpinWeight()-1);
  EdthBarA.SetBoostWeight(WeightError);

  return EdthBarA;
}

/// Geroch-Held-Penrose edth operator
GWFrames::Waveform GWFrames::Waveform::GHPEdth() const {
  /// This operator is the one defined by Geroch et al. (1973).  It
  /// raises the spin weight of any field on the sphere by 1, while
  /// leaving the boost weight unchanged.
  ///
  /// This operator is very similar to the basic Newman-Penrose edth
  /// operator, except that it preserves boost weights.  Its effect in
  /// this implementation is identical (up to a factor of
  /// \f$\sqrt{2}\f$) to the NP edth.  There is an additional term in
  /// the definition of the GHP operator, but its value is zero.  (It
  /// transforms nontrivially, though.)  In this context, we have
  /// `NPEdth() = sqrt(2)*GHPEdth()`.
  ///
  /// The complex shear \f$\sigma\f$ has spin weight +2 and boost
  /// weight +1.  The radial coordinate \f$r\f$ has boost weight -1,
  /// and the derivative with respect to time \f$d/du\f$ has boost
  /// weight -1.  The asymptotic metric shear \f$r\, h\f$ has spin
  /// weight -2 and boost weight -1.  In particular, it seems that
  /// \f$r\, h = r^2\, \bar{\sigma}\f$.
  ///
  /// The Newman-Penrose scalars \f$\Psi_i\f$ have spin weight and
  /// boost weight equal to \f$2-i\f$.  (E.g., \f$\Psi_4\f$ has \f$s =
  /// b = -2\f$.)  However, when these are multiplied by the
  /// appropriate factors of \f$r\f$ to find the leading-order terms,
  /// they acquire boost weights.  In particular, we need to multiply
  /// \f$\Psi_i\f$ by \f$r^{5-i}\f$ to get nonzero values at scri,
  /// which adds \f$i-5\f$ to the boost weight, so that the asymptotic
  /// NP scalars all have boost weight -3.
  ///
  /// \sa NPEdth
  /// \sa NPEdthBar
  /// \sa GHPEdthBar
  /// \sa IntegrateNPEdth
  /// \sa IntegrateNPEdthBar
  /// \sa IntegrateGHPEdth
  /// \sa IntegrateGHPEdthBar

  const Waveform& A = *this;
  const int NModes = A.NModes();
  const int NTimes = A.NTimes();
  const int s = A.SpinWeight();
  const int lMin = std::abs(s+1);
  Waveform EdthA(A);
  EdthA.history << "this->GHPEdth();" << endl;

  for(int i_m=0; i_m<NModes; ++i_m) {
    const int l = A.LM(i_m)[0];
    if(l<lMin) {
      for(int i_t=0; i_t<NTimes; ++i_t) {
        EdthA.SetData(i_m, i_t, 0.0);
      }
    } else {
      const double factor = std::sqrt((l-s)*(l+s+1.)/2.);
      for(int i_t=0; i_t<NTimes; ++i_t) {
        EdthA.SetData(i_m, i_t, EdthA.Data(i_m, i_t)*factor);
      }
    }
  }

  EdthA.SetSpinWeight(A.SpinWeight()+1);
  //EdthA.SetBoostWeight(A.BoostWeight()); // No change

  return EdthA;
}

/// Geroch-Held-Penrose edth operator conjugate
GWFrames::Waveform GWFrames::Waveform::GHPEdthBar() const {
  /// This operator is the one defined by Geroch et al. (1973).  It
  /// lowers the spin weight of any field on the sphere by 1, while
  /// leaving the boost weight unchanged.
  ///
  /// This operator is very similar to the basic Newman-Penrose edth
  /// operator, except that it preserves boost weights.  Its effect in
  /// this implementation is identical (up to a factor of
  /// \f$\sqrt{2}\f$) to the NP edth.  There is an additional term in
  /// the definition of the GHP operator, but its value is zero.  (It
  /// transforms nontrivially, though.)  In this context, we have
  /// `NPEdthBar() = sqrt(2)*GHPEdthBar()`.
  ///
  /// The complex shear \f$\sigma\f$ has spin weight +2 and boost
  /// weight +1.  The radial coordinate \f$r\f$ has boost weight -1,
  /// and the derivative with respect to time \f$d/du\f$ has boost
  /// weight -1.  The asymptotic metric shear \f$r\, h\f$ has spin
  /// weight -2 and boost weight -1.  In particular, it seems that
  /// \f$r\, h = r^2\, \bar{\sigma}\f$.
  ///
  /// The Newman-Penrose scalars \f$\Psi_i\f$ have spin weight and
  /// boost weight equal to \f$2-i\f$.  (E.g., \f$\Psi_4\f$ has \f$s =
  /// b = -2\f$.)  However, when these are multiplied by the
  /// appropriate factors of \f$r\f$ to find the leading-order terms,
  /// they acquire boost weights.  In particular, we need to multiply
  /// \f$\Psi_i\f$ by \f$r^{5-i}\f$ to get nonzero values at scri,
  /// which adds \f$i-5\f$ to the boost weight, so that the asymptotic
  /// NP scalars all have boost weight -3.
  ///
  /// \sa NPEdth
  /// \sa NPEdthBar
  /// \sa GHPEdth
  /// \sa IntegrateNPEdth
  /// \sa IntegrateNPEdthBar
  /// \sa IntegrateGHPEdth
  /// \sa IntegrateGHPEdthBar

  const Waveform& A = *this;
  const int NModes = A.NModes();
  const int NTimes = A.NTimes();
  const int s = A.SpinWeight();
  const int lMin = std::abs(s-1);
  Waveform EdthBarA(A);
  EdthBarA.history << "this->GHPEdthBar();" << endl;

  for(int i_m=0; i_m<NModes; ++i_m) {
    const int l = A.LM(i_m)[0];
    if(l<lMin) {
      for(int i_t=0; i_t<NTimes; ++i_t) {
        EdthBarA.SetData(i_m, i_t, 0.0);
      }
    } else {
      const double factor = -std::sqrt((l+s)*(l-s+1.)/2.);
      for(int i_t=0; i_t<NTimes; ++i_t) {
        EdthBarA.SetData(i_m, i_t, EdthBarA.Data(i_m, i_t)*factor);
      }
    }
  }

  EdthBarA.SetSpinWeight(A.SpinWeight()-1);
  //EdthBarA.SetBoostWeight(A.BoostWeight()); // No change

  return EdthBarA;
}




/// Integrate the Newman-Penrose edth operator
GWFrames::Waveform GWFrames::Waveform::IntegrateNPEdth() const {
  /// This operator inverts the action of the Newman-Penrose edth
  /// operator.  This is not a perfect inverse, because the l=s-1 term
  /// is set to zero.  To be precise, if `Waveform A` has spin weight
  /// \f$s\f$, then `A.NPEdth().IntegrateNPEdth()` has the effect of
  /// setting the \f$\ell=s\f$ term in `A` to zero.
  ///
  /// Note that the N-P edth operator does not preserve boost weights,
  /// so the boost weight is set to the value of `WeightError`, which
  /// is just meant to be large enough that it will give improbable
  /// values if used.  This is not fool-proof.  See the GHP edth
  /// operator for a weight-preserving version.
  ///
  /// \sa NPEdth
  /// \sa NPEdthBar
  /// \sa GHPEdth
  /// \sa GHPEdthBar
  /// \sa IntegrateNPEdthBar
  /// \sa IntegrateGHPEdth
  /// \sa IntegrateGHPEdthBar

  const Waveform& A = *this;
  const int NModes = A.NModes();
  const int NTimes = A.NTimes();
  const int s = A.SpinWeight()-1;
  const int lMin = std::abs(s);
  Waveform IntegralEdthA(A);
  IntegralEdthA.history << "this->IntegrateNPEdth();" << endl;

  for(int i_m=0; i_m<NModes; ++i_m) {
    const int l = A.LM(i_m)[0];
    const double factor = std::sqrt((l-s)*(l+s+1));
    if(l<=lMin || factor == 0.0) {
      for(int i_t=0; i_t<NTimes; ++i_t) {
        IntegralEdthA.SetData(i_m, i_t, 0.0);
      }
    } else {
      for(int i_t=0; i_t<NTimes; ++i_t) {
        IntegralEdthA.SetData(i_m, i_t, IntegralEdthA.Data(i_m, i_t)/factor);
      }
    }
  }

  IntegralEdthA.SetSpinWeight(A.SpinWeight()-1);
  IntegralEdthA.SetBoostWeight(WeightError);

  return IntegralEdthA;
}

/// Integrate the Newman-Penrose edth operator conjugate
GWFrames::Waveform GWFrames::Waveform::IntegrateNPEdthBar() const {
  /// This operator inverts the action of the conjugated
  /// Newman-Penrose edth operator.  This is not a perfect inverse,
  /// because the l=s-1 term is set to zero.  To be precise, if
  /// `Waveform A` has spin weight \f$s\f$, then
  /// `A.NPEdthBar().IntegrateNPEdthBar()` has the effect of setting
  /// the \f$\ell=s\f$ term in `A` to zero.
  ///
  /// Note that the N-P edth operator does not preserve boost weights,
  /// so the boost weight is set to the value of `WeightError`, which
  /// is just meant to be large enough that it will give improbable
  /// values if used.  This is not fool-proof.  See the GHP edth
  /// operator for a weight-preserving version.
  ///
  /// \sa NPEdth
  /// \sa NPEdthBar
  /// \sa GHPEdth
  /// \sa GHPEdthBar
  /// \sa IntegrateNPEdthBar
  /// \sa IntegrateGHPEdth
  /// \sa IntegrateGHPEdthBar

  const Waveform& A = *this;
  const int NModes = A.NModes();
  const int NTimes = A.NTimes();
  const int s = A.SpinWeight()+1;
  const int lMin = std::abs(s);
  Waveform IntegralEdthBarA(A);
  IntegralEdthBarA.history << "this->IntegrateNPEdthBar();" << endl;

  for(int i_m=0; i_m<NModes; ++i_m) {
    const int l = A.LM(i_m)[0];
    const double factor = -std::sqrt((l+s)*(l-s+1));
    if(l<=lMin || factor == 0.0) {
      for(int i_t=0; i_t<NTimes; ++i_t) {
        IntegralEdthBarA.SetData(i_m, i_t, 0.0);
      }
    } else {
      for(int i_t=0; i_t<NTimes; ++i_t) {
        IntegralEdthBarA.SetData(i_m, i_t, IntegralEdthBarA.Data(i_m, i_t)/factor);
      }
    }
  }

  IntegralEdthBarA.SetSpinWeight(A.SpinWeight()+1);
  IntegralEdthBarA.SetBoostWeight(WeightError);

  return IntegralEdthBarA;
}

/// Integrate the Geroch-Held-Penrose edth operator
GWFrames::Waveform GWFrames::Waveform::IntegrateGHPEdth() const {
  /// This operator inverts the action of the GHP edth operator.  This
  /// is not a perfect inverse, because the l=s-1 term is set to zero.
  /// To be precise, if `Waveform A` has spins weight \f$s\f$, then
  /// `A.GHPEdth().IntegrateGHPEdth()` has the effect of setting the
  /// \f$\ell=s\f$ term in `A` to zero.
  ///
  /// \sa NPEdth
  /// \sa NPEdthBar
  /// \sa GHPEdth
  /// \sa GHPEdthBar
  /// \sa IntegrateNPEdth
  /// \sa IntegrateNPEdthBar
  /// \sa IntegrateGHPEdthBar

  const Waveform& A = *this;
  const int NModes = A.NModes();
  const int NTimes = A.NTimes();
  const int s = A.SpinWeight()-1;
  const int lMin = std::abs(s);
  Waveform IntegralEdthA(A);
  IntegralEdthA.history << "this->IntegrateGHPEdth();" << endl;

  for(int i_m=0; i_m<NModes; ++i_m) {
    const int l = A.LM(i_m)[0];
    const double factor = std::sqrt((l-s)*(l+s+1.)/2.);
    if(l<=lMin || factor == 0.0) {
      for(int i_t=0; i_t<NTimes; ++i_t) {
        IntegralEdthA.SetData(i_m, i_t, 0.0);
      }
    } else {
      for(int i_t=0; i_t<NTimes; ++i_t) {
        IntegralEdthA.SetData(i_m, i_t, IntegralEdthA.Data(i_m, i_t)/factor);
      }
    }
  }

  IntegralEdthA.SetSpinWeight(A.SpinWeight()-1);
  //IntegralEdthA.SetBoostWeight(A.BoostWeight()); // No change

  return IntegralEdthA;
}

/// Integrate the Geroch-Held-Penrose edth operator conjugate
GWFrames::Waveform GWFrames::Waveform::IntegrateGHPEdthBar() const {
  /// This operator inverts the action of the GHP edth operator.  This
  /// is not a perfect inverse, because the l=s-1 term is set to zero.
  /// To be precise, if `Waveform A` has spins weight \f$s\f$, then
  /// `A.GHPEdth().IntegrateGHPEdth()` has the effect of setting the
  /// \f$\ell=s\f$ term in `A` to zero.
  ///
  /// \sa NPEdth
  /// \sa NPEdthBar
  /// \sa GHPEdth
  /// \sa GHPEdthBar
  /// \sa IntegrateNPEdth
  /// \sa IntegrateNPEdthBar
  /// \sa IntegrateGHPEdth

  const Waveform& A = *this;
  const int NModes = A.NModes();
  const int NTimes = A.NTimes();
  const int s = A.SpinWeight()+1;
  const int lMin = std::abs(s);
  Waveform IntegralEdthBarA(A);
  IntegralEdthBarA.history << "this->GHPEdthBar();" << endl;

  for(int i_m=0; i_m<NModes; ++i_m) {
    const int l = A.LM(i_m)[0];
    const double factor = -std::sqrt((l+s)*(l-s+1.)/2.);
    if(l<lMin || factor == 0.0) {
      for(int i_t=0; i_t<NTimes; ++i_t) {
        IntegralEdthBarA.SetData(i_m, i_t, 0.0);
      }
    } else {
      for(int i_t=0; i_t<NTimes; ++i_t) {
        IntegralEdthBarA.SetData(i_m, i_t, IntegralEdthBarA.Data(i_m, i_t)/factor);
      }
    }
  }

  IntegralEdthBarA.SetSpinWeight(A.SpinWeight()+1);
  //IntegralEdthBarA.SetBoostWeight(A.BoostWeight()); // No change

  return IntegralEdthBarA;
}


// #ifndef DOXYGEN
// #ifdef __restrict
// #define restrict __restrict
// #endif
// using namespace std;
// extern "C" {
//   #include <stdlib.h>
//   #include <stdio.h>
//   #include <math.h>
//   #include <complex.h>
//   #include "fftw3.h"
//   #include "alm.h"
//   #include "wigner_d_halfpi.h"
//   #include "spinsfast_forward.h"
//   #include "spinsfast_backward.h"
// }
// #endif // DOXYGEN

/// Re-interpolate data to new time slices given by this supertranslation
GWFrames::Waveform GWFrames::Waveform::ApplySupertranslation(std::vector<std::complex<double> >& gamma) const {
  /// This function takes the current data decomposed as spherical
  /// harmonics on a given slicing, transforms to physical space,
  /// re-interpolates the data at each point to a new set of time
  /// slices, and transforms back to spherical-harmonic coefficients.
  ///
  /// The supertranslation data input `gamma` is a vector of
  /// complex numbers representing the (scalar) spherical-harmonic
  /// components of the supertranslation, stored in the order (0,0),
  /// (1,-1), (1,0), (1,1), (2,-2), ...  The overall time translation
  /// is given by the first component; the spatial translation is
  /// given by the second through fourth componentes; higher
  /// components give the proper supertranslations.  In particular, a
  /// proper supertranslation will have its first four coefficients
  /// equal to 0.0.
  ///
  /// Note that, for general spin-weighted spherical-harmonic
  /// components \f${}_{s}a_{l,m}\f$, a real function results when
  /// \f${}_{-s}a_{l,-m} = {}_{s}a_{l,m}^\ast\f$.  In particular, the
  /// input `gamma` data are assumed to satisfy this formula with
  /// \f$s=0\f$.

  // Find the max ell present in the input gamma data
  int lMax=0;
  for(; lMax<=SphericalFunctions::ellMax; ++lMax) {
    if(N_lm(lMax)==int(gamma.size())) {
      break;
    }
  }

  // Make sure gamma doesn't have more ell modes than SphericalFunctions::ellMax
  if(lMax>=SphericalFunctions::ellMax) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nError: Input supertranslation data has length " << gamma.size() << "."
              << "\n       This is not a recognized length for spherical-harmonic data.\n"
              << std::endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  const unsigned int Nlm = N_lm(lMax);

  // Make sure the Waveform is in the inertial frame
  if(frameType != GWFrames::Inertial) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nError: Asking to supertranslate a Waveform in the " << GWFrames::WaveformFrameNames[frameType] << " frame."
              << "\n       This only makes sense with Waveforms in the " << GWFrames::WaveformFrameNames[GWFrames::Inertial] << " frame.\n"
              << std::endl;
    throw(GWFrames_WrongFrameType);
  }

  // Set up some constants
  const Waveform& A = *this;
  const unsigned int ntimes = A.NTimes();
  const unsigned int nmodes = A.NModes();
  const complex<double> zero(0.0,0.0);
  const complex<double> I(0.0,1.0);

  // Copy infrastructure to new Waveform
  Waveform B;
  B.spinweight = A.spinweight;
  B.boostweight = A.boostweight;
  B.history.str(A.history.str());
  B.history.clear();
  B.history.seekp(0, ios_base::end);
  B.history << "*this = B.ApplySupertranslation([" << gamma[0];
  for(unsigned int i=1; i<gamma.size(); ++i) {
    B.history << ", " << gamma[i];
  }
  B. history << "]);"<< std::endl;
  B.frame = A.frame;
  B.frameType = A.frameType;
  B.dataType = A.dataType;
  B.rIsScaledOut = A.rIsScaledOut;
  B.mIsScaledOut = A.mIsScaledOut;
  B.lm = A.lm;

  // These numbers determine the equi-angular grid on which we will do
  // the interpolation.  For best accuracy, have N_phi > 2*lMax and
  // N_theta > 2*lMax; but for speed, don't make them much greater.
  const unsigned int N_phi = 2*lMax + 1;
  const unsigned int N_theta = 2*lMax + 1;
  const unsigned int N_tot = N_phi*N_theta;

  // Transform times to physical space
  vector<complex<double> > DeltaT_complex(N_tot, zero);
  spinsfast_salm2map(reinterpret_cast<fftw_complex*>(&gamma[0]),
                     reinterpret_cast<fftw_complex*>(&DeltaT_complex[0]),
                     0, N_theta, N_phi, lMax);
  vector<double> DeltaT(DeltaT_complex.size());
  for(unsigned int i=0; i<DeltaT.size(); ++i) {
    DeltaT[i] = real(DeltaT_complex[i]);
  }

  // Find largest and smallest time excursions
  double MinDeltaT = 0.0;
  double MaxDeltaT = 0.0;
  for(unsigned int i_t=0; i_t<N_tot; ++i_t) {
    // cerr << "DeltaT[" << i_t << "] = " << DeltaT[i_t] << endl;
    if(DeltaT[i_t] < MinDeltaT) { MinDeltaT = DeltaT[i_t]; }
    if(DeltaT[i_t] > MaxDeltaT) { MaxDeltaT = DeltaT[i_t]; }
  }
  // cerr << "Min/MaxDeltaT: " << MinDeltaT << " " << MaxDeltaT << endl;

  // Set up new time slices, beginning with an offset of MinDeltaT
  // and ending with an offset of -MaxDeltaT (so that the
  // interpolation does not need to extrapolate)
  const double FirstT = t[0]-MinDeltaT;
  const double LastT = t.back()-MaxDeltaT;
  unsigned int i_Min, i_Max;
  for(i_Min=0; i_Min<ntimes; ++i_Min) {
    if(t[i_Min]>=FirstT) {
      break;
    }
  }
  for(i_Max=ntimes-1; i_Max>0; --i_Max) {
    if(t[i_Max]<=LastT) {
      break;
    }
  }
  B.t = std::vector<double>(t.begin()+i_Min, t.begin()+i_Max+1);
  const unsigned int ntimesB = B.NTimes();
  B.ResizeData(nmodes, ntimesB);
  // cerr << "First/LastT: " << FirstT << " " << LastT << endl;
  // cerr << "i_Min/Max: " << i_Min << " " << i_Max << endl;
  // cerr << "Times: " << B.t[0] << " " << B.t.back() << endl;

  // Storage arrays
  std::vector<std::complex<double> > alm(Nlm); // Work array
  std::vector<std::complex<double> > mapA(N_tot); // Work array
  std::vector<std::vector<double> > fARe(N_tot, std::vector<double>(ntimes));
  std::vector<std::vector<double> > fAIm(N_tot, std::vector<double>(ntimes));
  std::vector<std::vector<std::complex<double> > > fB(ntimesB, std::vector<std::complex<double> >(N_tot));

  // Transform to physical space (before interpolation)
  for(unsigned int i_t=0; i_t<ntimes; ++i_t) {
    // Get alm data at this time
    for(unsigned int i_m=0; i_m<nmodes; ++i_m) {
      alm[i_m] = A.Data(i_m,i_t);
    }
    // Transform
    spinsfast_salm2map(reinterpret_cast<fftw_complex*>(&alm[0]),
                       reinterpret_cast<fftw_complex*>(&mapA[0]),
                       0, N_theta, N_phi, lMax);
    // Save to main storage arrays
    for(unsigned int i_p=0; i_p<N_tot; ++i_p) {
      fARe[i_p][i_t] = real(mapA[i_p]);
      fAIm[i_p][i_t] = imag(mapA[i_p]);
    }
  }

  // Declare the GSL interpolators for the data
  gsl_interp_accel* accRe = gsl_interp_accel_alloc();
  gsl_interp_accel* accIm = gsl_interp_accel_alloc();
  gsl_spline* splineRe = gsl_spline_alloc(gsl_interp_cspline, ntimes);
  gsl_spline* splineIm = gsl_spline_alloc(gsl_interp_cspline, ntimes);

  // Loop over all points doing interpolation
  // cerr << "Begin\n=====" << endl;
  for(unsigned int i_p=0; i_p<N_tot; ++i_p) {
    const double t_p = DeltaT[i_p];
    // cerr << i_p << ", " << t_p << endl;
    // Initialize interpolators for this point
    gsl_spline_init(splineRe, &(A.t)[0], &fARe[i_p][0], ntimes);
    gsl_spline_init(splineIm, &(A.t)[0], &fAIm[i_p][0], ntimes);
    // Loop over all times at a given point
    for(unsigned int i_t=0; i_t<ntimesB; ++i_t) {
      // cerr << "\t" << B.t[i_t] << flush;
      // Interpolate data onto new time slices
      fB[i_t][i_p] = complex<double>( gsl_spline_eval(splineRe, B.t[i_t]+t_p, accRe), gsl_spline_eval(splineIm, B.t[i_t]+t_p, accIm) );
      // cerr << "+" << endl;
    }
    // cerr << "+" << endl;
  }
  // cerr << "End\n===" << endl;

  // Free the data and interpolators
  fARe.clear();
  fAIm.clear();
  gsl_interp_accel_free(accRe);
  gsl_interp_accel_free(accIm);
  gsl_spline_free(splineRe);
  gsl_spline_free(splineIm);

  // Transform back to spectral space (after interpolation)
  for(unsigned int i_t=0; i_t<ntimesB; ++i_t) {
    // Transform
    spinsfast_map2salm(reinterpret_cast<fftw_complex*>(&fB[i_t][0]),
                       reinterpret_cast<fftw_complex*>(&alm[0]),
                       0, N_theta, N_phi, lMax);
    // Save to main storage arrays
    for(unsigned int i_m=0; i_m<Nlm; ++i_m) {
      B.data[i_m][i_t] = alm[i_m];
    }
  }

  return B;
}


// /// Apply a boost to a boost-weighted function
// GWFrames::Waveform GWFrames::Waveform::Boost(const std::vector<double>& v) const {
//   /// This function does three things.  First, it evaluates the
//   /// Waveform on what will become an equi-angular grid after
//   /// transformation by the boost.  Second, it multiplies each of
//   /// those points by the appropriate conformal factor
//   /// \f$K^b(\vartheta, \varphi)\f$, where \f$b\f$ is the boost weight
//   /// stored with the Waveform.  Finally, it transforms back to
//   /// Fourier space using that new equi-angular grid.


//   std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Not properly implemented yet." << std::endl;
//   throw(GWFrames_NotYetImplemented);

//   //return Waveform();
//   // return B;
// }


/// Apply a boost to Psi4 data
GWFrames::Waveform& GWFrames::Waveform::BoostPsi4(const std::vector<std::vector<double> >& v) {
  /// This function does three things.  First, it evaluates the
  /// Waveform on what will become an equi-angular grid after
  /// transformation by the boost.  Second, at each point of that
  /// grid, it takes the appropriate combinations of the present value
  /// of Psi_4 and its conjugate to give the value of Psi_4 as
  /// observed in the boosted frame.  Finally, it transforms back to
  /// Fourier space using that new equi-angular grid.
  ///
  /// The input three-velocities are assumed to give the velocities of
  /// the boosted frame relative to the present frame.

  // Check the size of the input velocity
  if(v.size()!=NTimes()) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": (v.size()=" << v.size() << ") != (NTimes()=" << NTimes() << ")" << std::endl;
    throw(GWFrames_VectorSizeMismatch);
  }

  // Set up storage and calculate useful constants
  const int ellMax = this->EllMax();
  const int n_thetaRotated= 2*ellMax+1;
  const int n_phiRotated = 2*ellMax+1;
  const double dthetaRotated = M_PI/double(n_thetaRotated-1); // thetaRotated should return to M_PI
  const double dphiRotated = 2*M_PI/double(n_phiRotated); // phiRotated should not return to 2*M_PI
  vector<complex<double> > Grid(n_phiRotated*n_thetaRotated);
  SphericalFunctions::SWSH sYlm(SpinWeight());

  SpacetimeAlgebra::vector tPz;
  tPz.set_gamma_0(1./std::sqrt(2));
  tPz.set_gamma_3(1./std::sqrt(2));
  SpacetimeAlgebra::vector tMz;
  tMz.set_gamma_0(1./std::sqrt(2));
  tMz.set_gamma_3(-1./std::sqrt(2));
  SpacetimeAlgebra::vector xPiyRe;
  xPiyRe.set_gamma_1(1./std::sqrt(2));
  SpacetimeAlgebra::vector xPiyIm;
  xPiyIm.set_gamma_2(1./std::sqrt(2));
  SpacetimeAlgebra::vector xMiyRe;
  xMiyRe.set_gamma_1(1./std::sqrt(2));
  SpacetimeAlgebra::vector xMiyIm;
  xMiyIm.set_gamma_2(-1./std::sqrt(2));

  // Main loop over time steps
  for(unsigned int i_t=0; i_t<NTimes(); ++i_t) {
    const vector<double>& v_i = v[i_t];
    if(v_i.size()!=3) {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": v[" << i_t << "].size()=" << v_i.size()
                << ".  Input is assumed to be a vector of three-velocities." << std::endl;
      throw(GWFrames_VectorSizeMismatch);
    }

    const double beta = std::sqrt(v_i[0]*v_i[0] + v_i[1]*v_i[1] + v_i[2]*v_i[2]);
    if(beta<1.e-9) { continue; } // TODO: This may need to be adjusted, or other statements made smarter about using the value of gamma
    vector<double> vHat(3);
    vHat[0] = v_i[0]/beta;
    vHat[1] = v_i[1]/beta;
    vHat[2] = v_i[2]/beta;
    const double gamma = 1.0/std::sqrt(1.0-beta*beta);
    const double sqrtplus = std::sqrt((gamma+1)/2);
    const double sqrtminus = std::sqrt((gamma-1)/2);

    // Calculate the boost rotor
    SpacetimeAlgebra::spinor BoostRotor;
    BoostRotor.set_scalar(sqrtplus);
    BoostRotor.set_gamma_0_gamma_1(sqrtminus*vHat[0]);
    BoostRotor.set_gamma_0_gamma_2(sqrtminus*vHat[1]);
    BoostRotor.set_gamma_0_gamma_3(sqrtminus*vHat[2]);

    vector<complex<double> > Modes(N_lm(ellMax), 0.0);
    vector<complex<double> > Modes2(N_lm(ellMax), 0.0);

    // Fill the Modes data for this time step
    for(int i_mode=0; i_mode<N_lm(std::abs(SpinWeight())-1); ++i_mode) {
      // Explicitly zero the modes with ell<|s|
      Modes[i_mode] = 0.0;
    }
    for(int i_mode=N_lm(std::abs(SpinWeight())-1), ell=std::abs(SpinWeight()); ell<=ellMax; ++ell) {
      // Now, fill modes with ell>=|s| in the correct order
      for(int m=-ell; m<=ell; ++m, ++i_mode) {
        Modes[i_mode] = this->Data(FindModeIndex(ell,m), i_t);
      }
    }

    // Construct the data on the distorted grid
    for(int i_g=0, i_thetaRotated=0; i_thetaRotated<n_thetaRotated; ++i_thetaRotated) {
      for(int i_phiRotated=0; i_phiRotated<n_phiRotated; ++i_phiRotated, ++i_g) {
        const double thetaRotated = dthetaRotated*i_thetaRotated;
        const double phiRotated = dphiRotated*i_phiRotated;

        // Calculate the rotation rotors
        SpacetimeAlgebra::spinor Rotor_thetaRotated;
        Rotor_thetaRotated.set_scalar(std::cos(thetaRotated/2));
        Rotor_thetaRotated.set_gamma_1_gamma_3(std::sin(thetaRotated/2));
        SpacetimeAlgebra::spinor Rotor_phiRotated;
        Rotor_phiRotated.set_scalar(std::cos(phiRotated/2));
        Rotor_phiRotated.set_gamma_1_gamma_2(-std::sin(phiRotated/2));
        const SpacetimeAlgebra::spinor RotationRotorRotated(Rotor_phiRotated * Rotor_thetaRotated);

        // This is the complete transformation rotor for going from
        // (t,x,y,z) in the present frame to (t,theta,phi,r) in the
        // boosted frame:
        const SpacetimeAlgebra::spinor LorentzRotor(BoostRotor * RotationRotorRotated);

        // The following give the important tetrad elements in the boosted frame
        const int Filler=0; // Useless constant for Gaigen code
        const SpacetimeAlgebra::vector lRotated(LorentzRotor*tPz*SpacetimeAlgebra::reverse(LorentzRotor), Filler);
        const SpacetimeAlgebra::vector nRotated(LorentzRotor*tMz*SpacetimeAlgebra::reverse(LorentzRotor), Filler);
        const SpacetimeAlgebra::vector mBarReRotated(LorentzRotor*xMiyRe*SpacetimeAlgebra::reverse(LorentzRotor), Filler);
        const SpacetimeAlgebra::vector mBarImRotated(LorentzRotor*xMiyIm*SpacetimeAlgebra::reverse(LorentzRotor), Filler);

        // Figure out the coordinates in the present frame
        // corresponding to the given coordinates in the boosted frame
        vector<double> r(3);
        r[0] = lRotated.get_gamma_1();
        r[1] = lRotated.get_gamma_2();
        r[2] = lRotated.get_gamma_3();
        const double rMag = std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
        const double theta = std::acos(r[2]/rMag);
        const double phi = std::atan2(r[1],r[0]);

        if(i_g==0 && i_t==0) {
          std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ":\n"
                    << "    Note that the (theta,phi) coordinates produced here are not in the same range\n"
                    << "    as the (thetaRotated,phiRotated) coordinates because of (1) the range of atan2,\n"
                    << "    which is in (-pi,pi), rather than (0,2*pi); and (2) at (theta=0), the phi value\n"
                    << "    comes out as 0, even though phiRotated may not be.\n\n"
                    << "    Fortunately, I think both these problems are handled automatically by taking the\n"
                    << "    tetrad components as we do.  Of course, I may be missing something problematic...\n" << std::endl;
        }

        // This gives us the rotor to get from the z axis to the
        // spherical coordinates in the present frame
        SpacetimeAlgebra::spinor Rotor_theta;
        Rotor_theta.set_scalar(std::cos(theta/2));
        Rotor_theta.set_gamma_1_gamma_3(std::sin(theta/2));
        SpacetimeAlgebra::spinor Rotor_phi;
        Rotor_phi.set_scalar(std::cos(phi/2));
        Rotor_phi.set_gamma_1_gamma_2(-std::sin(phi/2));
        const SpacetimeAlgebra::spinor RotationRotor(Rotor_phi * Rotor_theta);

        // The following give the important tetrad elements in the present frame
        const SpacetimeAlgebra::vector l(RotationRotor*tPz*SpacetimeAlgebra::reverse(RotationRotor), Filler);
        const SpacetimeAlgebra::vector n(RotationRotor*tMz*SpacetimeAlgebra::reverse(RotationRotor), Filler);
        const SpacetimeAlgebra::vector mRe(RotationRotor*xPiyRe*SpacetimeAlgebra::reverse(RotationRotor), Filler);
        const SpacetimeAlgebra::vector mIm(RotationRotor*xPiyIm*SpacetimeAlgebra::reverse(RotationRotor), Filler);
        const SpacetimeAlgebra::vector mBarRe(RotationRotor*xMiyRe*SpacetimeAlgebra::reverse(RotationRotor), Filler);
        const SpacetimeAlgebra::vector mBarIm(RotationRotor*xMiyIm*SpacetimeAlgebra::reverse(RotationRotor), Filler);

        // Get the value of Psi4 in this frame at the appropriate
        // point of this frame
        const Quaternion Rp(theta, phi);
        sYlm.SetRotation(Rp);
        const complex<double> Psi_4 = sYlm.Evaluate(Modes);

        // Get the components of the other frame's tetrad in the basis
        // of this tetrad.  In particular, these are *not* the dot
        // products of the other frame's basis vectors with this
        // frame's basis vectors.  Instead, we expand, e.g., nRotated
        // in terms of this frame's (l,n,m,mbar) basis, and just take
        // the coefficients in that expansion.  [This distinction
        // matters because, e.g., n.n = 0 but n.l \neq 0.]  Also note
        // that we will not need any components involving the l vector
        // in either frame, because that will just give us terms
        // proportional to Psi3, etc., which are assumed to fall off
        // more quickly than we care to bother with.
        const complex<double> i_complex(0.,1.);
        const complex<double> nRotated_n = -SpacetimeAlgebra::sp(nRotated, l);
        const complex<double> nRotated_m = SpacetimeAlgebra::sp(nRotated, mBarRe) + i_complex*SpacetimeAlgebra::sp(nRotated, mBarIm);
        const complex<double> nRotated_mBar = SpacetimeAlgebra::sp(nRotated, mRe) + i_complex*SpacetimeAlgebra::sp(nRotated, mIm);
        const complex<double> mBarRotated_n =
          - ( SpacetimeAlgebra::sp(mBarReRotated, l) + i_complex*SpacetimeAlgebra::sp(mBarImRotated, l) );
        const complex<double> mBarRotated_m =
          SpacetimeAlgebra::sp(mBarReRotated, mBarRe) + i_complex*SpacetimeAlgebra::sp(mBarReRotated, mBarIm)
          + i_complex * ( SpacetimeAlgebra::sp(mBarImRotated, mBarRe) + i_complex*SpacetimeAlgebra::sp(mBarImRotated, mBarIm) );
        const complex<double> mBarRotated_mBar =
          SpacetimeAlgebra::sp(mBarReRotated, mRe) + i_complex*SpacetimeAlgebra::sp(mBarReRotated, mIm)
          + i_complex * ( SpacetimeAlgebra::sp(mBarImRotated, mRe) + i_complex*SpacetimeAlgebra::sp(mBarImRotated, mIm) );

        // Evaluate the data for the boosted frame at this point
        Grid[i_g] =
          (nRotated_n * mBarRotated_mBar * nRotated_n * mBarRotated_mBar
           - nRotated_mBar * mBarRotated_n * nRotated_n * mBarRotated_mBar
           - nRotated_n * mBarRotated_mBar * nRotated_mBar * mBarRotated_n
           + nRotated_mBar * mBarRotated_n * nRotated_mBar * mBarRotated_n) * Psi_4
          + (nRotated_n * mBarRotated_m * nRotated_n * mBarRotated_m
             - nRotated_m * mBarRotated_n * nRotated_n * mBarRotated_m
             - nRotated_n * mBarRotated_m * nRotated_m * mBarRotated_n
             + nRotated_m * mBarRotated_n * nRotated_m * mBarRotated_n) * std::conj(Psi_4);

        if(i_t%5000==0) {
          std::cerr << thetaRotated << "," << phiRotated << "; " << theta << "," << phi
                    << "; \t" << Grid[i_g] << "," << Psi_4 << std::endl;
        }
      }
    }

    // Decompose the data into modes
    spinsfast_map2salm(reinterpret_cast<fftw_complex*>(&Grid[0]),
                       reinterpret_cast<fftw_complex*>(&Modes2[0]),
                       SpinWeight(), n_thetaRotated, n_phiRotated, ellMax);

    // Set new data at this time step
    for(int i_mode=N_lm(std::abs(SpinWeight())-1), ell=std::abs(SpinWeight()); ell<=ellMax; ++ell) {
      for(int m=-ell; m<=ell; ++m, ++i_mode) {
        this->SetData(this->FindModeIndex(ell,m), i_t, Modes2[i_mode]);
      }
    }

  }

  return *this;
}

/// Translate the waveform data by some series of spatial translations
GWFrames::Waveform GWFrames::Waveform::Translate(const std::vector<std::vector<double> >& deltax) const {
  /// \param x Array of 3-vectors by which to translate
  if(deltax.size() != NTimes()) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nERROR: (deltax.size()=" << deltax.size() << ") != (NTimes()=" << NTimes() << ")"
              << "\n       If you're going to translate, you need to do it at each time step.\n" << std::endl;
    throw(GWFrames_VectorSizeMismatch);
  }

  const unsigned int ntimes = NTimes();

  for(unsigned int i=0; i<ntimes; ++i) {
    if(deltax[i].size() != 3) {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
                << "\nERROR: (deltax[" << i << "].size()=" << deltax[i].size() << ") != 3"
                << "\n       Each translation should be a 3-vector.\n" << std::endl;
      throw(GWFrames_VectorSizeMismatch);
    }
  }

  // Copy infrastructure to new Waveform
  const Waveform& A = *this;
  Waveform B;
  B.spinweight = A.spinweight;
  B.boostweight = A.boostweight;
  B.history.str(A.history.str());
  B.history.clear();
  B.history.seekp(0, ios_base::end);
  B.history << "*this = B.Translate(...);"<< std::endl;
  B.frame = std::vector<Quaternions::Quaternion>(0);
  B.frameType = A.frameType;
  B.dataType = A.dataType;
  B.rIsScaledOut = A.rIsScaledOut;
  B.mIsScaledOut = A.mIsScaledOut;
  B.lm = A.lm;
  B.t = A.t; // This will get reset later

  // These numbers determine the equi-angular grid on which we will do
  // the interpolation.  For best accuracy, have N_phi > 2*ellMax and
  // N_theta > 2*ellMax; but for speed, don't make them much greater.
  const int ellMax(EllMax());
  const unsigned int N_phi = 2*ellMax + 1;
  const unsigned int N_theta = 2*ellMax + 1;
  const double dtheta = M_PI/double(N_theta-1); // theta should return to M_PI
  const double dphi = 2*M_PI/double(N_phi); // phi should not return to 2*M_PI
  vector<complex<double> > Grid(N_phi*N_theta);

  // Find earliest and latest times we can use for our new data set
  unsigned int iEarliest = 0;
  unsigned int iLatest = ntimes-1;
  const double tEarliest = t[0];
  const double tLatest = t.back();
  { // Do the 0th point explicitly for earliest time only
    const unsigned int i=0;
    const double deltaxiMag = std::sqrt(deltax[i][0]*deltax[i][0] + deltax[i][1]*deltax[i][1] + deltax[i][2]*deltax[i][2]);
    if(t[i]-deltaxiMag<tEarliest) {
      iEarliest = std::max(iEarliest, i+1);
    }
  }
  for(unsigned int i=1; i<ntimes-1; ++i) { // Do all points in between
    const double deltaxiMag = std::sqrt(deltax[i][0]*deltax[i][0] + deltax[i][1]*deltax[i][1] + deltax[i][2]*deltax[i][2]);
    if(t[i]-deltaxiMag<tEarliest) {
      iEarliest = std::max(iEarliest, i+1);
    }
    if(t[i]+deltaxiMag>tLatest) {
      iLatest = std::min(iLatest, i-1);
    }
  }
  { // Do the last point explicitly for latest time
    const unsigned int i=ntimes-1;
    const double deltaxiMag = std::sqrt(deltax[i][0]*deltax[i][0] + deltax[i][1]*deltax[i][1] + deltax[i][2]*deltax[i][2]);
    if(t[i]+deltaxiMag>tLatest) {
      iLatest = std::min(iLatest, i-1);
    }
  }

  // Now check that we actually have something to work with here
  if(iLatest<iEarliest) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
              << "\nERROR: (iLatest=" << iLatest << ") < (iEarliest=" << iEarliest << ")"
              << "\n       The translation is too extreme, and there is not enough data for even one complete time step.\n" << std::endl;
    throw(GWFrames_EmptyIntersection);
  }

  // The new times will just be that subset of the old times for which
  // data exist in every direction after translation.
  B.t.erase(B.t.begin()+iLatest+1, B.t.end());
  B.t.erase(B.t.begin(), B.t.begin()+iEarliest);
  B.data.resize(NModes(), B.NTimes()); // Each row (first index, nn) corresponds to a mode

  // Main loop over time steps
  for(unsigned int i_t_B=0; i_t_B<B.NTimes(); ++i_t_B) {
    // These will hold the input and output data.  Note that spinsfast
    // has had some weird behavior in the past when the output data is
    // not initialized to zero, so just do this each time, even though
    // it's slow.
    vector<complex<double> > Modes(N_lm(ellMax), 0.0);

    // Construct the data on the translated grid
    for(int i_g=0, i_theta=0; i_theta<N_theta; ++i_theta) {
      for(int i_phi=0; i_phi<N_phi; ++i_phi, ++i_g) {
        const double theta = dtheta*i_theta;
        const double phi = dphi*i_phi;

        const double rHat_dot_deltax =
          deltax[i_t_B][0]*std::sin(theta)*std::cos(phi)
          + deltax[i_t_B][1]*std::sin(theta)*std::sin(phi)
          + deltax[i_t_B][2]*std::cos(theta);

        // Evaluate the data for the translated frame at this point
        Grid[i_g] = A.InterpolateToPoint(theta, phi, B.T(i_t_B)-rHat_dot_deltax);
      }
    }

    // Decompose the data into modes
    spinsfast_map2salm(reinterpret_cast<fftw_complex*>(&Grid[0]),
                       reinterpret_cast<fftw_complex*>(&Modes[0]),
                       SpinWeight(), N_theta, N_phi, ellMax);

    // Set new data at this time step
    for(int i_mode=N_lm(std::abs(SpinWeight())-1), ell=std::abs(SpinWeight()); ell<=ellMax; ++ell) {
      for(int m=-ell; m<=ell; ++m, ++i_mode) {
        B.SetData(B.FindModeIndex(ell,m), i_t_B, Modes[i_mode]);
      }
    }

  } // i_t_B loop

  return B;
}

















//////////////////////////
/// Waveforms (plural!) //
//////////////////////////


// Constructors

/// Empty constructor of N empty objects.
GWFrames::Waveforms::Waveforms(const int N) : Ws(N), CommonTimeSet(false) { }

/// Basic copy constructor.
GWFrames::Waveforms::Waveforms(const Waveforms& In) : Ws(In.Ws), CommonTimeSet(In.CommonTimeSet) { }

/// Basic copy constructor.
GWFrames::Waveforms::Waveforms(const std::vector<Waveform>& In) : Ws(In), CommonTimeSet(false) { }

/// Interpolate to a common set of times.
void GWFrames::Waveforms::SetCommonTime(std::vector<std::vector<double> >& Radii,
                                                        const double MinTimeStep, const double EarliestTime, const double LatestTime) {
  const unsigned int NWaveforms = Radii.size();
  vector<double> TLimits(2);
  TLimits[0] = EarliestTime;
  TLimits[1] = LatestTime;
  vector<double> T = GWFrames::Intersection(TLimits, Ws[0].T(), MinTimeStep, EarliestTime, LatestTime);
  for(unsigned int i_W=1; i_W<Ws.size(); ++i_W) {
    T = GWFrames::Intersection(T, Ws[i_W].T());
  }
  // Interpolate radii
  gsl_interp_accel* acc = gsl_interp_accel_alloc();
  gsl_spline* spline = gsl_spline_alloc(gsl_interp_cspline, Ws[0].NTimes());
  for(unsigned int i_W=0; i_W<NWaveforms; ++i_W) {
    vector<double> OriginalTime = Ws[i_W].T();
    vector<double> NewRadii(T.size());
    gsl_spline_init(spline, &OriginalTime[0], &Radii[i_W][0], Radii[i_W].size());
    for(unsigned int i_t=0; i_t<NewRadii.size(); ++i_t) {
      NewRadii[i_t] = gsl_spline_eval(spline, T[i_t], acc);
    }
    Radii[i_W] = NewRadii;
  }
  gsl_interp_accel_free(acc);
  gsl_spline_free(spline);
  // Interpolate Waveforms
  for(unsigned int i_W=0; i_W<NWaveforms; ++i_W) {
    Ws[i_W] = Ws[i_W].Interpolate(T);
  }
  CommonTimeSet = true;
  return;
}

template <typename T>
string NumberToString ( T Number ) {
  std::ostringstream ss;
  ss << Number;
  return ss.str();
}

/// Main extrapolation routine.
GWFrames::Waveforms GWFrames::Waveforms::Extrapolate(std::vector<std::vector<double> >& Radii,
                                                     const std::vector<int>& ExtrapolationOrders,
                                                     const std::vector<double>& Omegas)
{
  ///
  /// \param Radii Array of radii for each Waveform (first index) and each time (second index)
  /// \param ExtrapolationOrders List of integers denote extrapolation orders
  /// \param Omegas Optional list of angular frequencies for scaling extrapolation polynomial
  ///
  /// The input FiniteRadiusWaveforms are assumed to be properly
  /// scaled and time-retarded, and interpolated to a uniform set of
  /// retarded times.  This function simply steps through the indices,
  /// fitting those data to polynomials in 1/radius, and evaluating at
  /// 0 (for infinity).
  ///
  /// The extrapolation orders can be negative.  In this case, the
  /// scaled, time-retarded waveform at finite radius is given, where
  /// N=-1 is the outermost Waveform, N=-2 is the second to outermost,
  /// etc.
  ///
  /// Note that the fitting uses gsl_multifit_linear_usvd, which is
  /// GSL's fitting function that does NOT use column scaling
  /// (specified by the 'u' in front of 'svd' in the function name).
  /// The basic GSL fitting function uses column scaling "to improve
  /// the accuracy of the singular values".  However, for convergent
  /// series, this scaling can make all the coefficients roughly equal
  /// (just as the `Omegas` option does), which defeats the SVD.
  ///

  Waveforms& FiniteRadiusWaveforms = *this;
  if(! FiniteRadiusWaveforms.CommonTimeSet) { FiniteRadiusWaveforms.SetCommonTime(Radii); }

  // Get the various dimensions, etc.
  const int MaxN = *std::max_element(ExtrapolationOrders.begin(), ExtrapolationOrders.end());
  const int MinN = *std::min_element(ExtrapolationOrders.begin(), ExtrapolationOrders.end());
  const bool UseOmegas = (Omegas.size()!=0);
  const int NTimes = FiniteRadiusWaveforms[0].NTimes();
  const unsigned int NModes = FiniteRadiusWaveforms[0].NModes();
  const unsigned int NFiniteRadii = FiniteRadiusWaveforms.size();
  const unsigned int NExtrapolations = ExtrapolationOrders.size();
  const double SVDTol = 1.e-12; // Same as Numerical Recipes default in fitsvd.h

  // Make sure everyone is playing with a full deck
  if((unsigned int)(std::abs(MinN))>NFiniteRadii) {
      cerr << "\n\n" << __FILE__ << ":" << __LINE__
           << "\nERROR: Asking for finite-radius waveform " << MinN << ", but only got " << NFiniteRadii << " finite-radius Waveform objects."
           << "\n       Need at least as " << std::abs(MinN) << " finite-radius waveforms."
           << "\n\n" << endl;
      throw(GWFrames_IndexOutOfBounds);
  }
  if(MaxN>0 && (unsigned int)(MaxN+1)>=NFiniteRadii) {
      cerr << "\n\n" << __FILE__ << ":" << __LINE__
           << "\nERROR: Asking for extrapolation up to order " << MaxN << ", but only got " << NFiniteRadii << " finite-radius Waveform objects."
           << "\n       Need at least " << MaxN+1 << " waveforms."
           << "\n\n" << endl;
      throw(GWFrames_IndexOutOfBounds);
  }
  if(Radii.size()!=NFiniteRadii) {
      cerr << "\n\n" << __FILE__ << ":" << __LINE__
           << "\nERROR: Mismatch in data to be extrapolated; there are different numbers of waveforms and radius vectors."
           << "\n       FiniteRadiusWaveforms.size()=" << NFiniteRadii
           << "\n       Radii.size()=" << Radii.size()
           << "\n\n" << endl;
      throw(GWFrames_VectorSizeMismatch);
  }
  if(UseOmegas && int(Omegas.size())!=NTimes) {
      cerr << "\n\n" << __FILE__ << ":" << __LINE__
           << "\nERROR: NTimes mismatch in data to be extrapolated."
           << "\n       FiniteRadiusWaveforms[0].NTimes()=" << NTimes
           << "\n       Omegas.size()=" << Omegas.size()
           << "\n\n" << endl;
      throw(GWFrames_VectorSizeMismatch);
  }
  for(unsigned int i_W=1; i_W<NFiniteRadii; ++i_W) {
    if(int(FiniteRadiusWaveforms[i_W].NTimes()) != NTimes) {
      cerr << "\n\n" << __FILE__ << ":" << __LINE__
           << "\nERROR: NTimes mismatch in data to be extrapolated."
           << "\n       FiniteRadiusWaveforms[0].NTimes()=" << NTimes
           << "\n       FiniteRadiusWaveforms[" << i_W << "].NTimes()=" << FiniteRadiusWaveforms[i_W].NTimes()
           << "\n\n" << endl;
      throw(GWFrames_VectorSizeMismatch);
    }
    if(FiniteRadiusWaveforms[i_W].NModes() != NModes) {
      cerr << "\n\n" << __FILE__ << ":" << __LINE__
           << "\nERROR: NModes mismatch in data to be extrapolated."
           << "\n       FiniteRadiusWaveforms[0].NModes()=" << NModes
           << "\n       FiniteRadiusWaveforms[" << i_W << "].NModes()=" << FiniteRadiusWaveforms[i_W].NModes()
           << "\n\n" << endl;
      throw(GWFrames_VectorSizeMismatch);
    }
    if(int(Radii[i_W].size()) != NTimes) {
      cerr << "\n\n" << __FILE__ << ":" << __LINE__
           << "\nERROR: NTimes mismatch in data to be extrapolated."
           << "\n       FiniteRadiusWaveforms[0].NTimes()=" << NTimes
           << "\n       Radii[" << i_W << "].size()=" << Radii[i_W].size()
           << "\n\n" << endl;
      throw(GWFrames_VectorSizeMismatch);
    }
  }

  // Set up the output data, recording everything but the mode data
  GWFrames::Waveforms ExtrapolatedWaveforms(2*NExtrapolations);
  for(unsigned int i_N=0; i_N<NExtrapolations; ++i_N) {
    const int N = ExtrapolationOrders[i_N];
    if(N<0) {
      ExtrapolatedWaveforms[i_N] = FiniteRadiusWaveforms[NFiniteRadii+N];
    } else {
      // Do everything but set the data
      ExtrapolatedWaveforms[i_N].HistoryStream() << "### Extrapolating with N=" << N << "\n";
      ExtrapolatedWaveforms[i_N].AppendHistory("######## Begin old history ########\n");
      ExtrapolatedWaveforms[i_N].AppendHistory(FiniteRadiusWaveforms[NFiniteRadii-1].HistoryStr());
      ExtrapolatedWaveforms[i_N].AppendHistory("######## End old history ########\n");
      ExtrapolatedWaveforms[i_N].SetT(FiniteRadiusWaveforms[NFiniteRadii-1].T());
      ExtrapolatedWaveforms[i_N].SetFrame(FiniteRadiusWaveforms[NFiniteRadii-1].Frame());
      ExtrapolatedWaveforms[i_N].SetFrameType(WaveformFrameType(FiniteRadiusWaveforms[NFiniteRadii-1].FrameType()));
      ExtrapolatedWaveforms[i_N].SetDataType(WaveformDataType(FiniteRadiusWaveforms[NFiniteRadii-1].DataType()));
      ExtrapolatedWaveforms[i_N].SetRIsScaledOut(FiniteRadiusWaveforms[NFiniteRadii-1].RIsScaledOut());
      ExtrapolatedWaveforms[i_N].SetMIsScaledOut(FiniteRadiusWaveforms[NFiniteRadii-1].MIsScaledOut());
      ExtrapolatedWaveforms[i_N].SetLM(FiniteRadiusWaveforms[NFiniteRadii-1].LM());
      ExtrapolatedWaveforms[i_N].ResizeData(NModes, NTimes);
      // Set up the waveforms for sigma
      ExtrapolatedWaveforms[i_N+NExtrapolations].HistoryStream() << "### Extrapolating with N=" << N << "\n";
      ExtrapolatedWaveforms[i_N+NExtrapolations].AppendHistory("######## Begin old history ########\n");
      ExtrapolatedWaveforms[i_N+NExtrapolations].AppendHistory(FiniteRadiusWaveforms[NFiniteRadii-1].HistoryStr());
      ExtrapolatedWaveforms[i_N+NExtrapolations].AppendHistory("######## End old history ########\n");
      ExtrapolatedWaveforms[i_N+NExtrapolations].AppendHistory("### # This Waveform stores the uncertainty estimate for this extrapolation.\n");
      ExtrapolatedWaveforms[i_N+NExtrapolations].SetT(FiniteRadiusWaveforms[NFiniteRadii-1].T());
      ExtrapolatedWaveforms[i_N+NExtrapolations].SetFrame(FiniteRadiusWaveforms[NFiniteRadii-1].Frame());
      ExtrapolatedWaveforms[i_N+NExtrapolations].SetFrameType(WaveformFrameType(FiniteRadiusWaveforms[NFiniteRadii-1].FrameType()));
      ExtrapolatedWaveforms[i_N+NExtrapolations].SetDataType(WaveformDataType(FiniteRadiusWaveforms[NFiniteRadii-1].DataType()));
      ExtrapolatedWaveforms[i_N+NExtrapolations].SetRIsScaledOut(FiniteRadiusWaveforms[NFiniteRadii-1].RIsScaledOut());
      ExtrapolatedWaveforms[i_N+NExtrapolations].SetMIsScaledOut(FiniteRadiusWaveforms[NFiniteRadii-1].MIsScaledOut());
      ExtrapolatedWaveforms[i_N+NExtrapolations].SetLM(FiniteRadiusWaveforms[NFiniteRadii-1].LM());
      ExtrapolatedWaveforms[i_N+NExtrapolations].ResizeData(NModes, NTimes);
    }
  }

  if(MaxN<0) { return ExtrapolatedWaveforms; }
  const unsigned int MaxCoefficients = (unsigned int)(MaxN+1);

  // // Clear the extrapolation-coefficient files and write the headers
  // for(unsigned int i_N=0; i_N<NExtrapolations; ++i_N) {
  //   const int N = ExtrapolationOrders[i_N];
  //   if(N<0) { continue; }
  //   const string FileName = "Extrapolation_N"+NumberToString<int>(N)+"_re_l2_m2.dat";
  //   FILE* ofp = fopen(FileName.c_str(), "w");
  //   fprintf(ofp, "# [1] = t\n");
  //   for(int i_c=0; i_c<N+1; ++i_c) {
  //     fprintf(ofp, "# [%d] = i%d\n", 2*i_c+2, i_c);
  //     fprintf(ofp, "# [%d] = sigma%d\n", 2*i_c+3, i_c);
  //   }
  //   fclose(ofp);
  // }

  // Loop over time
  for(int i_t=0; i_t<NTimes; ++i_t) {

    // Set up the gsl storage variables
    size_t EffectiveRank;
    double ChiSquared;
    gsl_matrix *OneOverRadii, *Covariance;
    gsl_vector *Re, *Im, *FitCoefficients;
    OneOverRadii = gsl_matrix_alloc(NFiniteRadii, MaxCoefficients);
    Covariance = gsl_matrix_alloc(MaxCoefficients, MaxCoefficients);
    Re = gsl_vector_alloc(NFiniteRadii);
    Im = gsl_vector_alloc(NFiniteRadii);
    FitCoefficients = gsl_vector_alloc(MaxCoefficients);

    // Set up the radius data (if we're not using Omega)
    if(! UseOmegas) {
      for(unsigned int i_W=0; i_W<NFiniteRadii; ++i_W) {
        const double OneOverRadius = 1.0/Radii[i_W][i_t];
        gsl_matrix_set(OneOverRadii, i_W, 0, 1.0);
        for(unsigned int i_N=1; i_N<MaxCoefficients; ++i_N) {
          gsl_matrix_set(OneOverRadii, i_W, i_N, OneOverRadius*gsl_matrix_get(OneOverRadii, i_W, i_N-1));
        }
      }
    }

    // Loop over modes
    for(unsigned int i_m=0; i_m<NModes; ++i_m) {

      // Set up the radius data (if we are using Omega)
      const int M = FiniteRadiusWaveforms[0].LM(i_m)[1];
      if(UseOmegas) {
        for(unsigned int i_W=0; i_W<NFiniteRadii; ++i_W) {
          const double OneOverRadius = (M!=0 ? 1.0/(Radii[i_W][i_t]*M*Omegas[i_t]) : 1.0/(Radii[i_W][i_t]));
          gsl_matrix_set(OneOverRadii, i_W, 0, 1.0);
          for(unsigned int i_N=1; i_N<MaxCoefficients; ++i_N) {
            gsl_matrix_set(OneOverRadii, i_W, i_N, OneOverRadius*gsl_matrix_get(OneOverRadii, i_W, i_N-1));
          }
        }
      }

      // Set up the mode data
      for(unsigned int i_W=0; i_W<NFiniteRadii; ++i_W) {
        gsl_vector_set(Re, i_W, FiniteRadiusWaveforms[i_W].Re(i_m,i_t));
        gsl_vector_set(Im, i_W, FiniteRadiusWaveforms[i_W].Im(i_m,i_t));
      }

      // Loop over extrapolation orders
      for(unsigned int i_N=0; i_N<NExtrapolations; ++i_N) {
        const int N = ExtrapolationOrders[i_N];

        // If non-extrapolating, skip to the next one (the copying was
        // done when ExtrapolatedWaveforms[i_N] was constructed)
        if(N<0) {
          continue;
        }

        // Do the extrapolations
        gsl_multifit_linear_workspace* Workspace = gsl_multifit_linear_alloc(NFiniteRadii, N+1);
        gsl_matrix_view OneOverRadii_N = gsl_matrix_submatrix(OneOverRadii, 0, 0, NFiniteRadii, N+1);
        gsl_matrix_view Covariance_N = gsl_matrix_submatrix(Covariance, 0, 0, N+1, N+1);
        gsl_vector_view FitCoefficients_N = gsl_vector_subvector(FitCoefficients, 0, N+1);
        gsl_multifit_linear_usvd(&OneOverRadii_N.matrix, Re, SVDTol, &EffectiveRank, &FitCoefficients_N.vector, &Covariance_N.matrix, &ChiSquared, Workspace);
        const double re = gsl_vector_get(&FitCoefficients_N.vector, 0);
        const double re_err = std::sqrt(2*M_PI*(NFiniteRadii-EffectiveRank)*gsl_matrix_get(&Covariance_N.matrix, 0, 0));
        // if(i_m==0) {
        //   const string FileName = "Extrapolation_N"+NumberToString<int>(N)+"_re_l2_m2.dat";
        //   FILE* ofp = fopen(FileName.c_str(), "a");
        //   fprintf(ofp, "%g ", FiniteRadiusWaveforms[0].T(i_t));
        //   for(int i_c=0; i_c<N+1; ++i_c) {
        //     fprintf(ofp, "%g %g ", gsl_vector_get(&FitCoefficients_N.vector, i_c), std::sqrt(gsl_matrix_get(&Covariance_N.matrix, i_c, i_c)));
        //   }
        //   fprintf(ofp, "\n");
        //   fclose(ofp);
        // }
        gsl_multifit_linear_usvd(&OneOverRadii_N.matrix, Im, SVDTol, &EffectiveRank, &FitCoefficients_N.vector, &Covariance_N.matrix, &ChiSquared, Workspace);
        const double im = gsl_vector_get(&FitCoefficients_N.vector, 0);
        const double im_err = std::sqrt(2*M_PI*(NFiniteRadii-EffectiveRank)*gsl_matrix_get(&Covariance_N.matrix, 0, 0));
        gsl_multifit_linear_free(Workspace);
        ExtrapolatedWaveforms[i_N].SetData(i_m, i_t, std::complex<double>(re,im));
        ExtrapolatedWaveforms[i_N+NExtrapolations].SetData(i_m, i_t, std::complex<double>(re_err,im_err));

      } // Loop over extrapolation orders

    } // Loop over modes

    gsl_vector_free(FitCoefficients);
    gsl_vector_free(Im);
    gsl_vector_free(Re);
    gsl_matrix_free(Covariance);
    gsl_matrix_free(OneOverRadii);

  } // Loop over time

  return ExtrapolatedWaveforms;
}
