// Copyright (c) 2014, Michael Boyle
// See LICENSE file for details

#include <unistd.h>
#include <sys/param.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <climits>
#include <cmath>

#include <sys/time.h>
#include <ctime>

#include <functional>
#include <algorithm>
#include <complex>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multifit.h>
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#include "SpacetimeAlgebra.hpp"
#pragma clang diagnostic pop
#include "Waveforms.hpp"
#include "Utilities.hpp"
#include "Quaternions.hpp"
#include "Quaternions/QuaternionUtilities.hpp"
#include "IntegrateAngularVelocity.hpp"
#include "SphericalFunctions/SWSHs.hpp"
#include "Errors.hpp"

using Quaternions::Quaternion;
using Quaternions::QuaternionArray;
using GWFrames::Matrix;
using SphericalFunctions::LadderOperatorFactorSingleton;
using SphericalFunctions::Wigner3j;
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


// These macros are useful for debugging
#define INFOTOCERR std::cerr << __FILE__ << ":" << __LINE__ << ":" << __func__ << ": "
#define INFOTOCOUT std::cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << ": "

const complex<double> ImaginaryI(0.0,1.0);

const LadderOperatorFactorSingleton& LadderOperatorFactor = LadderOperatorFactorSingleton::Instance();

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
    char* result = getcwd(path, MAXPATHLEN);
    if(!result) {
      cerr << "\n" << __FILE__ << ":" << __LINE__ << ": getcwd error." << endl;
      throw(GWFrames_FailedSystemCall);
    }
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
    char* result = getcwd(path, MAXPATHLEN);
    if(!result) {
      cerr << "\n" << __FILE__ << ":" << __LINE__ << ": getcwd error." << endl;
      throw(GWFrames_FailedSystemCall);
    }
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

/// Copy the Waveform, except for the data (t, frame, lm, data)
GWFrames::Waveform GWFrames::Waveform::CopyWithoutData() const {
  Waveform that;
  that.spinweight = (*this).spinweight;
  that.boostweight = (*this).boostweight;
  that.history.str((*this).history.str());
  that.history.clear();
  that.history.seekp(0, ios_base::end);
  that.frameType = (*this).frameType;
  that.dataType = (*this).dataType;
  that.rIsScaledOut = (*this).rIsScaledOut;
  that.mIsScaledOut = (*this).mIsScaledOut;
  // that.t = (*this).t;
  // that.frame = (*this).frame;
  // that.lm = (*this).lm;
  // that.data = (*this).data;
  that.history << "this->CopyWithoutData();" << std::endl;
  return that;
}

/// Copy the Waveform between indices i_t_a and i_t_b
GWFrames::Waveform GWFrames::Waveform::SliceOfTimeIndices(const unsigned int i_t_a, unsigned int i_t_b) const {
  ///
  /// \param i_t_a Index of initial time
  /// \param i_t_b Index just beyond final time
  ///
  /// i_t_a and i_t_b should hold the indices pointing to the first
  /// time in `t` after `t_a`, and the first time in `t` after `t_b`
  /// (or one-past-the-end of `t` if necessary)
  if(i_t_b==0) {
    i_t_b = i_t_a+1;
  }
  if(i_t_a>i_t_b) {
    INFOTOCERR << ": Requesting impossible slice"
               << "\ni_t_a=" << i_t_a << "  >  i_t_b=" << i_t_b << std::endl;
    throw(GWFrames_EmptyIntersection);
  }
  if(i_t_b>NTimes()) {
    INFOTOCERR << ": Requesting impossible slice"
               << "\ni_t_b=" << i_t_b << "  >  NTimes()=" << NTimes() << std::endl;
    throw(GWFrames_IndexOutOfBounds);
  }
  Waveform Slice = this->CopyWithoutData();
  Slice.history << "this->SliceOfTimeIndices(" << i_t_a << ", " << i_t_b << ");" << std::endl;
  Slice.lm = lm;
  const unsigned int ntimes = i_t_b-i_t_a;
  const unsigned int nmodes = NModes();
  Slice.data.resize(nmodes, ntimes);
  for(unsigned int i_m=0; i_m<nmodes; ++i_m) {
    for(unsigned int i_t=0; i_t<ntimes; ++i_t) {
      Slice.data[i_m][i_t] = data[i_m][i_t+i_t_a];
    }
  }
  if(frame.size() == NTimes()) {
    Slice.frame = vector<Quaternion>(frame.begin()+i_t_a, frame.begin()+i_t_b);
  } else if(frame.size()==1) {
    Slice.frame = frame;
  } else if(frame.size()!=0) {
    INFOTOCERR << " I don't understand what to do with frame data of length " << frame.size() << " in a Waveform with " << NTimes() << " times." << std::endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  Slice.t = vector<double>(t.begin()+i_t_a, t.begin()+i_t_b);
  return Slice;
}

/// Copy the Waveform between t_a and t_b
GWFrames::Waveform GWFrames::Waveform::SliceOfTimes(const double t_a, const double t_b) const {
  const unsigned int i_t_a = Quaternions::hunt(t, t_a);
  const unsigned int i_t_b = Quaternions::huntRight(t, t_b, i_t_a);
  return this->SliceOfTimeIndices(i_t_a, i_t_b);
}

/// Copy of the Waveform between indices i_t_a and i_t_b, only ell=2 modes
GWFrames::Waveform GWFrames::Waveform::SliceOfTimeIndicesWithEll2(const unsigned int i_t_a, unsigned int i_t_b) const {
  /// i_t_a and i_t_b should hold the indices pointing to the first
  /// time in `t` after `t_a`, and the first time in `t` after `t_b`
  /// (or one-past-the-end of `t` if necessary)
  if(i_t_b==0) {
    i_t_b = i_t_a+1;
  }
  if(i_t_a>i_t_b) {
    INFOTOCERR << ": Requesting impossible slice"
               << "\ni_t_a=" << i_t_a << "  >  i_t_b=" << i_t_b << std::endl;
    throw(GWFrames_EmptyIntersection);
  }
  if(i_t_b>NTimes()) {
    INFOTOCERR << ": Requesting impossible slice"
               << "\ni_t_b=" << i_t_b << "  >  NTimes()=" << NTimes() << std::endl;
    throw(GWFrames_IndexOutOfBounds);
  }
  Waveform Slice = this->CopyWithoutData();
  Slice.history << "this->SliceOfTimeIndicesWithEll2(" << i_t_a << ", " << i_t_b << ");" << std::endl;
  Slice.lm = vector<vector<int> >(5, vector<int>(2));
  const unsigned int ntimes = i_t_b-i_t_a;
  const unsigned int nmodes = 5;
  Slice.data.resize(nmodes, ntimes);
  for(int m=-2; m<3; ++m) {
    Slice.lm[m+2][0] = 2;
    Slice.lm[m+2][1] = m;
    unsigned int i_m = FindModeIndex(2, m);
    for(unsigned int i_t=0; i_t<ntimes; ++i_t) {
      Slice.data[m+2][i_t] = data[i_m][i_t+i_t_a];
    }
  }
  if(frame.size() == NTimes()) {
    Slice.frame = vector<Quaternion>(frame.begin()+i_t_a, frame.begin()+i_t_b);
  } else if(frame.size()==1) {
    Slice.frame = frame;
  } else if(frame.size()!=0) {
    INFOTOCERR << " I don't understand what to do with frame data of length " << frame.size() << " in a Waveform with " << NTimes() << " times." << std::endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  Slice.t = vector<double>(t.begin()+i_t_a, t.begin()+i_t_b);
  return Slice;
}

/// Copy of the Waveform between t_a and t_b, only ell=2 modes
GWFrames::Waveform GWFrames::Waveform::SliceOfTimesWithEll2(const double t_a, const double t_b) const {
  const unsigned int i_t_a = Quaternions::hunt(t, t_a);
  const unsigned int i_t_b = Quaternions::huntRight(t, t_b, i_t_a);
  return this->SliceOfTimeIndicesWithEll2(i_t_a, i_t_b);
}

/// Copy of the Waveform between indices i_t_a and i_t_b without mode data
GWFrames::Waveform GWFrames::Waveform::SliceOfTimeIndicesWithoutModes(const unsigned int i_t_a, unsigned int i_t_b) const {
  /// i_t_a and i_t_b should hold the indices pointing to the first
  /// time in `t` after `t_a`, and the first time in `t` after `t_b`
  /// (or one-past-the-end of `t` if necessary)
  Waveform Slice = this->CopyWithoutData();
  Slice.history << "this->SliceOfTimeIndicesWithoutModes(" << i_t_a << ", " << i_t_b << ");" << std::endl;
  if(i_t_b==0) {
    i_t_b = i_t_a+1;
  }
  Slice.lm = vector<vector<int> >(0, vector<int>(2));
  const unsigned int ntimes = i_t_b-i_t_a;
  const unsigned int nmodes = 0;
  Slice.data.resize(nmodes, ntimes);
  if(frame.size() == NTimes()) {
    Slice.frame = vector<Quaternion>(frame.begin()+i_t_a, frame.begin()+i_t_b);
  }
  Slice.t = vector<double>(t.begin()+i_t_a, t.begin()+i_t_b);
  return Slice;
}

/// Copy of the Waveform between t_a and t_b without mode data
GWFrames::Waveform GWFrames::Waveform::SliceOfTimesWithoutModes(const double t_a, const double t_b) const {
  const unsigned int i_t_a = Quaternions::hunt(t, t_a);
  const unsigned int i_t_b = Quaternions::huntRight(t, t_b, i_t_a);
  return this->SliceOfTimeIndicesWithoutModes(i_t_a, i_t_b);
}

/// Remove all data relating to times outside of the given range
GWFrames::Waveform& GWFrames::Waveform::DropTimesOutside(const double t_a, const double t_b) {
  history << "this->DropTimesOutside(" << t_a << ", " << t_b << ");" << std::endl;
  const unsigned int i_t_a = Quaternions::huntRight(t, t_a);
  const unsigned int i_t_b = Quaternions::huntRight(t, t_b, i_t_a);
  vector<vector<complex<double> > > newdata(NModes(), vector<complex<double> >(t.size()));
  for(unsigned int mode=0; mode<NModes(); ++mode) {
    vector<complex<double> > ModeData = Data(mode);
    ModeData.erase(ModeData.begin()+i_t_b, ModeData.end());
    ModeData.erase(ModeData.begin(), ModeData.begin()+i_t_a);
    newdata[mode].swap(ModeData);
  }
  data = MatrixC(newdata);
  frame.erase(frame.begin()+i_t_b, frame.end());
  frame.erase(frame.begin(), frame.begin()+i_t_a);
  t.erase(t.begin()+i_t_b, t.end());
  t.erase(t.begin(), t.begin()+i_t_a);
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

/// Remove data relating to all but the ell=2 modes
GWFrames::Waveform& GWFrames::Waveform::KeepOnlyEll2() {
  vector<unsigned int> Two(1);
  Two[0] = 2;
  return this->KeepOnlyEllModes(Two);
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
    INFOTOCERR << "\nError: Asking to Differentiate a Waveform in the " << GWFrames::WaveformFrameNames[frameType] << " frame."
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
  const unsigned int i_expected = l*(l+1) + m - SpinWeight()*SpinWeight();
  if (lm.size() > i_expected) {
    if(lm[i_expected][0]==l && lm[i_expected][1]==m) { return i_expected; }
  }
  //ORIENTATION!!! following loop
  for(unsigned int i=0; i<NModes(); ++i) {
    if(lm[i][0]==l && lm[i][1]==m) { return i; }
  }
  INFOTOCERR << " Can't find (ell,m)=(" << l << ", " << m << ")" << endl;
  throw(GWFrames_WaveformMissingLMIndex);
}

/// Find index of mode with given (l,m) data without the chance of throwing an exception.
unsigned int GWFrames::Waveform::FindModeIndexWithoutError(const int l, const int m) const {
  /// If the requested mode is not present, the returned index is 1
  /// beyond the end of the mode vector.
  const unsigned int i_expected = l*(l+1) + m - SpinWeight()*SpinWeight();
  if(lm[i_expected][0]==l && lm[i_expected][1]==m) { return i_expected; }
  // ORIENTATION!!! following loop
  unsigned int i=0;
  for(; i<NModes(); ++i) {
    if(lm[i][0]==l && lm[i][1]==m) { return i; }
  }
  ++i;
  return i;
}

/// Return the normalized asymmetry as a function of time
std::vector<double> GWFrames::Waveform::NormalizedAntisymmetry(std::vector<int> LModesForAsymmetry) const {
  /// \param LModesForAsymmetry \f$\ell\f$ modes to use when calculating numerator
  ///
  /// This function just returns the value of the normalized asymmetry
  /// \f$a\f$ defined by Boyle et al. (2014), which is the difference
  /// between the field at a point and its conjugate at the antipodal
  /// point, the amplitude squared and integrated over the sphere.
  /// The normalization is just the `Norm` of the waveform---its
  /// overall power at each instant.
  ///
  /// By default, all ell modes in the data are used for both the
  /// asymmetry and the normalization factor.  If an argument is
  /// input, only modes with ell values in that argument will be used
  /// to calculate the asymmetry.  All modes will always be used to
  /// calculate the normalization.

  const unsigned int ntimes = NTimes();
  const int ellMax = EllMax();
  std::vector<double> asymmetry(ntimes);
  double diff,norm;
  for(unsigned int i_t=0; i_t<ntimes; ++i_t) {
    diff = 0.;
    norm = 0.;
    for(int ell=2; ell<=ellMax; ++ell) {
      if(LModesForAsymmetry.size()==0 || GWFrames::xINy(ell,LModesForAsymmetry)) {
        for(int m=-ell; m<=ell; ++m) {
          const complex<double> h_ell_m = Data(FindModeIndexWithoutError(ell,m), i_t);
          const complex<double> hbar_ell_mm = std::conj(Data(FindModeIndexWithoutError(ell,-m), i_t));
          if(((ell+m)%2)==0) {
            diff += std::norm(h_ell_m - hbar_ell_mm);
          } else {
            diff += std::norm(h_ell_m + hbar_ell_mm);
          }
          norm += std::norm(h_ell_m);
        }
      } else{
        for(int m=-ell; m<=ell; ++m) {
          const complex<double> h_ell_m = Data(FindModeIndexWithoutError(ell,m), i_t);
          norm += std::norm(h_ell_m);
        }
      }
    }
    asymmetry[i_t] = std::sqrt(diff/(4*norm));
  }
  return asymmetry;
}

/// Evaluate the dipole moment of the waveform
std::vector<std::vector<double> > GWFrames::Waveform::DipoleMoment(int ellMax) const {
  /// \param ellMax Maximum ell mode to include [default: all]
  ///
  /// This function evaluates the dipole moment of the waveform's
  /// magnitude, defined as \f$\vec{d} = \int \hat{n} \lvert f
  /// \rvert^2 d\Omega\f$.  Up to a geometric factor, this function
  /// applied to \f$\dot{h}\f$ is the rate of emission of momentum in
  /// gravitational waves.

  if(ellMax==0) {
    ellMax = EllMax();
  }

  vector<vector<double> > D(NTimes(), vector<double>(3));
  vector<complex<double> > d(3);

  for(int i_t=0; i_t<NTimes(); ++i_t) {
    d[0] = 0.0;
    d[1] = 0.0;
    d[2] = 0.0;
    for(int ell=2; ell<=ellMax; ++ell) {
      for(int m=-ell; m<=ell; ++m) {
        for(int ellPrime=std::max(ell-1,2); ellPrime<=std::min(ell+1,ellMax); ++ellPrime) {
          const double sqrtFactor = std::sqrt((2*ell+1)*(2*ellPrime+1)/2.);
          const double Wigner3j_A = Wigner3j(ell, ellPrime, 1, 2, -2, 0);
          for(int mPrime=std::max(m-1,-ellPrime); mPrime<=std::min(m+1,ellPrime); ++mPrime) {
            // This is the whole thing, except for the n_j modes
            if(mPrime%2 == 0) {
              const complex<double> BasicFactor = Data(FindModeIndex(ell,m), i_t) * std::conj(Data(FindModeIndex(ellPrime,mPrime), i_t))
                * sqrtFactor * Wigner3j(ell, ellPrime, 1, m, -mPrime, mPrime-m) * Wigner3j_A;
              if(mPrime==m) { // This will only affect the z component
                d[2] += std::sqrt(2) * BasicFactor;
              } else { // This will only affect the x and y components
                d[0] += (mPrime-m==1 ? -1. : 1.) * BasicFactor;
                d[1] += ImaginaryI * BasicFactor;
              }
            } else {
              const complex<double> BasicFactor = -Data(FindModeIndex(ell,m), i_t) * std::conj(Data(FindModeIndex(ellPrime,mPrime), i_t))
                * sqrtFactor * Wigner3j(ell, ellPrime, 1, m, -mPrime, mPrime-m) * Wigner3j_A;
              if(mPrime==m) { // This will only affect the z component
                d[2] += std::sqrt(2) * BasicFactor;
              } else { // This will only affect the x and y components
                d[0] += (mPrime-m==1 ? -1. : 1.) * BasicFactor;
                d[1] += ImaginaryI * BasicFactor;
              }
            }
            // if(d[0]!=d[0] || d[1]!=d[1] || d[2]!=d[2]) { // DEBUG NANs
            //   INFOTOCERR << i_t << ": (" << ell << "," << m << "); (" << ellPrime << "," << mPrime
            //              << ")\n" << Data(FindModeIndex(ell,m), i_t) << " * " << std::conj(Data(FindModeIndex(ellPrime,mPrime), i_t))
            //              << " * " << std::pow(-1,mPrime) << "*" << sqrtFactor << " * " << Wigner3j(ell, ellPrime, 1, m, -mPrime, mPrime-m)
            //              << " * " << Wigner3j_A
            //              << "=" << BasicFactor << "\t[" << d[0] << "," << d[1] << "," << d[2] << "]" << std::endl;
            //   throw(GWFrames_ValueError);
            // }
          }
        }
      }
    }
    D[i_t][0] = d[0].real();
    D[i_t][1] = d[1].real();
    D[i_t][2] = d[2].real();
  }

  return D;
}


// Local utility function for evaluating involution violation
std::vector<double> GWFrames::Waveform::InvolutionViolationSquared(GWFrames::Waveform::WaveformInvolutionFunction Invol, std::vector<int> Lmodes) const {
  vector<double> violation(NTimes(),0.0);

  for(int ell=std::abs(SpinWeight()); ell<=EllMax(); ++ell) {
    if(Lmodes.size()==0 || GWFrames::xINy(ell,Lmodes)) {
      for(int m=-ell; m<=ell; ++m) {
        const unsigned int i_m = FindModeIndex(ell,m);
        const unsigned int i_mm = FindModeIndex(ell,-m);
        for(unsigned int i_t=0; i_t<NTimes(); ++i_t) {
          violation[i_t] += 0.25 * std::norm( Data(i_m,i_t) - (this->*Invol)(i_m,i_mm,i_t) );
        }
      }
    }
  }

  return violation;
}


// Private utility function for taking involution of a Waveform
GWFrames::Waveform GWFrames::Waveform::Involution(WaveformInvolutionFunction WInvol, QuaternionInvolutionFunction QInvol) const {
  const Waveform& A = *this;
  Waveform B = A.CopyWithoutData();
  const int ellMax = A.EllMax();
  const unsigned int nTimes = A.NTimes();
  B.t = A.t;
  B.frame.resize(A.frame.size());
  for(unsigned int i=0; i<A.frame.size(); ++i) {
    B.frame[i] = QInvol(A.frame[i]);
  }
  B.lm = A.lm;
  B.data.resize(A.NModes(), A.NTimes());
  for(int ell=std::abs(SpinWeight()); ell<=ellMax; ++ell) {
    for(int m=-ell; m<=ell; ++m) {
      const unsigned int i_m = FindModeIndex(ell,m);
      const unsigned int i_mm = FindModeIndex(ell,-m);
      for(unsigned int i_t=0; i_t<nTimes; ++i_t) {
        B.data[i_m][i_t] = (A.*WInvol)(i_m,i_mm,i_t);
      }
    }
  }
  return B;
}

// Private utility function for taking involution of a Waveform
GWFrames::Waveform GWFrames::Waveform::InvolutionSymmetricPart(WaveformInvolutionFunction WInvol, QuaternionInvolutionFunction QInvol) const {
  const Waveform& A = *this;
  Waveform B = A.CopyWithoutData();
  const int ellMax = A.EllMax();
  const unsigned int nTimes = A.NTimes();
  B.t = A.t;
  B.frame.resize(A.frame.size());
  for(unsigned int i=0; i<A.frame.size(); ++i) {
    B.frame[i] = 0.5 * ( A.frame[i] +  QInvol(A.frame[i]) );
  }
  B.lm = A.lm;
  B.data.resize(A.NModes(), A.NTimes());
  for(int ell=std::abs(SpinWeight()); ell<=ellMax; ++ell) {
    for(int m=-ell; m<=ell; ++m) {
      const unsigned int i_m = FindModeIndex(ell,m);
      const unsigned int i_mm = FindModeIndex(ell,-m);
      for(unsigned int i_t=0; i_t<nTimes; ++i_t) {
        B.data[i_m][i_t] = 0.5 * ( A.data[i_m][i_t] + (A.*WInvol)(i_m,i_mm,i_t) );
      }
    }
  }
  return B;
}

// Private utility function for taking involution of a Waveform
GWFrames::Waveform GWFrames::Waveform::InvolutionAntisymmetricPart(WaveformInvolutionFunction WInvol, QuaternionInvolutionFunction QInvol) const {
  const Waveform& A = *this;
  Waveform B = A.CopyWithoutData();
  const int ellMax = A.EllMax();
  const unsigned int nTimes = A.NTimes();
  B.t = A.t;
  B.frame.resize(A.frame.size());
  for(unsigned int i=0; i<A.frame.size(); ++i) {
    B.frame[i] = 0.5 * ( A.frame[i] -  QInvol(A.frame[i]) );
  }
  B.lm = A.lm;
  B.data.resize(A.NModes(), A.NTimes());
  for(int ell=std::abs(SpinWeight()); ell<=ellMax; ++ell) {
    for(int m=-ell; m<=ell; ++m) {
      const unsigned int i_m = FindModeIndex(ell,m);
      const unsigned int i_mm = FindModeIndex(ell,-m);
      for(unsigned int i_t=0; i_t<nTimes; ++i_t) {
        B.data[i_m][i_t] = 0.5 * ( A.data[i_m][i_t] - (A.*WInvol)(i_m,i_mm,i_t) );
      }
    }
  }
  return B;
}



// The following are local objects used by `MinimalParityViolation`
#ifndef DOXYGEN
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
void myGSLErrorHandler2 (const char * reason, const char * file, int line, int gsl_errno) {
  std::cerr << "\n\n" << file << ":" << line << ": " << reason << std::endl;
  throw(gsl_errno);
}
gsl_error_handler_t* defaultGSLErrorHandler3 = gsl_set_error_handler((gsl_error_handler_t*) &myGSLErrorHandler2);
double minfunc_MinimalParityViolation (const gsl_vector* delta, void* params);
class MinimalParityViolationMinimizer {
public:
  const GWFrames::Waveform& Win;
  GWFrames::Waveform W;
  Quaternions::Quaternion R_last;
  const std::vector<int> EllEqualsTwo;
  gsl_multimin_fminimizer* s;
  gsl_vector* ss;
  gsl_vector* x;
  gsl_multimin_function min_func;
public:
  MinimalParityViolationMinimizer(const GWFrames::Waveform& WIN)
    : Win(WIN), W(Win.SliceOfTimeIndicesWithEll2(Win.NTimes()/2)), R_last(), EllEqualsTwo(1,2)
  {
    min_func.n = 2;
    min_func.f = &minfunc_MinimalParityViolation;
    min_func.params = (void*) this;
    R_last = Quaternions::sqrtOfRotor(-Quaternions::zHat*Quaternions::Quaternion(W.LLDominantEigenvector()[0]));
    s = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, min_func.n);
    ss = gsl_vector_alloc(min_func.n);
    x = gsl_vector_alloc(min_func.n);
  }
  ~MinimalParityViolationMinimizer() {
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);
  }
  double EvaluateMinimizationQuantity(const double deltax, const double deltay) const {
    const GWFrames::Waveform W2 = GWFrames::Waveform(W).RotateDecompositionBasis(Quaternions::exp(Quaternions::Quaternion(0.0,deltax,deltay,0.0)));
    return W2.ZParityViolationSquared()[0];
  }
  double Minimize(const unsigned int i, const int direction) {
    switch (direction) {
    case 0:
      R_last = Quaternions::Quaternion(0.707106781187, 0, 0.707106781187, 0);
      break;
    case 1:
      R_last = Quaternions::Quaternion(0.707106781187, -0.707106781187, 0, 0);
      break;
    case 2:
      R_last = Quaternions::One;
      break;
    default:
      throw(GWFrames_ValueError);
    }
    W = Win.SliceOfTimeIndices(i);
    W.RotateDecompositionBasis(R_last);
    const unsigned int MaxIterations = 2000;
    const double MinSimplexSize = 1.0e-8;
    const double InitialTrialAngleStep = M_PI/8.0;
    size_t iter = 0;
    int status = GSL_CONTINUE;
    double size = 0.0;

    // Explicit optimization
    gsl_vector_set(x, 0, 0.);
    gsl_vector_set(x, 1, 0.);
    gsl_vector_set(ss, 0, InitialTrialAngleStep);
    gsl_vector_set(ss, 1, InitialTrialAngleStep);
    gsl_multimin_fminimizer_set(s, &min_func, x, ss);
    // Run the minimization
    while(status == GSL_CONTINUE && iter < MaxIterations) {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      if(status==GSL_EBADFUNC) {
        INFOTOCERR << ":\nThe iteration encountered a singular point where the function evaluated to Inf or NaN"
                   << "\nwhile minimizing at (" << gsl_vector_get(s->x, 0) << ", " << gsl_vector_get(s->x, 1)
                   << ")." << std::endl;
      }
      if(status==GSL_FAILURE) {
        INFOTOCERR << ":\nThe algorithm could not improve the current best approximation or bounding interval." << std::endl;
      }
      if(status==GSL_ENOPROG) {
        INFOTOCERR << ":\nThe minimizer is unable to improve on its current estimate, either due to"
                   << "\nnumerical difficulty or because a genuine local minimum has been reached." << std::endl;
      }
      if(status) break;
      size = gsl_multimin_fminimizer_size(s);
      status = gsl_multimin_test_size(size, MinSimplexSize);
    }
    if(iter>=MaxIterations) {
      INFOTOCERR << "Warning: stopped minimizing due to too many iterations." << std::endl;
    }
    // Get rotation and normalized value
    const Quaternions::Quaternion R_delta = Quaternions::exp(Quaternions::Quaternion(0.0, gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1), 0.0));
    R_last = R_delta*R_last;
    W.RotateDecompositionBasis(R_delta);
    return W.ZParityViolationNormalized()[0];
  }
};
double minfunc_MinimalParityViolation (const gsl_vector* delta, void* params) {
  MinimalParityViolationMinimizer* Minimizer = (MinimalParityViolationMinimizer*) params;
  return Minimizer->EvaluateMinimizationQuantity(gsl_vector_get(delta,0),
                                                 gsl_vector_get(delta,1));
}
#endif // DOXYGEN

/// Measure the relative magnitude of the violation of parity in the z direction
std::vector<double> GWFrames::Waveform::MinimalParityViolation() const {
  /// This function measures the violation of invariance under
  /// z-parity (reflection across the x-y plane).  Nonprecessing
  /// systems in a suitable frame should have zero violation.
  /// Precessing systems in any frame and nonprecessing systems in the
  /// wrong frame will show violations.  The quantity is normalized by
  /// the overall norm of the data at each instant, and the
  /// square-root of that ratio is taken.
  ///
  /// This is performed iteratively at each time step, as the system
  /// is rotated, and the parity violation is minimized.

  const unsigned int ntimes = NTimes();
  vector<double> violations(ntimes);
  MinimalParityViolationMinimizer Minimizer(*this);

  for(unsigned int i_t=0; i_t<ntimes; ++i_t) {
    INFOTOCERR << i_t+1 << "/" << ntimes << std::endl;
    violations[i_t] = std::min(std::min(Minimizer.Minimize(i_t,0), Minimizer.Minimize(i_t,1)), Minimizer.Minimize(i_t,2));
  }

  return violations;
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
  /// Note that this function does not change the `frameType`; that is
  /// left to the calling function.
  ///

  history << "this->RotatePhysicalSystem(R_phys); // R_phys=[";
  if(R_phys.size()>0) {
    history << std::setprecision(16) << R_phys[0];
  }
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
    INFOTOCERR << "\nWarning: You have passed an empty set of rotors to `RotateDecompositionBasis`.  I will"
               << "\n         do nothing, and simply return the Waveform as is.  But this is most unusual!"
               << std::endl;
    return *this;
  }
  if(R_frame.size()!=NTimes()) {
    INFOTOCERR << "\nError: (R_frame.size()=" << R_frame.size() << ") != (NTimes()=" << NTimes() << ")" << std::endl;
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

/// Calculate the \f$<L \partial_t>\f$ quantity defined in the paper.
vector<vector<double> > GWFrames::Waveform::LdtVector(vector<int> Lmodes) const {
  ///
  /// \param Lmodes L modes to evaluate
  ///
  /// If Lmodes is empty (default), all L modes are used.  Setting
  /// Lmodes to [2] or [2,3,4], for example, restricts the range of
  /// the sum.  The vector is given with respect to the (possibly
  /// rotating) mode frame (X,Y,Z), rather than the inertial frame
  /// (x,y,z).
  ///
  /// \f$<L \partial_t>^a = \sum_{\ell,m,m'} \Im [ \bar{f}^{\ell,m'} < \ell,m' | L_a | \ell,m > \dot{f}^{\ell,m} ]\f$

  // L+ = Lx + i Ly      Lx =    (L+ + L-) / 2     Im(Lx) =  ( Im(L+) + Im(L-) ) / 2
  // L- = Lx - i Ly      Ly = -i (L+ - L-) / 2     Im(Ly) = -( Re(L+) - Re(L-) ) / 2
  // Lz = Lz             Lz = Lz                   Im(Lz) = Im(Lz)

  const LadderOperatorFactorSingleton& LadderOperatorFactor = LadderOperatorFactorSingleton::Instance();
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
  /// the sum.  The matrix is given in the (possibly rotating) mode
  /// frame (X,Y,Z), rather than the inertial frame (x,y,z).
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
  const complex<double> ImaginaryI(0.0,1.0);
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
        const complex<double> LxLy = -0.25 * ImaginaryI * (LpLp - LmLm + LmLp - LpLm);
        const complex<double> LxLz = 0.5 * (LpLz + LmLz);
        const complex<double> LyLx = -0.25 * ImaginaryI * (LpLp - LmLp + LpLm - LmLm);
        const complex<double> LyLy = -0.25 * (LpLp - LmLp - LpLm + LmLm);
        const complex<double> LyLz = -0.5 * ImaginaryI * (LpLz - LmLz);
        const complex<double> LzLx = 0.5 * (LzLp + LzLm);
        const complex<double> LzLy = -0.5 * ImaginaryI * (LzLp - LzLm);
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
std::vector<std::vector<double> > GWFrames::Waveform::LLDominantEigenvector(const std::vector<int>& Lmodes,
                                                                            const Quaternions::Quaternion& RoughInitialEllDirection) const {
  ///
  /// \param Lmodes L modes to evaluate (optional)
  /// \param RoughInitialEllDirection Vague guess about the preferred initial (optional)
  ///
  /// If Lmodes is empty (default), all L modes are used.  Setting
  /// Lmodes to [2] or [2,3,4], for example, restricts the range of
  /// the sum.
  ///
  /// Ell is the direction of the angular velocity for a PN system, so
  /// some rough guess about that direction allows us to choose the
  /// direction of the eigenvectors output by this function to be more
  /// parallel than anti-parallel to that direction.  The default is
  /// to simply choose the z axis, since this is most often the
  /// correct choice anyway.
  ///
  /// The vector is given in the (possibly rotating) mode frame
  /// (X,Y,Z), rather than the inertial frame (x,y,z).

  // Calculate the LL matrix at each instant
  vector<Matrix> ll = LLMatrix(Lmodes);

  // Calculate the dominant principal axis (dpa) of LL at each instant
  vector<vector<double> > dpa(NTimes(), vector<double>(3));
  for(unsigned int i=0; i<ll.size(); ++i) {
    dpa[i] = GWFrames::DominantPrincipalAxis(ll[i]);
  }

  // Make the initial direction closer to RoughInitialEllDirection than not
  if(RoughInitialEllDirection.dot(dpa[0])<0.) {
    dpa[0][0] = -dpa[0][0];
    dpa[0][1] = -dpa[0][1];
    dpa[0][2] = -dpa[0][2];
  }

  // Now, go through and make the vectors reasonably continuous.
  if(dpa.size()>0) {
    double LastNorm = std::sqrt(dpa[0][0]*dpa[0][0]+dpa[0][1]*dpa[0][1]+dpa[0][2]*dpa[0][2]);
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
      // While we're here, let's just normalize that last one
      if(LastNorm!=0.0) {
        dpa[i-1][0] = dpa[i-1][0] / LastNorm;
        dpa[i-1][1] = dpa[i-1][1] / LastNorm;
        dpa[i-1][2] = dpa[i-1][2] / LastNorm;
      }
      LastNorm = Norm;
    }
    if(LastNorm!=0.0) {
      const unsigned int i=dpa.size();
      dpa[i-1][0] = dpa[i-1][0] / LastNorm;
      dpa[i-1][1] = dpa[i-1][1] / LastNorm;
      dpa[i-1][2] = dpa[i-1][2] / LastNorm;
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
  /// frame"</a>.  The vector is given in the (possibly rotating) mode
  /// frame (X,Y,Z), rather than the inertial frame (x,y,z).
  ///
  /// If Lmodes is empty (default), all L modes are used.  Setting
  /// Lmodes to [2] or [2,3,4], for example, restricts the range of
  /// the sum.
  ///

  // Calculate the L vector and LL matrix at each instant
  vector<vector<double> > l = LdtVector(Lmodes);
  vector<Matrix> ll = LLMatrix(Lmodes);

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
  }

  // Free the memory
  gsl_permutation_free(p);
  gsl_vector_free(x);

  return omega;
}

/// Calculate the angular velocity of the Waveform.
vector<vector<double> > GWFrames::Waveform::AngularVelocityVectorRelativeToInertial(const vector<int>& Lmodes) const {
  ///
  /// \param Lmodes L modes to evaluate
  ///
  /// This returns the angular velocity of the Waveform, as defined in
  /// Sec. II of <a href="http://arxiv.org/abs/1302.2919">"Angular
  /// velocity of gravitational radiation and the corotating
  /// frame"</a>.  The vector is given in the inertial frame (x,y,z),
  /// rather than the (possibly rotating) mode frame (X,Y,Z), which is
  /// unlike most other methods for `Waveform`.
  ///
  /// If Lmodes is empty (default), all L modes are used.  Setting
  /// Lmodes to [2] or [2,3,4], for example, restricts the range of
  /// the sum.
  ///

  vector<vector<double> > omega = this->AngularVelocityVector(Lmodes);

  // If the frame is nontrivial, include its contribution
  const bool TimeDependentFrame = (frame.size()>1);
  const vector<Quaternion> Rdot = (TimeDependentFrame
                                   ? Quaternions::QuaternionDerivative(frame, t)
                                   : vector<Quaternion>(0));
  const bool ConstantNontrivialFrame = (frame.size()==1);
  const Quaternion R0 = (ConstantNontrivialFrame
                         ? frame[0]
                         : Quaternion(1,0,0,0));

  if(!TimeDependentFrame && !ConstantNontrivialFrame) {
    return omega;
  }

  // Loop through time steps
  if(TimeDependentFrame) { // Include frame-rotation effects
    for(unsigned int iTime=0; iTime<omega.size(); ++iTime) {
      const Quaternion& R = frame[iTime];
      omega[iTime] = (R*Quaternion(omega[iTime])*R.conjugate() + 2*Rdot[iTime]*R.conjugate()).vec();
    }
  } else if(ConstantNontrivialFrame) { // Just rotate the resul
    for(unsigned int iTime=0; iTime<omega.size(); ++iTime) {
      omega[iTime] = (R0*Quaternion(omega[iTime])*R0.conjugate()).vec();
    }
  }

  return omega;
}

/// Frame in which the rotation is minimal.
std::vector<Quaternions::Quaternion> GWFrames::Waveform::CorotatingFrame(const std::vector<int>& Lmodes) const {
  ///
  /// \param Lmodes L modes to evaluate
  ///
  /// This function combines the steps required to obtain the
  /// corotating frame.  The result is given relative to the (possibly
  /// rotating) mode frame (X,Y,Z), rather than the inertial frame
  /// (x,y,z).
  ///
  /// If Lmodes is empty (default), all L modes are used.  Setting
  /// Lmodes to [2] or [2,3,4], for example, restricts the range of
  /// the sum.
  ///

  return Quaternions::FrameFromAngularVelocity(QuaternionArray(AngularVelocityVector(Lmodes)), T());
}

/// Transform Waveform to co-precessing frame.
GWFrames::Waveform& GWFrames::Waveform::TransformToCoprecessingFrame(const std::vector<int>& Lmodes) {
  ///
  /// \param Lmodes L modes to evaluate
  ///
  /// This function combines the steps required to obtain the Waveform
  /// in a frame with its `z` axis aligned with the vector V_h defined
  /// by O'Shaughnessy et al. [2011], and its alignment about the `z`
  /// axis satisfying the minimal-rotation condition, as given by
  /// Boyle, Owen, Pfeiffer [2011].
  ///
  /// If Lmodes is empty (default), all L modes are used.  Setting
  /// Lmodes to [2] or [2,3,4], for example, restricts the range of
  /// the sum.
  ///
  Quaternions::Quaternion RoughInitialEllDirection;
  const unsigned int NPointsForDeriv = 7;
  if(NTimes()<=NPointsForDeriv || FrameType()==GWFrames::Coorbital || FrameType()==GWFrames::Corotating) {
    // In either of these frame types, the angular velocity is likely
    // to be nonsense.  So we might as well revert to the default,
    // which is at least more likely to be correct than the ~random
    // output of GSL.
    RoughInitialEllDirection = Quaternions::zHat;
  } else {
    const Waveform Segment = SliceOfTimeIndicesWithEll2(0, NPointsForDeriv);
    RoughInitialEllDirection = Quaternions::Quaternion(Segment.AngularVelocityVector()[NPointsForDeriv/2]); // Using integer division
  }
  history << "this->TransformToCoprecessingFrame(" << StringForm(Lmodes) << ")\n#";
  SetFrameType(GWFrames::Coprecessing);
  return this->RotateDecompositionBasis(FrameFromZ(normalized(QuaternionArray(this->LLDominantEigenvector(Lmodes, RoughInitialEllDirection))), T()));
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
  /// constant unset.  To set it, the modes can be rotated so that
  /// they are aligned with the frame using, for example, the function
  /// `AlignDecompositionFrameToModes`.
  ///
  /// If Lmodes is empty (default), all L modes are used.  Setting
  /// Lmodes to [2] or [2,3,4], for example, restricts the range of
  /// the sum.
  ///

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
    INFOTOCERR << "\nWarning: Waveform is already in the " << GWFrames::WaveformFrameNames[GWFrames::Inertial] << " frame;"
               << " this function will have no effect."
               << "\n         If it's not really in that frame, tell the Waveform first.\n"
               << std::endl;
    return *this;
  }

  if(frameType == GWFrames::UnknownFrameType) {
    INFOTOCERR << "\nWarning: Waveform frame type is " << GWFrames::WaveformFrameNames[GWFrames::UnknownFrameType] << "."
               << "\n         I assume you know what you're doing...\n"
               << std::endl;
    return *this;
  }

  if(frame.size() == 0) {
    INFOTOCERR << "\nWarning: frame.size()=" << frame.size() << ".  I will assume that this *is* the"
               << "\n         inertial frame, so this function will have no effect.\n"
               << "\n         And I will continue to assume you know what you're doing...\n"
               << std::endl;
    return *this;
  }

  if(frame.size() != NTimes()) {
    INFOTOCERR << "\nError: (frame.size()=" << frame.size() << ") != (NTimes()=" << NTimes() << ")."
               << "\n       I don't know what to do with this data, or what the inertial frame is."
               << std::endl;
    throw(GWFrames_VectorSizeMismatch);
  }

  history << "this->TransformToInertialFrame();\n#";
  this->frameType = GWFrames::Inertial; // Must come first
  this->RotateDecompositionBasis(Quaternions::conjugate(frame));
  this->SetFrame(vector<Quaternion>(0));
  return *this;
}

/// Interpolate the Waveform to a new set of time instants.
GWFrames::Waveform GWFrames::Waveform::Interpolate(const std::vector<double>& NewTime, const bool AllowTimesOutsideCurrentDomain) const {
  /// \param NewTime New vector of times to which this interpolates
  /// \param AllowTimesOutsideCurrentDomain [Default: false]
  ///
  /// If `AllowTimesOutsideCurrentDomain` is true, the values of all
  /// modes will be set to 0.0 for times outside the current set of
  /// time data.  If false, and such times are requested, an error
  /// will be thrown.
  ///
  if(NewTime.size()==0) {
    INFOTOCERR << ": Asking for empty Waveform." << std::endl;
    throw(GWFrames_EmptyIntersection);
  }
  unsigned int i0=0, i1=NewTime.size()-1;
  const unsigned int i2 = NewTime.size();
  vector<double> NewTimesInsideCurrentDomain;
  if(AllowTimesOutsideCurrentDomain) {
    // Set the indices in the dumbest way possible
    while(NewTime[i0]<t[0]) { ++i0; }
    while(NewTime[i1]>t.back() && i1>0) { --i1; }
    ++i1;
    // Now, i0 is the first index in NewTime for which a current time
    // exists, and i1 is 1 beyond the last index in NewTime for which
    // a current time exists.
    NewTimesInsideCurrentDomain.resize(i1-i0);
    std::copy(NewTime.begin()+i0, NewTime.begin()+i1, NewTimesInsideCurrentDomain.begin());
  } else {
    if(NewTime[0]<t[0]) {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Asking for extrapolation; we only do interpolation.\n"
                << "NewTime[0]=" << NewTime[0] << "\tt[0]=" << t[0]
                << "\nMaybe you meant to pass the `AllowTimesOutsideCurrentDomain=true` flag..." << std::endl;
      throw(GWFrames_EmptyIntersection);
    }
    if(NewTime.back()>t.back()) {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Asking for extrapolation; we only do interpolation.\n"
                << "NewTime.back()=" << NewTime.back() << "\tt.back()=" << t.back()
                << "\nMaybe you meant to pass the `AllowTimesOutsideCurrentDomain=true` flag..."  << std::endl;
      throw(GWFrames_EmptyIntersection);
    }
  }

  Waveform C;
  C.history << HistoryStr()
            << "### *this = this->Interpolate(NewTime," << AllowTimesOutsideCurrentDomain << ");" << std::endl;
  C.t = NewTime;
  if(frame.size()==1) { // Assume we have just a constant non-trivial frame
    C.frame = frame;
  } else if(frame.size()>1) { // Assume we have frame data for each time step
    if(AllowTimesOutsideCurrentDomain) {
      C.frame.resize(NewTime.size());
      const std::vector<Quaternion> NewFrame = Squad(frame, t, NewTimesInsideCurrentDomain);
      for(unsigned int i=0; i<i0; ++i) {
        C.frame[i] = NewFrame[0];
      }
      for(unsigned int i=i0; i<i1; ++i) {
        C.frame[i] = NewFrame[i-i0];
      }
      for(unsigned int i=i1; i<i2; ++i) {
        C.frame[i] = NewFrame.back();
      }
    } else {
      C.frame = Squad(frame, t, NewTime);
    }
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
    if(AllowTimesOutsideCurrentDomain) {
      for(unsigned int i_t=0; i_t<i0; ++i_t) {
        C.data[i_m][i_t] = complex<double>( 0., 0. );
      }
      for(unsigned int i_t=i0; i_t<i1; ++i_t) {
        C.data[i_m][i_t] = complex<double>( gsl_spline_eval(splineRe, NewTimesInsideCurrentDomain[i_t], accRe),
                                            gsl_spline_eval(splineIm, NewTimesInsideCurrentDomain[i_t], accIm) );
      }
      for(unsigned int i_t=i1; i_t<i2; ++i_t) {
        C.data[i_m][i_t] = complex<double>( 0., 0. );
      }
    } else {
      for(unsigned int i_t=0; i_t<C.t.size(); ++i_t) {
        C.data[i_m][i_t] = complex<double>( gsl_spline_eval(splineRe, C.t[i_t], accRe), gsl_spline_eval(splineIm, C.t[i_t], accIm) );
      }
    }
  }
  // Free the interpolators
  gsl_interp_accel_free(accRe);
  gsl_interp_accel_free(accIm);
  gsl_spline_free(splineRe);
  gsl_spline_free(splineIm);

  return C;
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

/// Find the appropriate rotations to fix the attitude of the corotating frame.
std::vector<Quaternions::Quaternion> GWFrames::Waveform::GetAlignmentsOfDecompositionFrameToModes(const std::vector<int>& Lmodes) const {
  ///
  /// \param Lmodes Lmodes to use in computing \f$<LL>\f$
  ///
  /// This function finds the appropriate pre-multiplied rotation
  /// \f$R_{\varepsilon}\f$ so that the decomposition frame is aligned
  /// to the waveform.  This particular version finds the appropriate
  /// \f$R_{\varepsilon}\f$ at each time in the input Waveform.  This
  /// is useful in cases where we need to try many such alignments,
  /// because the setup for interpolation is very slow.
  ///
  /// Note that this function has no option to choose the direction of
  /// X based on some nHat vector, as other similar functions have.
  /// That issue is assumed to be handled elsewhere.

  if(frameType!=GWFrames::Coprecessing && frameType!=GWFrames::Coorbital && frameType!=GWFrames::Corotating) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ":"
              << "\nError: GetAlignmentOfDecompositionFrameToModes only takes Waveforms in the "
              << GWFrames::WaveformFrameNames[GWFrames::Coprecessing] << ", "
              << GWFrames::WaveformFrameNames[GWFrames::Coorbital] << ", or "
              << GWFrames::WaveformFrameNames[GWFrames::Corotating] << " frames."
              << "\n       This Waveform is in the " << GWFrames::WaveformFrameNames[frameType] << " frame." << std::endl;
    throw(GWFrames_WrongFrameType);
  }

  if(frame.size()!=NTimes()) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ":"
              << "\nError: GetAlignmentOfDecompositionFrameToModes requires full information about the Waveform's frame."
              << "\n       This Waveform has " << NTimes() << " time steps, but " << frame.size() << " rotors in its frame data." << std::endl;
    throw(GWFrames_VectorSizeMismatch);
  }

  // This will be the returned quantity
  std::vector<Quaternion> R_eps(NTimes());

  // Get direction of angular-velocity vector at each time step, in this frame
  const vector<Quaternion> omegaHat
    = Quaternions::inverse(frame)
    * Quaternions::normalized(Quaternions::QuaternionArray(this->SliceOfTimesWithEll2().AngularVelocityVectorRelativeToInertial())) * frame;

  const vector<vector<double> > V_h = this->LLDominantEigenvector(Lmodes);

  const unsigned int i_22 = FindModeIndex(2,2);
  const unsigned int i_2m2 = FindModeIndex(2,-2);

  for(unsigned int i_t=0; i_t<NTimes(); ++i_t) {
    // Choose the normalized eigenvector more parallel to omegaHat than anti-parallel
    const Quaternion V_hi = (omegaHat[i_t].dot(V_h[i_t]) < 0 ? -Quaternions::normalized(V_h[i_t]) : Quaternions::normalized(V_h[i_t]));

    // R_V_hi is the rotor taking the Z axis onto V_hi
    const Quaternion R_V_hi = Quaternions::sqrtOfRotor(-V_hi*Quaternions::zHat);

    // Now rotate Instant so that its z axis is aligned with V_f
    Waveform Instant = this->SliceOfTimeIndicesWithEll2(i_t);
    Instant.RotateDecompositionBasis(R_V_hi);

    // Get the phase of the (2,+/-2) modes after rotation
    const double phase_22 = std::atan2(Instant.Im(i_22,0),Instant.Re(i_22,0));
    const double phase_2m2 = std::atan2(Instant.Im(i_2m2,0),Instant.Re(i_2m2,0));

    // R_eps is the rotation we will be applying on the right-hand side
    R_eps[i_t] = R_V_hi * Quaternions::exp(Quaternions::Quaternion(0,0,0,(-(phase_22+phase_2m2)/8.)));
  }

  return UnflipRotors(R_eps);
}

/// Find the appropriate rotation to fix the attitude of the corotating frame.
Quaternions::Quaternion GWFrames::Waveform::GetAlignmentOfDecompositionFrameToModes(const double t_fid, const Quaternions::Quaternion& nHat_t_fid,
                                                                                    const std::vector<int>& Lmodes) const {
  ///
  /// \param t_fid Fiducial time at which the alignment should happen
  /// \param nHat_t_fid The approximate direction of nHat at t_fid
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

  if(frame.size()!=NTimes()) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ":"
              << "\nError: GetAlignmentOfDecompositionFrameToModes requires full information about the Waveform's frame."
              << "\n       This Waveform has " << NTimes() << " time steps, but only " << frame.size() << " rotors in its frame." << std::endl;
  }

  if(t_fid<t[0] || t_fid>t.back()) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ":"
              << "\nError: The requested alignment time t_fid=" << t_fid << " is outside the range of times in this waveform ("
              << t[0] << ", " << t.back() << ")." << std::endl;
    throw(GWFrames_EmptyIntersection);
  }

  Quaternions::Quaternion R_eps;

  // Get direction of angular-velocity vector near t_fid
  int i_t_fid = Quaternions::huntRight(t, t_fid);
  unsigned int i1 = (i_t_fid-5<0 ? 0 : i_t_fid-5);
  unsigned int i2 = (i1+11>int(t.size()) ? t.size() : i1+11);
  const Waveform Region = (this->SliceOfTimeIndicesWithEll2(i1,i2)).TransformToInertialFrame();
  Quaternion omegaHat = Quaternion(Region.AngularVelocityVector()[i_t_fid-i1]).normalized();
  // omegaHat contains the components of that vector relative to the
  // inertial frame.  To get its components in this Waveform's
  // (possibly rotating) frame, we need to rotate it by the inverse
  // of this Waveform's `frame` data:
  if(Frame().size()>1) {
    const Quaternion& R = Frame(i_t_fid);
    // INFOTOCERR << "Rotating omegaHat by " << R.inverse() << endl;
    omegaHat = R.inverse() * omegaHat * R;
  } else if(Frame().size()==1) {
    const Quaternion& R = Frame(0);
    // INFOTOCERR << "Rotating omegaHat by " << R.inverse() << endl;
    omegaHat = R.inverse() * omegaHat * R;
  }

  // Interpolate the Waveform to t_fid
  Waveform Instant = this->SliceOfTimeIndices(i1,i2).Interpolate(vector<double>(1,t_fid));
  const Quaternion R_f0 = Instant.Frame(0);

  // V_f is the dominant eigenvector of <LL>, suggested by O'Shaughnessy et al.
  const Quaternion V_f = Quaternions::Quaternion(Instant.LLDominantEigenvector(Lmodes)[0]).normalized();
  const Quaternion V_f_aligned = (omegaHat.dot(V_f) < 0 ? -V_f : V_f);

  // R_V_f is the rotor taking the Z axis onto V_f
  const Quaternion R_V_f = Quaternions::sqrtOfRotor(-V_f_aligned*Quaternions::zHat);
  // INFOTOCERR << omegaHat << "\n"
  //            << V_f << "\n"
  //            << V_f_aligned << "\n"
  //            << R_V_f * Quaternions::zHat * R_V_f.conjugate() << "\n" << std::endl;

  // Now rotate Instant so that its z axis is aligned with V_f
  Instant.RotateDecompositionBasis(R_V_f);

  // Get the phase of the (2,+/-2) modes after rotation
  const unsigned int i_22 = Instant.FindModeIndex(2,2);
  const unsigned int i_2m2 = Instant.FindModeIndex(2,-2);
  const double phase_22 = std::atan2(Instant.Im(i_22,0),Instant.Re(i_22,0));
  const double phase_2m2 = std::atan2(Instant.Im(i_2m2,0),Instant.Re(i_2m2,0));

  // R_eps is the rotation we will be applying on the right-hand side
  R_eps = R_V_f * Quaternions::exp(Quaternions::Quaternion(0,0,0,(-(phase_22-phase_2m2)/8.)));

  // Without changing anything else (the direction of V_f or the
  // phase), make sure that the rotating frame's XHat axis is more
  // parallel to the input nHat_t_fid than anti-parallel.
  if(nHat_t_fid.dot(R_f0*R_eps*Quaternions::xHat*R_eps.inverse()*R_f0.inverse()) < 0) {
    R_eps = R_eps * Quaternions::exp((M_PI/2.)*Quaternions::zHat);
    // INFOTOCERR << "Rotating by pi about the z axis initially.\n"
    //            << nHat_t_fid << "\n"
    //            << R_f0*R_eps*Quaternions::xHat*R_eps.inverse()*R_f0.inverse() << "\n"
    //            << R_f0 << "\n"
    //            << R_eps << "\n" << std::endl;
  }

  return R_eps;
}

/// Fix the attitude of the corotating frame.
GWFrames::Waveform& GWFrames::Waveform::AlignDecompositionFrameToModes(const double t_fid, const Quaternions::Quaternion& nHat_t_fid, const std::vector<int>& Lmodes) {
  ///
  /// \param t_fid Fiducial time at which the alignment should happen
  /// \param nHat_t_fid The approximate direction of nHat at t_fid
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
  const Quaternion R_eps = GetAlignmentOfDecompositionFrameToModes(t_fid, nHat_t_fid, Lmodes);

  // Record what happened
  history << "this->AlignDecompositionFrameToModes(" << std::setprecision(16) << t_fid << ", " << nHat_t_fid << ", " << Lmodes << ");  # R_eps=" << R_eps << std::endl;

  // Now, apply the rotation
  this->RotateDecompositionBasis(R_eps);

  return *this;
}

#ifndef DOXYGEN
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
// This is a local object used by `AlignWaveforms`
class WaveformAligner {
public:
  std::vector<Quaternions::Quaternion> R_fA;
  std::vector<double> t_A;
  const GWFrames::Waveform& W_B;
  const double t_mid;
  std::vector<Quaternion> R_epsB;
  bool R_epsB_is_set;
  mutable unsigned int nHat_B_i; // Just a guess to speed up hunting for the index
  mutable unsigned int Rbar_epsB_i; // Just a guess to speed up hunting for the index
  mutable bool Flip;
  const bool Debug;
public:
  WaveformAligner(const GWFrames::Waveform& W_A, const GWFrames::Waveform& iW_B,
                  const double t_1, const double t_2, const bool iDebug)
    : R_fA(W_A.Frame()), t_A(W_A.T()), W_B(iW_B), t_mid((t_1+t_2)/2.),
      R_epsB(0), R_epsB_is_set(false), nHat_B_i(0), Rbar_epsB_i(0), Flip(false), Debug(iDebug)
  {
    // Check to make sure we have sufficient times before any offset.
    // (This is necessary but not sufficient for the method to work.)
    if(t_1<t_A[0]) {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Alignment time t_1=" << t_1
                << " does not occur in t_A (which has t_A[0]=" << t_A[0] << ")." << std::endl;
      throw(GWFrames_IndexOutOfBounds);
    }
    if(t_2>t_A.back()) {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Alignment time t_2=" << t_2
                << " does not occur in t_A (which has t_A.back()=" << t_A.back() << ")." << std::endl;
      throw(GWFrames_IndexOutOfBounds);
    }
    if(W_B.NTimes()>0) {
      if(t_1<W_B.T(0)) {
        std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Alignment time t_1=" << t_1
                  << " does not occur in W_B (which has W_B.T(0)=" << W_B.T(0) << ")." << std::endl;
        throw(GWFrames_IndexOutOfBounds);
      }
      if(t_2>W_B.T(W_B.NTimes()-1)) {
        std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Alignment time t_2=" << t_2
                  << " does not occur in W_B (which has W_B.T(-1)=" << W_B.T(W_B.NTimes()-1) << ")." << std::endl;
        throw(GWFrames_IndexOutOfBounds);
      }
    }
    // Trim the fixed frame (R_fA) and its set of times, to which we
    // will interpolate.
    unsigned int i=t_A.size()-1;
    while(t_A[i]>t_2 && i>0) { --i; }
    t_A.erase(t_A.begin()+i, t_A.end());
    R_fA.erase(R_fA.begin()+i, R_fA.end());
    i=0;
    while(i<t_A.size() && t_A[i]<t_1) { ++i; }
    t_A.erase(t_A.begin(), t_A.begin()+i);
    R_fA.erase(R_fA.begin(), R_fA.begin()+i);

    if(Debug) {
      INFOTOCERR << "\tOutput to /tmp/XiIntegral.dat" << std::endl;
      ofstream myfile;
      myfile.open ("/tmp/XiIntegralPoints.dat");
      myfile.close();
    }
  }

  void SetR_epsB(const std::vector<Quaternion>& iR_epsB) {
    if(iR_epsB.size()!=W_B.NTimes()) {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": R_epsB is assumed to have the same size as W_B.T()"
                << "\n    R_epsB.size()=" << R_epsB.size() << "\tW_B.NTimes()=" << W_B.NTimes() << std::endl;
      throw(GWFrames_VectorSizeMismatch);
    }
    R_epsB = iR_epsB;
    R_epsB_is_set = true;
    return;
  }

  std::vector<Quaternion> Rbar_fB(const std::vector<double>& t) const {
    return Quaternions::conjugate(Quaternions::Squad(W_B.Frame(), W_B.T(), t));
  }

  Quaternion Rbar_epsB(const double t) const {
    if(!R_epsB_is_set) {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": R_epsB has not yet been set." << std::endl;
      throw(GWFrames_ValueError);
    }
    Rbar_epsB_i = Quaternions::hunt(W_B.T(), t, Rbar_epsB_i);
    return Quaternions::conjugate(R_epsB[Rbar_epsB_i]);
  }

  double EvaluateMinimizationQuantity(const double deltat, const double deltax, const double deltay, const double deltaz) const {
    using namespace Quaternions; // Allow me to subtract a double from a vector<double> below
    const Quaternions::Quaternion R_eps = W_B.GetAlignmentOfDecompositionFrameToModes(t_mid+deltat, Quaternions::xHat);
    const Quaternions::Quaternion R_delta = Quaternions::exp(Quaternions::Quaternion(0, deltax, deltay, deltaz));
    const std::vector<Quaternions::Quaternion> R_Bprime = Quaternions::Squad(R_delta * W_B.Frame() * R_eps, W_B.T(), t_A+deltat);
    const unsigned int Size=R_Bprime.size();
    double f1 = 0.0;
    double f2 = 0.0;
    double fdot_last1 = 4 * Quaternions::normsquared( Quaternions::logRotor( R_fA[0] * Quaternions::inverse(R_Bprime[0]) ) );
    double fdot_last2 = 4 * Quaternions::normsquared( Quaternions::logRotor( R_fA[0] * Quaternions::inverse(R_Bprime[0]*Quaternions::zHat) ) );
    for(unsigned int i=1; i<Size; ++i) {
      const double fdot1 = 4 * Quaternions::normsquared( Quaternions::logRotor( R_fA[i] * Quaternions::inverse(R_Bprime[i]) ) );
      const double fdot2 = 4 * Quaternions::normsquared( Quaternions::logRotor( R_fA[i] * Quaternions::inverse(R_Bprime[i]*Quaternions::zHat) ) );
      f1 += (t_A[i]-t_A[i-1])*(fdot1+fdot_last1)/2.0;
      f2 += (t_A[i]-t_A[i-1])*(fdot2+fdot_last2)/2.0;
      fdot_last1 = fdot1;
      fdot_last2 = fdot2;
    }
    if(Debug) {
      // INFOTOCERR << "\tOutput to /tmp/XiIntegralPoints.dat" << std::endl;
      ofstream myfile;
      myfile.open("/tmp/XiIntegralPoints.dat", std::ofstream::app);
      myfile << std::setprecision(15);
      myfile << deltat << " " << f1 << " " << f2 << std::endl;
      myfile.close();
    }
    Flip = (f2<f1);
    return std::min(f1,f2);
  }
};
double minfunc (const gsl_vector* delta, void* params) {
  WaveformAligner* Aligner = (WaveformAligner*) params;
  return Aligner->EvaluateMinimizationQuantity(gsl_vector_get(delta,0),
                                               gsl_vector_get(delta,1),
                                               gsl_vector_get(delta,2),
                                               gsl_vector_get(delta,3));
}
#endif // DOXYGEN

/// Do everything necessary to align two waveform objects
void GWFrames::AlignWaveforms(GWFrames::Waveform& W_A, GWFrames::Waveform& W_B,
                              const double t_1, const double t_2, unsigned int InitialEvaluations, std::vector<double> nHat_A, const bool Debug)
{
  /// \param W_A Fixed waveform (though modes are re-aligned)
  /// \param W_B Adjusted waveform (modes are re-aligned and frame and time are offset)
  /// \param t_1 Beginning of alignment interval
  /// \param t_2 End of alignment interval
  /// \param InitialEvaluations Number of evaluations for dumb initial optimization
  /// \param nHat_A Approximate nHat vector at (t_1+t_2)/2. [optional]
  ///
  /// This function aligns the frame to the waveform modes for both
  /// input Waveform objects at time t_mid = (t_1+t_2)/2.  It also
  /// optimizes the alignment of `W_B` by adjusting its time and
  /// overall attitudes to align with `W_A` as well as possible.
  /// While doing so, it re-adjusts the frame alignment to the modes
  /// for `W_B` to account for the changing meaning of t_mid.
  ///
  /// Note that `t_1` and `t_2` refer to fixed times with respect to
  /// the time axis of `W_A`.
  ///
  /// The input waveforms are transformed to their co-rotating frames
  /// if they are in the inertial frame.  Otherwise, they must already
  /// be in the co-rotating frame.  (E.g., the co-orbital frame is an
  /// error.)
  ///
  /// The `nHat` quantity is just the approximate direction for that
  /// vector (pointing from black hole A to black hole B) in the
  /// systems, used to set the direction of the x axis for the
  /// rotating frame.  Only the value at t_mid for `W_A` is needed.
  ///
  /// The alignment algorithm assumes that the waveforms are already
  /// reasonably well aligned in time.  In particular, the final value
  /// of t_mid+deltat for `W_B` must lie somewhere in the interval
  /// (t_1, t_2) at least, and after the time shift, `W_B` must have
  /// data over all of that interval.
  ///
  /// As long as this last condition is satisfied, and the waveforms
  /// are even remotely well sampled, and assuming the code's logic
  /// does not fail due to otherwise weird or unexpected inputs, I
  /// (Mike Boyle) do hereby guarantee that this algorithm will find
  /// the optimal alignment in both time and attitude.  Or your money
  /// back.

  if(nHat_A.size()==0) {
    nHat_A = Quaternions::xHat.vec();
  }

  // Make sure W_A and W_B are in their co-rotating frames
  if(W_A.FrameType()==GWFrames::Inertial) {
    INFOTOCOUT << "\tTransforming input Waveform A (in place) to co-rotating frame." << std::endl;
    W_A.TransformToCorotatingFrame();
  }
  if(W_B.FrameType()==GWFrames::Inertial) {
    INFOTOCOUT << "\tTransforming input Waveform B (in place) to co-rotating frame." << std::endl;
    W_B.TransformToCorotatingFrame();
  }
  if(W_A.FrameType()!=GWFrames::Corotating || W_B.FrameType()!=GWFrames::Corotating) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ":"
              << "\nError: `AlignWaveforms` takes Waveforms in the "
              << GWFrames::WaveformFrameNames[GWFrames::Inertial] << " or "
              << GWFrames::WaveformFrameNames[GWFrames::Corotating] << " frames only."
              << "\n       These Waveforms are in the " << W_A.FrameTypeString() << " and " << W_B.FrameTypeString() << " frames." << std::endl;
    throw(GWFrames_WrongFrameType);
  }

  // Make sure the various times fit together
  if(t_1>=t_2 || t_1<W_A.T(0) || t_2>W_A.T(W_A.NTimes()-1) || t_1<W_B.T(0) || t_2>W_B.T(W_B.NTimes()-1)) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ":"
              << "\nError: Incompatible input times:"
              << "\n       t_1 = " << t_1
              << "\n       t_2 = " << t_2
              << "\n       W_A.T(0) = " << W_A.T(0)
              << "\n       W_A.T(" << W_A.NTimes()-1 << ") = " << W_A.T(W_A.NTimes()-1)
              << "\n       W_B.T(0) = " << W_B.T(0)
              << "\n       W_B.T(" << W_B.NTimes()-1 << ") = " << W_B.T(W_B.NTimes()-1)
              << std::endl;
    throw(GWFrames_EmptyIntersection);
  }

  Quaternion R_delta;
  const double t_mid = (t_1+t_2)/2.;

  // We have two time offsets: deltat_1 being the most negative
  // number; deltat_2 being the most positive number.  These are the
  // offsets given to the time window on which we will evaluate W_B.
  // That is, W_B is left in place, but we evaluate it on a shifting
  // window of time (equivalent to shifting W_B by the opposite
  // amount); W_A is always left in place, and we only evaluate it on
  // the original window.  Thus, we set the bounds to ensure that
  // W_B.T(0)<t_1+deltat_1 and W_B.T(-1)>t_2+deltat_2 -- which
  // translate into deltat_1>W_B.T(0)-t_1 and deltat_2<W_B.T(-1)-t_2.
  // Just for good measure, let's move those W_B.T indices in one.
  // Also, we don't search more than (t2-t1)/2.0 to either left or
  // right.
  const double deltat_1 = std::max(W_B.T(1)-t_1, -(t_2-t_1)/2.);
  const double deltat_2 = std::min(W_B.NTimes()-2-t_2, (t_2-t_1)/2.);


  // Align W_A forever, and align W_B initially as a first guess
  W_A.AlignDecompositionFrameToModes(t_mid, nHat_A);
  W_B.AlignDecompositionFrameToModes(t_mid, Quaternions::xHat);

  WaveformAligner Aligner(W_A, W_B, t_1, t_2, Debug);
  const std::vector<double>& t_A = Aligner.t_A;
  const std::vector<Quaternions::Quaternion>& R_fA = Aligner.R_fA;

  struct timeval now;
  gettimeofday(&now, NULL); unsigned long long tNow = now.tv_usec + (unsigned long long)now.tv_sec * 1000000;

  // First, minimize the dumb way, by just evaluating at every deltat
  // in W_B so that we don't have to interpolate (which takes a *lot*
  // of time).  This should get us a very good estimate of the true
  // minimum.
  {
    // R_epsB is an array of R_eps rotors for waveform B, assuming the
    // various deltat values
    Aligner.SetR_epsB(W_B.GetAlignmentsOfDecompositionFrameToModes());

    // Evaluate Xi_c for every deltat that won't require interpolating
    // W_B to find R_eps_B (because interpolation is really slow)
    const GWFrames::Waveform W_B_Interval = W_B.SliceOfTimesWithoutModes(t_mid+deltat_1, t_mid+deltat_2);
    using namespace GWFrames; // To subtract double from vector<double> below
    vector<double> deltats = W_B_Interval.T()-t_mid;
    if(InitialEvaluations>0 && InitialEvaluations<deltats.size()) { // make sure deltats is small enough
      vector<double> deltats_tmp;
      const unsigned int step = deltats.size()/InitialEvaluations + 1;
      deltats_tmp.reserve(InitialEvaluations);
      for(unsigned int i=0; i<deltats.size(); i+=step) {
        deltats_tmp.push_back(deltats[i]);
      }
      deltats_tmp.swap(deltats);
    }
    vector<Quaternion> XiIntegral1(deltats.size());
    vector<Quaternion> XiIntegral2(deltats.size());
    for(unsigned int i=0; i<XiIntegral1.size(); ++i) {
      XiIntegral1[i] = Quaternions::DefiniteIntegral(R_fA*Aligner.Rbar_epsB(t_mid+deltats[i])*Aligner.Rbar_fB(t_A+deltats[i]), t_A);
      XiIntegral2[i] = Quaternions::DefiniteIntegral(R_fA*(-Quaternions::zHat)*Aligner.Rbar_epsB(t_mid+deltats[i])*Aligner.Rbar_fB(t_A+deltats[i]), t_A);
    }

    if(Debug) {
      INFOTOCERR << "\tOutput to /tmp/XiIntegral.dat" << std::endl;
      ofstream myfile;
      myfile.open ("/tmp/XiIntegral.dat");
      myfile << std::setprecision(15);
      for(unsigned int i=0; i<XiIntegral1.size(); ++i) {
        myfile << deltats[i] << " "
               << 2*(t_2 - t_1 - Quaternions::abs(XiIntegral1[i])) << " "
               << 2*(t_2 - t_1 - Quaternions::abs(XiIntegral2[i])) << std::endl;
      }
      myfile.close();
    }

    // Find the best value
    double Xi_c_min = 1e300;
    unsigned int i_Xi_c_min = 0;
    bool Flip = false;
    for(unsigned int i=0; i<XiIntegral1.size(); ++i) {
      // Note that, presumably because of the `abs` here, the sign of
      // the argmin will be ambiguous.  This needs to be detected
      // below, when we set the initial conditions for the minimizer.
      const double Xi_c_i1 = 2*(t_2 - t_1 - Quaternions::abs(XiIntegral1[i]));
      const double Xi_c_i2 = 2*(t_2 - t_1 - Quaternions::abs(XiIntegral2[i]));
      if(Xi_c_i1<Xi_c_min) {
        Xi_c_min = Xi_c_i1;
        i_Xi_c_min = i;
        Flip = false;
      }
      if(Xi_c_i2<Xi_c_min) {
        Xi_c_min = Xi_c_i2;
        i_Xi_c_min = i;
        Flip = true;
      }
    }
    const double deltat = deltats[i_Xi_c_min];
    // W_B.RotateDecompositionBasis(Aligner.Rbar_epsB(t_mid+deltat).conjugate() * (Flip ? Quaternions::zHat : Quaternions::One));
    W_B.SetTime(W_B.T()-deltat);
    R_delta = (Flip ? XiIntegral2[i_Xi_c_min].normalized() : XiIntegral1[i_Xi_c_min].normalized());

    INFOTOCOUT << "Objective function=" << Xi_c_min << " at " << deltat << " with" << (Flip ? " " : " no ") << "flip." << std::endl;
  }

  gettimeofday(&now, NULL); unsigned long long tThen = now.tv_usec + (unsigned long long)now.tv_sec * 1000000;
  INFOTOCOUT << "\tFirst stage took " << (tThen-tNow)/1000000.0L << " seconds." << std::endl;

  // Next, minimize algorithmically, in four dimensions, accounting
  // for all adjustments in generality.  This is very slow, but we've
  // gotten a very good initial guess from the dumb way above.
  {
    const unsigned int NDimensions = 4;
    const unsigned int MaxIterations = 2000;
    const double MinSimplexSize = 2.0e-9; // This can be less than sqrt(machine precision) because of the integral nature of our objective function
    // const double MinSimplexSize = 2.0e-13; // This can be less than sqrt(machine precision) because of the integral nature of our objective function
    double deltat=0.0;

    const double InitialTrialTimeStep = std::max(W_A.T(1)-W_A.T(0), W_B.T(1)-W_B.T(0))/2.;
    const double InitialTrialAngleStep = 1.0/(t_2-t_1);

    // Use Nelder-Mead simplex minimization
    const gsl_multimin_fminimizer_type* T =
      gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer* s = NULL;
    gsl_vector* ss;
    gsl_vector* x;
    gsl_multimin_function min_func;
    size_t iter = 0;
    int status = GSL_CONTINUE;
    double size = 0.0;

    // Set initial values
    x = gsl_vector_alloc(NDimensions);
    const Quaternions::Quaternion R_delta_log = Quaternions::logRotor(R_delta);
    const Quaternions::Quaternion NegativeR_delta_log = Quaternions::logRotor(-R_delta);
    gsl_vector_set(x, 0, deltat);
    if(Aligner.EvaluateMinimizationQuantity(deltat, R_delta_log[1], R_delta_log[2], R_delta_log[3])
       <= Aligner.EvaluateMinimizationQuantity(deltat, NegativeR_delta_log[1], NegativeR_delta_log[2], NegativeR_delta_log[3])) {
      gsl_vector_set(x, 1, R_delta_log[1]);
      gsl_vector_set(x, 2, R_delta_log[2]);
      gsl_vector_set(x, 3, R_delta_log[3]);
    } else {
      gsl_vector_set(x, 1, NegativeR_delta_log[1]);
      gsl_vector_set(x, 2, NegativeR_delta_log[2]);
      gsl_vector_set(x, 3, NegativeR_delta_log[3]);
    }

    // Set initial step sizes
    ss = gsl_vector_alloc(NDimensions);
    gsl_vector_set(ss, 0, InitialTrialTimeStep);
    gsl_vector_set(ss, 1, InitialTrialAngleStep);
    gsl_vector_set(ss, 2, InitialTrialAngleStep);
    gsl_vector_set(ss, 3, InitialTrialAngleStep);

    min_func.n = NDimensions;
    min_func.f = &minfunc;
    min_func.params = (void*) &Aligner;

    s = gsl_multimin_fminimizer_alloc(T, NDimensions);
    gsl_multimin_fminimizer_set(s, &min_func, x, ss);

    // Run the minimization
    while(status == GSL_CONTINUE && iter < MaxIterations) {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);

      if(status==GSL_EBADFUNC) {
        INFOTOCERR << ":\nThe iteration encountered a singular point where the function evaluated to Inf or NaN"
                   << "\nwhile minimizing at (" << gsl_vector_get(s->x, 0) << ", " << gsl_vector_get(s->x, 1)
                   << ", " << gsl_vector_get(s->x, 2) << ", " << gsl_vector_get(s->x, 3) << ")." << std::endl;
      }

      if(status==GSL_FAILURE) {
        INFOTOCERR << ":\nThe algorithm could not improve the current best approximation or bounding interval." << std::endl;
      }

      if(status==GSL_ENOPROG) {
        INFOTOCERR << ":\nThe minimizer is unable to improve on its current estimate, either due to"
                   << "\nnumerical difficulty or because a genuine local minimum has been reached." << std::endl;
      }

      if(status) break;
      size = gsl_multimin_fminimizer_size(s);
      status = gsl_multimin_test_size(size, MinSimplexSize);
    }

    if(iter==MaxIterations) {
      INFOTOCERR << "\nWarning: Minimization ended because it went through " << MaxIterations << " iterations."
                 << "\n         This may indicate failure.  You may want to try with a better initial guess." << std::endl;
    }

    // Get time shift and rotation
    deltat = gsl_vector_get(s->x, 0);
    R_delta = Quaternions::exp(Quaternions::Quaternion(0.0, gsl_vector_get(s->x, 1), gsl_vector_get(s->x, 2), gsl_vector_get(s->x, 3)));
    Aligner.EvaluateMinimizationQuantity(gsl_vector_get(s->x,0), gsl_vector_get(s->x,1), gsl_vector_get(s->x,2), gsl_vector_get(s->x,3));
    const bool Flip = Aligner.Flip;

    INFOTOCOUT << "Objective function=" << s->fval << " at " << deltat << std::endl;

    // Free allocated memory
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);

    // Now, apply the transformations
    W_B.AlignDecompositionFrameToModes(t_mid+deltat, (Flip ? -Quaternions::xHat : Quaternions::xHat));
    W_B.SetTime(W_B.T()-deltat);
    W_B.SetFrame(R_delta*W_B.Frame());

    gettimeofday(&now, NULL); unsigned long long tWhen = now.tv_usec + (unsigned long long)now.tv_sec * 1000000;
    INFOTOCOUT << "\tSecond stage took " << (tWhen-tThen)/1000000.0L << " seconds with " << iter << " iterations." << std::endl;
  }

  return;
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
    INFOTOCERR << "\nWarning:"
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

  // If the average frame rotor is closer to -1 than to 1, flip the sign
  if(C.frame.size()==C.NTimes()) {
    const Quaternions::Quaternion R_m = Quaternions::ApproximateMeanRotor(C.frame, C.t);
    if( Quaternions::ChordalDistance(R_m, -Quaternions::One) < Quaternions::ChordalDistance(R_m, Quaternions::One) ) {
      for(unsigned int i_t=0; i_t<C.NTimes(); ++i_t) {
        C.frame[i_t] = -C.frame[i_t];
      }
    }
  } else if(C.frame.size()==1) {
    if( Quaternions::ChordalDistance(C.frame[0], -Quaternions::One) < Quaternions::ChordalDistance(C.frame[0], Quaternions::One) ) {
      C.frame[0] = -C.frame[0];
    }
  }

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
    gsl_interp_accel_reset(accReA);
    gsl_interp_accel_reset(accImA);
    gsl_interp_accel_reset(accReB);
    gsl_interp_accel_reset(accImB);
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
    INFOTOCERR << "\nWarning:"
               << "\n       This Waveform has spin weight " << A.spinweight << "."
               << "\n       The Waveform in the argument has spin weight " << B.spinweight << "."
               << "\n       Hybridizing them probably does not make sense.\n"
               << std::endl;
  }
  if(A.boostweight != B.boostweight) {
    INFOTOCERR << "\nWarning:"
               << "\n       This Waveform has boost weight " << A.boostweight << "."
               << "\n       The Waveform in the argument has boost weight " << B.boostweight << "."
               << "\n       Hybridizing them probably does not make sense.\n"
               << std::endl;
  }
  if(A.frameType != B.frameType) {
    INFOTOCERR << "\nWarning:"
               << "\n       This Waveform is in the " << GWFrames::WaveformFrameNames[A.frameType] << " frame."
               << "\n       The Waveform in the argument is in the " << GWFrames::WaveformFrameNames[B.frameType] << " frame."
               << "\n       Hybridizing them probably does not make sense.\n"
               << std::endl;
  }
  if(A.dataType != B.dataType) {
    INFOTOCERR << "\nWarning:"
               << "\n       This Waveform has data type " << GWFrames::WaveformDataNames[A.dataType] << "."
               << "\n       The Waveform in the argument has data type " << GWFrames::WaveformDataNames[B.dataType] << "."
               << "\n       Hybridizing them probably does not make sense.\n"
               << std::endl;
  }
  if(A.rIsScaledOut != B.rIsScaledOut) {
    INFOTOCERR << "\nWarning:"
               << "\n       This Waveform claims radius is " << (A.rIsScaledOut ? "" : "not ") << "scaled out."
               << "\n       The Waveform in the argument claims radius is " << (B.rIsScaledOut ? "" : "not ") << "scaled out."
               << "\n       Hybridizing them probably does not make sense.\n"
               << std::endl;
  }
  if(A.mIsScaledOut != B.mIsScaledOut) {
    INFOTOCERR << "\nWarning:"
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
      INFOTOCERR << ": These Waveforms do not overlap on the requested range."
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
std::vector<std::complex<double> > GWFrames::Waveform::EvaluateAtPoint(const double vartheta, const double varphi, const unsigned int i_0, int i_1) const {
  ///
  /// \param vartheta Polar angle of detector
  /// \param varphi Azimuthal angle of detector
  /// \param i_0 Optional initial index to evaluate
  /// \param i_1 Optional one-past-final index to evaluate
  ///
  /// Note that the input angle parameters are measured relative to
  /// the inertial coordinate system.  If the modes are decomposed in
  /// a rotating frame, the angles will be adjusted appropriately at
  /// each time step (and the spin-weight will be compensated for
  /// automatically).
  ///
  /// Basically, that means that this function works correctly for
  /// Waveforms in rotating frames, without first rotating the
  /// Waveform into the inertial frame.  This saves significant
  /// computational cost.
  ///

  if(frameType == GWFrames::UnknownFrameType) {
    INFOTOCERR << "\nWarning: Asking for a Waveform in the " << GWFrames::WaveformFrameNames[GWFrames::UnknownFrameType] << " frame to be evaluated at a point."
               << "\n         This assumes that the Waveform::frame member data is correct...\n"
               << std::endl;
  }
  if(i_1==-1) {
    i_1 = NTimes();
  }
  if(i_0>=i_1) {
    INFOTOCERR << "\nError: Asking to EvaluateAtPoint on indices (i_0=" << i_0 << ") >= (i_1=" << i_1 << ")."
               << "\n       This is impossible; i_1 should be at least 1 more than i_0." << std::endl;
    throw(GWFrames_IndexOutOfBounds);
  }
  if(i_1>NTimes()) {
    INFOTOCERR << "\nError: Asking to EvaluateAtPoint on indices [i_0,i_1)=[" << i_0 << "," << i_1 << ") in a Waveform with " << NTimes() << " time steps." << std::endl;
    throw(GWFrames_IndexOutOfBounds);
  }

  const int NM = NModes();

  vector<complex<double> > d(i_1-i_0, complex<double>(0.,0.)); // Be sure to initialize to 0.0
  const Quaternions::Quaternion R_thetaphi(vartheta, varphi);
  SphericalFunctions::SWSH Y(SpinWeight()); // Y can be evaluated in terms of a unit quaternion

  if(frame.size()<2) {
    if(frame.size()==0) {
      Y.SetRotation(R_thetaphi);
    } else { // frame.size()==1
      Y.SetRotation(frame[0].inverse()*R_thetaphi);
    }
    for(int i_m=0; i_m<NM; ++i_m) {
      const int ell = LM(i_m)[0];
      const int m   = LM(i_m)[1];
      const complex<double> Ylm = Y(ell,m);
      for(int i_t=i_0; i_t<i_1; ++i_t) {
        d[i_t-i_0] += Data(i_m, i_t) * Ylm;
      }
    }
  } else {
    for(unsigned int i_t=i_0; i_t<i_1; ++i_t) {
      Y.SetRotation(frame[i_t].inverse()*R_thetaphi);
      for(int i_m=0; i_m<NM; ++i_m) {
        const int ell = LM(i_m)[0];
        const int m   = LM(i_m)[1];
        d[i_t-i_0] += Data(i_m, i_t) * Y(ell,m);
      }
    }
  }

  return d;
}

/// Evaluate Waveform at a particular sky location and an instant of time
std::complex<double> GWFrames::Waveform::InterpolateToPoint(const double vartheta, const double varphi, const double t_i,
                                                            gsl_interp_accel* accRe, gsl_interp_accel* accIm, gsl_spline* splineRe, gsl_spline* splineIm) const {
  ///
  /// \param vartheta Polar angle of complex detector
  /// \param varphi Azimuthal angle of complex detector
  /// \param t_i New time to interpolate to
  ///
  /// Note that the input angle parameters are measured relative to
  /// the binary's coordinate system.  In particular, this will make
  /// no sense if the frame type is something other than inertial, and
  /// will fail if the `FrameType` is neither `UnknownFrameType` nor
  /// `Inertial`.
  ///
  /// Pointers to GSL interpolation objects can be passed in, which
  /// eliminates the need to re-allocate them for each interpolation.

  bool ThisFunctionOwnsThePointers = (accRe==0);

  if(NTimes()<4) {
    INFOTOCERR << "\nError: " << NTimes() << " is not enough points to interpolate.\n"
               << std::endl;
    throw(GWFrames_BadWaveformInformation);
  }

  if(t_i<T(0)) {
    INFOTOCERR << "\nError: (t_i=" << t_i << ") is earlier than the earliest time in the data (T(0)=" << T(0) << ").\n"
               << std::endl;
    throw(GWFrames_ValueError);
  }

  if(t_i>T(NTimes()-1)) {
    INFOTOCERR << "\nError: (t_i=" << t_i << ") is later than the latest time in the data (T(" << NTimes()-1 << ")=" << T(NTimes()-1) << ").\n"
               << std::endl;
    throw(GWFrames_ValueError);
  }

  vector<double> dRe(4), dIm(dRe.size()); // N.B.: If `4` is changed, any function calling this with nonzero pointers must be changed.

  // Find a series of 4 time steps around the requested time
  unsigned int i_t_0(std::max(int(Quaternions::hunt(t, t_i))-(int(dRe.size()/2)-1), 0));
  if(i_t_0+(dRe.size()-1)>=NTimes()) {
    i_t_0 = NTimes()-dRe.size();
  }

  // Evaluate at this point at a series of times around the requested time
  const std::vector<complex<double> > val = this->EvaluateAtPoint(vartheta, varphi, i_t_0, i_t_0+dRe.size());

  // Copy that result into separate real, imaginary arrays for interpolation
  for(unsigned int i_t=0; i_t<dRe.size(); ++i_t) {
    dRe[i_t] = std::real(val[i_t]);
    dIm[i_t] = std::imag(val[i_t]);
  }

  // Now interpolate in time
  if(ThisFunctionOwnsThePointers) {
    accRe = gsl_interp_accel_alloc();
    accIm = gsl_interp_accel_alloc();
    splineRe = gsl_spline_alloc(gsl_interp_cspline, dRe.size());
    splineIm = gsl_spline_alloc(gsl_interp_cspline, dIm.size());
  }
  gsl_spline_init(splineRe, &(t)[i_t_0], &dRe[0], dRe.size());
  gsl_spline_init(splineIm, &(t)[i_t_0], &dIm[0], dIm.size());
  const complex<double> value( gsl_spline_eval(splineRe, t_i, accRe), gsl_spline_eval(splineIm, t_i, accIm) );
  if(ThisFunctionOwnsThePointers) {
    gsl_interp_accel_free(accRe);
    gsl_interp_accel_free(accIm);
    gsl_spline_free(splineRe);
    gsl_spline_free(splineIm);
  }

  return value;
}

/// Translate the waveform data by some series of spatial translations
GWFrames::Waveform GWFrames::Waveform::Translate(const std::vector<std::vector<double> >& deltax) const {
  /// \param deltax Array of 3-vectors by which to translate (function of time)
  ///
  /// The `deltax` parameter is assumed to be given relative to the
  /// inertial frame; if this Waveform has `frame` data, it will be
  /// used to rotate the `deltax` vector appropriately at each time
  /// step.
  ///
  /// The output Waveform is always in the inertial frame.  The input
  /// Waveform can be in any frame.  And though this routine is
  /// fastest if the Waveform is in the inertial frame, it would be
  /// more expensive to transform it first.  (Basically, try not to
  /// bother transforming the Waveform before calling this function.)
  ///
  /// This function is very slow because it has to interpolate to
  /// NTimes*(2*ellMax+1)^2 different completely unique points.  Since
  /// interpolation is so slow, this takes a lot of time.  The
  /// function also has to transform back from physical space to
  /// modes, which is another significant chunk of time.

  if(frameType == GWFrames::UnknownFrameType) {
    INFOTOCERR << "\nWarning: Asking to Translate a Waveform in an `" << GWFrames::WaveformFrameNames[frameType] << "` frame."
               << "\n         This assumes that the Waveform::frame member data is correct...."
               << std::endl;
  }

  if(deltax.size() != NTimes()) {
    INFOTOCERR << "\nERROR: (deltax.size()=" << deltax.size() << ") != (NTimes()=" << NTimes() << ")"
               << "\n       If you're going to translate, you need to do it at each time step.\n" << std::endl;
    throw(GWFrames_VectorSizeMismatch);
  }

  const unsigned int ntimes = NTimes();

  for(unsigned int i=0; i<ntimes; ++i) {
    if(deltax[i].size() != 3) {
      INFOTOCERR << "\nERROR: (deltax[" << i << "].size()=" << deltax[i].size() << ") != 3"
                 << "\n       Each translation should be a 3-vector.\n" << std::endl;
      throw(GWFrames_VectorSizeMismatch);
    }
  }

  // Copy infrastructure to new Waveform
  const Waveform& A = *this;
  Waveform B = A.CopyWithoutData();
  B.history << "*this = this->.Translate(...);"<< std::endl;
  B.frame = std::vector<Quaternions::Quaternion>(0);
  B.frameType = GWFrames::Inertial;
  B.lm = A.lm;
  B.t = A.t; // B.t will get reset later

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
    INFOTOCERR << "\nERROR: (iLatest=" << iLatest << ") < (iEarliest=" << iEarliest << ")"
               << "\n       The translation is too extreme, and there is not enough data for even one complete time step.\n" << std::endl;
    throw(GWFrames_EmptyIntersection);
  }

  // The new times will just be that subset of the old times for which
  // data exist in every direction after translation.
  B.t.erase(B.t.begin()+iLatest+1, B.t.end());
  B.t.erase(B.t.begin(), B.t.begin()+iEarliest);
  B.data.resize(NModes(), B.NTimes()); // Each row (first index, nn) corresponds to a mode

  // We allocate just once, for speed, and pass these to InterpolateToPoint
  gsl_interp_accel* accRe = gsl_interp_accel_alloc();
  gsl_interp_accel* accIm = gsl_interp_accel_alloc();
  gsl_spline* splineRe = gsl_spline_alloc(gsl_interp_cspline, 4);
  gsl_spline* splineIm = gsl_spline_alloc(gsl_interp_cspline, 4);

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
        Grid[i_g] = A.InterpolateToPoint(theta, phi, B.T(i_t_B)-rHat_dot_deltax, accRe, accIm, splineRe, splineIm);
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

  gsl_interp_accel_free(accRe);
  gsl_interp_accel_free(accIm);
  gsl_spline_free(splineRe);
  gsl_spline_free(splineIm);

  return B;
}

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

        // if(i_t%5000==0) {
        //   std::cerr << thetaRotated << "," << phiRotated << "; " << theta << "," << phi
        //             << "; \t" << Grid[i_g] << "," << Psi_4 << std::endl;
        // }
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


/// Apply a boost to h data, with nontrivial assumptions
GWFrames::Waveform& GWFrames::Waveform::BoostHFaked(const std::vector<std::vector<double> >& v) {
  /// This function does three things.  First, it evaluates the
  /// Waveform on what will become an equi-angular grid after
  /// transformation by the boost.  Second, at each point of that
  /// grid, it takes the appropriate combinations of the present value
  /// of h and its conjugate to give the value of h as observed in the
  /// boosted frame.  Finally, it transforms back to Fourier space
  /// using that new equi-angular grid.
  ///
  /// The input three-velocities are assumed to give the velocities of
  /// the boosted frame relative to the present frame.
  ///
  /// Note that this function simply uses the correct transformation
  /// of Psi_4, then multiplies by the appropriate power of gamma (-2)
  /// at each point.  This, of course, assumes that \f$\ddot{h} =
  /// \Psi_4\f$ in both frames.  That need not be the case, which is
  /// why "Faked" is in the name of this function.

  // Check the size of the input velocity
  if(v.size()!=NTimes()) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": (v.size()=" << v.size() << ") != (NTimes()=" << NTimes() << ")" << std::endl;
    throw(GWFrames_VectorSizeMismatch);
  }

  INFOTOCERR << "\nCAUTION!!!  This function relies on an imperfect formula."
             << "\nIt assumes that the second time derivative of h equals"
             << "\n(plus or minus) Psi_4, which need not be exactly true.\n" << std::endl;

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

        // Get the value of h in this frame at the appropriate point
        // of this frame
        const Quaternion Rp(theta, phi);
        sYlm.SetRotation(Rp);
        const complex<double> h = sYlm.Evaluate(Modes);

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
          ( (nRotated_n * mBarRotated_mBar * nRotated_n * mBarRotated_mBar
             - nRotated_mBar * mBarRotated_n * nRotated_n * mBarRotated_mBar
             - nRotated_n * mBarRotated_mBar * nRotated_mBar * mBarRotated_n
             + nRotated_mBar * mBarRotated_n * nRotated_mBar * mBarRotated_n) * h
            + (nRotated_n * mBarRotated_m * nRotated_n * mBarRotated_m
               - nRotated_m * mBarRotated_n * nRotated_n * mBarRotated_m
               - nRotated_n * mBarRotated_m * nRotated_m * mBarRotated_n
               + nRotated_m * mBarRotated_n * nRotated_m * mBarRotated_n) * std::conj(h) ) / (gamma*gamma);

        // if(i_t%5000==0) {
        //   std::cerr << thetaRotated << "," << phiRotated << "; " << theta << "," << phi
        //             << "; \t" << Grid[i_g] << "," << h << std::endl;
        // }
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

GWFrames::Waveform GWFrames::Waveform::operator+(const GWFrames::Waveform& B) const {
  const Waveform& A = *this;

  if(A.spinweight != B.spinweight) {
    INFOTOCERR << "\nError: Asking for the sum of two Waveform objects with different spin weights."
               << "\n       A.SpinWeight()=" << A.SpinWeight() << "\tB.SpinWeight()=" << B.SpinWeight()
               << std::endl;
    throw(GWFrames_BadWaveformInformation);
  }

  if(A.frameType != GWFrames::Inertial || B.frameType != GWFrames::Inertial) {
    if(A.frameType != B.frameType) {
      INFOTOCERR << "\nError: Asking for the sum of Waveforms in " << GWFrames::WaveformFrameNames[A.frameType]
                 << " and " << GWFrames::WaveformFrameNames[B.frameType] << " frames."
                 << "\n       This should only be applied to Waveforms in the same frame.\n"
                 << std::endl;
      throw(GWFrames_WrongFrameType);
    } else if(A.frame.size() != B.frame.size()) {
      INFOTOCERR << "\nError: Asking for the sum of Waveforms with " << A.frame.size() << " and " << B.frame.size() << " frame data points."
                 << "\n       This should only be applied to Waveforms in the same frame.\n"
                 << std::endl;
      throw(GWFrames_WrongFrameType);
    }
  }

  if(A.NTimes() != B.NTimes()) {
    INFOTOCERR << "\nError: Asking for the sum of two Waveform objects with different time data."
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
    INFOTOCERR << "\nError: Asking for the difference of two Waveform objects with different spin weights."
               << "\n       A.SpinWeight()=" << A.SpinWeight() << "\tB.SpinWeight()=" << B.SpinWeight()
               << std::endl;
    throw(GWFrames_BadWaveformInformation);
  }

  if(A.frameType != GWFrames::Inertial || B.frameType != GWFrames::Inertial) {
    if(A.frameType != B.frameType) {
      INFOTOCERR << "\nError: Asking for the difference of Waveforms in " << GWFrames::WaveformFrameNames[A.frameType]
                 << " and " << GWFrames::WaveformFrameNames[B.frameType] << " frames."
                 << "\n       This should only be applied to Waveforms in the same frame.\n"
                 << std::endl;
      throw(GWFrames_WrongFrameType);
    } else if(A.frame.size() != B.frame.size()) {
      INFOTOCERR << "\nError: Asking for the difference of Waveforms with " << A.frame.size() << " and " << B.frame.size() << " frame data points."
                 << "\n       This should only be applied to Waveforms in the same frame.\n"
                 << std::endl;
      throw(GWFrames_WrongFrameType);
    }
  }

  if(A.NTimes() != B.NTimes()) {
    INFOTOCERR << "\nError: Asking for the difference of two Waveform objects with different time data."
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
