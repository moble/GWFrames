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
#include "IntegrateAngularVelocity.hpp"
#include "SphericalFunctions/SWSHs.hpp"
#include "Errors.hpp"

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
    history << "### Code revision (`git rev-parse HEAD` or arXiv version) = " << CodeRevision << endl
	    << "### pwd = " << pwd << endl
	    << "### hostname = " << hostname << endl
	    << "### date = " << date // comes with a newline
	    << "### Waveform(); // empty constructor" << endl;
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
    history << "### Code revision (`git rev-parse HEAD` or arXiv version) = " << CodeRevision << endl
	    << "### pwd = " << pwd << endl
	    << "### hostname = " << hostname << endl
	    << "### date = " << date // comes with a newline
	    << "### Waveform(" << FileName << ", " << DataFormat << "); // Constructor from data file" << endl;
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
  const std::complex<double>* D = this->operator()(Mode);
  return ComplexDerivative(vector<std::complex<double> >(D,D+NTimes()), T());
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

/// Rotate the physical content of the Waveform by a constant rotor.
GWFrames::Waveform& GWFrames::Waveform::RotatePhysicalSystem(const Quaternions::Quaternion& R_phys) {
  history << "### this->RotatePhysicalSystem(" << R_phys << ");" << endl;

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

  history << "### this->RotatePhysicalSystem(R_phys); // R_phys=["
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
  history << "### this->RotateDecompositionBasis(" << R_frame << ");" << endl;

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

  history << "### this->RotateDecompositionBasis(R_frame); // R_frame=["
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

  history << "### this->RotateDecompositionBasisOfUncertainties(R_frame); // R_frame=["
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


double SchmidtEtAl_f(const gsl_vector *v, void *params) {
  double *h = (double *)params;
  const double hRe2Neg2 = h[0]; const double hIm2Neg2 = h[1];
  const double hRe2Neg1 = h[2]; const double hIm2Neg1 = h[3];
  const double hRe20 = h[4];    const double hIm20 = h[5];
  const double hRe21 = h[6];    const double hIm21 = h[7];
  const double hRe22 = h[8];    const double hIm22 = h[9];
  const double x0 = gsl_vector_get(v, 0); const double x1 = gsl_vector_get(v, 1);
  const double cbetao2 = cos(x1/2.0);
  const double sbetao2 = sin(x1/2.0);
  const double c1alpha = cos(x0);
  const double s1alpha = sin(x0);
  const double c2alpha = c1alpha*c1alpha-s1alpha*s1alpha;
  const double s2alpha = 2*s1alpha*c1alpha;
  const double c3alpha = (4*c1alpha*c1alpha-3)*c1alpha;
  const double s3alpha = s1alpha*(3-4*s1alpha*s1alpha);
  const double c4alpha = c2alpha*c2alpha-s2alpha*s2alpha;
  const double s4alpha = 2*s2alpha*c2alpha;
  return -1.*pow(cbetao2,8)*(pow(hIm22,2) + pow(hRe22,2) + pow(hRe2Neg2,2)) - 4.*pow(cbetao2,7)*(hIm22*hRe21 - 1.*hIm21*hRe22 + hIm2Neg2*hRe2Neg1 - 1.*hIm2Neg1*hRe2Neg2)*s1alpha*sbetao2 - 2.*pow(cbetao2,6)*(2.*(pow(hIm21,2) + pow(hRe21,2) + pow(hRe2Neg1,2)) + 2.4494897427831781*(c2alpha*(hIm20*(hIm22 + hIm2Neg2) + hRe20*(hRe22 + hRe2Neg2)) + ((hIm22 - 1.*hIm2Neg2)*hRe20 + hIm20*(-1.*hRe22 + hRe2Neg2))*s2alpha))*pow(sbetao2,2) - 4.*pow(cbetao2,5)*(c3alpha*(hIm22*hIm2Neg1 - 1.*hIm21*hIm2Neg2 + hRe22*hRe2Neg1 - 1.*hRe21*hRe2Neg2) + 2.4494897427831781*((hIm21 + hIm2Neg1)*hRe20 - 1.*hIm20*(hRe21 + hRe2Neg1))*s1alpha + (hIm2Neg2*hRe21 - 1.*hIm2Neg1*hRe22 + hIm22*hRe2Neg1 - 1.*hIm21*hRe2Neg2)*s3alpha)*pow(sbetao2,3) - 4.*pow(cbetao2,4)*(3.*(pow(hIm20,2) + pow(hRe20,2)) + 4.*c2alpha*(hIm21*hIm2Neg1 + hRe21*hRe2Neg1) + c4alpha*(hIm22*hIm2Neg2 + hRe22*hRe2Neg2) + 4.*(-1.*hIm2Neg1*hRe21 + hIm21*hRe2Neg1)*s2alpha + (-1.*hIm2Neg2*hRe22 + hIm22*hRe2Neg2)*s4alpha)*pow(sbetao2,4) - 4.*pow(cbetao2,3)*(c3alpha*(-1.*hIm22*hIm2Neg1 + hIm21*hIm2Neg2 - 1.*hRe22*hRe2Neg1 + hRe21*hRe2Neg2) + 2.4494897427831781*(-1.*(hIm21 + hIm2Neg1)*hRe20 + hIm20*(hRe21 + hRe2Neg1))*s1alpha + (-1.*hIm2Neg2*hRe21 + hIm2Neg1*hRe22 - 1.*hIm22*hRe2Neg1 + hIm21*hRe2Neg2)*s3alpha)*pow(sbetao2,5) - 2.*pow(cbetao2,2)*(2.*(pow(hIm21,2) + pow(hRe21,2) + pow(hRe2Neg1,2)) + 2.4494897427831781*(c2alpha*(hIm20*(hIm22 + hIm2Neg2) + hRe20*(hRe22 + hRe2Neg2)) + ((hIm22 - 1.*hIm2Neg2)*hRe20 + hIm20*(-1.*hRe22 + hRe2Neg2))*s2alpha))*pow(sbetao2,6) - 4.*cbetao2*(-1.*hIm22*hRe21 + hIm21*hRe22 - 1.*hIm2Neg2*hRe2Neg1 + hIm2Neg1*hRe2Neg2)*s1alpha*pow(sbetao2,7) - 1.*(pow(hIm22,2) + pow(hRe22,2) + pow(hRe2Neg2,2))*pow(sbetao2,8) - 4.*pow(cbetao2,2)*pow(hIm2Neg1,2)*pow(sbetao2,2)*(pow(cbetao2,4) + pow(sbetao2,4)) + 4.*c1alpha*cbetao2*(cbetao2 - 1.*sbetao2)*sbetao2*(cbetao2 + sbetao2)*(pow(cbetao2,4)*(-1.*hIm21*hIm22 + hIm2Neg1*hIm2Neg2 - 1.*hRe21*hRe22 + hRe2Neg1*hRe2Neg2) + pow(cbetao2,2)*(-1.*hIm21*(2.4494897427831781*hIm20 + hIm22) + hIm2Neg1*(2.4494897427831781*hIm20 + hIm2Neg2) - 1.*hRe21*(2.4494897427831781*hRe20 + hRe22) + hRe2Neg1*(2.4494897427831781*hRe20 + hRe2Neg2))*pow(sbetao2,2) + (-1.*hIm21*hIm22 + hIm2Neg1*hIm2Neg2 - 1.*hRe21*hRe22 + hRe2Neg1*hRe2Neg2)*pow(sbetao2,4)) - 1.*pow(hIm2Neg2,2)*(pow(cbetao2,8) + pow(sbetao2,8));
}
void SchmidtEtAl_df(const gsl_vector *v, void *params, gsl_vector *df) {
  double *h = (double *)params;
  const double hRe2Neg2 = h[0]; const double hIm2Neg2 = h[1];
  const double hRe2Neg1 = h[2]; const double hIm2Neg1 = h[3];
  const double hRe20 = h[4];    const double hIm20 = h[5];
  const double hRe21 = h[6];    const double hIm21 = h[7];
  const double hRe22 = h[8];    const double hIm22 = h[9];
  const double x0 = gsl_vector_get(v, 0); const double x1 = gsl_vector_get(v, 1);
  const double cbetao2 = cos(x1/2.0);
  const double sbetao2 = sin(x1/2.0);
  const double c1alpha = cos(x0);
  const double s1alpha = sin(x0);
  const double c2alpha = c1alpha*c1alpha-s1alpha*s1alpha;
  const double s2alpha = 2*s1alpha*c1alpha;
  const double c3alpha = (4*c1alpha*c1alpha-3)*c1alpha;
  const double s3alpha = s1alpha*(3-4*s1alpha*s1alpha);
  const double c4alpha = c2alpha*c2alpha-s2alpha*s2alpha;
  const double s4alpha = 2*s2alpha*c2alpha;
  gsl_vector_set(df, 0, -4.*cbetao2*sbetao2*(4.*c4alpha*pow(cbetao2,3)*(-1.*hIm2Neg2*hRe22 + hIm22*hRe2Neg2)*pow(sbetao2,3) - 4.*pow(cbetao2,3)*(hIm22*hIm2Neg2 + hRe22*hRe2Neg2)*s4alpha*pow(sbetao2,3) + 3.*c3alpha*pow(cbetao2,2)*(hIm2Neg2*hRe21 - 1.*hIm2Neg1*hRe22 + hIm22*hRe2Neg1 - 1.*hIm21*hRe2Neg2)*(cbetao2 - 1.*sbetao2)*pow(sbetao2,2)*(cbetao2 + sbetao2) + 3.*pow(cbetao2,2)*(-1.*hIm22*hIm2Neg1 + hIm21*hIm2Neg2 - 1.*hRe22*hRe2Neg1 + hRe21*hRe2Neg2)*s3alpha*(cbetao2 - 1.*sbetao2)*pow(sbetao2,2)*(cbetao2 + sbetao2) - 1.*c1alpha*(cbetao2 - 1.*sbetao2)*(cbetao2 + sbetao2)*(pow(cbetao2,4)*(-1.*hIm22*hRe21 + hIm21*hRe22 - 1.*hIm2Neg2*hRe2Neg1 + hIm2Neg1*hRe2Neg2) + pow(cbetao2,2)*(-1.*hIm22*hRe21 + hIm21*(-2.4494897427831781*hRe20 + hRe22) - 1.*hIm2Neg2*hRe2Neg1 + 2.4494897427831781*hIm20*(hRe21 + hRe2Neg1) + hIm2Neg1*(-2.4494897427831781*hRe20 + hRe2Neg2))*pow(sbetao2,2) + (-1.*hIm22*hRe21 + hIm21*hRe22 - 1.*hIm2Neg2*hRe2Neg1 + hIm2Neg1*hRe2Neg2)*pow(sbetao2,4)) + s1alpha*(cbetao2 - 1.*sbetao2)*(cbetao2 + sbetao2)*(pow(cbetao2,4)*(-1.*hIm21*hIm22 + hIm2Neg1*hIm2Neg2 - 1.*hRe21*hRe22 + hRe2Neg1*hRe2Neg2) + pow(cbetao2,2)*(-1.*hIm21*(2.4494897427831781*hIm20 + hIm22) + hIm2Neg1*(2.4494897427831781*hIm20 + hIm2Neg2) - 1.*hRe21*(2.4494897427831781*hRe20 + hRe22) + hRe2Neg1*(2.4494897427831781*hRe20 + hRe2Neg2))*pow(sbetao2,2) + (-1.*hIm21*hIm22 + hIm2Neg1*hIm2Neg2 - 1.*hRe21*hRe22 + hRe2Neg1*hRe2Neg2)*pow(sbetao2,4)) + c2alpha*cbetao2*sbetao2*(2.4494897427831781*pow(cbetao2,4)*((hIm22 - 1.*hIm2Neg2)*hRe20 + hIm20*(-1.*hRe22 + hRe2Neg2)) + 8.*pow(cbetao2,2)*(-1.*hIm2Neg1*hRe21 + hIm21*hRe2Neg1)*pow(sbetao2,2) + 2.4494897427831781*((hIm22 - 1.*hIm2Neg2)*hRe20 + hIm20*(-1.*hRe22 + hRe2Neg2))*pow(sbetao2,4)) + cbetao2*s2alpha*sbetao2*(-2.4494897427831781*pow(cbetao2,4)*(hIm20*(hIm22 + hIm2Neg2) + hRe20*(hRe22 + hRe2Neg2)) - 8.*pow(cbetao2,2)*(hIm21*hIm2Neg1 + hRe21*hRe2Neg1)*pow(sbetao2,2) - 2.4494897427831781*(hIm20*(hIm22 + hIm2Neg2) + hRe20*(hRe22 + hRe2Neg2))*pow(sbetao2,4))));
  gsl_vector_set(df, 1, -2.*(4.*c4alpha*pow(cbetao2,3)*(hIm22*hIm2Neg2 + hRe22*hRe2Neg2)*(cbetao2 - 1.*sbetao2)*pow(sbetao2,3)*(cbetao2 + sbetao2) + 4.*pow(cbetao2,3)*(-1.*hIm2Neg2*hRe22 + hIm22*hRe2Neg2)*s4alpha*(cbetao2 - 1.*sbetao2)*pow(sbetao2,3)*(cbetao2 + sbetao2) - 1.*c3alpha*pow(cbetao2,2)*(-1.*hIm22*hIm2Neg1 + hIm21*hIm2Neg2 - 1.*hRe22*hRe2Neg1 + hRe21*hRe2Neg2)*pow(sbetao2,2)*(3.*pow(cbetao2,4) - 10.*pow(cbetao2,2)*pow(sbetao2,2) + 3.*pow(sbetao2,4)) + pow(cbetao2,2)*(hIm2Neg2*hRe21 - 1.*hIm2Neg1*hRe22 + hIm22*hRe2Neg1 - 1.*hIm21*hRe2Neg2)*s3alpha*pow(sbetao2,2)*(3.*pow(cbetao2,4) - 10.*pow(cbetao2,2)*pow(sbetao2,2) + 3.*pow(sbetao2,4)) + 2.*sbetao2*(-1.*pow(cbetao2,3) + cbetao2*pow(sbetao2,2))*(pow(cbetao2,4)*(-1.*pow(hIm21,2) + pow(hIm22,2) - 1.*pow(hIm2Neg1,2) + pow(hIm2Neg2,2) - 1.*pow(hRe21,2) + pow(hRe22,2) - 1.*pow(hRe2Neg1,2) + pow(hRe2Neg2,2)) + pow(cbetao2,2)*(-6.*pow(hIm20,2) + 2.*pow(hIm21,2) + pow(hIm22,2) + 2.*pow(hIm2Neg1,2) + pow(hIm2Neg2,2) + pow(hRe22,2) + 2.*(-3.*pow(hRe20,2) + pow(hRe21,2) + pow(hRe2Neg1,2)) + pow(hRe2Neg2,2))*pow(sbetao2,2) + (-1.*pow(hIm21,2) + pow(hIm22,2) - 1.*pow(hIm2Neg1,2) + pow(hIm2Neg2,2) - 1.*pow(hRe21,2) + pow(hRe22,2) - 1.*pow(hRe2Neg1,2) + pow(hRe2Neg2,2))*pow(sbetao2,4)) + cbetao2*s2alpha*(cbetao2 - 1.*sbetao2)*sbetao2*(cbetao2 + sbetao2)*(2.4494897427831781*pow(cbetao2,4)*((hIm22 - 1.*hIm2Neg2)*hRe20 + hIm20*(-1.*hRe22 + hRe2Neg2)) - 2.*pow(cbetao2,2)*(2.4494897427831781*(hIm22 - 1.*hIm2Neg2)*hRe20 + 8.*hIm2Neg1*hRe21 - 8.*hIm21*hRe2Neg1 + 2.4494897427831781*hIm20*(-1.*hRe22 + hRe2Neg2))*pow(sbetao2,2) + 2.4494897427831781*((hIm22 - 1.*hIm2Neg2)*hRe20 + hIm20*(-1.*hRe22 + hRe2Neg2))*pow(sbetao2,4)) + c2alpha*cbetao2*(cbetao2 - 1.*sbetao2)*sbetao2*(cbetao2 + sbetao2)*(2.4494897427831781*pow(cbetao2,4)*(hIm20*(hIm22 + hIm2Neg2) + hRe20*(hRe22 + hRe2Neg2)) - 2.*pow(cbetao2,2)*(-8.*hIm21*hIm2Neg1 + 2.4494897427831781*hIm20*(hIm22 + hIm2Neg2) - 8.*hRe21*hRe2Neg1 + 2.4494897427831781*hRe20*(hRe22 + hRe2Neg2))*pow(sbetao2,2) + 2.4494897427831781*(hIm20*(hIm22 + hIm2Neg2) + hRe20*(hRe22 + hRe2Neg2))*pow(sbetao2,4)) + s1alpha*(pow(cbetao2,8)*(hIm22*hRe21 - 1.*hIm21*hRe22 + hIm2Neg2*hRe2Neg1 - 1.*hIm2Neg1*hRe2Neg2) + pow(cbetao2,6)*(-7.*hIm22*hRe21 + hIm21*(7.34846922834953429*hRe20 + 7.*hRe22) - 7.*hIm2Neg2*hRe2Neg1 - 7.34846922834953429*hIm20*(hRe21 + hRe2Neg1) + hIm2Neg1*(7.34846922834953429*hRe20 + 7.*hRe2Neg2))*pow(sbetao2,2) + 24.494897427831781*pow(cbetao2,4)*(-1.*(hIm21 + hIm2Neg1)*hRe20 + hIm20*(hRe21 + hRe2Neg1))*pow(sbetao2,4) + pow(cbetao2,2)*(-7.*hIm22*hRe21 + hIm21*(7.34846922834953429*hRe20 + 7.*hRe22) - 7.*hIm2Neg2*hRe2Neg1 - 7.34846922834953429*hIm20*(hRe21 + hRe2Neg1) + hIm2Neg1*(7.34846922834953429*hRe20 + 7.*hRe2Neg2))*pow(sbetao2,6) + (hIm22*hRe21 - 1.*hIm21*hRe22 + hIm2Neg2*hRe2Neg1 - 1.*hIm2Neg1*hRe2Neg2)*pow(sbetao2,8)) + c1alpha*(pow(cbetao2,8)*(hIm21*hIm22 - 1.*hIm2Neg1*hIm2Neg2 + hRe21*hRe22 - 1.*hRe2Neg1*hRe2Neg2) + pow(cbetao2,6)*(hIm21*(7.34846922834953429*hIm20 - 7.*hIm22) + hIm2Neg1*(-7.34846922834953429*hIm20 + 7.*hIm2Neg2) + hRe21*(7.34846922834953429*hRe20 - 7.*hRe22) + hRe2Neg1*(-7.34846922834953429*hRe20 + 7.*hRe2Neg2))*pow(sbetao2,2) + 24.494897427831781*pow(cbetao2,4)*(hIm20*(-1.*hIm21 + hIm2Neg1) + hRe20*(-1.*hRe21 + hRe2Neg1))*pow(sbetao2,4) + pow(cbetao2,2)*(hIm21*(7.34846922834953429*hIm20 - 7.*hIm22) + hIm2Neg1*(-7.34846922834953429*hIm20 + 7.*hIm2Neg2) + hRe21*(7.34846922834953429*hRe20 - 7.*hRe22) + hRe2Neg1*(-7.34846922834953429*hRe20 + 7.*hRe2Neg2))*pow(sbetao2,6) + (hIm21*hIm22 - 1.*hIm2Neg1*hIm2Neg2 + hRe21*hRe22 - 1.*hRe2Neg1*hRe2Neg2)*pow(sbetao2,8))));
  return;
}
void SchmidtEtAl_fdf (const gsl_vector* v, void* params, double* f, gsl_vector* df) {
  double* h = (double*)params;
  const double hRe2Neg2 = h[0]; const double hIm2Neg2 = h[1];
  const double hRe2Neg1 = h[2]; const double hIm2Neg1 = h[3];
  const double hRe20 = h[4];    const double hIm20 = h[5];
  const double hRe21 = h[6];    const double hIm21 = h[7];
  const double hRe22 = h[8];    const double hIm22 = h[9];
  const double x0 = gsl_vector_get(v, 0); const double x1 = gsl_vector_get(v, 1);
  const double cbetao2 = cos(x1/2.0);
  const double sbetao2 = sin(x1/2.0);
  const double c1alpha = cos(x0);
  const double s1alpha = sin(x0);
  const double c2alpha = c1alpha*c1alpha-s1alpha*s1alpha;
  const double s2alpha = 2*s1alpha*c1alpha;
  const double c3alpha = (4*c1alpha*c1alpha-3)*c1alpha;
  const double s3alpha = s1alpha*(3-4*s1alpha*s1alpha);
  const double c4alpha = c2alpha*c2alpha-s2alpha*s2alpha;
  const double s4alpha = 2*s2alpha*c2alpha;
  *f = -1.*pow(cbetao2,8)*(pow(hIm22,2) + pow(hRe22,2) + pow(hRe2Neg2,2)) - 4.*pow(cbetao2,7)*(hIm22*hRe21 - 1.*hIm21*hRe22 + hIm2Neg2*hRe2Neg1 - 1.*hIm2Neg1*hRe2Neg2)*s1alpha*sbetao2 - 2.*pow(cbetao2,6)*(2.*(pow(hIm21,2) + pow(hRe21,2) + pow(hRe2Neg1,2)) + 2.4494897427831781*(c2alpha*(hIm20*(hIm22 + hIm2Neg2) + hRe20*(hRe22 + hRe2Neg2)) + ((hIm22 - 1.*hIm2Neg2)*hRe20 + hIm20*(-1.*hRe22 + hRe2Neg2))*s2alpha))*pow(sbetao2,2) - 4.*pow(cbetao2,5)*(c3alpha*(hIm22*hIm2Neg1 - 1.*hIm21*hIm2Neg2 + hRe22*hRe2Neg1 - 1.*hRe21*hRe2Neg2) + 2.4494897427831781*((hIm21 + hIm2Neg1)*hRe20 - 1.*hIm20*(hRe21 + hRe2Neg1))*s1alpha + (hIm2Neg2*hRe21 - 1.*hIm2Neg1*hRe22 + hIm22*hRe2Neg1 - 1.*hIm21*hRe2Neg2)*s3alpha)*pow(sbetao2,3) - 4.*pow(cbetao2,4)*(3.*(pow(hIm20,2) + pow(hRe20,2)) + 4.*c2alpha*(hIm21*hIm2Neg1 + hRe21*hRe2Neg1) + c4alpha*(hIm22*hIm2Neg2 + hRe22*hRe2Neg2) + 4.*(-1.*hIm2Neg1*hRe21 + hIm21*hRe2Neg1)*s2alpha + (-1.*hIm2Neg2*hRe22 + hIm22*hRe2Neg2)*s4alpha)*pow(sbetao2,4) - 4.*pow(cbetao2,3)*(c3alpha*(-1.*hIm22*hIm2Neg1 + hIm21*hIm2Neg2 - 1.*hRe22*hRe2Neg1 + hRe21*hRe2Neg2) + 2.4494897427831781*(-1.*(hIm21 + hIm2Neg1)*hRe20 + hIm20*(hRe21 + hRe2Neg1))*s1alpha + (-1.*hIm2Neg2*hRe21 + hIm2Neg1*hRe22 - 1.*hIm22*hRe2Neg1 + hIm21*hRe2Neg2)*s3alpha)*pow(sbetao2,5) - 2.*pow(cbetao2,2)*(2.*(pow(hIm21,2) + pow(hRe21,2) + pow(hRe2Neg1,2)) + 2.4494897427831781*(c2alpha*(hIm20*(hIm22 + hIm2Neg2) + hRe20*(hRe22 + hRe2Neg2)) + ((hIm22 - 1.*hIm2Neg2)*hRe20 + hIm20*(-1.*hRe22 + hRe2Neg2))*s2alpha))*pow(sbetao2,6) - 4.*cbetao2*(-1.*hIm22*hRe21 + hIm21*hRe22 - 1.*hIm2Neg2*hRe2Neg1 + hIm2Neg1*hRe2Neg2)*s1alpha*pow(sbetao2,7) - 1.*(pow(hIm22,2) + pow(hRe22,2) + pow(hRe2Neg2,2))*pow(sbetao2,8) - 4.*pow(cbetao2,2)*pow(hIm2Neg1,2)*pow(sbetao2,2)*(pow(cbetao2,4) + pow(sbetao2,4)) + 4.*c1alpha*cbetao2*(cbetao2 - 1.*sbetao2)*sbetao2*(cbetao2 + sbetao2)*(pow(cbetao2,4)*(-1.*hIm21*hIm22 + hIm2Neg1*hIm2Neg2 - 1.*hRe21*hRe22 + hRe2Neg1*hRe2Neg2) + pow(cbetao2,2)*(-1.*hIm21*(2.4494897427831781*hIm20 + hIm22) + hIm2Neg1*(2.4494897427831781*hIm20 + hIm2Neg2) - 1.*hRe21*(2.4494897427831781*hRe20 + hRe22) + hRe2Neg1*(2.4494897427831781*hRe20 + hRe2Neg2))*pow(sbetao2,2) + (-1.*hIm21*hIm22 + hIm2Neg1*hIm2Neg2 - 1.*hRe21*hRe22 + hRe2Neg1*hRe2Neg2)*pow(sbetao2,4)) - 1.*pow(hIm2Neg2,2)*(pow(cbetao2,8) + pow(sbetao2,8));
  gsl_vector_set(df, 0, -4.*cbetao2*sbetao2*(4.*c4alpha*pow(cbetao2,3)*(-1.*hIm2Neg2*hRe22 + hIm22*hRe2Neg2)*pow(sbetao2,3) - 4.*pow(cbetao2,3)*(hIm22*hIm2Neg2 + hRe22*hRe2Neg2)*s4alpha*pow(sbetao2,3) + 3.*c3alpha*pow(cbetao2,2)*(hIm2Neg2*hRe21 - 1.*hIm2Neg1*hRe22 + hIm22*hRe2Neg1 - 1.*hIm21*hRe2Neg2)*(cbetao2 - 1.*sbetao2)*pow(sbetao2,2)*(cbetao2 + sbetao2) + 3.*pow(cbetao2,2)*(-1.*hIm22*hIm2Neg1 + hIm21*hIm2Neg2 - 1.*hRe22*hRe2Neg1 + hRe21*hRe2Neg2)*s3alpha*(cbetao2 - 1.*sbetao2)*pow(sbetao2,2)*(cbetao2 + sbetao2) - 1.*c1alpha*(cbetao2 - 1.*sbetao2)*(cbetao2 + sbetao2)*(pow(cbetao2,4)*(-1.*hIm22*hRe21 + hIm21*hRe22 - 1.*hIm2Neg2*hRe2Neg1 + hIm2Neg1*hRe2Neg2) + pow(cbetao2,2)*(-1.*hIm22*hRe21 + hIm21*(-2.4494897427831781*hRe20 + hRe22) - 1.*hIm2Neg2*hRe2Neg1 + 2.4494897427831781*hIm20*(hRe21 + hRe2Neg1) + hIm2Neg1*(-2.4494897427831781*hRe20 + hRe2Neg2))*pow(sbetao2,2) + (-1.*hIm22*hRe21 + hIm21*hRe22 - 1.*hIm2Neg2*hRe2Neg1 + hIm2Neg1*hRe2Neg2)*pow(sbetao2,4)) + s1alpha*(cbetao2 - 1.*sbetao2)*(cbetao2 + sbetao2)*(pow(cbetao2,4)*(-1.*hIm21*hIm22 + hIm2Neg1*hIm2Neg2 - 1.*hRe21*hRe22 + hRe2Neg1*hRe2Neg2) + pow(cbetao2,2)*(-1.*hIm21*(2.4494897427831781*hIm20 + hIm22) + hIm2Neg1*(2.4494897427831781*hIm20 + hIm2Neg2) - 1.*hRe21*(2.4494897427831781*hRe20 + hRe22) + hRe2Neg1*(2.4494897427831781*hRe20 + hRe2Neg2))*pow(sbetao2,2) + (-1.*hIm21*hIm22 + hIm2Neg1*hIm2Neg2 - 1.*hRe21*hRe22 + hRe2Neg1*hRe2Neg2)*pow(sbetao2,4)) + c2alpha*cbetao2*sbetao2*(2.4494897427831781*pow(cbetao2,4)*((hIm22 - 1.*hIm2Neg2)*hRe20 + hIm20*(-1.*hRe22 + hRe2Neg2)) + 8.*pow(cbetao2,2)*(-1.*hIm2Neg1*hRe21 + hIm21*hRe2Neg1)*pow(sbetao2,2) + 2.4494897427831781*((hIm22 - 1.*hIm2Neg2)*hRe20 + hIm20*(-1.*hRe22 + hRe2Neg2))*pow(sbetao2,4)) + cbetao2*s2alpha*sbetao2*(-2.4494897427831781*pow(cbetao2,4)*(hIm20*(hIm22 + hIm2Neg2) + hRe20*(hRe22 + hRe2Neg2)) - 8.*pow(cbetao2,2)*(hIm21*hIm2Neg1 + hRe21*hRe2Neg1)*pow(sbetao2,2) - 2.4494897427831781*(hIm20*(hIm22 + hIm2Neg2) + hRe20*(hRe22 + hRe2Neg2))*pow(sbetao2,4))));
  gsl_vector_set(df, 1, -2.*(4.*c4alpha*pow(cbetao2,3)*(hIm22*hIm2Neg2 + hRe22*hRe2Neg2)*(cbetao2 - 1.*sbetao2)*pow(sbetao2,3)*(cbetao2 + sbetao2) + 4.*pow(cbetao2,3)*(-1.*hIm2Neg2*hRe22 + hIm22*hRe2Neg2)*s4alpha*(cbetao2 - 1.*sbetao2)*pow(sbetao2,3)*(cbetao2 + sbetao2) - 1.*c3alpha*pow(cbetao2,2)*(-1.*hIm22*hIm2Neg1 + hIm21*hIm2Neg2 - 1.*hRe22*hRe2Neg1 + hRe21*hRe2Neg2)*pow(sbetao2,2)*(3.*pow(cbetao2,4) - 10.*pow(cbetao2,2)*pow(sbetao2,2) + 3.*pow(sbetao2,4)) + pow(cbetao2,2)*(hIm2Neg2*hRe21 - 1.*hIm2Neg1*hRe22 + hIm22*hRe2Neg1 - 1.*hIm21*hRe2Neg2)*s3alpha*pow(sbetao2,2)*(3.*pow(cbetao2,4) - 10.*pow(cbetao2,2)*pow(sbetao2,2) + 3.*pow(sbetao2,4)) + 2.*sbetao2*(-1.*pow(cbetao2,3) + cbetao2*pow(sbetao2,2))*(pow(cbetao2,4)*(-1.*pow(hIm21,2) + pow(hIm22,2) - 1.*pow(hIm2Neg1,2) + pow(hIm2Neg2,2) - 1.*pow(hRe21,2) + pow(hRe22,2) - 1.*pow(hRe2Neg1,2) + pow(hRe2Neg2,2)) + pow(cbetao2,2)*(-6.*pow(hIm20,2) + 2.*pow(hIm21,2) + pow(hIm22,2) + 2.*pow(hIm2Neg1,2) + pow(hIm2Neg2,2) + pow(hRe22,2) + 2.*(-3.*pow(hRe20,2) + pow(hRe21,2) + pow(hRe2Neg1,2)) + pow(hRe2Neg2,2))*pow(sbetao2,2) + (-1.*pow(hIm21,2) + pow(hIm22,2) - 1.*pow(hIm2Neg1,2) + pow(hIm2Neg2,2) - 1.*pow(hRe21,2) + pow(hRe22,2) - 1.*pow(hRe2Neg1,2) + pow(hRe2Neg2,2))*pow(sbetao2,4)) + cbetao2*s2alpha*(cbetao2 - 1.*sbetao2)*sbetao2*(cbetao2 + sbetao2)*(2.4494897427831781*pow(cbetao2,4)*((hIm22 - 1.*hIm2Neg2)*hRe20 + hIm20*(-1.*hRe22 + hRe2Neg2)) - 2.*pow(cbetao2,2)*(2.4494897427831781*(hIm22 - 1.*hIm2Neg2)*hRe20 + 8.*hIm2Neg1*hRe21 - 8.*hIm21*hRe2Neg1 + 2.4494897427831781*hIm20*(-1.*hRe22 + hRe2Neg2))*pow(sbetao2,2) + 2.4494897427831781*((hIm22 - 1.*hIm2Neg2)*hRe20 + hIm20*(-1.*hRe22 + hRe2Neg2))*pow(sbetao2,4)) + c2alpha*cbetao2*(cbetao2 - 1.*sbetao2)*sbetao2*(cbetao2 + sbetao2)*(2.4494897427831781*pow(cbetao2,4)*(hIm20*(hIm22 + hIm2Neg2) + hRe20*(hRe22 + hRe2Neg2)) - 2.*pow(cbetao2,2)*(-8.*hIm21*hIm2Neg1 + 2.4494897427831781*hIm20*(hIm22 + hIm2Neg2) - 8.*hRe21*hRe2Neg1 + 2.4494897427831781*hRe20*(hRe22 + hRe2Neg2))*pow(sbetao2,2) + 2.4494897427831781*(hIm20*(hIm22 + hIm2Neg2) + hRe20*(hRe22 + hRe2Neg2))*pow(sbetao2,4)) + s1alpha*(pow(cbetao2,8)*(hIm22*hRe21 - 1.*hIm21*hRe22 + hIm2Neg2*hRe2Neg1 - 1.*hIm2Neg1*hRe2Neg2) + pow(cbetao2,6)*(-7.*hIm22*hRe21 + hIm21*(7.34846922834953429*hRe20 + 7.*hRe22) - 7.*hIm2Neg2*hRe2Neg1 - 7.34846922834953429*hIm20*(hRe21 + hRe2Neg1) + hIm2Neg1*(7.34846922834953429*hRe20 + 7.*hRe2Neg2))*pow(sbetao2,2) + 24.494897427831781*pow(cbetao2,4)*(-1.*(hIm21 + hIm2Neg1)*hRe20 + hIm20*(hRe21 + hRe2Neg1))*pow(sbetao2,4) + pow(cbetao2,2)*(-7.*hIm22*hRe21 + hIm21*(7.34846922834953429*hRe20 + 7.*hRe22) - 7.*hIm2Neg2*hRe2Neg1 - 7.34846922834953429*hIm20*(hRe21 + hRe2Neg1) + hIm2Neg1*(7.34846922834953429*hRe20 + 7.*hRe2Neg2))*pow(sbetao2,6) + (hIm22*hRe21 - 1.*hIm21*hRe22 + hIm2Neg2*hRe2Neg1 - 1.*hIm2Neg1*hRe2Neg2)*pow(sbetao2,8)) + c1alpha*(pow(cbetao2,8)*(hIm21*hIm22 - 1.*hIm2Neg1*hIm2Neg2 + hRe21*hRe22 - 1.*hRe2Neg1*hRe2Neg2) + pow(cbetao2,6)*(hIm21*(7.34846922834953429*hIm20 - 7.*hIm22) + hIm2Neg1*(-7.34846922834953429*hIm20 + 7.*hIm2Neg2) + hRe21*(7.34846922834953429*hRe20 - 7.*hRe22) + hRe2Neg1*(-7.34846922834953429*hRe20 + 7.*hRe2Neg2))*pow(sbetao2,2) + 24.494897427831781*pow(cbetao2,4)*(hIm20*(-1.*hIm21 + hIm2Neg1) + hRe20*(-1.*hRe21 + hRe2Neg1))*pow(sbetao2,4) + pow(cbetao2,2)*(hIm21*(7.34846922834953429*hIm20 - 7.*hIm22) + hIm2Neg1*(-7.34846922834953429*hIm20 + 7.*hIm2Neg2) + hRe21*(7.34846922834953429*hRe20 - 7.*hRe22) + hRe2Neg1*(-7.34846922834953429*hRe20 + 7.*hRe2Neg2))*pow(sbetao2,6) + (hIm21*hIm22 - 1.*hIm2Neg1*hIm2Neg2 + hRe21*hRe22 - 1.*hRe2Neg1*hRe2Neg2)*pow(sbetao2,8))));
  return;
}
vector<vector<double> > GWFrames::Waveform::SchmidtEtAlVector(const double alpha0Guess, const double beta0Guess) const {
  std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Not properly implemented yet.  Need correct function (and gradients) to minimize." << std::endl;
  throw(GWFrames_NotYetImplemented);
  const double InitialStepSize = 0.01;
  const double LineMinimizationParameter = 1.e-4;
  const double NormOfGradientOnStop = 1.e-5;
  const unsigned int i2Neg2 = FindModeIndex(2,-2);
  const unsigned int i2Neg1 = FindModeIndex(2,-1);
  const unsigned int i20    = FindModeIndex(2, 0);
  const unsigned int i2Pos1 = FindModeIndex(2, 1);
  const unsigned int i2Pos2 = FindModeIndex(2, 2);
  size_t iterMax = 20;
  size_t iter = 0;
  int status = GSL_CONTINUE;
  const gsl_multimin_fdfminimizer_type* T;
  gsl_multimin_fdfminimizer* s;
  double h[10] = {this->Re(i2Neg2, 0), this->Im(i2Neg2, 0),
		  this->Re(i2Neg1, 0), this->Im(i2Neg1, 0),
		  this->Re(i20,    0), this->Im(i20,    0),
		  this->Re(i2Pos1, 0), this->Im(i2Pos1, 0),
		  this->Re(i2Pos2, 0), this->Im(i2Pos2, 0)};
  gsl_vector* x;
  gsl_multimin_function_fdf SchmidtEtAl_func;
  SchmidtEtAl_func.n = 2;
  SchmidtEtAl_func.f = SchmidtEtAl_f;
  SchmidtEtAl_func.df = SchmidtEtAl_df;
  SchmidtEtAl_func.fdf = SchmidtEtAl_fdf;
  SchmidtEtAl_func.params = h;
  x = gsl_vector_alloc(2);
  gsl_vector_set(x, 0, alpha0Guess);
  gsl_vector_set(x, 1, beta0Guess);
  T = gsl_multimin_fdfminimizer_conjugate_fr;
  s = gsl_multimin_fdfminimizer_alloc (T, 2);

  vector<vector<double> > v(NTimes(), vector<double>(3));
  for(unsigned int iTime=0; iTime<NTimes(); ++iTime) {
    iter = 0;
    double h[10] = {this->Re(i2Neg2, iTime), this->Im(i2Neg2, iTime),
		    this->Re(i2Neg1, iTime), this->Im(i2Neg1, iTime),
		    this->Re(i20,    iTime), this->Im(i20,    iTime),
		    this->Re(i2Pos1, iTime), this->Im(i2Pos1, iTime),
		    this->Re(i2Pos2, iTime), this->Im(i2Pos2, iTime)};
    SchmidtEtAl_func.params = h;
    gsl_multimin_fdfminimizer_set (s, &SchmidtEtAl_func, x, InitialStepSize, LineMinimizationParameter);
    do
      {
	iter++;
	status = gsl_multimin_fdfminimizer_iterate (s);
	if (status) {
	  cerr << "STATUS: " << status << endl;
	  break;
	}
	status = gsl_multimin_test_gradient (s->gradient, NormOfGradientOnStop);
	if (status == GSL_SUCCESS)
	  printf ("Minimum found at:\n");
	printf ("%5d %.16f %.16f %21.16f\n", int(iter),
		gsl_vector_get (s->x, 0),
		gsl_vector_get (s->x, 1),
		s->f);
      }
    while (status == GSL_CONTINUE && iter < iterMax);
    gsl_vector_set(x, 0, gsl_vector_get (s->x, 0)+0e-14);
    gsl_vector_set(x, 1, gsl_vector_get (s->x, 1)+0e-15);
    const double alpha = gsl_vector_get (s->x, 0);
    const double beta = gsl_vector_get (s->x, 1);
    v[iTime][0] = std::sin(beta)*std::cos(alpha);
    v[iTime][1] = std::sin(beta)*std::sin(alpha);
    v[iTime][2] = std::cos(beta);
  }
  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);
  return v;
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


/// Transform Waveform to Schmidt et al. frame.
GWFrames::Waveform& GWFrames::Waveform::TransformToSchmidtEtAlFrame(const double alpha0Guess, const double beta0Guess) {
  ///
  /// \param alpha0Guess Initial guess for optimal direction alpha
  /// \param beta0Guess  Initial guess for optimal direction beta
  ///
  /// This function combines the steps required to obtain the Waveform
  /// in the Schmidt et al. frame.
  history << "### this->TransformToSchmidtEtAlFrame(" << alpha0Guess << ", " << beta0Guess << ")\n#";
  this->frameType = GWFrames::Aligned;
  return this->RotateDecompositionBasis(FrameFromZ(normalized(QuaternionArray(this->SchmidtEtAlVector(alpha0Guess, beta0Guess))), T()));
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
  history << "### this->TransformToOShaughnessyEtAlFrame(" << StringForm(Lmodes) << ")\n#";
  this->frameType = GWFrames::Aligned;
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
  history << "### this->TransformToAngularVelocityFrame(" << StringForm(Lmodes) << ")\n#";
  vector<Quaternion> R_AV = normalized(QuaternionArray(this->AngularVelocityVector(Lmodes)));
  this->frameType = GWFrames::Aligned;
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
  history << "### this->TransformToCorotatingFrame(" << StringForm(Lmodes) << ")\n#";
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
  history << "### this->TransformToInertialFrame();\n#";
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
	      << "\nWarning: Asking a Waveform in the " << GWFrames::WaveformFrameNames[frameType] << " frame  into the corotating frame."
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
  history << "### this->TransformUncertaintiesToCorotatingFrame(R_frame); // R_frame=[" << setprecision(16) << R_frame[0] << ", " << R_frame[1] << ", ...]\n#";
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
  history << "### this->TransformUncertaintiesToInertialFrame();\n#";
  this->frameType = GWFrames::Inertial; // Must come first
  this->RotateDecompositionBasisOfUncertainties(Quaternions::conjugate(frame));
  this->SetFrame(vector<Quaternion>(0));
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
	    << "### *this = this->Interpolate(NewTime);" << std::endl;
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
	    << "### this->Segment(" << i1 << ", " << i2 << ");" << std::endl;
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
  history << "### this->AlignTime(A, " << std::setprecision(16) << t_fid << ");  # deltat=" << deltat << std::endl;

  return *this;
}

/// Find the appropriate rotation to fix the orientation of the corotating frame.
void GWFrames::Waveform::GetAlignmentOfDecompositionFrameToModes(const double t_fid, Quaternion& R_delta, const std::vector<int>& Lmodes) const {
  ///
  /// \param t_fid Fiducial time at which the alignment should happen
  /// \param R_delta Returned rotor
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

  if(frameType != GWFrames::Corotating) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ":"
	      << "\nError: GetAlignmentOfDecompositionFrameToModes only takes Waveforms in the " << GWFrames::WaveformFrameNames[GWFrames::Corotating] << " frame."
	      << "\n       This Waveform is in the " << GWFrames::WaveformFrameNames[frameType] << " frame." << std::endl;
    throw(GWFrames_WrongFrameType);
  }

  // Get direction of angular-velocity vector near t_fid
  Quaternion omegaHat;
  {
    int i_t_fid = 0;
    while(t[i_t_fid]<t_fid && i_t_fid<int(t.size())) { ++i_t_fid; }
    unsigned int i1 = (i_t_fid-10<0 ? 0 : i_t_fid-10);
    unsigned int i2 = (i_t_fid+11>int(t.size()) ? t.size() : i_t_fid+11);
    vector<double> tRegion(&t[i1], &t[i2]);
    Waveform Region = this->Interpolate(tRegion);
    omegaHat = Quaternion(Region.AngularVelocityVector()[i_t_fid-i1]).normalized();
    // Rotate omegaHat back into this frame, if necessary
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

  // R_delta is the rotation we will be applying on the right-hand side
  R_delta = R_V_f * Quaternions::exp(Quaternions::Quaternion(0,0,0,(-phase_22/4)));
}

/// Fix the orientation of the corotating frame.
GWFrames::Waveform& GWFrames::Waveform::AlignDecompositionFrameToModes(const double t_fid, const std::vector<int>& Lmodes) {
  ///
  /// \param t_fid Fiducial time at which the alignment should happen
  /// \param Lmodes Lmodes to use in computing \f$<LL>\f$
  ///
  /// The corotating frame is only defined up to some constant rotor
  /// R_const; if R_corot is corotating, then so is R_corot*R_const.
  /// This function uses that freedom to ensure that the frame is
  /// aligned with the Waveform modes at the fiducial time.  In
  /// particular, it ensures that the Z axis of the frame in which the
  /// decomposition is done is along the dominant eigenvector of
  /// \f$<LL>\f$ (suggested by O'Shaughnessy et al.), and the phase of
  /// the (2,2) mode is zero.
  ///
  /// If Lmodes is empty (default), all L modes are used.  Setting
  /// Lmodes to [2] or [2,3,4], for example, restricts the range of
  /// the sum.
  ///

  // Find the appropriate rotation
  Quaternion R_c;
  GetAlignmentOfDecompositionFrameToModes(t_fid, R_c, Lmodes);

  // Record what happened
  history << "### this->AlignDecompositionFrameToModes(" << std::setprecision(16) << t_fid << ", ...);  # R_c=" << R_c << std::endl;

  // Now, apply the rotation
  this->RotateDecompositionBasis(R_c);

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

  if(frameType != A.frameType && (frameType!=GWFrames::Corotating || frameType!=GWFrames::Aligned)) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
	      << "\nError: GetAlignmentOfFrame assumes that the two Waveforms are in the same type of frame, and that"
	      << "\n       that frame is physically meaningful (either " << GWFrames::WaveformFrameNames[GWFrames::Corotating]
	      << " or " << GWFrames::WaveformFrameNames[GWFrames::Aligned] << "),"
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
  history << "### this->AlignFrame(A, " << std::setprecision(16) << t_fid << ");  # R_delta=" << R_delta << std::endl;

  return *this;
}

#ifndef DOXYGEN
// This is a local object used by AlignTimeAndFrame
class WaveformAligner {
private:
  const GWFrames::Waveform& a, b;
  std::vector<double> t;
  std::vector<Quaternions::Quaternion> aFrame;
public:
  double Multiplier;
  WaveformAligner(const GWFrames::Waveform& A, const GWFrames::Waveform& B, const double t1, const double t2, const double M=1.0)
    : a(A), b(B), t(a.T()), aFrame(a.Frame()), Multiplier(M)
  {
    // Check to make sure we have sufficient times before any offset.
    // (This is necessary but not sufficient for the method to work.)
    if(t1<a.T(0)) {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Alignment time t1=" << t1
		<< " does not occur in Waveform A (which has A.T(0)=" << A.T(0) << ")." << std::endl;
      throw(GWFrames_IndexOutOfBounds);
    }
    if(t1<b.T(0)) {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Alignment time t1=" << t1
		<< " does not occur in Waveform B (which has B.T(0)=" << B.T(0) << ")." << std::endl;
      throw(GWFrames_IndexOutOfBounds);
    }
    if(t2>a.T().back()) {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Alignment time t2=" << t2
		<< " does not occur in Waveform A (which has A.T(-1)=" << A.T().back() << ")." << std::endl;
      throw(GWFrames_IndexOutOfBounds);
    }
    if(t2>B.T().back()) {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Alignment time t2=" << t2
		<< " does not occur in Waveform B (which has B.T(-1)=" << B.T().back() << ")." << std::endl;
      throw(GWFrames_IndexOutOfBounds);
    }
    // Store the fixed frame (A) and its set of times, to which we
    // will interpolate.
    unsigned int i=t.size()-1;
    while(t[i]>t2 && i>0) { --i; }
    t.erase(t.begin()+i, t.end());
    aFrame.erase(aFrame.begin()+i, aFrame.end());
    i=0;
    while(i<a.T().size() && a.T(i)<t1) { ++i; }
    t.erase(t.begin(), t.begin()+i);
    aFrame.erase(aFrame.begin(), aFrame.begin()+i);
  }

  GWFrames::Quaternion AvgLogRatio(const double deltat, const double deltax, const double deltay, const double deltaz) const {
    using namespace GWFrames; // Allow me to subtract a double from a vector<double> below
    const GWFrames::Quaternion R_delta = GWFrames::exp(GWFrames::Quaternion(0, deltax, deltay, deltaz));
    const std::vector<GWFrames::Quaternion> bFrame = GWFrames::Squad(R_delta * b.Frame(), (b.T())+Multiplier*deltat, t);
    const unsigned int Size=bFrame.size();
    GWFrames::Quaternion avg(0.0, 0.0, 0.0, 0.0);
    GWFrames::Quaternion avgdot_last = GWFrames::log( aFrame[0] * GWFrames::inverse(bFrame[0]) );
    for(unsigned int i=1; i<Size; ++i) {
      const GWFrames::Quaternion avgdot = GWFrames::log( aFrame[i] * GWFrames::inverse(bFrame[i]) );
      avg = avg + ((t[i]-t[i-1])/2.0)*(avgdot+avgdot_last);
      avgdot_last = avgdot;
    }
    return avg;
  }

  double EvaluateMinimizationQuantity(const double deltat, const double deltax, const double deltay, const double deltaz) const {
    using namespace GWFrames; // Allow me to subtract a double from a vector<double> below
    const Quaternions::Quaternion R_delta = Quaternions::exp(Quaternions::Quaternion(0, deltax, deltay, deltaz));
    const std::vector<Quaternions::Quaternion> bFrame = Quaternions::Squad(R_delta * b.Frame(), (b.T())+Multiplier*deltat, t);
    const unsigned int Size=bFrame.size();
    double f = 0.0;
    double fdot_last = 4 * Quaternions::normsquared( Quaternions::log( aFrame[0] * Quaternions::inverse(bFrame[0]) ) );
    for(unsigned int i=1; i<Size; ++i) {
      const double fdot = 4 * Quaternions::normsquared( Quaternions::log( aFrame[i] * Quaternions::inverse(bFrame[i]) ) );
      f += (t[i]-t[i-1])*(fdot+fdot_last)/2.0;
      fdot_last = fdot;
    }
    return f;
  }
};
// Local utility function used by GSL below
double minfunc (const gsl_vector* delta, void* params) {
  WaveformAligner* Aligner = (WaveformAligner*) params;
  return Aligner->EvaluateMinimizationQuantity(gsl_vector_get(delta,0),
					       gsl_vector_get(delta,1),
					       gsl_vector_get(delta,2),
					       gsl_vector_get(delta,3));
}
#endif // DOXYGEN

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
  /// \sa AlignTimeAndFrame
  ///

  if(frameType != A.frameType || (frameType!=GWFrames::Corotating && frameType!=GWFrames::Aligned)) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
	      << "\nError: GetAlignmentOfTimeAndFrame assumes that the two Waveforms are in the type of same frame, and that"
	      << "\n       that frame is physically meaningful (either " << GWFrames::WaveformFrameNames[GWFrames::Corotating]
	      << " or " << GWFrames::WaveformFrameNames[GWFrames::Aligned] << "),"
	      << "\n       so that it makes sense to align the frames."
	      << "\n       This Waveform is in the " << GWFrames::WaveformFrameNames[frameType] << " frame,"
	      << "\n       The Waveform in the argument is in the " << GWFrames::WaveformFrameNames[A.frameType] << " frame.\n"
	      << "\n       Also note that each waveform should have its decomposition frame aligned to its modes (see `AlignDecompositionFrameToModes`).\n"
	      << std::endl;
    throw(GWFrames_WrongFrameType);
  }

  const Waveform& B = *this;

  const unsigned int MaxIterations = 2000;
  const double MinSimplexSize = 1.0e-10;
  const double InitialTrialTimeStep = 10.0;
  const double InitialTrialAngleStep = 1.0;

  WaveformAligner Aligner(A, B, t1, t2);
  const unsigned int NDimensions = 4;

  // Use Nelder-Mead simplex minimization
  const gsl_multimin_fminimizer_type* T =
    gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer* s = NULL;
  gsl_vector* ss;
  gsl_vector* x;
  gsl_multimin_function minex_func;
  size_t iter = 0;
  int status = GSL_CONTINUE;
  double size = 0.0;

  // Set initial values
  x = gsl_vector_alloc(NDimensions);
  if(Aligner.EvaluateMinimizationQuantity(0.,0.,0.,0.)/std::abs(t2-t1)<M_PI) {
    // Rotors are nearly aligned
    gsl_vector_set(x, 0, 0.0);
    gsl_vector_set(x, 1, 0.0);
    gsl_vector_set(x, 2, 0.0);
    gsl_vector_set(x, 3, 0.0);
  } else {
    // Rotors are off by too much; start from a better spot -- the avg log of the "difference" ratio
    const GWFrames::Quaternion avglog = Aligner.AvgLogRatio(0., 0., 0., 0.);
    gsl_vector_set(x, 0, 0.0);
    gsl_vector_set(x, 1, avglog[1]);
    gsl_vector_set(x, 2, avglog[2]);
    gsl_vector_set(x, 3, avglog[3]);
  }

  // Set initial step sizes
  ss = gsl_vector_alloc(NDimensions);
  gsl_vector_set(ss, 0, InitialTrialTimeStep);
  gsl_vector_set(ss, 1, InitialTrialAngleStep);
  gsl_vector_set(ss, 2, InitialTrialAngleStep);
  gsl_vector_set(ss, 3, InitialTrialAngleStep);

  minex_func.n = NDimensions;
  minex_func.f = &minfunc;
  minex_func.params = (void*) &Aligner;

  s = gsl_multimin_fminimizer_alloc(T, NDimensions);
  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

  // Run the minimization
  while(status == GSL_CONTINUE && iter < MaxIterations) {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
    // std::cout << iter
    if(status) break;
    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, MinSimplexSize);
  }

  // Apply time shift and rotation
  deltat = gsl_vector_get(s->x, 0) * Aligner.Multiplier;
  R_delta = Quaternions::exp(Quaternions::Quaternion(0.0, gsl_vector_get(s->x, 1), gsl_vector_get(s->x, 2), gsl_vector_get(s->x, 3)));

  // Free allocated memory
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);

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
  /// As noted above, it is implicitly assumed that both Waveforms are
  /// in their corotating frames, with the modes appropriately aligned
  /// to the frames at t_fid.  The assumption is that the frames
  /// actually represent something physically meaningful, so that it
  /// is meaningful to insist that they be the same.
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
  history << "### this->AlignTimeAndFrame(A, " << std::setprecision(16) << t1 << ", " << t2 << ");  # deltat=" << deltat << "; R_delta=" << R_delta << std::endl;

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
  C.history << "### B.Compare(A)\n"
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
  /// is assumed to have been done already.  The transition function is a smooth
  ///
  /// Note that this function does NOT operate in place; a new
  /// Waveform object is constructed and returned.

  // Make A a convenient alias
  const GWFrames::Waveform& A = *this;

  // Check to see if the frameTypes are the same
  if(A.frameType != B.frameType) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
	      << "\nWarning:"
	      << "\n       This Waveform is in the " << GWFrames::WaveformFrameNames[A.frameType] << " frame,"
	      << "\n       The Waveform in the argument is in the " << GWFrames::WaveformFrameNames[B.frameType] << " frame."
	      << "\n       Comparing them probably does not make sense.\n"
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
  // Store both old histories in C's
  C.history << "### A.Hybridize(B, " << t1 << ", " << t2 << ", " << tMinStep << ")\n"
	    << "#### A.history.str():\n" << A.history.str()
	    << "#### B.history.str():\n" << B.history.str()
	    << "#### End of old histories from `Hybridize`" << std::endl;
  // The new time axis will be the union of the two old ones
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


/// Output Waveform object to data file.
const GWFrames::Waveform& GWFrames::Waveform::Output(const std::string& FileName, const unsigned int precision) const {
  const std::string Descriptor = DescriptorString();
  ofstream ofs(FileName.c_str(), ofstream::out);
  ofs << setprecision(precision) << flush;
  ofs << history.str() << "### this->Output(" << FileName << ", " << precision << ")" << endl;
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
// 	for(unsigned int i_t=0; i_t<N_t; ++i_t) {
// 	  data[iPositive][i_t].imag() *= -1;
// 	  data[iNegative][i_t].imag() *= -1;
// 	}
//       } else {
// 	for(unsigned int i_t=0; i_t<N_t; ++i_t) {
// 	  data[iPositive][i_t].real() *= -1;
// 	  data[iNegative][i_t].real() *= -1;
// 	}
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
  C.history << "### *this = *this+B\n"
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
  C.history << "### *this = *this-B\n"
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
  C.history << "### *this = (*this) * " << b << std::endl;

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
  C.history << "### *this = (*this) / " << b << std::endl;

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
  EdthA.history << "### this->NPEdth();" << endl;

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
  EdthBarA.history << "### this->NPEdthBar();" << endl;

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
  EdthA.history << "### this->GHPEdth();" << endl;

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
  EdthBarA.history << "### this->GHPEdthBar();" << endl;

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
  IntegralEdthA.history << "### this->IntegrateNPEdth();" << endl;

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
  IntegralEdthBarA.history << "### this->IntegrateNPEdthBar();" << endl;

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
  IntegralEdthA.history << "### this->IntegrateGHPEdth();" << endl;

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
  IntegralEdthBarA.history << "### this->GHPEdthBar();" << endl;

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
  B.history << "### *this = B.ApplySupertranslation([" << gamma[0];
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


/// Apply a boost to a boost-weighted function
GWFrames::Waveform GWFrames::Waveform::Boost(const std::vector<double>& v) const {
  /// This function does three things.  First, it evaluates the
  /// Waveform on what will become an equi-angular grid after
  /// transformation by the boost.  Second, it multiplies each of
  /// those points by the appropriate conformal factor
  /// \f$K^b(\vartheta, \varphi)\f$, where \f$b\f$ is the boost weight
  /// stored with the Waveform.  Finally, it transforms back to
  /// Fourier space using that new equi-angular grid.


  std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Not properly implemented yet." << std::endl;
  throw(GWFrames_NotYetImplemented);

  //return Waveform();
  // return B;
}


/// Apply a boost to a boost- and spin-weighted function
GWFrames::Waveform& GWFrames::Waveform::Boost(const std::vector<std::vector<double> >& v) {
  /// This function does three things.  First, it evaluates the
  /// Waveform on what will become an equi-angular grid after
  /// transformation by the boost.  Second, it multiplies each of
  /// those points by the appropriate conformal factor
  /// \f$K^b(\vartheta, \varphi)\f$, where \f$b\f$ is the boost weight
  /// stored with the Waveform.  Finally, it transforms back to
  /// Fourier space using that new equi-angular grid.

  // Check the size of the input velocity
  if(v.size()!=NTimes()) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": (v.size()=" << v.size() << ") != (NTimes()=" << NTimes() << ")" << std::endl;
    throw(GWFrames_VectorSizeMismatch);
  }

  // Set up storage and calculate useful constants
  const int ellMax = this->EllMax();
  const int n_theta= 2*ellMax+1;
  const int n_phi = 2*ellMax+1;
  const double dtheta = M_PI/double(n_theta-1); // theta should return to M_PI
  const double dphi = 2*M_PI/double(n_phi); // phi should not return to 2*M_PI
  vector<complex<double> > Grid(n_phi*n_theta);
  SphericalFunctions::SWSH sYlm(SpinWeight());

  // Main loop over time steps
  for(unsigned int i_t=0; i_t<NTimes(); ++i_t) {
    if(v[i_t].size()!=3) {
      std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": v[" << i_t << "].size()=" << v[i_t].size()
		<< ".  Input is assumed to be a vector of three-velocities." << std::endl;
      throw(GWFrames_VectorSizeMismatch);
    }

    vector<complex<double> > Modes(N_lm(ellMax), 0.0);
    vector<complex<double> > Modes2(N_lm(ellMax), 0.0);

    const double gamma = 1.0/std::sqrt(1.0-v[i_t][0]*v[i_t][0]-v[i_t][1]*v[i_t][1]-v[i_t][2]*v[i_t][2]);

    // Fill the Modes data for this time step
    for(int i_mode=0; i_mode<N_lm(std::abs(SpinWeight())-1); ++i_mode) {
      Modes[i_mode] = 0.0;
    }
    for(int i_mode=N_lm(std::abs(SpinWeight())-1), ell=std::abs(SpinWeight()); ell<=ellMax; ++ell) {
      for(int m=-ell; m<=ell; ++m, ++i_mode) {
	Modes[i_mode] = this->Data(FindModeIndex(ell,m), i_t);
      }
    }

    // Construct the data on the distorted grid
    for(int i_g=0, i_theta=0; i_theta<n_theta; ++i_theta) {
      for(int i_phi=0; i_phi<n_phi; ++i_phi, ++i_g) {
	const Quaternion Rp(dtheta*i_theta, dphi*i_phi);
	const Quaternion nHat = Rp*Quaternions::zHat*Rp.conjugate();
	const double ConformalFactor = std::pow(gamma*(1-nHat.dot(Quaternions::Quaternion(v[i_t]))), -BoostWeight());
	// const Quaternion R_b = Quaternions::BoostRotor(-v[i_t], nHat.vec());
	const Quaternion R_b = Quaternions::BoostRotor(v[i_t], nHat.vec());
	if(i_t==0 && i_theta==0 && i_phi==0) {
	  std::cerr << __FILE__ << ":" << __LINE__ << ": Previous line is given wrong sign for testing purposes!\n\n" << std::endl;
	}
	sYlm.SetRotation(R_b*Rp);

	// Evaluate the data at this point, incorporating the spin due
	// to the boost, multiplying by the conformal factor at this
	// point
	Grid[i_g] = ConformalFactor*sYlm.Evaluate(Modes);
      }
    }

    // Decompose the data into modes
    spinsfast_map2salm(reinterpret_cast<fftw_complex*>(&Grid[0]),
		       reinterpret_cast<fftw_complex*>(&Modes2[0]),
		       SpinWeight(), n_theta, n_phi, ellMax);

    // Set new data at this time step
    for(int i_mode=N_lm(std::abs(SpinWeight())-1), ell=std::abs(SpinWeight()); ell<=ellMax; ++ell) {
      for(int m=-ell; m<=ell; ++m, ++i_mode) {
	this->SetData(this->FindModeIndex(ell,m), i_t, Modes2[i_mode]);
	// if(i_t%1000==0 && ell==2 && m==2) {
	//   std::cerr << i_t << "/" << NTimes()
	// 	    << "\tModes[" << i_mode << "] = " << Modes[i_mode]
	// 	    << "\tModes2[" << i_mode << "] = " << Modes2[i_mode]
	// 	    << "\tData((" << ell << ", " << m << "), " << i_t << ") = " << this->Data(this->FindModeIndex(ell,m),i_t) << std::endl;
	// }
      }
    }

  }

  return *this;
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
