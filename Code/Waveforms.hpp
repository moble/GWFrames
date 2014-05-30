// Copyright (c) 2014, Michael Boyle
// See LICENSE file for details

#ifndef WAVEFORMS_HPP
#define WAVEFORMS_HPP

#include <complex>
#ifndef DOXYGEN
#ifdef __restrict
#define restrict __restrict
#endif
extern "C" {
  #include "fftw3.h"
  #include "alm.h"
  #include "wigner_d_halfpi.h"
  #include "spinsfast_forward.h"
  #include "spinsfast_backward.h"
}
#undef complex
#endif // DOXYGEN
#include <cmath>
#include <string>
#include <sstream>

#include "Quaternions.hpp"
#include "Utilities.hpp"
#include "Errors.hpp"

namespace GWFrames {

  // Note on Waveform Types:
  // In any system, h -- being strain -- should be dimensionless.
  // When G=c=1, the dimensionless quantities are rMPsi4, rhdot, and
  // rhOverM; as are rOverM and tOverM.  When G and c are
  // dimensionful, the dimensionless quantities are
  //     -  (r/c) * (M*G/c^3) * Psi4
  //     -  (r/c) * hdot
  //     -  (r/c) * h / (M*G/c^3)
  //     -  (r/c) / (M*G/c^3)
  //     -  t / (M*G/c^3)
  // To regain the dimensionful quantities, we simply need to remove
  // the relevant dimensionful elements (e.g., the r and M factors).

  enum WaveformFrameType { UnknownFrameType, Inertial, Coprecessing, Coorbital, Corotating };
  static const std::string WaveformFrameNames[5] = { "UnknownFrameType", "Inertial", "Coprecessing", "Coorbital", "Corotating" };
  enum WaveformDataType { UnknownDataType, h, hdot, Psi4 };
  static const std::string WaveformDataNames[4] = { "UnknownDataType", "h", "hdot", "Psi4" };
  static const std::string WaveformDataNamesLaTeX[4] = { "\\mathrm{unknown data type}", "h", "\\dot{h}", "\\Psi_4" };
  const int WeightError = 1000;

  /// Object storing data and other information for a single waveform
  class Waveform {

  protected:  // Member data
    int spinweight;
    int boostweight;
    std::stringstream history;
    std::vector<double> t;
    std::vector<Quaternions::Quaternion> frame;
    WaveformFrameType frameType;
    WaveformDataType dataType;
    bool rIsScaledOut;
    bool mIsScaledOut;
    std::vector<std::vector<int> > lm;
    MatrixC data; // Each row (first index, nn) corresponds to a mode

  public:  // Constructors and Destructor
    Waveform();
    Waveform(const Waveform& W);
    Waveform(const std::string& FileName, const std::string& DataFormat);
    Waveform(const std::vector<double>& T, const std::vector<std::vector<int> >& LM,
             const std::vector<std::vector<std::complex<double> > >& Data);
    ~Waveform() { }
    Waveform& operator=(const Waveform&);

  public:  // Copy-ish constructoroids
    Waveform CopyWithoutData() const;
    Waveform SliceOfTimeIndices(const unsigned int i_t_a, unsigned int t_b=0) const;
    Waveform SliceOfTimeIndicesWithEll2(const unsigned int i_t_a, unsigned int i_t_b=0) const;
    Waveform SliceOfTimeIndicesWithoutModes(const unsigned int i_t_a, unsigned int i_t_b=0) const;
    Waveform SliceOfTimes(const double t_a=-1e300, const double t_b=1e300) const;
    Waveform SliceOfTimesWithEll2(const double t_a=-1e300, const double t_b=1e300) const;
    Waveform SliceOfTimesWithoutModes(const double t_a=-1e300, const double t_b=1e300) const;
    Waveform Segment(const unsigned int i1, const unsigned int i2) const;
    Waveform Interpolate(const std::vector<double>& NewTime, const bool AllowTimesOutsideCurrentDomain=false) const;
    Waveform& InterpolateInPlace(const std::vector<double>& NewTime);

  public: // Data alteration functions -- USE AT YOUR OWN RISK!
    Waveform& DropTimesOutside(const double ta, const double tb);
    Waveform& DropEllModes(const std::vector<unsigned int>& EllModesToDrop);
    Waveform& KeepOnlyEllModes(const std::vector<unsigned int>& EllModesToKeep);
    Waveform& KeepOnlyEll2();
    inline Waveform& SetSpinWeight(const int NewSpinWeight) { spinweight=NewSpinWeight; return *this; }
    inline Waveform& SetBoostWeight(const int NewBoostWeight) { boostweight=NewBoostWeight; return *this; }
    inline Waveform& AppendHistory(const std::string& Hist) { history << Hist; return *this; }
    inline Waveform& SetHistory(const std::string& Hist) { history.str(Hist); history.seekp(0, std::ios_base::end); return *this; }
    inline Waveform& SetT(const std::vector<double>& a) { t = a; return *this; }
    inline Waveform& SetTime(const std::vector<double>& a) { t = a; return *this; }
    inline Waveform& SetFrame(const std::vector<Quaternions::Quaternion>& a) { frame = a; return *this; }
    inline Waveform& SetFrameType(const WaveformFrameType Type) { frameType = Type; return *this; }
    inline Waveform& SetDataType(const WaveformDataType Type) { dataType = Type; return *this; }
    inline Waveform& SetRIsScaledOut(const bool Scaled) { rIsScaledOut = Scaled; return *this; }
    inline Waveform& SetMIsScaledOut(const bool Scaled) { mIsScaledOut = Scaled; return *this; }
    inline Waveform& SetLM(const std::vector<std::vector<int> >& a) { lm = a; return *this; }
    inline Waveform& SetData(const std::vector<std::vector<std::complex<double> > >& a) { data = MatrixC(a); return *this; }
    inline Waveform& SetData(const unsigned int i_Mode, const unsigned int i_Time, const std::complex<double>& a) { data[i_Mode][i_Time] = a; return *this; }
    inline Waveform& ResizeData(const unsigned int NModes, const unsigned int NTimes) { data.resize(NModes, NTimes); return *this; }
    void swap(Waveform& b);

  public:  // Data access functions
    inline unsigned int NTimes() const { return t.size(); }
    inline unsigned int NModes() const { return data.nrows(); }
    inline int SpinWeight() const { return spinweight; }
    inline int BoostWeight() const { return boostweight; }
    inline std::string HistoryStr() const { return history.str(); }
    inline std::stringstream& HistoryStream() { return history; }
    inline int FrameType() const { return frameType; }
    inline int DataType() const { return dataType; }
    inline std::string FrameTypeString() const { return WaveformFrameNames[frameType]; }
    inline std::string DataTypeString() const { return WaveformDataNames[dataType]; }
    inline std::string DataTypeLaTeXString() const { return WaveformDataNamesLaTeX[dataType]; }
    std::string DescriptorString() const;
    inline bool RIsScaledOut() const { return rIsScaledOut; }
    inline bool MIsScaledOut() const { return mIsScaledOut; }
    inline double T(const unsigned int TimeIndex) const { return t[TimeIndex]; }
    inline Quaternions::Quaternion Frame(const unsigned int TimeIndex) const { return (frame.size()>1 ? frame[TimeIndex] : frame[0]); }
    inline double Re(const unsigned int Mode, const unsigned int TimeIndex) const { return std::real(data[Mode][TimeIndex]); }
    inline double Im(const unsigned int Mode, const unsigned int TimeIndex) const { return std::imag(data[Mode][TimeIndex]); }
    inline double Abs(const unsigned int Mode, const unsigned int TimeIndex) const { return std::abs(data[Mode][TimeIndex]); }
    inline double Arg(const unsigned int Mode, const unsigned int TimeIndex) const { return std::arg(data[Mode][TimeIndex]); }
    inline std::complex<double> Data(const unsigned int Mode, const unsigned int TimeIndex) const { return data[Mode][TimeIndex]; }
    inline std::complex<double> operator()(const unsigned int Mode, const unsigned int TimeIndex) const { return data[Mode][TimeIndex]; }
    inline const std::vector<int>& LM(const unsigned int Mode) const { return lm[Mode]; }
    std::vector<double> Re(const unsigned int Mode) const;
    std::vector<double> Im(const unsigned int Mode) const;
    std::vector<double> Abs(const unsigned int Mode) const;
    std::vector<double> Arg(const unsigned int Mode) const;
    inline std::vector<double> ArgUnwrapped(const unsigned int Mode) const { return Unwrap(Arg(Mode)); }
    std::vector<std::complex<double> > Data(const unsigned int Mode) const;
    inline const std::complex<double>* operator()(const unsigned int Mode) const { return data[Mode]; }
    inline const std::vector<double>& T() const { return t; }
    inline const std::vector<Quaternions::Quaternion>& Frame() const { return frame; }
    inline const std::vector<std::vector<int> >& LM() const { return lm; }
    std::vector<std::vector<double> > Re() const;
    std::vector<std::vector<double> > Im() const;
    std::vector<std::vector<double> > Abs() const;
    std::vector<std::vector<double> > Arg() const;
    std::vector<std::vector<double> > ArgUnwrapped() const;
    std::vector<std::vector<std::complex<double> > > Data() const;

  public: // Data characterization
    int EllMax() const;
    unsigned int FindModeIndex(const int L, const int M) const;
    unsigned int FindModeIndexWithoutError(const int L, const int M) const;
    std::vector<std::complex<double> > DataDot(const unsigned int Mode) const;
    Waveform& Differentiate();
    std::vector<double> Norm(const bool TakeSquareRoot=false) const;
    unsigned int MaxNormIndex(const unsigned int SkipFraction=4) const;
    inline double MaxNormTime(const unsigned int SkipFraction=4) const { return T(MaxNormIndex(SkipFraction)); }
    std::vector<double> NormalizedAntisymmetry(std::vector<int> LModesForAsymmetry=std::vector<int>(0)) const;
    std::vector<std::vector<double> > DipoleMoment(int ellMax=0) const;
    std::vector<double> ZParityViolationSquared(std::vector<int> Lmodes=std::vector<int>(0)) const;
    std::vector<double> ZParityViolationNormalized(std::vector<int> Lmodes=std::vector<int>(0)) const;
    std::vector<double> ZParityViolationMinimized() const;

  private: // Private function for use in public rotation methods
    Waveform& TransformModesToRotatedFrame(const std::vector<Quaternions::Quaternion>& R_frame);

  public:  // Member functions
    // Rotate by the given Quaternion(s)
    Waveform& RotatePhysicalSystem(const Quaternions::Quaternion& R_phys);
    Waveform& RotatePhysicalSystem(std::vector<Quaternions::Quaternion> R_phys);
    Waveform& RotateDecompositionBasis(const Quaternions::Quaternion& R_frame);
    Waveform& RotateDecompositionBasis(const std::vector<Quaternions::Quaternion>& R_frame);

    std::vector<std::vector<double> > LdtVector(std::vector<int> Lmodes=std::vector<int>(0)) const;
    std::vector<Matrix> LLMatrix(std::vector<int> Lmodes=std::vector<int>(0)) const;
    std::vector<std::vector<double> > LLDominantEigenvector(const std::vector<int>& Lmodes=std::vector<int>(0)) const;
    std::vector<std::vector<double> > AngularVelocityVector(const std::vector<int>& Lmodes=std::vector<int>(0)) const;
    std::vector<Quaternions::Quaternion> CorotatingFrame(const std::vector<int>& Lmodes=std::vector<int>(0)) const;

    // Convenient transformations
    Waveform& TransformToCoprecessingFrame(const std::vector<int>& Lmodes=std::vector<int>(0));
    Waveform& TransformToAngularVelocityFrame(const std::vector<int>& Lmodes=std::vector<int>(0));
    Waveform& TransformToCorotatingFrame(const std::vector<int>& Lmodes=std::vector<int>(0));
    Waveform& TransformToInertialFrame();

    // Alignment, comparison, and hybridization
    void GetAlignmentOfTime(const Waveform& A, const double t_fid, double& deltat) const;
    Waveform& AlignTime(const Waveform& A, const double t_fid);
    std::vector<Quaternions::Quaternion> GetAlignmentsOfDecompositionFrameToModes(const std::vector<Quaternions::Quaternion>& nHat_t_fid,
                                                                                  const std::vector<int>& Lmodes=std::vector<int>(0)) const;
    // void GetAlignmentOfDecompositionFrameToModes(const double t_fid, Quaternions::Quaternion& R_eps, const std::vector<int>& Lmodes=std::vector<int>(0)) const;
    // Waveform& AlignDecompositionFrameToModes(const double t_fid, const std::vector<int>& Lmodes=std::vector<int>(0));
    Quaternions::Quaternion GetAlignmentOfDecompositionFrameToModes(const double t_fid, const Quaternions::Quaternion& nHat_t_fid,
                                                                    const std::vector<int>& Lmodes=std::vector<int>(0)) const {
      Quaternions::Quaternion R_eps;
      GetAlignmentOfDecompositionFrameToModes(t_fid, nHat_t_fid, R_eps, Lmodes);
      return R_eps;
    }
    void GetAlignmentOfDecompositionFrameToModes(const double t_fid, const Quaternions::Quaternion& nHat_t_fid, Quaternions::Quaternion& R_eps,
                                                 const std::vector<int>& Lmodes=std::vector<int>(0)) const;
    Waveform& AlignDecompositionFrameToModes(const double t_fid, const Quaternions::Quaternion& nHat_t_fid, const std::vector<int>& Lmodes=std::vector<int>(0));
    void GetAlignmentOfDecompositionFrameToModes(const double t1, const double t2, const Quaternions::Quaternion& nHat_t1,
                                                 Quaternions::Quaternion& R_eps, const std::vector<int>& Lmodes=std::vector<int>(0)) const;
    Waveform& AlignDecompositionFrameToModes(const double t1, const double t2, const Quaternions::Quaternion& nHat_t1, const std::vector<int>& Lmodes=std::vector<int>(0));
    void GetAlignmentOfFrame(const Waveform& A, const double t_fid, Quaternions::Quaternion& R_delta) const;
    Waveform& AlignFrame(const Waveform& A, const double t_fid);
    void GetAlignmentOfTimeAndFrame(const Waveform& A, const double t1, const double t2, double& deltat, Quaternions::Quaternion& R_delta) const;
    Waveform& AlignTimeAndFrame(const Waveform& A, const double t1, const double t2);
    Waveform Compare(const Waveform& B, const double MinTimeStep=0.005, const double MinTime=-3.0e300) const;
    Waveform Hybridize(const Waveform& B, const double t1, const double t2, const double tMinStep=0.005) const;
    // See also GWFrames::AlignWaveforms below; that is not a member
    // of this class, to make clear that both waveforms will be
    // altered inside that function.

    // Pointwise operations and spin-weight operators
    std::vector<std::complex<double> > EvaluateAtPoint(const double vartheta, const double varphi) const;
    std::complex<double> EvaluateAtPoint(const double vartheta, const double varphi, const unsigned int i_t) const;
    std::complex<double> InterpolateToPoint(const double vartheta, const double varphi, const double t_i) const;
    template <typename Op> Waveform BinaryOp(const Waveform& b) const;
    Waveform operator+(const Waveform& B) const;
    Waveform operator-(const Waveform& B) const;
    Waveform operator*(const Waveform& B) const;
    Waveform operator/(const Waveform& B) const;
    // Waveform operator+(const double b) const;
    // Waveform operator-(const double b) const;
    Waveform operator*(const double b) const;
    Waveform operator/(const double b) const;
    Waveform NPEdth() const;
    Waveform NPEdthBar() const;
    Waveform GHPEdth() const;
    Waveform GHPEdthBar() const;
    Waveform IntegrateNPEdth() const;
    Waveform IntegrateNPEdthBar() const;
    Waveform IntegrateGHPEdth() const;
    Waveform IntegrateGHPEdthBar() const;
    Waveform ApplySupertranslation(std::vector<std::complex<double> >& gamma) const;
    // Waveform Boost(const std::vector<double>& v) const;
    Waveform& BoostPsi4(const std::vector<std::vector<double> >& v);
    Waveform Translate(const std::vector<std::vector<double> >& x) const;

    // Output to data file
    const Waveform& Output(const std::string& FileName, const unsigned int precision=14) const;

    // // Correct the error in older RWZ data from SpEC
    // Waveform& HackSpECSignError();

  }; // class Waveform
  inline Waveform operator*(const double b, const Waveform& A) { return A*b; }
  #include "Waveforms_BinaryOp.ipp"

  void AlignWaveforms(Waveform& A, Waveform& B, const std::vector<double>& nHat_A,
                      const std::vector<std::vector<double> >& nHat_B, const std::vector<double>& t_B, const double t_1, const double t_2);



  /// Object storing a collection of Waveform objects to be operated on uniformly
  class Waveforms { // (plural!)

  private:  // Member data
    std::vector<Waveform> Ws;
    bool CommonTimeSet;

  public:  // Constructors and Destructor
    Waveforms(const int N=0);
    Waveforms(const Waveforms& In);
    Waveforms(const std::vector<Waveform>& In);
    ~Waveforms() { }

  public:  // Operators
    inline const Waveform& operator[](const int i) const { return Ws[i]; }
    inline Waveform& operator[](const int i) { return Ws[i]; }

  public:  // Members
    void clear() { Ws.clear(); }
    inline unsigned int size() const { return Ws.size(); }
    void SetCommonTime(std::vector<std::vector<double> >& Radii,
                       const double MinTimeStep=0.005, const double EarliestTime=-3e300, const double LatestTime=3e300);
    Waveforms Extrapolate(std::vector<std::vector<double> >& Radii,
                          const std::vector<int>& ExtrapolationOrders,
                          const std::vector<double>& Omegas=std::vector<double>(0));

  }; // class Waveforms (plural!)

} // namespace GWFrames

#endif // WAVEFORMS_HPP
