// Copyright (c) 2026, Michael Boyle
// See LICENSE file for details

#include "bindings.hpp"

#include <iomanip>
#include <sstream>

#include "PNWaveforms.hpp"
#include "Waveforms.hpp"
#include "WaveformsAtAPointFT.hpp"

namespace GWFramesBindings {

using GWFrames::PNWaveform;
using GWFrames::Waveform;
using GWFrames::WaveformAtAPointFT;
using Quaternions::Quaternion;
using DVec = std::vector<double>;
using IVec = std::vector<int>;
using QVec = std::vector<Quaternion>;
using CVec = std::vector<std::complex<double>>;

namespace {

std::string WaveformStr(const Waveform& W) {
  std::ostringstream S;
  S << W.HistoryStr() << "#\n"
    << "# <in_python>\n"
    << "import GWFrames\n"
    << "print(this)\n"
    << "# </in_python>" << std::endl
    << std::setprecision(14);
  for (unsigned int t = 0; t < W.NTimes(); ++t) {
    S << W.T(t) << " ";
    for (unsigned int mode = 0; mode < W.NModes(); ++mode) {
      S << W.Re(mode, t) << " " << W.Im(mode, t) << " ";
    }
    S << std::endl;
  }
  return S.str();
}

}  // namespace

void BindWaveforms(py::module_& m) {
  // Frame types and data types, as plain integers
  m.attr("UnknownFrameType") = int(GWFrames::UnknownFrameType);
  m.attr("Inertial") = int(GWFrames::Inertial);
  m.attr("Coprecessing") = int(GWFrames::Coprecessing);
  m.attr("Coorbital") = int(GWFrames::Coorbital);
  m.attr("Corotating") = int(GWFrames::Corotating);
  m.attr("UnknownDataType") = int(GWFrames::UnknownDataType);
  m.attr("h") = int(GWFrames::h);
  m.attr("hdot") = int(GWFrames::hdot);
  m.attr("Psi4") = int(GWFrames::Psi4);
  m.attr("Psi3") = int(GWFrames::Psi3);
  m.attr("Psi2") = int(GWFrames::Psi2);
  m.attr("Psi1") = int(GWFrames::Psi1);
  m.attr("Psi0") = int(GWFrames::Psi0);

  const auto noL = IVec();
  py::class_<Waveform> W(m, "Waveform", Doc("Waveform"));
  W.def(py::init<>())
      .def(py::init<const Waveform&>())
      .def(py::init<const std::string&, const std::string&>(), py::arg("FileName"), py::arg("DataFormat"))
      .def(py::init<const DVec&, const std::vector<IVec>&, const std::vector<CVec>&>(), py::arg("T"),
           py::arg("LM"), py::arg("Data"))
      .def("__str__", &WaveformStr)
      .def("__repr__", &Waveform::HistoryStr)
      .def("__copy__", [](const Waveform& w) { return Waveform(w); })
      .def("__deepcopy__", [](const Waveform& w, py::dict) { return Waveform(w); }, py::arg("memo"))

      // Copy-ish constructoroids
      .def("CopyWithoutData", &Waveform::CopyWithoutData, Doc("Waveform::CopyWithoutData"))
      .def("SliceOfTimeIndices", &Waveform::SliceOfTimeIndices, py::arg("i_t_a"), py::arg("i_t_b") = 0,
           Doc("Waveform::SliceOfTimeIndices"))
      .def("SliceOfTimeIndicesWithEll2", &Waveform::SliceOfTimeIndicesWithEll2, py::arg("i_t_a"),
           py::arg("i_t_b") = 0, Doc("Waveform::SliceOfTimeIndicesWithEll2"))
      .def("SliceOfTimeIndicesWithoutModes", &Waveform::SliceOfTimeIndicesWithoutModes, py::arg("i_t_a"),
           py::arg("i_t_b") = 0, Doc("Waveform::SliceOfTimeIndicesWithoutModes"))
      .def("SliceOfTimes", &Waveform::SliceOfTimes, py::arg("t_a") = -1e300, py::arg("t_b") = 1e300,
           Doc("Waveform::SliceOfTimes"))
      .def("SliceOfTimesWithEll2", &Waveform::SliceOfTimesWithEll2, py::arg("t_a") = -1e300,
           py::arg("t_b") = 1e300, Doc("Waveform::SliceOfTimesWithEll2"))
      .def("SliceOfTimesWithoutModes", &Waveform::SliceOfTimesWithoutModes, py::arg("t_a") = -1e300,
           py::arg("t_b") = 1e300, Doc("Waveform::SliceOfTimesWithoutModes"))
      .def("Interpolate", &Waveform::Interpolate, py::arg("NewTime"),
           py::arg("AllowTimesOutsideCurrentDomain") = false, Doc("Waveform::Interpolate"))
      .def("InterpolateInPlace", ReturnSelf(&Waveform::InterpolateInPlace), py::arg("NewTime"),
           Doc("Waveform::InterpolateInPlace"))

      // Data alteration
      .def("DropTimesOutside", ReturnSelf(&Waveform::DropTimesOutside), py::arg("ta"), py::arg("tb"),
           Doc("Waveform::DropTimesOutside"))
      .def("DropEllModes", ReturnSelf(&Waveform::DropEllModes), py::arg("EllModesToDrop"),
           Doc("Waveform::DropEllModes"))
      .def("KeepOnlyEllModes", ReturnSelf(&Waveform::KeepOnlyEllModes), py::arg("EllModesToKeep"),
           Doc("Waveform::KeepOnlyEllModes"))
      .def("KeepOnlyEll2", ReturnSelf(&Waveform::KeepOnlyEll2), Doc("Waveform::KeepOnlyEll2"))
      .def("SetSpinWeight", ReturnSelf(&Waveform::SetSpinWeight), py::arg("NewSpinWeight"))
      .def("SetBoostWeight", ReturnSelf(&Waveform::SetBoostWeight), py::arg("NewBoostWeight"))
      .def("AppendHistory", ReturnSelf(&Waveform::AppendHistory), py::arg("Hist"))
      .def("SetHistory", ReturnSelf(&Waveform::SetHistory), py::arg("Hist"))
      .def("SetVersionHist", ReturnSelf(&Waveform::SetVersionHist), py::arg("VerHist"))
      .def("SetT", ReturnSelf(&Waveform::SetT), py::arg("a"))
      .def("SetTime", ReturnSelf(&Waveform::SetTime), py::arg("a"))
      .def("SetFrame", ReturnSelf(&Waveform::SetFrame), py::arg("a"))
      .def("SetFrameType", [](py::object self, const int type) {
        self.cast<Waveform&>().SetFrameType(static_cast<GWFrames::WaveformFrameType>(type));
        return self;
      }, py::arg("Type"))
      .def("SetDataType", [](py::object self, const int type) {
        self.cast<Waveform&>().SetDataType(static_cast<GWFrames::WaveformDataType>(type));
        return self;
      }, py::arg("Type"))
      .def("SetRIsScaledOut", ReturnSelf(&Waveform::SetRIsScaledOut), py::arg("Scaled"))
      .def("SetMIsScaledOut", ReturnSelf(&Waveform::SetMIsScaledOut), py::arg("Scaled"))
      .def("SetLM", ReturnSelf(&Waveform::SetLM), py::arg("a"))
      .def("SetData", ReturnSelf(py::overload_cast<const std::vector<CVec>&>(&Waveform::SetData)), py::arg("a"))
      .def("SetData",
           ReturnSelf(py::overload_cast<const unsigned int, const unsigned int, const std::complex<double>&>(
               &Waveform::SetData)),
           py::arg("i_Mode"), py::arg("i_Time"), py::arg("a"))
      .def("ResizeData", ReturnSelf(&Waveform::ResizeData), py::arg("NModes"), py::arg("NTimes"))
      .def("swap", &Waveform::swap, py::arg("b"))

      // Data access
      .def("NTimes", &Waveform::NTimes)
      .def("NModes", &Waveform::NModes)
      .def("SpinWeight", &Waveform::SpinWeight)
      .def("BoostWeight", &Waveform::BoostWeight)
      .def("HistoryStr", &Waveform::HistoryStr)
      .def("VersionHist", &Waveform::VersionHist)
      .def("FrameType", &Waveform::FrameType)
      .def("DataType", &Waveform::DataType)
      .def("FrameTypeString", &Waveform::FrameTypeString)
      .def("DataTypeString", &Waveform::DataTypeString)
      .def("DataTypeLaTeXString", &Waveform::DataTypeLaTeXString)
      .def("DescriptorString", &Waveform::DescriptorString, Doc("Waveform::DescriptorString"))
      .def("RIsScaledOut", &Waveform::RIsScaledOut)
      .def("MIsScaledOut", &Waveform::MIsScaledOut)
      .def("T", py::overload_cast<>(&Waveform::T, py::const_))
      .def("T", py::overload_cast<const unsigned int>(&Waveform::T, py::const_), py::arg("TimeIndex"))
      .def("Frame", py::overload_cast<>(&Waveform::Frame, py::const_))
      .def("Frame", py::overload_cast<const unsigned int>(&Waveform::Frame, py::const_), py::arg("TimeIndex"))
      .def("LM", py::overload_cast<>(&Waveform::LM, py::const_))
      .def("LM", py::overload_cast<const unsigned int>(&Waveform::LM, py::const_), py::arg("Mode"))
      .def("Re", py::overload_cast<>(&Waveform::Re, py::const_), Doc("Waveform::Re"))
      .def("Re", py::overload_cast<const unsigned int>(&Waveform::Re, py::const_), py::arg("Mode"))
      .def("Re", py::overload_cast<const unsigned int, const unsigned int>(&Waveform::Re, py::const_),
           py::arg("Mode"), py::arg("TimeIndex"))
      .def("Im", py::overload_cast<>(&Waveform::Im, py::const_), Doc("Waveform::Im"))
      .def("Im", py::overload_cast<const unsigned int>(&Waveform::Im, py::const_), py::arg("Mode"))
      .def("Im", py::overload_cast<const unsigned int, const unsigned int>(&Waveform::Im, py::const_),
           py::arg("Mode"), py::arg("TimeIndex"))
      .def("Abs", py::overload_cast<>(&Waveform::Abs, py::const_), Doc("Waveform::Abs"))
      .def("Abs", py::overload_cast<const unsigned int>(&Waveform::Abs, py::const_), py::arg("Mode"))
      .def("Abs", py::overload_cast<const unsigned int, const unsigned int>(&Waveform::Abs, py::const_),
           py::arg("Mode"), py::arg("TimeIndex"))
      .def("Arg", py::overload_cast<>(&Waveform::Arg, py::const_), Doc("Waveform::Arg"))
      .def("Arg", py::overload_cast<const unsigned int>(&Waveform::Arg, py::const_), py::arg("Mode"))
      .def("Arg", py::overload_cast<const unsigned int, const unsigned int>(&Waveform::Arg, py::const_),
           py::arg("Mode"), py::arg("TimeIndex"))
      .def("ArgUnwrapped", py::overload_cast<>(&Waveform::ArgUnwrapped, py::const_),
           Doc("Waveform::ArgUnwrapped"))
      .def("ArgUnwrapped", py::overload_cast<const unsigned int>(&Waveform::ArgUnwrapped, py::const_),
           py::arg("Mode"))
      .def("Data", py::overload_cast<>(&Waveform::Data, py::const_), Doc("Waveform::Data"))
      .def("Data", py::overload_cast<const unsigned int>(&Waveform::Data, py::const_), py::arg("Mode"))
      .def("Data", py::overload_cast<const unsigned int, const unsigned int>(&Waveform::Data, py::const_),
           py::arg("Mode"), py::arg("TimeIndex"))

      // Data characterization
      .def("EllMax", &Waveform::EllMax, Doc("Waveform::EllMax"))
      .def("FindModeIndex", &Waveform::FindModeIndex, py::arg("L"), py::arg("M"), Doc("Waveform::FindModeIndex"))
      .def("FindModeIndexWithoutError", &Waveform::FindModeIndexWithoutError, py::arg("L"), py::arg("M"),
           Doc("Waveform::FindModeIndexWithoutError"))
      .def("DataDot", &Waveform::DataDot, py::arg("Mode"), Doc("Waveform::DataDot"))
      .def("Differentiate", ReturnSelf(&Waveform::Differentiate), Doc("Waveform::Differentiate"))
      .def("Norm", &Waveform::Norm, py::arg("TakeSquareRoot") = false, Doc("Waveform::Norm"))
      .def("MaxNormIndex", &Waveform::MaxNormIndex, py::arg("SkipFraction") = 4, Doc("Waveform::MaxNormIndex"))
      .def("MaxNormTime", &Waveform::MaxNormTime, py::arg("SkipFraction") = 4, Doc("Waveform::MaxNormTime"))

      // Involutions
      .def("NormalizedAntisymmetry", &Waveform::NormalizedAntisymmetry, py::arg("LModesForAsymmetry") = noL,
           Doc("Waveform::NormalizedAntisymmetry"))
      .def("DipoleMoment", &Waveform::DipoleMoment, py::arg("ellMax") = 0, Doc("Waveform::DipoleMoment"))
      .def("MinimalParityViolation", &Waveform::MinimalParityViolation, Doc("Waveform::MinimalParityViolation"))
      .def("XParityInvolution", &Waveform::XParityInvolution)
      .def("XParitySymmetricPart", &Waveform::XParitySymmetricPart)
      .def("XParityAntisymmetricPart", &Waveform::XParityAntisymmetricPart)
      .def("XParityViolationSquared", &Waveform::XParityViolationSquared, py::arg("Lmodes") = noL)
      .def("XParityViolationNormalized", &Waveform::XParityViolationNormalized, py::arg("Lmodes") = noL)
      .def("YParityInvolution", &Waveform::YParityInvolution)
      .def("YParitySymmetricPart", &Waveform::YParitySymmetricPart)
      .def("YParityAntisymmetricPart", &Waveform::YParityAntisymmetricPart)
      .def("YParityViolationSquared", &Waveform::YParityViolationSquared, py::arg("Lmodes") = noL)
      .def("YParityViolationNormalized", &Waveform::YParityViolationNormalized, py::arg("Lmodes") = noL)
      .def("ZParityInvolution", &Waveform::ZParityInvolution)
      .def("ZParitySymmetricPart", &Waveform::ZParitySymmetricPart)
      .def("ZParityAntisymmetricPart", &Waveform::ZParityAntisymmetricPart)
      .def("ZParityViolationSquared", &Waveform::ZParityViolationSquared, py::arg("Lmodes") = noL)
      .def("ZParityViolationNormalized", &Waveform::ZParityViolationNormalized, py::arg("Lmodes") = noL)
      .def("ParityInvolution", &Waveform::ParityInvolution)
      .def("ParitySymmetricPart", &Waveform::ParitySymmetricPart)
      .def("ParityAntisymmetricPart", &Waveform::ParityAntisymmetricPart)
      .def("ParityViolationSquared", &Waveform::ParityViolationSquared, py::arg("Lmodes") = noL)
      .def("ParityViolationNormalized", &Waveform::ParityViolationNormalized, py::arg("Lmodes") = noL)
      .def("ConjugateAntipodalEvaluation", py::overload_cast<>(&Waveform::ConjugateAntipodalEvaluation, py::const_))
      .def("ConjugateAntipodalSymmetricPart", &Waveform::ConjugateAntipodalSymmetricPart)
      .def("ConjugateAntipodalAntisymmetricPart", &Waveform::ConjugateAntipodalAntisymmetricPart)

      // Rotations and frames
      .def("RotatePhysicalSystem", ReturnSelf(py::overload_cast<const Quaternion&>(&Waveform::RotatePhysicalSystem)),
           py::arg("R_phys"), Doc("Waveform::RotatePhysicalSystem"))
      .def("RotatePhysicalSystem", ReturnSelf(py::overload_cast<QVec>(&Waveform::RotatePhysicalSystem)), py::arg("R_phys"))
      .def("RotateDecompositionBasis", ReturnSelf(py::overload_cast<const Quaternion&>(&Waveform::RotateDecompositionBasis)),
           py::arg("R_frame"), Doc("Waveform::RotateDecompositionBasis"))
      .def("RotateDecompositionBasis", ReturnSelf(py::overload_cast<const QVec&>(&Waveform::RotateDecompositionBasis)),
           py::arg("R_frame"))
      .def("LdtVector", &Waveform::LdtVector, py::arg("Lmodes") = noL, Doc("Waveform::LdtVector"))
      .def("LLMatrix", &Waveform::LLMatrix, py::arg("Lmodes") = noL, Doc("Waveform::LLMatrix"))
      .def("LLDominantEigenvector", &Waveform::LLDominantEigenvector, py::arg("Lmodes") = noL,
           py::arg("RoughInitialEllDirection") = Quaternions::zHat, Doc("Waveform::LLDominantEigenvector"))
      .def("AngularVelocityVector", &Waveform::AngularVelocityVector, py::arg("Lmodes") = noL,
           Doc("Waveform::AngularVelocityVector"))
      .def("AngularVelocityVectorRelativeToInertial", &Waveform::AngularVelocityVectorRelativeToInertial,
           py::arg("Lmodes") = noL, Doc("Waveform::AngularVelocityVectorRelativeToInertial"))
      .def("CorotatingFrame", &Waveform::CorotatingFrame, py::arg("Lmodes") = noL, Doc("Waveform::CorotatingFrame"))
      .def("TransformToCoprecessingFrame", ReturnSelf(&Waveform::TransformToCoprecessingFrame), py::arg("Lmodes") = noL, Doc("Waveform::TransformToCoprecessingFrame"))
      .def("TransformToAngularVelocityFrame", ReturnSelf(&Waveform::TransformToAngularVelocityFrame), py::arg("Lmodes") = noL, Doc("Waveform::TransformToAngularVelocityFrame"))
      .def("TransformToCorotatingFrame", ReturnSelf(&Waveform::TransformToCorotatingFrame), py::arg("Lmodes") = noL, Doc("Waveform::TransformToCorotatingFrame"))
      .def("TransformToInertialFrame", ReturnSelf(&Waveform::TransformToInertialFrame),
           Doc("Waveform::TransformToInertialFrame"))

      // Alignment, comparison, and hybridization
      .def("GetAlignmentsOfDecompositionFrameToModes", &Waveform::GetAlignmentsOfDecompositionFrameToModes,
           py::arg("Lmodes") = noL, Doc("Waveform::GetAlignmentsOfDecompositionFrameToModes"))
      .def("GetAlignmentOfDecompositionFrameToModes", &Waveform::GetAlignmentOfDecompositionFrameToModes,
           py::arg("t_fid"), py::arg("nHat_t_fid") = Quaternions::xHat, py::arg("Lmodes") = noL,
           Doc("Waveform::GetAlignmentOfDecompositionFrameToModes"))
      .def("AlignDecompositionFrameToModes", ReturnSelf(&Waveform::AlignDecompositionFrameToModes), py::arg("t_fid"),
           py::arg("nHat_t_fid") = Quaternions::xHat, py::arg("Lmodes") = noL,
           Doc("Waveform::AlignDecompositionFrameToModes"))
      .def("Compare", &Waveform::Compare, py::arg("B"), py::arg("MinTimeStep") = 0.005,
           py::arg("MinTime") = -3.0e300, Doc("Waveform::Compare"))
      .def("Hybridize", &Waveform::Hybridize, py::arg("B"), py::arg("t1"), py::arg("t2"),
           py::arg("tMinStep") = 0.005, Doc("Waveform::Hybridize"))

      // Pointwise operations
      .def("EvaluateAtPoint", &Waveform::EvaluateAtPoint, py::arg("vartheta"), py::arg("varphi"),
           py::arg("i_0") = 0, py::arg("i_1") = -1, Doc("Waveform::EvaluateAtPoint"))
      .def("InterpolateToPoint", [](const Waveform& w, const double vartheta, const double varphi,
                                    const double t_i) { return w.InterpolateToPoint(vartheta, varphi, t_i); },
           py::arg("vartheta"), py::arg("varphi"), py::arg("t_i"), Doc("Waveform::InterpolateToPoint"))
      .def(py::self + py::self)
      .def(py::self - py::self)
      .def(py::self * py::self)
      .def(py::self / py::self)
      .def(py::self * double())
      .def(py::self / double())
      .def(double() * py::self)
      .def("Translate", &Waveform::Translate, py::arg("deltax"), Doc("Waveform::Translate"))
      .def("BoostPsi4", ReturnSelf(&Waveform::BoostPsi4), py::arg("v"), Doc("Waveform::BoostPsi4"))
      .def("BoostHFaked", ReturnSelf(&Waveform::BoostHFaked), py::arg("v"), Doc("Waveform::BoostHFaked"))
      .def("Output", ReturnSelf(&Waveform::Output), py::arg("FileName"), py::arg("precision") = 14,
           Doc("Waveform::Output"));

  m.def("AlignWaveforms", &GWFrames::AlignWaveforms, py::arg("A"), py::arg("B"), py::arg("t_1"), py::arg("t_2"),
        py::arg("InitialEvaluations") = 0, py::arg("nHat_A") = DVec(), py::arg("Debug") = false,
        Doc("AlignWaveforms"));

  py::class_<PNWaveform, Waveform>(m, "PNWaveform", Doc("PNWaveform"))
      .def(py::init<>())
      .def(py::init<const PNWaveform&>())
      .def(py::init<const std::string&, const double, const DVec&, const DVec&, const double, double,
                    const Quaternion&, const unsigned int, const double, const double>(),
           py::arg("Approximant"), py::arg("delta"), py::arg("chi1_i"), py::arg("chi2_i"), py::arg("Omega_orb_i"),
           py::arg("Omega_orb_0") = -1.0, py::arg("R_frame_i") = Quaternion(1, 0, 0, 0),
           py::arg("MinStepsPerOrbit") = 32, py::arg("PNWaveformModeOrder") = 3.5,
           py::arg("PNOrbitalEvolutionOrder") = 4.0)
      .def("__copy__", [](const PNWaveform& w) { return PNWaveform(w); })
      .def("__deepcopy__", [](const PNWaveform& w, py::dict) { return PNWaveform(w); }, py::arg("memo"))
#define GWFRAMES_PN_VECTOR(name)                                                                     \
  .def(#name, py::overload_cast<>(&PNWaveform::name, py::const_), Doc("PNWaveform::" #name))         \
      .def(#name, py::overload_cast<const unsigned int>(&PNWaveform::name, py::const_), py::arg("iTime"))
      GWFRAMES_PN_VECTOR(chi1)
      GWFRAMES_PN_VECTOR(chi2)
      GWFRAMES_PN_VECTOR(Omega_orb)
      GWFRAMES_PN_VECTOR(Omega_prec)
      GWFRAMES_PN_VECTOR(Omega_tot)
      GWFRAMES_PN_VECTOR(L)
      GWFRAMES_PN_VECTOR(Omega_orbMag)
      GWFRAMES_PN_VECTOR(Omega_precMag)
      GWFRAMES_PN_VECTOR(Omega_totMag)
      GWFRAMES_PN_VECTOR(LMag)
      GWFRAMES_PN_VECTOR(chiHat1)
      GWFRAMES_PN_VECTOR(chiHat2)
      GWFRAMES_PN_VECTOR(OmegaHat_orb)
      GWFRAMES_PN_VECTOR(OmegaHat_prec)
      GWFRAMES_PN_VECTOR(OmegaHat_tot)
      GWFRAMES_PN_VECTOR(LHat)
      GWFRAMES_PN_VECTOR(Phi_orb)
#undef GWFRAMES_PN_VECTOR
      .def("chi1Mag", &PNWaveform::chi1Mag, py::arg("iTime"))
      .def("chi2Mag", &PNWaveform::chi2Mag, py::arg("iTime"));

  py::class_<WaveformAtAPointFT>(m, "WaveformAtAPointFT", Doc("WaveformAtAPointFT"))
      .def(py::init<const Waveform&, const double, const double, const double, const double, const unsigned int,
                    const double, const double, const unsigned int>(),
           py::arg("W"), py::arg("Dt"), py::arg("Vartheta"), py::arg("Varphi"), py::arg("TotalMass"),
           py::arg("WindowNCycles") = 1, py::arg("DetectorResponseAmp") = 1.0,
           py::arg("DetectorResponsePhase") = 0.0, py::arg("ExtraZeroPadPowers") = 0)
      .def("F", py::overload_cast<>(&WaveformAtAPointFT::F, py::const_), Doc("WaveformAtAPointFT::F"))
      .def("F", py::overload_cast<const unsigned int>(&WaveformAtAPointFT::F, py::const_), py::arg("f"))
      .def("NFreq", &WaveformAtAPointFT::NFreq)
      .def("IsNormalized", &WaveformAtAPointFT::IsNormalized)
      .def("Dt", &WaveformAtAPointFT::Dt)
      .def("Vartheta", &WaveformAtAPointFT::Vartheta)
      .def("Varphi", &WaveformAtAPointFT::Varphi)
      .def("Re", py::overload_cast<>(&WaveformAtAPointFT::Re, py::const_))
      .def("Re", py::overload_cast<const unsigned int>(&WaveformAtAPointFT::Re, py::const_), py::arg("f"))
      .def("Im", py::overload_cast<>(&WaveformAtAPointFT::Im, py::const_))
      .def("Im", py::overload_cast<const unsigned int>(&WaveformAtAPointFT::Im, py::const_), py::arg("f"))
      .def("InversePSD", &WaveformAtAPointFT::InversePSD, py::arg("Detector") = "AdvLIGO_ZeroDet_HighP",
           Doc("WaveformAtAPointFT::InversePSD"))
      .def("SNR", py::overload_cast<const DVec&>(&WaveformAtAPointFT::SNR, py::const_), py::arg("InversePSD"),
           Doc("WaveformAtAPointFT::SNR"))
      .def("SNR", py::overload_cast<const std::string&>(&WaveformAtAPointFT::SNR, py::const_),
           py::arg("Detector") = "AdvLIGO_ZeroDet_HighP")
      .def("InnerProduct", &WaveformAtAPointFT::InnerProduct, py::arg("B"), py::arg("InversePSD"),
           Doc("WaveformAtAPointFT::InnerProduct"))
      .def("Match",
           py::overload_cast<const WaveformAtAPointFT&, const DVec&>(&WaveformAtAPointFT::Match, py::const_),
           py::arg("B"), py::arg("InversePSD"), Doc("WaveformAtAPointFT::Match"))
      .def("Match",
           py::overload_cast<const WaveformAtAPointFT&, const std::string&>(&WaveformAtAPointFT::Match,
                                                                           py::const_),
           py::arg("B"), py::arg("Detector") = "AdvLIGO_ZeroDet_HighP")
      .def("Match",
           [](const WaveformAtAPointFT& A, const WaveformAtAPointFT& B, const DVec& InversePSD, double timeOffset,
              double phaseOffset, double match) {
             A.Match(B, InversePSD, timeOffset, phaseOffset, match);
             return py::make_tuple(timeOffset, phaseOffset, match);
           },
           py::arg("B"), py::arg("InversePSD"), py::arg("timeOffset"), py::arg("phaseOffset"), py::arg("match"))
      .def("Match",
           [](const WaveformAtAPointFT& A, const WaveformAtAPointFT& B, double timeOffset, double phaseOffset,
              double match, const std::string& Detector) {
             A.Match(B, timeOffset, phaseOffset, match, Detector);
             return py::make_tuple(timeOffset, phaseOffset, match);
           },
           py::arg("B"), py::arg("timeOffset"), py::arg("phaseOffset"), py::arg("match"),
           py::arg("Detector") = "AdvLIGO_ZeroDet_HighP")
      .def("Normalize", ReturnSelf(py::overload_cast<const DVec&>(&WaveformAtAPointFT::Normalize)), py::arg("InversePSD"), Doc("WaveformAtAPointFT::Normalize"))
      .def("Normalize", ReturnSelf(py::overload_cast<const std::string&>(&WaveformAtAPointFT::Normalize)),
           py::arg("Detector") = "AdvLIGO_ZeroDet_HighP")
      .def("ZeroAbove", ReturnSelf(&WaveformAtAPointFT::ZeroAbove), py::arg("Frequency"),
           Doc("WaveformAtAPointFT::ZeroAbove"));
}

}  // namespace GWFramesBindings
