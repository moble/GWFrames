// Copyright (c) 2026, Michael Boyle
// See LICENSE file for details

#include "bindings.hpp"

#include <iomanip>
#include <sstream>

#include "IntegrateAngularVelocity.hpp"
#include "QuaternionUtilities.hpp"
#include "Quaternions.hpp"

namespace GWFramesBindings {

using Quaternions::Quaternion;
using QVec = std::vector<Quaternion>;
using DVec = std::vector<double>;

namespace {

std::string Components(const Quaternion& q) {
  std::ostringstream s;
  s << std::setprecision(15) << q[0] << ", " << q[1] << ", " << q[2] << ", " << q[3];
  return s.str();
}

unsigned int ComponentIndex(long i) {
  if (i < 0) {
    i += 4;
  }
  if (i < 0 || i > 3) {
    throw py::index_error("Quaternion index out of range");
  }
  return static_cast<unsigned int>(i);
}

}  // namespace

void BindQuaternions(py::module_& m) {
  py::class_<Quaternion>(m, "Quaternion", Doc("Quaternion"))
      .def(py::init<>())
      .def(py::init<const Quaternion&>(), py::arg("Q"))
      .def(py::init<const double, const double>(), py::arg("vartheta"), py::arg("varphi"))
      .def(py::init<const double, const double, const double>(), py::arg("alpha"), py::arg("beta"),
           py::arg("gamma"))
      .def(py::init<const double, const double, const double, const double>(), py::arg("w"), py::arg("x"),
           py::arg("y"), py::arg("z"))
      .def(py::init<const DVec&>(), py::arg("q"))
      .def(py::init<const double, const DVec&>(), py::arg("angle"), py::arg("axis"))

      // Sequence-like access
      .def("__len__", [](const Quaternion&) { return 4; })
      .def("__getitem__", [](const Quaternion& q, long i) { return q[ComponentIndex(i)]; }, py::arg("i"))
      .def("__setitem__", [](Quaternion& q, long i, double a) { q[ComponentIndex(i)] = a; }, py::arg("i"),
           py::arg("value"))
      .def("__str__", [](const Quaternion& q) { return "[" + Components(q) + "]"; })
      .def("__repr__", [](const Quaternion& q) { return "Quaternion(" + Components(q) + ")"; })
      .def("__copy__", [](const Quaternion& q) { return Quaternion(q); })
      .def("__deepcopy__", [](const Quaternion& q, py::dict) { return Quaternion(q); }, py::arg("memo"))
      .def(py::pickle([](const Quaternion& q) { return py::make_tuple(q[0], q[1], q[2], q[3]); },
                      [](py::tuple t) {
                        return Quaternion(t[0].cast<double>(), t[1].cast<double>(), t[2].cast<double>(),
                                          t[3].cast<double>());
                      }))

      // Arithmetic
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def(-py::self)
      .def(py::self + py::self)
      .def(py::self - py::self)
      .def(py::self * py::self)
      .def(py::self / py::self)
      .def(py::self + double())
      .def(py::self - double())
      .def(py::self * double())
      .def(py::self / double())
      .def(double() + py::self)
      .def(double() - py::self)
      .def(double() * py::self)
      .def(double() / py::self)
      .def("__pow__", [](const Quaternion& q, const double t) { return q.pow(t); }, py::arg("t"))
      .def("__pow__", [](const Quaternion& q, const Quaternion& p) { return q.pow(p); }, py::arg("Q"))

      // Methods
      .def("pow", py::overload_cast<const double>(&Quaternion::pow, py::const_), py::arg("t"))
      .def("pow", py::overload_cast<const Quaternion&>(&Quaternion::pow, py::const_), py::arg("Q"),
           Doc("Quaternion::pow"))
      .def("cross", &Quaternion::cross, py::arg("Q"), Doc("Quaternion::cross"))
      .def("commutator", &Quaternion::commutator, py::arg("Q"), Doc("Quaternion::commutator"))
      .def("dot", &Quaternion::dot, py::arg("Q"), Doc("Quaternion::dot"))
      .def("abs", &Quaternion::abs, Doc("Quaternion::abs"))
      .def("log", &Quaternion::log, Doc("Quaternion::log"))
      .def("logRotor", &Quaternion::logRotor, Doc("Quaternion::logRotor"))
      .def("exp", &Quaternion::exp, Doc("Quaternion::exp"))
      .def("sqrt", &Quaternion::sqrt, Doc("Quaternion::sqrt"))
      .def("angle", &Quaternion::angle, Doc("Quaternion::angle"))
      .def("inverse", &Quaternion::inverse, Doc("Quaternion::inverse"))
      .def("conjugate", &Quaternion::conjugate, Doc("Quaternion::conjugate"))
      .def("normalized", &Quaternion::normalized, Doc("Quaternion::normalized"))
      .def("normsquared", &Quaternion::normsquared, Doc("Quaternion::normsquared"))
      .def("sqrtOfRotor", &Quaternion::sqrtOfRotor, Doc("Quaternion::sqrtOfRotor"))
      .def("vec", &Quaternion::vec, Doc("Quaternion::vec"))
      .def("XParityConjugateSpinor", &Quaternion::XParityConjugateSpinor, Doc("Quaternion::XParityConjugateSpinor"))
      .def("YParityConjugateSpinor", &Quaternion::YParityConjugateSpinor, Doc("Quaternion::YParityConjugateSpinor"))
      .def("ZParityConjugateSpinor", &Quaternion::ZParityConjugateSpinor, Doc("Quaternion::ZParityConjugateSpinor"))
      .def("ParityConjugateSpinor", &Quaternion::ParityConjugateSpinor, Doc("Quaternion::ParityConjugateSpinor"))
      .def("IntrinsicDistance", &Quaternion::IntrinsicDistance, py::arg("Q"), Doc("Quaternion::IntrinsicDistance"))
      .def("ChordalDistance", &Quaternion::ChordalDistance, py::arg("Q"), Doc("Quaternion::ChordalDistance"))
      .def("str", &Quaternion::str, Doc("Quaternion::str"));

  // Lists and tuples of numbers are accepted wherever a Quaternion is.
  py::implicitly_convertible<DVec, Quaternion>();

  m.attr("Zero") = Quaternions::Zero;
  m.attr("One") = Quaternions::One;
  m.attr("xHat") = Quaternions::xHat;
  m.attr("yHat") = Quaternions::yHat;
  m.attr("zHat") = Quaternions::zHat;

  // Functions of single quaternions, each with a version for arrays
  m.def("pow", py::overload_cast<const Quaternion&, const double>(&Quaternions::pow), py::arg("Q"), py::arg("x"));
  m.def("pow", py::overload_cast<const Quaternion&, const Quaternion&>(&Quaternions::pow), py::arg("Q"), py::arg("P"));
  m.def("pow", py::overload_cast<const QVec&, const double>(&Quaternions::pow), py::arg("Q"), py::arg("x"));
  m.def("pow", py::overload_cast<const QVec&, const Quaternion&>(&Quaternions::pow), py::arg("Q"), py::arg("P"));
  m.def("pow", py::overload_cast<const Quaternion&, const DVec&>(&Quaternions::pow), py::arg("Q"), py::arg("x"));
  m.def("pow", py::overload_cast<const Quaternion&, const QVec&>(&Quaternions::pow), py::arg("Q"), py::arg("P"));
  m.def("pow", py::overload_cast<const QVec&, const DVec&>(&Quaternions::pow), py::arg("Q"), py::arg("x"));
  m.def("pow", py::overload_cast<const QVec&, const QVec&>(&Quaternions::pow), py::arg("Q"), py::arg("P"), Doc("pow"));
  m.def("cross", py::overload_cast<const Quaternion&, const Quaternion&>(&Quaternions::cross), py::arg("Q"),
        py::arg("P"), Doc("cross"));
  m.def("commutator", py::overload_cast<const Quaternion&, const Quaternion&>(&Quaternions::commutator), py::arg("Q"),
        py::arg("P"), Doc("commutator"));
  m.def("dot", py::overload_cast<const Quaternion&, const Quaternion&>(&Quaternions::dot), py::arg("Q"), py::arg("P"),
        Doc("dot"));
  m.def("abs", py::overload_cast<const Quaternion&>(&Quaternions::abs), py::arg("Q"));
  m.def("abs", py::overload_cast<const QVec&>(&Quaternions::abs), py::arg("Q"), Doc("abs"));
  m.def("log", py::overload_cast<const Quaternion&>(&Quaternions::log), py::arg("Q"));
  m.def("log", py::overload_cast<const QVec&>(&Quaternions::log), py::arg("Q"), Doc("log"));
  m.def("logRotor", py::overload_cast<const Quaternion&>(&Quaternions::logRotor), py::arg("Q"));
  m.def("logRotor", py::overload_cast<const QVec&>(&Quaternions::logRotor), py::arg("Q"), Doc("logRotor"));
  m.def("exp", py::overload_cast<const Quaternion&>(&Quaternions::exp), py::arg("Q"));
  m.def("exp", py::overload_cast<const QVec&>(&Quaternions::exp), py::arg("Q"), Doc("exp"));
  m.def("sqrt", py::overload_cast<const Quaternion&>(&Quaternions::sqrt), py::arg("Q"));
  m.def("sqrt", py::overload_cast<const QVec&>(&Quaternions::sqrt), py::arg("Q"), Doc("sqrt"));
  m.def("angle", py::overload_cast<const Quaternion&>(&Quaternions::angle), py::arg("Q"));
  m.def("angle", py::overload_cast<const QVec&>(&Quaternions::angle), py::arg("Q"), Doc("angle"));
  m.def("inverse", py::overload_cast<const Quaternion&>(&Quaternions::inverse), py::arg("Q"));
  m.def("inverse", py::overload_cast<const QVec&>(&Quaternions::inverse), py::arg("Q"), Doc("inverse"));
  m.def("conjugate", py::overload_cast<const Quaternion&>(&Quaternions::conjugate), py::arg("Q"));
  m.def("conjugate", py::overload_cast<const QVec&>(&Quaternions::conjugate), py::arg("Q"), Doc("conjugate"));
  m.def("normalized", py::overload_cast<const Quaternion&>(&Quaternions::normalized), py::arg("Q"));
  m.def("normalized", py::overload_cast<const QVec&>(&Quaternions::normalized), py::arg("Q"), Doc("normalized"));
  m.def("normsquared", py::overload_cast<const Quaternion&>(&Quaternions::normsquared), py::arg("Q"));
  m.def("normsquared", py::overload_cast<const QVec&>(&Quaternions::normsquared), py::arg("Q"), Doc("normsquared"));
  m.def("sqrtOfRotor", py::overload_cast<const Quaternion&>(&Quaternions::sqrtOfRotor), py::arg("Q"));
  m.def("sqrtOfRotor", py::overload_cast<const QVec&>(&Quaternions::sqrtOfRotor), py::arg("Q"), Doc("sqrtOfRotor"));
  m.def("vec", py::overload_cast<const Quaternion&>(&Quaternions::vec), py::arg("Q"));
  m.def("vec", py::overload_cast<const QVec&>(&Quaternions::vec), py::arg("Q"), Doc("vec"));
  m.def("XParityConjugateSpinor", py::overload_cast<const Quaternion&>(&Quaternions::XParityConjugateSpinor),
        py::arg("Q"), Doc("XParityConjugateSpinor"));
  m.def("YParityConjugateSpinor", py::overload_cast<const Quaternion&>(&Quaternions::YParityConjugateSpinor),
        py::arg("Q"), Doc("YParityConjugateSpinor"));
  m.def("ZParityConjugateSpinor", py::overload_cast<const Quaternion&>(&Quaternions::ZParityConjugateSpinor),
        py::arg("Q"), Doc("ZParityConjugateSpinor"));
  m.def("ParityConjugateSpinor", py::overload_cast<const Quaternion&>(&Quaternions::ParityConjugateSpinor),
        py::arg("Q"), Doc("ParityConjugateSpinor"));
  m.def("IntrinsicDistance", py::overload_cast<const Quaternion&, const Quaternion&>(&Quaternions::IntrinsicDistance),
        py::arg("R_A"), py::arg("R_B"), Doc("IntrinsicDistance"));
  m.def("ChordalDistance", py::overload_cast<const Quaternion&, const Quaternion&>(&Quaternions::ChordalDistance),
        py::arg("R_A"), py::arg("R_B"), Doc("ChordalDistance"));
  m.def("Slerp", &Quaternions::Slerp, py::arg("tau"), py::arg("Qa"), py::arg("Qb"), Doc("Slerp"));

  // Construction of arrays
  m.def("QuaternionArray", py::overload_cast<const DVec&, const DVec&>(&Quaternions::QuaternionArray), py::arg("vartheta"),
        py::arg("varphi"));
  m.def("QuaternionArray", py::overload_cast<const DVec&, const DVec&, const DVec&>(&Quaternions::QuaternionArray),
        py::arg("alpha"), py::arg("beta"), py::arg("gamma"));
  m.def("QuaternionArray", py::overload_cast<const DVec&, const DVec&, const DVec&, const DVec&>(&Quaternions::QuaternionArray),
        py::arg("w"), py::arg("x"), py::arg("y"), py::arg("z"));
  m.def("QuaternionArray", py::overload_cast<const std::vector<DVec>&>(&Quaternions::QuaternionArray), py::arg("q"));
  m.def("QuaternionArray", py::overload_cast<const DVec&, const std::vector<DVec>&>(&Quaternions::QuaternionArray),
        py::arg("angle"), py::arg("axis"), Doc("QuaternionArray"));
  m.def("Component", &Quaternions::Component, py::arg("Q"), py::arg("i"), Doc("Component"));
  m.def("Component0", &Quaternions::Component0, py::arg("Q"), Doc("Component0"));
  m.def("Component1", &Quaternions::Component1, py::arg("Q"), Doc("Component1"));
  m.def("Component2", &Quaternions::Component2, py::arg("Q"), Doc("Component2"));
  m.def("Component3", &Quaternions::Component3, py::arg("Q"), Doc("Component3"));

  // Calculus and frames
  m.def("IndefiniteIntegral", &Quaternions::IndefiniteIntegral, py::arg("RIn"), py::arg("tIn"), Doc("IndefiniteIntegral"));
  m.def("DefiniteIntegral", &Quaternions::DefiniteIntegral, py::arg("RIn"), py::arg("tIn"), Doc("DefiniteIntegral"));
  m.def("DifferentiateRotorByLogarithm", &Quaternions::DifferentiateRotorByLogarithm, py::arg("RIn"), py::arg("tIn"),
        Doc("DifferentiateRotorByLogarithm"));
  m.def("QuaternionDerivative", &Quaternions::QuaternionDerivative, py::arg("f"), py::arg("t"), Doc("QuaternionDerivative"));
  m.def("MinimalRotation", &Quaternions::MinimalRotation, py::arg("R"), py::arg("T"), py::arg("NIterations") = 5,
        Doc("MinimalRotation"));
  m.def("PrescribedRotation", &Quaternions::PrescribedRotation, py::arg("RotationRateAboutZ"), py::arg("R"), py::arg("T"),
        py::arg("NIterations") = 5, Doc("PrescribedRotation"));
  m.def("FrameFromXY", &Quaternions::FrameFromXY, py::arg("X"), py::arg("Y"), Doc("FrameFromXY"));
  m.def("FrameFromZ", &Quaternions::FrameFromZ, py::arg("Z"), py::arg("T"), py::arg("NIterations") = 5, Doc("FrameFromZ"));
  m.def("FrameFromPrescribedRotation", &Quaternions::FrameFromPrescribedRotation, py::arg("omega"), py::arg("T"),
        py::arg("NIterations") = 5, Doc("FrameFromPrescribedRotation"));
  m.def("FrameFromAngularVelocity",
        py::overload_cast<const QVec&, const DVec&>(&Quaternions::FrameFromAngularVelocity), py::arg("Omega"), py::arg("T"),
        Doc("FrameFromAngularVelocity"));
  m.def("FrameFromAngularVelocity_2D",
        py::overload_cast<const QVec&, const DVec&>(&Quaternions::FrameFromAngularVelocity_2D), py::arg("Omega"),
        py::arg("T"), Doc("FrameFromAngularVelocity_2D"));
  m.def("FrameFromAngularVelocity_Integrand", &Quaternions::FrameFromAngularVelocity_Integrand, py::arg("r"), py::arg("Omega"),
        Doc("FrameFromAngularVelocity_Integrand"));
  m.def(
      "FrameFromAngularVelocity_2D_Integrand",
      [](const double r_x, const double r_y, const DVec& Omega) {
        double rdot_x = 0.0, rdot_y = 0.0;
        Quaternions::FrameFromAngularVelocity_2D_Integrand(r_x, r_y, Omega, rdot_x, rdot_y);
        return py::make_tuple(rdot_x, rdot_y);
      },
      py::arg("r_x"), py::arg("r_y"), py::arg("Omega"),
      "Return the time derivative of a two-dimensional rotor logarithm, as the tuple ``(rdot_x, rdot_y)``.");
  m.def("UnflipRotors", &Quaternions::UnflipRotors, py::arg("R"), py::arg("discont") = 1.4142135623730951, Doc("UnflipRotors"));
  m.def("RDelta", &Quaternions::RDelta, py::arg("R1"), py::arg("R2"), py::arg("IndexOfFiducialTime") = -1, Doc("RDelta"));
  m.def("Squad", &Quaternions::Squad, py::arg("RIn"), py::arg("tIn"), py::arg("tOut"), Doc("Squad"));
  m.def("FrameAngularVelocity", &Quaternions::FrameAngularVelocity, py::arg("f"), py::arg("t"), Doc("FrameAngularVelocity"));
  m.def("BoostRotor", &Quaternions::BoostRotor, py::arg("ThreeVelocity"), py::arg("ThreeVectorToBeBoosted"), Doc("BoostRotor"));

  // Mean rotors and alignment
  m.def("ApproximateMeanRotor", py::overload_cast<const QVec&>(&Quaternions::ApproximateMeanRotor), py::arg("R"));
  m.def("ApproximateMeanRotor",
        py::overload_cast<const QVec&, const DVec&, const double, const double>(&Quaternions::ApproximateMeanRotor),
        py::arg("R"), py::arg("t"), py::arg("t1") = -1e300, py::arg("t2") = 1e300, Doc("ApproximateMeanRotor"));
  m.def("ApproximateOptimalAlignmentRotor", &Quaternions::ApproximateOptimalAlignmentRotor, py::arg("Ra"), py::arg("Rb"),
        py::arg("t"), py::arg("t1") = -1e300, py::arg("t2") = 1e300, Doc("ApproximateOptimalAlignmentRotor"));
  m.def(
      "ApproximateOptimalAlignment",
      [](const double t1, const double t2, const QVec& Ra, const DVec& ta, const QVec& Rb, const DVec& tb) {
        double deltat = 0.0;
        Quaternion R_delta;
        Quaternions::ApproximateOptimalAlignment(t1, t2, Ra, ta, Rb, tb, deltat, R_delta);
        return py::make_tuple(deltat, R_delta);
      },
      py::arg("t1"), py::arg("t2"), py::arg("Ra"), py::arg("ta"), py::arg("Rb"), py::arg("tb"),
      Doc("ApproximateOptimalAlignment"));
  m.def(
      "OptimalAlignment",
      [](const double t1, const double t2, const QVec& Ra, const DVec& ta, const QVec& Rb, const DVec& tb) {
        double deltat = 0.0;
        Quaternion R_delta;
        Quaternions::OptimalAlignment(t1, t2, Ra, ta, Rb, tb, deltat, R_delta);
        return py::make_tuple(deltat, R_delta);
      },
      py::arg("t1"), py::arg("t2"), py::arg("Ra"), py::arg("ta"), py::arg("Rb"), py::arg("tb"),
      Doc("OptimalAlignment"));
  m.def("MeanRotor", py::overload_cast<const QVec&>(&Quaternions::MeanRotor), py::arg("R"));
  m.def("MeanRotor", py::overload_cast<const QVec&, const DVec&, double, double>(&Quaternions::MeanRotor), py::arg("R"),
        py::arg("t"), py::arg("t1") = -1e300, py::arg("t2") = 1e300, Doc("MeanRotor"));

  // Searching sorted arrays
  m.def("hunt", &Quaternions::hunt, py::arg("x"), py::arg("x_i"), py::arg("i") = 0, Doc("hunt"));
  m.def("huntRight", &Quaternions::huntRight, py::arg("x"), py::arg("x_i"), py::arg("i") = 0, Doc("huntRight"));
}

}  // namespace GWFramesBindings
