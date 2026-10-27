// Copyright (c) 2014, Michael Boyle
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
      .def(py::init<const Quaternion&>())
      .def(py::init<const double, const double>(), py::arg("vartheta"), py::arg("varphi"))
      .def(py::init<const double, const double, const double>(), py::arg("alpha"), py::arg("beta"),
           py::arg("gamma"))
      .def(py::init<const double, const double, const double, const double>(), py::arg("w"), py::arg("x"),
           py::arg("y"), py::arg("z"))
      .def(py::init<const DVec&>(), py::arg("q"))
      .def(py::init<const double, const DVec&>(), py::arg("angle"), py::arg("axis"))

      // Sequence-like access
      .def("__len__", [](const Quaternion&) { return 4; })
      .def("__getitem__", [](const Quaternion& q, long i) { return q[ComponentIndex(i)]; })
      .def("__setitem__", [](Quaternion& q, long i, double a) { q[ComponentIndex(i)] = a; })
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
      .def("__pow__", [](const Quaternion& q, const double t) { return q.pow(t); })
      .def("__pow__", [](const Quaternion& q, const Quaternion& p) { return q.pow(p); })

      // Methods
      .def("pow", py::overload_cast<const double>(&Quaternion::pow, py::const_), Doc("Quaternion::pow"))
      .def("pow", py::overload_cast<const Quaternion&>(&Quaternion::pow, py::const_))
      .def("cross", &Quaternion::cross, Doc("Quaternion::cross"))
      .def("commutator", &Quaternion::commutator, Doc("Quaternion::commutator"))
      .def("dot", &Quaternion::dot, Doc("Quaternion::dot"))
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
      .def("XParityConjugateSpinor", &Quaternion::XParityConjugateSpinor)
      .def("YParityConjugateSpinor", &Quaternion::YParityConjugateSpinor)
      .def("ZParityConjugateSpinor", &Quaternion::ZParityConjugateSpinor)
      .def("ParityConjugateSpinor", &Quaternion::ParityConjugateSpinor)
      .def("IntrinsicDistance", &Quaternion::IntrinsicDistance, Doc("Quaternion::IntrinsicDistance"))
      .def("ChordalDistance", &Quaternion::ChordalDistance, Doc("Quaternion::ChordalDistance"))
      .def("str", &Quaternion::str);

  // Lists and tuples of numbers are accepted wherever a Quaternion is.
  py::implicitly_convertible<DVec, Quaternion>();

  m.attr("Zero") = Quaternions::Zero;
  m.attr("One") = Quaternions::One;
  m.attr("xHat") = Quaternions::xHat;
  m.attr("yHat") = Quaternions::yHat;
  m.attr("zHat") = Quaternions::zHat;

  // Functions of single quaternions, each with a version for arrays
  m.def("pow", py::overload_cast<const Quaternion&, const double>(&Quaternions::pow), Doc("pow"));
  m.def("pow", py::overload_cast<const Quaternion&, const Quaternion&>(&Quaternions::pow));
  m.def("pow", py::overload_cast<const QVec&, const double>(&Quaternions::pow));
  m.def("pow", py::overload_cast<const QVec&, const Quaternion&>(&Quaternions::pow));
  m.def("pow", py::overload_cast<const Quaternion&, const DVec&>(&Quaternions::pow));
  m.def("pow", py::overload_cast<const Quaternion&, const QVec&>(&Quaternions::pow));
  m.def("pow", py::overload_cast<const QVec&, const DVec&>(&Quaternions::pow));
  m.def("pow", py::overload_cast<const QVec&, const QVec&>(&Quaternions::pow));
  m.def("cross", py::overload_cast<const Quaternion&, const Quaternion&>(&Quaternions::cross), Doc("cross"));
  m.def("commutator", py::overload_cast<const Quaternion&, const Quaternion&>(&Quaternions::commutator), Doc("commutator"));
  m.def("dot", py::overload_cast<const Quaternion&, const Quaternion&>(&Quaternions::dot), Doc("dot"));
  m.def("abs", py::overload_cast<const Quaternion&>(&Quaternions::abs), Doc("abs"));
  m.def("abs", py::overload_cast<const QVec&>(&Quaternions::abs));
  m.def("log", py::overload_cast<const Quaternion&>(&Quaternions::log), Doc("log"));
  m.def("log", py::overload_cast<const QVec&>(&Quaternions::log));
  m.def("logRotor", py::overload_cast<const Quaternion&>(&Quaternions::logRotor), Doc("logRotor"));
  m.def("logRotor", py::overload_cast<const QVec&>(&Quaternions::logRotor));
  m.def("exp", py::overload_cast<const Quaternion&>(&Quaternions::exp), Doc("exp"));
  m.def("exp", py::overload_cast<const QVec&>(&Quaternions::exp));
  m.def("sqrt", py::overload_cast<const Quaternion&>(&Quaternions::sqrt), Doc("sqrt"));
  m.def("sqrt", py::overload_cast<const QVec&>(&Quaternions::sqrt));
  m.def("angle", py::overload_cast<const Quaternion&>(&Quaternions::angle), Doc("angle"));
  m.def("angle", py::overload_cast<const QVec&>(&Quaternions::angle));
  m.def("inverse", py::overload_cast<const Quaternion&>(&Quaternions::inverse), Doc("inverse"));
  m.def("inverse", py::overload_cast<const QVec&>(&Quaternions::inverse));
  m.def("conjugate", py::overload_cast<const Quaternion&>(&Quaternions::conjugate), Doc("conjugate"));
  m.def("conjugate", py::overload_cast<const QVec&>(&Quaternions::conjugate));
  m.def("normalized", py::overload_cast<const Quaternion&>(&Quaternions::normalized), Doc("normalized"));
  m.def("normalized", py::overload_cast<const QVec&>(&Quaternions::normalized));
  m.def("normsquared", py::overload_cast<const Quaternion&>(&Quaternions::normsquared), Doc("normsquared"));
  m.def("normsquared", py::overload_cast<const QVec&>(&Quaternions::normsquared));
  m.def("sqrtOfRotor", py::overload_cast<const Quaternion&>(&Quaternions::sqrtOfRotor), Doc("sqrtOfRotor"));
  m.def("sqrtOfRotor", py::overload_cast<const QVec&>(&Quaternions::sqrtOfRotor));
  m.def("vec", py::overload_cast<const Quaternion&>(&Quaternions::vec), Doc("vec"));
  m.def("vec", py::overload_cast<const QVec&>(&Quaternions::vec));
  m.def("XParityConjugateSpinor", py::overload_cast<const Quaternion&>(&Quaternions::XParityConjugateSpinor));
  m.def("YParityConjugateSpinor", py::overload_cast<const Quaternion&>(&Quaternions::YParityConjugateSpinor));
  m.def("ZParityConjugateSpinor", py::overload_cast<const Quaternion&>(&Quaternions::ZParityConjugateSpinor));
  m.def("ParityConjugateSpinor", py::overload_cast<const Quaternion&>(&Quaternions::ParityConjugateSpinor));
  m.def("IntrinsicDistance", py::overload_cast<const Quaternion&, const Quaternion&>(&Quaternions::IntrinsicDistance),
        Doc("IntrinsicDistance"));
  m.def("ChordalDistance", py::overload_cast<const Quaternion&, const Quaternion&>(&Quaternions::ChordalDistance),
        Doc("ChordalDistance"));
  m.def("Slerp", &Quaternions::Slerp, py::arg("tau"), py::arg("Qa"), py::arg("Qb"), Doc("Slerp"));

  // Construction of arrays
  m.def("QuaternionArray", py::overload_cast<const DVec&, const DVec&>(&Quaternions::QuaternionArray), py::arg("vartheta"),
        py::arg("varphi"), Doc("QuaternionArray"));
  m.def("QuaternionArray", py::overload_cast<const DVec&, const DVec&, const DVec&>(&Quaternions::QuaternionArray),
        py::arg("alpha"), py::arg("beta"), py::arg("gamma"));
  m.def("QuaternionArray", py::overload_cast<const DVec&, const DVec&, const DVec&, const DVec&>(&Quaternions::QuaternionArray),
        py::arg("w"), py::arg("x"), py::arg("y"), py::arg("z"));
  m.def("QuaternionArray", py::overload_cast<const std::vector<DVec>&>(&Quaternions::QuaternionArray), py::arg("q"));
  m.def("QuaternionArray", py::overload_cast<const DVec&, const std::vector<DVec>&>(&Quaternions::QuaternionArray),
        py::arg("angle"), py::arg("axis"));
  m.def("Component", &Quaternions::Component, py::arg("Q"), py::arg("i"), Doc("Component"));
  m.def("Component0", &Quaternions::Component0, Doc("Component0"));
  m.def("Component1", &Quaternions::Component1, Doc("Component1"));
  m.def("Component2", &Quaternions::Component2, Doc("Component2"));
  m.def("Component3", &Quaternions::Component3, Doc("Component3"));

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
      "Return the tuple (rdot_x, rdot_y) for the two-dimensional integrand.");
  m.def("UnflipRotors", &Quaternions::UnflipRotors, py::arg("R"), py::arg("discont") = 1.4142135623730951, Doc("UnflipRotors"));
  m.def("RDelta", &Quaternions::RDelta, py::arg("R1"), py::arg("R2"), py::arg("IndexOfFiducialTime") = -1, Doc("RDelta"));
  m.def("Squad", &Quaternions::Squad, py::arg("RIn"), py::arg("tIn"), py::arg("tOut"), Doc("Squad"));
  m.def("FrameAngularVelocity", &Quaternions::FrameAngularVelocity, py::arg("f"), py::arg("t"), Doc("FrameAngularVelocity"));
  m.def("BoostRotor", &Quaternions::BoostRotor, py::arg("ThreeVelocity"), py::arg("ThreeVectorToBeBoosted"), Doc("BoostRotor"));

  // Mean rotors and alignment
  m.def("ApproximateMeanRotor", py::overload_cast<const QVec&>(&Quaternions::ApproximateMeanRotor), py::arg("R"),
        Doc("ApproximateMeanRotor"));
  m.def("ApproximateMeanRotor",
        py::overload_cast<const QVec&, const DVec&, const double, const double>(&Quaternions::ApproximateMeanRotor),
        py::arg("R"), py::arg("t"), py::arg("t1") = -1e300, py::arg("t2") = 1e300);
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
  m.def("MeanRotor", py::overload_cast<const QVec&>(&Quaternions::MeanRotor), py::arg("R"), Doc("MeanRotor"));
  m.def("MeanRotor", py::overload_cast<const QVec&, const DVec&, double, double>(&Quaternions::MeanRotor), py::arg("R"),
        py::arg("t"), py::arg("t1") = -1e300, py::arg("t2") = 1e300);

  // Searching sorted arrays
  m.def("hunt", &Quaternions::hunt, py::arg("x"), py::arg("x_i"), py::arg("i") = 0, Doc("hunt"));
  m.def("huntRight", &Quaternions::huntRight, py::arg("x"), py::arg("x_i"), py::arg("i") = 0, Doc("huntRight"));
}

}  // namespace GWFramesBindings
