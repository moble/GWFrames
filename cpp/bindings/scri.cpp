// Copyright (c) 2026, Michael Boyle
// See LICENSE file for details

#include "bindings.hpp"

#include "Scri.hpp"

namespace GWFramesBindings {

using namespace GWFrames;
using Quaternions::Quaternion;
using CVec = std::vector<std::complex<double>>;

namespace {

// Allows Python subclasses of ScriFunctor to define `__call__`.
class PyScriFunctor : public ScriFunctor {
 public:
  double operator()(const Quaternion& R) const override {
    PYBIND11_OVERRIDE_NAME(double, ScriFunctor, "__call__", operator(), R);
  }
};

void CheckSliceIndex(const unsigned int i) {
  if (i > 6) {
    throw py::index_error("SliceOfScri index out of range");
  }
}

// The Python class name doubles as the docstring key, so that the two classes
// bound from this template get different docstrings.
template <class D>
void BindSliceOfScri(py::module_& m, const char* name) {
  using S = SliceOfScri<D>;
  py::class_<S>(m, name, Doc(name))
      .def(py::init<const int>(), py::arg("size") = 0)
      .def(py::init<const S&>(), py::arg("S"))
      .def_readwrite("psi0", &S::psi0)
      .def_readwrite("psi1", &S::psi1)
      .def_readwrite("psi2", &S::psi2)
      .def_readwrite("psi3", &S::psi3)
      .def_readwrite("psi4", &S::psi4)
      .def_readwrite("sigma", &S::sigma)
      .def_readwrite("sigmadot", &S::sigmadot)
      .def("__len__", [](const S&) { return 7; })
      .def("__getitem__", [](const S& s, const unsigned int i) -> const D& {
        CheckSliceIndex(i);
        return s[i];
      }, py::arg("i"), py::return_value_policy::reference_internal)
      .def("__setitem__", [](S& s, const unsigned int i, const D& a) {
        CheckSliceIndex(i);
        s[i] = a;
      }, py::arg("i"), py::arg("value"));
}

}  // namespace

void BindScri(py::module_& m) {
  py::class_<ScriFunctor, PyScriFunctor>(m, "ScriFunctor", Doc("ScriFunctor"))
      .def(py::init<>())
      .def("__call__", &ScriFunctor::operator(), py::arg("R"));

  m.def("Boost", &Boost, py::arg("v"), py::arg("n"), Doc("Boost"));

  py::class_<DataGrid>(m, "DataGrid", Doc("DataGrid"))
      .def(py::init<const int>(), py::arg("size") = 0)
      .def(py::init<const DataGrid&>(), py::arg("A"))
      .def(py::init<const int, const int, const int, const CVec&>(), py::arg("Spin"), py::arg("N_theta"),
           py::arg("N_phi"), py::arg("D"))
      .def(py::init<Modes, const int, const int>(), py::arg("M"), py::arg("N_theta") = 0, py::arg("N_phi") = 0)
      .def(py::init<const Modes&, const ThreeVector&, const int, const int>(), py::arg("M"), py::arg("v"),
           py::arg("N_theta") = 0, py::arg("N_phi") = 0)
      .def(py::init<const int, const int, const int, const ThreeVector&, const ScriFunctor&>(), py::arg("Spin"),
           py::arg("N_theta"), py::arg("N_phi"), py::arg("v"), py::arg("f"))
      .def("SetSpin", &DataGrid::SetSpin, py::arg("ess"), py::return_value_policy::reference_internal,
           Doc("DataGrid::SetSpin"))
      .def("SetNTheta", &DataGrid::SetNTheta, py::arg("N_theta"), py::return_value_policy::reference_internal,
           Doc("DataGrid::SetNTheta"))
      .def("SetNPhi", &DataGrid::SetNPhi, py::arg("N_phi"), py::return_value_policy::reference_internal,
           Doc("DataGrid::SetNPhi"))
      .def("size", &DataGrid::size, Doc("DataGrid::size"))
      .def("__len__", &DataGrid::size)
      .def("Spin", &DataGrid::Spin, Doc("DataGrid::Spin"))
      .def("N_theta", &DataGrid::N_theta, Doc("DataGrid::N_theta"))
      .def("N_phi", &DataGrid::N_phi, Doc("DataGrid::N_phi"))
      .def("Data", &DataGrid::Data, Doc("DataGrid::Data"))
      .def("__getitem__", [](const DataGrid& d, const unsigned int i) {
        if (i >= d.size()) {
          throw py::index_error("DataGrid index out of range");
        }
        return d[i];
      }, py::arg("i"))
      .def("__setitem__", [](DataGrid& d, const unsigned int i, const std::complex<double>& a) {
        if (i >= d.size()) {
          throw py::index_error("DataGrid index out of range");
        }
        d[i] = a;
      }, py::arg("i"), py::arg("value"))
      .def(py::self * py::self, py::arg("other"))
      .def(py::self / py::self, py::arg("other"))
      .def(py::self + py::self, py::arg("other"))
      .def(py::self - py::self, py::arg("other"))
      .def(double() * py::self, py::arg("other"))
      .def(double() / py::self, py::arg("other"))
      .def(double() + py::self, py::arg("other"))
      .def(double() - py::self, py::arg("other"))
      .def(py::self * double(), py::arg("other"))
      .def(py::self / double(), py::arg("other"))
      .def(py::self + double(), py::arg("other"))
      .def(py::self - double(), py::arg("other"))
      .def("pow", &DataGrid::pow, py::arg("p"), Doc("DataGrid::pow"));
  m.def("ConformalFactorGrid", &ConformalFactorGrid, py::arg("v"), py::arg("n_theta"), py::arg("n_phi"),
        Doc("ConformalFactorGrid"));
  m.def("InverseConformalFactorGrid", &InverseConformalFactorGrid, py::arg("v"), py::arg("n_theta"),
        py::arg("n_phi"), Doc("InverseConformalFactorGrid"));
  m.def("InverseConformalFactorBoostedGrid", &InverseConformalFactorBoostedGrid, py::arg("v"), py::arg("n_theta"),
        py::arg("n_phi"), Doc("InverseConformalFactorBoostedGrid"));

  py::class_<Modes>(m, "Modes", Doc("Modes"))
      .def(py::init<const int>(), py::arg("size") = 0)
      .def(py::init<const Modes&>(), py::arg("A"))
      .def(py::init<const int, const CVec&>(), py::arg("spin"), py::arg("Data"))
      .def(py::init<DataGrid, const int>(), py::arg("D"), py::arg("L") = -1)
      .def("SetSpin", &Modes::SetSpin, py::arg("ess"), py::return_value_policy::reference_internal,
           Doc("Modes::SetSpin"))
      .def("SetEllMax", &Modes::SetEllMax, py::arg("ell"), py::return_value_policy::reference_internal,
           Doc("Modes::SetEllMax"))
      .def("size", &Modes::size, Doc("Modes::size"))
      .def("__len__", &Modes::size)
      .def("Spin", &Modes::Spin, Doc("Modes::Spin"))
      .def("EllMax", &Modes::EllMax, Doc("Modes::EllMax"))
      .def("Data", &Modes::Data, Doc("Modes::Data"))
      .def("__getitem__", [](const Modes& M, const unsigned int i) {
        if (i >= M.size()) {
          throw py::index_error("Modes index out of range");
        }
        return M[i];
      }, py::arg("i"))
      .def("__setitem__", [](Modes& M, const unsigned int i, const std::complex<double>& a) {
        if (i >= M.size()) {
          throw py::index_error("Modes index out of range");
        }
        M[i] = a;
      }, py::arg("i"), py::arg("value"))
      .def("pow", &Modes::pow, py::arg("p"), Doc("Modes::pow"))
      .def("bar", &Modes::bar, Doc("Modes::bar"))
      .def(py::self * py::self, py::arg("other"))
      .def(py::self / py::self, py::arg("other"))
      .def(py::self + py::self, py::arg("other"))
      .def(py::self - py::self, py::arg("other"))
      .def("edth", &Modes::edth, Doc("Modes::edth"))
      .def("edthbar", &Modes::edthbar, Doc("Modes::edthbar"))
      .def("edth2edthbar2", &Modes::edth2edthbar2, Doc("Modes::edth2edthbar2"))
      .def("EvaluateAtPoint", py::overload_cast<const double, const double>(&Modes::EvaluateAtPoint, py::const_),
           py::arg("vartheta"), py::arg("varphi"))
      .def("EvaluateAtPoint", py::overload_cast<const Quaternion&>(&Modes::EvaluateAtPoint, py::const_),
           py::arg("R"), Doc("Modes::EvaluateAtPoint"));
  m.def("vFromOneOverK", &vFromOneOverK, py::arg("OneOverK"), Doc("vFromOneOverK"));

  BindSliceOfScri<DataGrid>(m, "SliceOfScriGrid");
  BindSliceOfScri<Modes>(m, "SliceOfScriModes");

  py::class_<SliceModes, SliceOfScri<Modes>>(m, "SliceModes", Doc("SliceModes"))
      .def(py::init<const int>(), py::arg("ellMax") = 0)
      .def(py::init<const SliceModes&>(), py::arg("S"))
      .def("EllMax", &SliceModes::EllMax, Doc("SliceModes::EllMax"))
      .def("Mass", &SliceModes::Mass, Doc("SliceModes::Mass"))
      .def("FourMomentum", &SliceModes::FourMomentum, Doc("SliceModes::FourMomentum"))
      .def("SuperMomentum", &SliceModes::SuperMomentum, Doc("SliceModes::SuperMomentum"))
      .def("BMSTransformationOnSlice", &SliceModes::BMSTransformationOnSlice, py::arg("u"), py::arg("v"),
           py::arg("delta"), Doc("SliceModes::BMSTransformationOnSlice"));

  py::class_<Scri>(m, "Scri", Doc("Scri"))
      .def(py::init<const Waveform&, const Waveform&, const Waveform&, const Waveform&, const Waveform&,
                    const Waveform&>(),
           py::arg("psi0"), py::arg("psi1"), py::arg("psi2"), py::arg("psi3"), py::arg("psi4"), py::arg("sigma"))
      .def(py::init<const Scri&>(), py::arg("S"))
      .def("BMSTransformation", &Scri::BMSTransformation, py::arg("u0"), py::arg("v"), py::arg("delta"),
           Doc("Scri::BMSTransformation"))
      .def("NTimes", &Scri::NTimes, Doc("Scri::NTimes"))
      .def("T", &Scri::T, Doc("Scri::T"))
      .def("__len__", &Scri::NTimes)
      .def("__getitem__", [](const Scri& s, const int i) {
        if (i < 0 || i >= s.NTimes()) {
          throw py::index_error("Scri index out of range");
        }
        return s[i];
      }, py::arg("i"))
      .def("__setitem__", [](Scri& s, const int i, const SliceModes& a) {
        if (i < 0 || i >= s.NTimes()) {
          throw py::index_error("Scri index out of range");
        }
        s[i] = a;
      }, py::arg("i"), py::arg("value"));

  py::class_<SuperMomenta>(m, "SuperMomenta", Doc("SuperMomenta"))
      .def(py::init<const unsigned int>(), py::arg("size"))
      .def(py::init<const SuperMomenta&>(), py::arg("S"))
      .def(py::init<const std::vector<double>&, const std::vector<Modes>&>(), py::arg("T"), py::arg("psi"))
      .def(py::init<const Scri&>(), py::arg("scri"))
      .def("NTimes", &SuperMomenta::NTimes, Doc("SuperMomenta::NTimes"))
      .def("T", &SuperMomenta::T, Doc("SuperMomenta::T"))
      .def("__len__", &SuperMomenta::NTimes)
      .def("__getitem__", [](const SuperMomenta& s, const int i) {
        if (i < 0 || i >= s.NTimes()) {
          throw py::index_error("SuperMomenta index out of range");
        }
        return s[i];
      }, py::arg("i"))
      .def("__setitem__", [](SuperMomenta& s, const int i, const Modes& a) {
        if (i < 0 || i >= s.NTimes()) {
          throw py::index_error("SuperMomenta index out of range");
        }
        s[i] = a;
      }, py::arg("i"), py::arg("value"))
      .def("BMSTransform", &SuperMomenta::BMSTransform, py::arg("OneOverK"), py::arg("delta"),
           Doc("SuperMomenta::BMSTransform"))
      .def("MoreschiIteration", &SuperMomenta::MoreschiIteration, py::arg("OneOverK"), py::arg("delta"),
           Doc("SuperMomenta::MoreschiIteration"));
}

}  // namespace GWFramesBindings
