// Copyright (c) 2014, Michael Boyle
// See LICENSE file for details

#include "bindings.hpp"

#include <sstream>

#include "Utilities.hpp"

namespace GWFramesBindings {

using GWFrames::Matrix;
using DVec = std::vector<double>;

void BindUtilities(py::module_& m) {
  py::class_<Matrix>(m, "Matrix", Doc("Matrix"))
      .def(py::init<>())
      .def(py::init<unsigned int, unsigned int>(), py::arg("rows"), py::arg("cols"))
      .def(py::init<unsigned int, unsigned int, const double>(), py::arg("rows"), py::arg("cols"), py::arg("a"))
      .def(py::init<const Matrix&>())
      .def(py::init<const std::vector<DVec>&>(), py::arg("DataIn"))
      .def("__call__", [](const Matrix& M, unsigned int r, unsigned int c) {
        if (r >= M.nrows() || c >= M.ncols()) {
          throw py::index_error("Matrix index out of range");
        }
        return M(r, c);
      })
      .def("set", [](Matrix& M, unsigned int r, unsigned int c, double v) -> Matrix& {
        if (r >= M.nrows() || c >= M.ncols()) {
          throw py::index_error("Matrix index out of range");
        }
        return M.set(r, c, v);
      }, py::return_value_policy::reference_internal)
      .def("nrows", &Matrix::nrows)
      .def("ncols", &Matrix::ncols)
      .def("resize", &Matrix::resize, py::arg("newNRows"), py::arg("newNCols"), py::arg("value") = 0.0)
      .def("clear", &Matrix::clear)
      .def("__sub__", &Matrix::operator-)
      .def("__mul__", &Matrix::operator*)
      .def("numpy", [](const Matrix& M) {
        py::array_t<double> a({static_cast<py::ssize_t>(M.nrows()), static_cast<py::ssize_t>(M.ncols())});
        for (unsigned int r = 0; r < M.nrows(); ++r) {
          for (unsigned int c = 0; c < M.ncols(); ++c) {
            a.mutable_at(r, c) = M(r, c);
          }
        }
        return a;
      }, "Return a copy of this Matrix as a two-dimensional numpy array.")
      .def("__repr__", [](const Matrix& M) {
        std::ostringstream s;
        s << "Matrix([";
        for (unsigned int r = 0; r < M.nrows(); ++r) {
          s << (r ? ",\n        [" : "[");
          for (unsigned int c = 0; c < M.ncols(); ++c) {
            s << (c ? ", " : "") << M(r, c);
          }
          s << "]";
        }
        s << "])";
        return s.str();
      });

  m.def("Unwrap", &GWFrames::Unwrap, py::arg("In"), Doc("Unwrap"));
  m.def("ScalarIntegral", &GWFrames::ScalarIntegral, py::arg("fdot"), py::arg("t"), Doc("ScalarIntegral"));
  m.def("CumulativeScalarIntegral", &GWFrames::CumulativeScalarIntegral, py::arg("fdot"), py::arg("t"),
        Doc("CumulativeScalarIntegral"));
  m.def("ScalarDerivative", &GWFrames::ScalarDerivative, py::arg("f"), py::arg("t"), Doc("ScalarDerivative"));
  m.def("ComplexDerivative", &GWFrames::ComplexDerivative, py::arg("f"), py::arg("t"), Doc("ComplexDerivative"));
  m.def("VectorIntegral", &GWFrames::VectorIntegral, py::arg("fdot"), py::arg("t"), Doc("VectorIntegral"));
  m.def("CumulativeVectorIntegral", &GWFrames::CumulativeVectorIntegral, py::arg("fdot"), py::arg("t"),
        Doc("CumulativeVectorIntegral"));
  m.def("Intersection", &GWFrames::Intersection, py::arg("t1"), py::arg("t2"), py::arg("MinStep") = 0.005,
        py::arg("MinTime") = -1e300, py::arg("MaxTime") = 1e300, Doc("Intersection"));
  m.def("Union", &GWFrames::Union, py::arg("t1"), py::arg("t2"), py::arg("MinStep") = 0.005, Doc("Union"));
  m.def("DominantPrincipalValue", &GWFrames::DominantPrincipalValue, py::arg("M"), Doc("DominantPrincipalValue"));
  m.def("SubordinatePrincipalValue", &GWFrames::SubordinatePrincipalValue, py::arg("M"),
        Doc("SubordinatePrincipalValue"));
  m.def("DominantPrincipalAxis", &GWFrames::DominantPrincipalAxis, py::arg("M"), Doc("DominantPrincipalAxis"));
  m.def("SubordinatePrincipalAxis", &GWFrames::SubordinatePrincipalAxis, py::arg("M"),
        Doc("SubordinatePrincipalAxis"));
  m.def("Eigenvalues", &GWFrames::Eigenvalues, py::arg("M"), Doc("Eigenvalues"));
  m.def("Eigenvectors", &GWFrames::Eigenvectors, py::arg("M"), Doc("Eigenvectors"));
  m.def("Eigensystem", &GWFrames::Eigensystem, py::arg("M"), Doc("Eigensystem"));
  m.def("Determinant", &GWFrames::Determinant, py::arg("M"), Doc("Determinant"));
}

}  // namespace GWFramesBindings
