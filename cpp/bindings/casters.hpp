// Copyright (c) 2014, Michael Boyle
// See LICENSE file for details

// Type casters reproducing the conversions of the old SWIG typemaps.
//
// Vectors of numbers (and vectors of such vectors) are accepted from numpy
// arrays or from any Python sequence, and are returned as numpy arrays.
// Vectors of quaternions are accepted from sequences of `Quaternion` objects
// or of length-4 sequences, or from an (N, 4) array of floats, and are
// returned as one-dimensional numpy arrays of `Quaternion` objects.
//
// These specializations must be visible before any binding code uses the
// types, so every binding source includes this header first, through
// `bindings.hpp`.

#ifndef GWFRAMES_BINDINGS_CASTERS_HPP
#define GWFRAMES_BINDINGS_CASTERS_HPP

#include <complex>
#include <type_traits>
#include <vector>

#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Quaternions.hpp"

namespace GWFramesBindings {

namespace py = pybind11;

template <typename T>
struct is_complex : std::false_type {};
template <typename T>
struct is_complex<std::complex<T>> : std::true_type {};

inline bool IsStringLike(py::handle src) {
  return py::isinstance<py::str>(src) || py::isinstance<py::bytes>(src);
}

// Convert `src` to a contiguous numpy array of `T`, or return a null array.
// Without `convert`, only arrays that already have the right dtype, and
// lists or tuples, are considered, so that overload resolution prefers exact
// matches.  Complex data are never silently converted to real data.
template <typename T>
py::array_t<T, py::array::c_style | py::array::forcecast> AsArray(py::handle src, bool convert) {
  using Array = py::array_t<T, py::array::c_style | py::array::forcecast>;
  // A default-constructed array_t is an empty array, not a null handle.
  const auto null = [] { return py::reinterpret_steal<Array>(py::handle()); };
  if (!src || IsStringLike(src)) {
    return null();
  }
  const bool is_array = py::isinstance<py::array>(src);
  if (!convert && !py::isinstance<py::array_t<T>>(src) && !py::isinstance<py::list>(src)
      && !py::isinstance<py::tuple>(src)) {
    return null();
  }
  if (is_array) {
    const char kind = py::reinterpret_borrow<py::array>(src).dtype().kind();
    if (kind == 'O' || kind == 'U' || kind == 'S' || (kind == 'c' && !is_complex<T>::value)) {
      return null();
    }
  }
  return Array::ensure(src);
}

template <typename T>
struct NumpyVectorCaster {
  using Vector = std::vector<T>;
  PYBIND11_TYPE_CASTER(Vector, py::detail::const_name("numpy.ndarray"));

  bool load(py::handle src, bool convert) {
    auto array = AsArray<T>(src, convert);
    if (!array) {
      PyErr_Clear();
      return false;
    }
    if (array.ndim() != 1) {
      return false;
    }
    value.assign(array.data(), array.data() + array.size());
    return true;
  }

  static py::handle cast(const Vector& src, py::return_value_policy, py::handle) {
    py::array_t<T> array(static_cast<py::ssize_t>(src.size()));
    std::copy(src.begin(), src.end(), array.mutable_data());
    return array.release();
  }
};

template <typename T>
struct NumpyVectorVectorCaster {
  using Vector = std::vector<std::vector<T>>;
  PYBIND11_TYPE_CASTER(Vector, py::detail::const_name("numpy.ndarray"));

  bool load(py::handle src, bool convert) {
    if (!src || IsStringLike(src)) {
      return false;
    }
    auto array = AsArray<T>(src, convert);
    if (array) {
      if (array.ndim() == 2) {
        const auto rows = array.shape(0);
        const auto cols = array.shape(1);
        value.assign(rows, std::vector<T>(cols));
        for (py::ssize_t i = 0; i < rows; ++i) {
          std::copy(array.data(i, 0), array.data(i, 0) + cols, value[i].begin());
        }
        return true;
      }
      if (array.ndim() == 1 && array.size() == 0) {
        value.clear();
        return true;
      }
      return false;
    }
    PyErr_Clear();
    // Ragged input is still acceptable as a sequence of sequences.
    if (!convert || !py::isinstance<py::sequence>(src)) {
      return false;
    }
    auto sequence = py::reinterpret_borrow<py::sequence>(src);
    value.clear();
    value.reserve(sequence.size());
    for (auto item : sequence) {
      NumpyVectorCaster<T> row;
      if (!row.load(item, convert)) {
        return false;
      }
      value.push_back(std::move(*row));
    }
    return true;
  }

  static py::handle cast(const Vector& src, py::return_value_policy policy, py::handle parent) {
    const std::size_t cols = src.empty() ? 0 : src[0].size();
    bool rectangular = true;
    for (const auto& row : src) {
      rectangular = rectangular && row.size() == cols;
    }
    if (!rectangular) {
      py::list list;
      for (const auto& row : src) {
        list.append(py::reinterpret_steal<py::object>(NumpyVectorCaster<T>::cast(row, policy, parent)));
      }
      return list.release();
    }
    py::array_t<T> array({static_cast<py::ssize_t>(src.size()), static_cast<py::ssize_t>(cols)});
    for (std::size_t i = 0; i < src.size(); ++i) {
      std::copy(src[i].begin(), src[i].end(), array.mutable_data(i, 0));
    }
    return array.release();
  }
};

struct QuaternionVectorCaster {
  using Vector = std::vector<Quaternions::Quaternion>;
  PYBIND11_TYPE_CASTER(Vector, py::detail::const_name("numpy.ndarray[Quaternion]"));

  bool load(py::handle src, bool convert) {
    if (!src || IsStringLike(src)) {
      return false;
    }
    // An (N, 4) array of floats
    if (py::isinstance<py::array>(src)
        && py::reinterpret_borrow<py::array>(src).dtype().kind() != 'O') {
      auto array = AsArray<double>(src, convert);
      if (!array) {
        PyErr_Clear();
        return false;
      }
      if (array.ndim() == 1 && array.size() == 0) {
        value.clear();
        return true;
      }
      if (array.ndim() != 2 || array.shape(1) != 4) {
        return false;
      }
      value.resize(array.shape(0));
      for (py::ssize_t i = 0; i < array.shape(0); ++i) {
        value[i] = Quaternions::Quaternion(array.at(i, 0), array.at(i, 1), array.at(i, 2), array.at(i, 3));
      }
      return true;
    }
    // Any other sequence of things that convert to Quaternion
    if (!py::isinstance<py::sequence>(src)) {
      return false;
    }
    auto sequence = py::reinterpret_borrow<py::sequence>(src);
    value.clear();
    value.reserve(sequence.size());
    for (auto item : sequence) {
      py::detail::make_caster<Quaternions::Quaternion> element;
      if (!element.load(item, convert)) {
        return false;
      }
      value.push_back(py::detail::cast_op<const Quaternions::Quaternion&>(element));
    }
    return true;
  }

  static py::handle cast(const Vector& src, py::return_value_policy, py::handle) {
    py::array array(py::dtype("O"), py::array::ShapeContainer{static_cast<py::ssize_t>(src.size())},
                    py::array::StridesContainer{});
    auto data = static_cast<PyObject**>(array.mutable_data());
    for (std::size_t i = 0; i < src.size(); ++i) {
      PyObject* q = py::cast(src[i]).release().ptr();
      Py_XDECREF(data[i]);
      data[i] = q;
    }
    return array.release();
  }
};

}  // namespace GWFramesBindings

namespace pybind11 {
namespace detail {

template <>
struct type_caster<std::vector<double>> : GWFramesBindings::NumpyVectorCaster<double> {};
template <>
struct type_caster<std::vector<int>> : GWFramesBindings::NumpyVectorCaster<int> {};
template <>
struct type_caster<std::vector<unsigned int>> : GWFramesBindings::NumpyVectorCaster<unsigned int> {};
template <>
struct type_caster<std::vector<std::complex<double>>>
    : GWFramesBindings::NumpyVectorCaster<std::complex<double>> {};
template <>
struct type_caster<std::vector<std::vector<double>>> : GWFramesBindings::NumpyVectorVectorCaster<double> {};
template <>
struct type_caster<std::vector<std::vector<int>>> : GWFramesBindings::NumpyVectorVectorCaster<int> {};
template <>
struct type_caster<std::vector<std::vector<std::complex<double>>>>
    : GWFramesBindings::NumpyVectorVectorCaster<std::complex<double>> {};
template <>
struct type_caster<std::vector<Quaternions::Quaternion>> : GWFramesBindings::QuaternionVectorCaster {};

}  // namespace detail
}  // namespace pybind11

#endif  // GWFRAMES_BINDINGS_CASTERS_HPP
