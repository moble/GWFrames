// Copyright (c) 2014, Michael Boyle
// See LICENSE file for details

#ifndef GWFRAMES_BINDINGS_HPP
#define GWFRAMES_BINDINGS_HPP

// The casters must come before anything else that uses pybind11.
#include "casters.hpp"

#include <pybind11/operators.h>

#include <string>
#include <unordered_set>

#include "docstrings.hpp"

namespace GWFramesBindings {

namespace py = pybind11;

// Return the docstring for `name` (without the leading `GWFrames::` or
// `Quaternions::`), or an empty string if there is none.  pybind11 joins the
// docstrings of all overloads of a function, so each docstring is returned
// only the first time it is requested.
inline const char* Doc(const std::string& name) {
  static std::unordered_set<std::string> used;
  const auto& docstrings = Docstrings();
  const auto it = docstrings.find(name);
  if (it == docstrings.end() || !used.insert(name).second) {
    return "";
  }
  return it->second;
}

// Wrap a method returning `*this` so that the Python method returns the very
// object it was called on.  (pybind11 would otherwise return a new Python
// object when called on an instance of a derived class, like `PNWaveform`.)
template <typename C, typename R, typename... Args>
auto ReturnSelf(R& (C::*f)(Args...)) {
  return [f](py::object self, Args... args) {
    (self.cast<C&>().*f)(std::forward<Args>(args)...);
    return self;
  };
}
template <typename C, typename R, typename... Args>
auto ReturnSelf(const R& (C::*f)(Args...) const) {
  return [f](py::object self, Args... args) {
    (self.cast<const C&>().*f)(std::forward<Args>(args)...);
    return self;
  };
}

void BindQuaternions(py::module_& m);
void BindUtilities(py::module_& m);
void BindWaveforms(py::module_& m);
void BindScri(py::module_& m);

}  // namespace GWFramesBindings

#endif  // GWFRAMES_BINDINGS_HPP
