#pragma once
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Implemented in src/gfoil_ad_bindings.cpp (separate TU — includes real_type.hpp,
// must NOT include real_type.h to avoid norm2 redefinition).
py::dict run_AD_py(py::dict input, py::dict jacobian);
