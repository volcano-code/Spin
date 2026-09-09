#pragma once
#include <Eigen/Core>
#include <pybind11/pybind11.h>

// Reuse the established lossless ingestion boundary without NumPy inference.
void registerSourceSpinMaterialBindings(pybind11::module_&,
    Eigen::MatrixXd (*vertices)(const pybind11::handle&),
    Eigen::MatrixXi (*faces)(const pybind11::handle&));
