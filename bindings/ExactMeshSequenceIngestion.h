#pragma once

#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS

#include <Eigen/Core>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <cmath>
#include <limits>
#include <string>


namespace exact_mesh_ingestion
{

namespace py = pybind11;


inline py::object numpyScalarType()
{
    return py::module_::import("numpy").attr("generic");
}


inline bool isNumpyScalar(
    const py::handle& value,
    const py::handle& scalarType
)
{
    return !py::isinstance<py::array>(value) &&
        py::isinstance(value, scalarType);
}


inline py::object integerObject(const py::handle& value, const char* message)
{
    PyObject* indexed = PyNumber_Index(value.ptr());
    if (indexed == nullptr)
    {
        PyErr_Clear();
        throw py::type_error(message);
    }
    return py::reinterpret_steal<py::object>(indexed);
}


inline double exactIntegerAsDouble(const py::handle& value)
{
    const py::object integer = integerObject(
        value,
        "vertex coordinates must be real Python or NumPy numeric scalars"
    );
    const double converted = PyLong_AsDouble(integer.ptr());
    if (PyErr_Occurred() != nullptr)
    {
        PyErr_Clear();
        throw py::value_error(
            "integer coordinate is not exactly representable as binary64"
        );
    }
    py::object roundTrip = py::reinterpret_steal<py::object>(
        PyLong_FromDouble(converted)
    );
    if (!roundTrip || !integer.equal(roundTrip))
    {
        throw py::value_error(
            "integer coordinate is not exactly representable as binary64"
        );
    }
    return converted;
}


inline double exactVertexScalar(
    const py::handle& value,
    const py::handle& numpyScalar
)
{
    if (PyBool_Check(value.ptr()))
    {
        throw py::type_error("boolean vertex coordinates are invalid");
    }
    if (PyLong_Check(value.ptr()))
    {
        return exactIntegerAsDouble(value);
    }
    if (PyFloat_Check(value.ptr()))
    {
        const double converted = PyFloat_AsDouble(value.ptr());
        if (!std::isfinite(converted))
        {
            throw py::value_error("vertex coordinates must be finite");
        }
        return converted;
    }
    if (py::isinstance<py::array>(value))
    {
        throw py::type_error("array-valued vertex coordinates are invalid");
    }
    if (!isNumpyScalar(value, numpyScalar))
    {
        throw py::type_error(
            "vertex coordinates must be real Python or NumPy numeric scalars"
        );
    }

    const py::array scalar = py::array::ensure(value);
    if (!scalar || scalar.ndim() != 0)
    {
        throw py::type_error(
            "vertex coordinates must be real Python or NumPy numeric scalars"
        );
    }
    const char kind = scalar.dtype().kind();
    if (kind == 'b')
    {
        throw py::type_error("boolean vertex coordinates are invalid");
    }
    if (kind == 'i' || kind == 'u')
    {
        return exactIntegerAsDouble(value);
    }
    if (kind != 'f' || scalar.itemsize() > 8)
    {
        throw py::type_error(
            "exact boundary accepts at most binary64 real coordinates"
        );
    }
    const double converted = PyFloat_AsDouble(value.ptr());
    if (PyErr_Occurred() != nullptr)
    {
        PyErr_Clear();
        throw py::type_error("vertex coordinate cannot be represented as binary64");
    }
    if (!std::isfinite(converted))
    {
        throw py::value_error("vertex coordinates must be finite");
    }
    return converted;
}


inline int exactFaceScalar(
    const py::handle& value,
    const py::handle& numpyScalar
)
{
    if (PyBool_Check(value.ptr()))
    {
        throw py::type_error("boolean face indices are invalid");
    }
    if (py::isinstance<py::array>(value))
    {
        throw py::type_error("array-valued face indices are invalid");
    }
    bool integer = PyLong_Check(value.ptr());
    if (!integer && isNumpyScalar(value, numpyScalar))
    {
        const py::array scalar = py::array::ensure(value);
        integer = scalar && scalar.ndim() == 0 &&
            (scalar.dtype().kind() == 'i' || scalar.dtype().kind() == 'u');
    }
    if (!integer)
    {
        throw py::type_error(
            "face indices must be nonboolean Python or NumPy integers"
        );
    }

    const py::object indexed = integerObject(
        value,
        "face indices must be nonboolean Python or NumPy integers"
    );
    int overflow = 0;
    const long long converted = PyLong_AsLongLongAndOverflow(
        indexed.ptr(), &overflow
    );
    if (PyErr_Occurred() != nullptr)
    {
        PyErr_Clear();
        overflow = 1;
    }
    if (
        overflow != 0 ||
        converted < static_cast<long long>(std::numeric_limits<int>::min()) ||
        converted > static_cast<long long>(std::numeric_limits<int>::max())
    )
    {
        throw py::value_error(
            "each face index must be representable as a C++ int"
        );
    }
    return static_cast<int>(converted);
}


inline py::sequence matrixRows(const py::handle& value, const char* name)
{
    if (!PyList_Check(value.ptr()) && !PyTuple_Check(value.ptr()))
    {
        throw py::type_error(
            std::string(name) + " must be a NumPy array, list, or tuple"
        );
    }
    py::sequence rows = py::reinterpret_borrow<py::sequence>(value);
    if (py::len(rows) == 0)
    {
        throw py::value_error(std::string(name) + " must have shape (N, 3)");
    }
    if (
        py::len(rows) >
        static_cast<py::ssize_t>(std::numeric_limits<Eigen::Index>::max())
    )
    {
        throw py::value_error(std::string(name) + " has too many rows");
    }
    return rows;
}


inline py::sequence matrixRow(
    const py::handle& row,
    const char* name,
    const bool vertices
)
{
    if (py::isinstance<py::array>(row))
    {
        const py::array array = py::reinterpret_borrow<py::array>(row);
        if (array.ndim() != 1 || array.shape(0) != 3)
        {
            throw py::value_error(std::string(name) + " must have shape (N, 3)");
        }
        const char kind = array.dtype().kind();
        if (vertices)
        {
            if ((kind != 'f' && kind != 'i' && kind != 'u') || array.itemsize() > 8)
            {
                throw py::type_error(
                    "vertex NumPy rows must use at most binary64 real numeric dtypes"
                );
            }
        }
        else if ((kind != 'i' && kind != 'u') || array.itemsize() > 8)
        {
            throw py::type_error("face NumPy rows must use integer dtypes");
        }
    }
    else if (!PyList_Check(row.ptr()) && !PyTuple_Check(row.ptr()))
    {
        throw py::type_error(
            std::string(name) + " rows must be lists, tuples, or one-dimensional NumPy arrays"
        );
    }

    py::sequence sequence = py::reinterpret_borrow<py::sequence>(row);
    if (py::len(sequence) != 3)
    {
        throw py::value_error(std::string(name) + " must have shape (N, 3)");
    }
    return sequence;
}


inline Eigen::MatrixXd verticesFromSequence(const py::handle& value)
{
    const py::sequence rows = matrixRows(value, "vertices");
    const py::object numpyScalar = numpyScalarType();
    Eigen::MatrixXd result(static_cast<Eigen::Index>(py::len(rows)), 3);
    for (py::ssize_t row = 0; row < py::len(rows); ++row)
    {
        const py::object originalRow = rows[row];
        const py::sequence columns = matrixRow(originalRow, "vertices", true);
        for (py::ssize_t column = 0; column < 3; ++column)
        {
            const py::object scalar = columns[column];
            result(
                static_cast<Eigen::Index>(row),
                static_cast<Eigen::Index>(column)
            ) = exactVertexScalar(scalar, numpyScalar);
        }
    }
    return result;
}


inline Eigen::MatrixXi facesFromSequence(const py::handle& value)
{
    const py::sequence rows = matrixRows(value, "faces");
    const py::object numpyScalar = numpyScalarType();
    Eigen::MatrixXi result(static_cast<Eigen::Index>(py::len(rows)), 3);
    for (py::ssize_t row = 0; row < py::len(rows); ++row)
    {
        const py::object originalRow = rows[row];
        const py::sequence columns = matrixRow(originalRow, "faces", false);
        for (py::ssize_t column = 0; column < 3; ++column)
        {
            const py::object scalar = columns[column];
            result(
                static_cast<Eigen::Index>(row),
                static_cast<Eigen::Index>(column)
            ) = exactFaceScalar(scalar, numpyScalar);
        }
    }
    return result;
}

} // namespace exact_mesh_ingestion

#endif
