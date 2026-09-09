// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/SparseCore>
#include <cstddef>
#include <cstdint>
#include <fmt/format.h>
#include <nanobind/eigen/sparse.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "gelex/bayes/draws.h"
#include "gelex/exception.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/csc_reader.h"
#include "gelex/io/csc_writer.h"
#include "gelex/io/dense_reader.h"
#include "gelex/io/dense_writer.h"

#include "gelex_py/register.h"

namespace nb = nanobind;

namespace gelex_py
{

namespace
{

// ---- format

auto register_format(nb::module_& m) -> void
{
    nb::enum_<gelex::BinaryType>(m, "BinaryType")
        .value("float64", gelex::BinaryType::float64)
        .value("float32", gelex::BinaryType::float32)
        .value("uint8", gelex::BinaryType::uint8);

    nb::class_<gelex::MatrixHeader>(m, "MatrixHeader")
        .def_prop_ro(
            "identifier",
            [](const gelex::MatrixHeader& header) { return header.identifier; })
        .def_prop_ro(
            "type",
            [](const gelex::MatrixHeader& header) { return header.type; })
        .def_prop_ro(
            "shape",
            [](const gelex::MatrixHeader& header)
            { return std::pair{header.shape[0], header.shape[1]}; })
        .def(
            "__repr__",
            [](const gelex::MatrixHeader& header)
            {
                return fmt::format(
                    "MatrixHeader(identifier='{}', shape=({}, {}))",
                    header.identifier,
                    header.shape[0],
                    header.shape[1]);
            });
}

// ---- writer

// Inputs are not converted: a dtype or layout mismatch is a TypeError rather
// than a silent copy. Column-major (rows, columns) is the on-disk layout, so
// a whole matrix is taken as an F-contiguous array.
template <typename T>
using ColumnArray
    = nb::ndarray<const T, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
template <typename T>
using MatrixArray
    = nb::ndarray<const T, nb::ndim<2>, nb::f_contig, nb::device::cpu>;

template <gelex::detail::SupportedDtype T>
auto register_dense_stream(nb::module_& m, const char* name) -> void
{
    using Stream = gelex::DenseStream<T>;
    nb::class_<Stream>(
        m,
        name,
        "Handle to one reserved matrix; append() adds one column (one draw) "
        "and write() stores the whole matrix at once.")
        .def(
            "append",
            [](Stream& self, const ColumnArray<T>& column)
            { self << std::span<const T>{column.data(), column.shape(0)}; },
            nb::arg("column").noconvert())
        .def(
            "write",
            [](Stream& self, const MatrixArray<T>& values)
            {
                self << std::span<const T>{
                    values.data(), values.shape(0) * values.shape(1)};
            },
            nb::arg("values").noconvert())
        .def_prop_ro("identifier", &Stream::identifier);
}

using AnyDenseStream = std::variant<
    gelex::DenseStream<double>,
    gelex::DenseStream<float>,
    gelex::DenseStream<std::uint8_t>>;

auto reserve(
    gelex::DenseWriter& writer,
    std::string_view identifier,
    gelex::BinaryType type,
    gelex::BinaryShape shape) -> AnyDenseStream
{
    switch (type)
    {
        case gelex::BinaryType::float64:
            return writer.reserve<double>(identifier, shape);
        case gelex::BinaryType::float32:
            return writer.reserve<float>(identifier, shape);
        case gelex::BinaryType::uint8:
            return writer.reserve<std::uint8_t>(identifier, shape);
    }
    throw gelex::GelexException("unsupported matrix dtype");
}

template <gelex::detail::SupportedDtype T>
auto register_csc_stream(nb::module_& m, const char* name) -> void
{
    using Stream = gelex::CscStream<T>;
    nb::class_<Stream>(
        m,
        name,
        "Handle to one reserved sparse matrix; append() adds one dense column "
        "(one draw) and stores its non-zero entries.")
        .def(
            "append",
            [](Stream& self, const ColumnArray<T>& column)
            { self << std::span<const T>{column.data(), column.shape(0)}; },
            nb::arg("column").noconvert());
}

using AnyCscStream = std::variant<
    gelex::CscStream<double>,
    gelex::CscStream<float>,
    gelex::CscStream<std::uint8_t>>;

auto reserve_csc(
    gelex::CscWriter& writer,
    std::string_view identifier,
    gelex::BinaryType type,
    gelex::BinaryShape shape) -> AnyCscStream
{
    switch (type)
    {
        case gelex::BinaryType::float64:
            return writer.reserve<double>(identifier, shape);
        case gelex::BinaryType::float32:
            return writer.reserve<float>(identifier, shape);
        case gelex::BinaryType::uint8:
            return writer.reserve<std::uint8_t>(identifier, shape);
    }
    throw gelex::GelexException("unsupported matrix dtype");
}

auto register_writer(nb::module_& m) -> void
{
    register_dense_stream<double>(m, "DenseStreamF64");
    register_dense_stream<float>(m, "DenseStreamF32");
    register_dense_stream<std::uint8_t>(m, "DenseStreamU8");

    nb::class_<gelex::DenseWriter>(
        m,
        "DenseWriter",
        "Writer for gelex dense containers. Reserve matrices with a dtype and "
        "(rows, columns) shape, fill them column by column, then close() (or "
        "leave the with-block) to publish the file; every matrix must be "
        "complete. An unclosed writer discards its output.")
        .def(nb::init<std::string_view>(), nb::arg("path"))
        .def(
            "reserve",
            &reserve,
            nb::arg("identifier"),
            nb::arg("type"),
            nb::arg("shape"),
            nb::keep_alive<0, 1>())
        .def("close", &gelex::DenseWriter::close)
        .def_prop_ro("is_open", &gelex::DenseWriter::is_open)
        .def(
            "__enter__",
            [](nb::handle_t<gelex::DenseWriter> self) { return self; })
        .def(
            "__exit__",
            [](gelex::DenseWriter& self, nb::args) { self.close(); },
            nb::arg("args"));

    register_csc_stream<double>(m, "CscStreamF64");
    register_csc_stream<float>(m, "CscStreamF32");
    register_csc_stream<std::uint8_t>(m, "CscStreamU8");

    nb::class_<gelex::CscWriter>(
        m,
        "CscWriter",
        "Writer for gelex CSC containers such as the MCMC .draws.csc output. "
        "Reserve matrices with a dtype and (rows, columns) shape, append "
        "dense columns whose zeros are dropped, then close() (or leave the "
        "with-block) to publish the file; every matrix must be complete. An "
        "unclosed writer discards its output.")
        .def(nb::init<std::string_view>(), nb::arg("path"))
        .def(
            "reserve",
            &reserve_csc,
            nb::arg("identifier"),
            nb::arg("type"),
            nb::arg("shape"),
            nb::keep_alive<0, 1>())
        .def("close", &gelex::CscWriter::close)
        .def_prop_ro("is_open", &gelex::CscWriter::is_open)
        .def(
            "__enter__",
            [](nb::handle_t<gelex::CscWriter> self) { return self; })
        .def(
            "__exit__",
            [](gelex::CscWriter& self, nb::args) { self.close(); },
            nb::arg("args"));
}

// ---- reader

using PayloadArray = nb::ndarray<nb::numpy, nb::ro>;

// The array aliases the reader's memory map: the reader object is registered
// as the array owner so Python keeps it alive for as long as any view exists.
template <gelex::detail::SupportedDtype T>
auto payload_view(
    const gelex::DenseReader& reader,
    std::string_view identifier,
    nb::handle owner) -> PayloadArray
{
    const auto map = reader.to_map<T>(identifier);
    const std::size_t shape[2]{
        static_cast<std::size_t>(map.rows()),
        static_cast<std::size_t>(map.cols())};
    const std::int64_t strides[2]{1, static_cast<std::int64_t>(map.rows())};
    return PayloadArray{map.data(), 2, shape, owner, strides, nb::dtype<T>()};
}

auto payload(nb::handle_t<gelex::DenseReader> self, std::string_view identifier)
    -> PayloadArray
{
    const auto& reader = nb::cast<const gelex::DenseReader&>(self);
    switch (reader.info(identifier).type)
    {
        case gelex::BinaryType::float64:
            return payload_view<double>(reader, identifier, self);
        case gelex::BinaryType::float32:
            return payload_view<float>(reader, identifier, self);
        case gelex::BinaryType::uint8:
            return payload_view<std::uint8_t>(reader, identifier, self);
    }
    throw gelex::GelexException("unsupported payload dtype");
}

auto register_reader(nb::module_& m) -> void
{
    nb::class_<gelex::DenseReader>(
        m,
        "DenseReader",
        "Memory-mapped reader for gelex binary containers such as the MCMC "
        ".draws output. Payloads are exposed as read-only, column-major "
        "(rows, columns) NumPy views that alias the mapped file.")
        .def(nb::init<std::string_view>(), nb::arg("path"))
        .def("__len__", &gelex::DenseReader::size)
        .def(
            "__contains__",
            &gelex::DenseReader::contains,
            nb::arg("identifier"))
        .def("__getitem__", &payload, nb::arg("identifier"))
        .def(
            "info",
            [](const gelex::DenseReader& reader, std::string_view identifier)
            { return reader.info(identifier); },
            nb::arg("identifier"))
        .def("payloads", &gelex::DenseReader::payloads)
        .def(
            "keys",
            [](const gelex::DenseReader& reader)
            {
                std::vector<std::string> keys;
                for (auto& info : reader.payloads())
                {
                    keys.push_back(std::move(info.identifier));
                }
                return keys;
            });
}

// ---- CSC reader

using AnySparseMatrix = std::variant<
    gelex::CscReader::sparse_matrix_type<double>,
    gelex::CscReader::sparse_matrix_type<float>,
    gelex::CscReader::sparse_matrix_type<std::uint8_t>>;

// The mapped arrays cannot be handed to SciPy without a copy, so the matrix
// is materialised; nanobind converts it to scipy.sparse.csc_matrix.
auto sparse_payload(const gelex::CscReader& reader, std::string_view identifier)
    -> AnySparseMatrix
{
    switch (reader.info(identifier).type)
    {
        case gelex::BinaryType::float64:
            return reader.to_mat<double>(identifier);
        case gelex::BinaryType::float32:
            return reader.to_mat<float>(identifier);
        case gelex::BinaryType::uint8:
            return reader.to_mat<std::uint8_t>(identifier);
    }
    throw gelex::GelexException("unsupported payload dtype");
}

auto register_csc_reader(nb::module_& m) -> void
{
    nb::class_<gelex::CscReader>(
        m,
        "CscReader",
        "Memory-mapped reader for gelex CSC containers such as the MCMC "
        ".draws.csc output. Matrices are returned as scipy.sparse.csc_matrix "
        "copies of shape (rows, columns).")
        .def(nb::init<std::string_view>(), nb::arg("path"))
        .def("__len__", &gelex::CscReader::size)
        .def("__contains__", &gelex::CscReader::contains, nb::arg("identifier"))
        .def("__getitem__", &sparse_payload, nb::arg("identifier"))
        .def(
            "info",
            [](const gelex::CscReader& reader, std::string_view identifier)
            { return reader.info(identifier); },
            nb::arg("identifier"))
        .def("nnz", &gelex::CscReader::nnz, nb::arg("identifier"))
        .def("payloads", &gelex::CscReader::payloads)
        .def(
            "keys",
            [](const gelex::CscReader& reader)
            {
                std::vector<std::string> keys;
                for (auto& info : reader.payloads())
                {
                    keys.push_back(std::move(info.identifier));
                }
                return keys;
            });

    m.def(
        "sparse_draws_path",
        &gelex::sparse_draws_path,
        nb::arg("draws_path"),
        "Path of the CSC companion written next to a .draws file.");
}

}  // namespace

void register_io(nb::module_& m)
{
    register_format(m);
    register_writer(m);
    register_reader(m);
    register_csc_reader(m);
}

}  // namespace gelex_py
