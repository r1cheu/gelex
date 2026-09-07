/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <cstddef>
#include <cstdint>
#include <fmt/format.h>
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

#include "gelex/exception.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/binary_reader.h"
#include "gelex/io/binary_writer.h"

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
        .value("int32", gelex::BinaryType::int32)
        .value("uint8", gelex::BinaryType::uint8);

    nb::class_<gelex::PayloadInfo>(m, "PayloadInfo")
        .def_prop_ro(
            "identifier",
            [](const gelex::PayloadInfo& info) { return info.identifier; })
        .def_prop_ro(
            "type",
            [](const gelex::PayloadInfo& info) { return info.descriptor.type; })
        .def_prop_ro(
            "shape",
            [](const gelex::PayloadInfo& info)
            {
                return std::pair{
                    info.descriptor.shape[0], info.descriptor.shape[1]};
            })
        .def(
            "__repr__",
            [](const gelex::PayloadInfo& info)
            {
                return fmt::format(
                    "PayloadInfo(identifier='{}', shape=({}, {}))",
                    info.identifier,
                    info.descriptor.shape[0],
                    info.descriptor.shape[1]);
            });
}

// ---- writer

// Inputs are not converted: a dtype or layout mismatch is a TypeError rather
// than a silent copy. Column-major (rows, columns) is the on-disk layout, so
// a whole payload is taken as an F-contiguous matrix.
template <typename T>
using ColumnArray
    = nb::ndarray<const T, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
template <typename T>
using MatrixArray
    = nb::ndarray<const T, nb::ndim<2>, nb::f_contig, nb::device::cpu>;

template <gelex::detail::SupportedDtype T>
auto register_payload_writer(nb::module_& m, const char* name) -> void
{
    using Writer = gelex::PayloadWriter<T>;
    nb::class_<Writer>(
        m,
        name,
        "Handle to one reserved payload; append() adds one column (one draw) "
        "and write() stores the whole payload at once.")
        .def(
            "append",
            [](Writer& self, const ColumnArray<T>& column)
            {
                if (column.shape(0) != self.rows())
                {
                    throw gelex::GelexException(
                        "append expects one value per payload row");
                }
                self.append(std::span<const T>{column.data(), column.shape(0)});
            },
            nb::arg("column").noconvert())
        .def(
            "write",
            [](Writer& self, const MatrixArray<T>& values)
            {
                if (values.shape(0) != self.rows())
                {
                    throw gelex::GelexException(
                        "write expects a (rows, columns) array with rows "
                        "matching the payload");
                }
                self.write(
                    std::span<const T>{
                        values.data(), values.shape(0) * values.shape(1)});
            },
            nb::arg("values").noconvert())
        .def_prop_ro("identifier", &Writer::identifier)
        .def_prop_ro("rows", &Writer::rows);
}

using AnyPayloadWriter = std::variant<
    gelex::PayloadWriter<double>,
    gelex::PayloadWriter<float>,
    gelex::PayloadWriter<std::int32_t>,
    gelex::PayloadWriter<std::uint8_t>>;

auto reserve(
    gelex::BinaryWriter& writer,
    std::string_view identifier,
    gelex::BinaryType type,
    gelex::BinaryShape shape) -> AnyPayloadWriter
{
    switch (type)
    {
        case gelex::BinaryType::float64:
            return writer.reserve<double>(identifier, shape);
        case gelex::BinaryType::float32:
            return writer.reserve<float>(identifier, shape);
        case gelex::BinaryType::int32:
            return writer.reserve<std::int32_t>(identifier, shape);
        case gelex::BinaryType::uint8:
            return writer.reserve<std::uint8_t>(identifier, shape);
    }
    throw gelex::GelexException("unsupported payload dtype");
}

auto register_writer(nb::module_& m) -> void
{
    register_payload_writer<double>(m, "PayloadWriterF64");
    register_payload_writer<float>(m, "PayloadWriterF32");
    register_payload_writer<std::int32_t>(m, "PayloadWriterI32");
    register_payload_writer<std::uint8_t>(m, "PayloadWriterU8");

    nb::class_<gelex::BinaryWriter>(
        m,
        "BinaryWriter",
        "Writer for gelex binary containers. Reserve payloads with a dtype and "
        "(rows, columns) shape, fill them column by column, then close() (or "
        "leave the with-block) to finalise the file.")
        .def(nb::init<std::string_view>(), nb::arg("path"))
        .def(
            "reserve",
            &reserve,
            nb::arg("identifier"),
            nb::arg("type"),
            nb::arg("shape"),
            nb::keep_alive<0, 1>())
        .def("close", &gelex::BinaryWriter::close)
        .def_prop_ro("is_open", &gelex::BinaryWriter::is_open)
        .def(
            "__enter__",
            [](nb::handle_t<gelex::BinaryWriter> self) { return self; })
        .def(
            "__exit__",
            [](gelex::BinaryWriter& self, nb::args) { self.close(); },
            nb::arg("args"));
}

// ---- reader

using PayloadArray = nb::ndarray<nb::numpy, nb::ro>;

// The array aliases the reader's memory map: the reader object is registered
// as the array owner so Python keeps it alive for as long as any view exists.
template <gelex::detail::SupportedDtype T>
auto payload_view(
    const gelex::BinaryReader& reader,
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

auto payload(
    nb::handle_t<gelex::BinaryReader> self,
    std::string_view identifier) -> PayloadArray
{
    const auto& reader = nb::cast<const gelex::BinaryReader&>(self);
    switch (reader.info(identifier).descriptor.type)
    {
        case gelex::BinaryType::float64:
            return payload_view<double>(reader, identifier, self);
        case gelex::BinaryType::float32:
            return payload_view<float>(reader, identifier, self);
        case gelex::BinaryType::int32:
            return payload_view<std::int32_t>(reader, identifier, self);
        case gelex::BinaryType::uint8:
            return payload_view<std::uint8_t>(reader, identifier, self);
    }
    throw gelex::GelexException("unsupported payload dtype");
}

auto register_reader(nb::module_& m) -> void
{
    nb::class_<gelex::BinaryReader>(
        m,
        "BinaryReader",
        "Memory-mapped reader for gelex binary containers such as the MCMC "
        ".draws output. Payloads are exposed as read-only, column-major "
        "(rows, columns) NumPy views that alias the mapped file.")
        .def(nb::init<std::string_view>(), nb::arg("path"))
        .def("__len__", &gelex::BinaryReader::size)
        .def(
            "__contains__",
            &gelex::BinaryReader::contains,
            nb::arg("identifier"))
        .def("__getitem__", &payload, nb::arg("identifier"))
        .def(
            "info",
            [](const gelex::BinaryReader& reader, std::string_view identifier)
            { return reader.info(identifier); },
            nb::arg("identifier"))
        .def("payloads", &gelex::BinaryReader::payloads)
        .def(
            "keys",
            [](const gelex::BinaryReader& reader)
            {
                std::vector<std::string> keys;
                for (auto& info : reader.payloads())
                {
                    keys.push_back(std::move(info.identifier));
                }
                return keys;
            });
}

}  // namespace

void register_io(nb::module_& m)
{
    register_format(m);
    register_writer(m);
    register_reader(m);
}

}  // namespace gelex_py
