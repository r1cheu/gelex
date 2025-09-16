#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <spdlog/logger.h>
#include <memory>

#include "../src/logger/logger_utils.h"

namespace bind
{
namespace nb = nanobind;
using nb::literals::operator""_a;
namespace gx = gelex;

void logger_bindings(nb::module_& module)
{
    nb::class_<spdlog::logger>(module, "Logger")
        .def(
            "info",
            [](spdlog::logger& self, const std::string& msg)
            { self.info(msg); },
            "msg"_a)
        .def(
            "warn",
            [](spdlog::logger& self, const std::string& msg)
            { self.warn(msg); },
            "msg"_a)
        .def(
            "error",
            [](spdlog::logger& self, const std::string& msg)
            { self.error(msg); },
            "msg"_a)
        .def(
            "debug",
            [](spdlog::logger& self, const std::string& msg)
            { self.debug(msg); },
            "msg"_a)
        .def(
            "critical",
            [](spdlog::logger& self, const std::string& msg)
            { self.critical(msg); },
            "msg"_a);

    module.def(
        "get_logger",
        []() -> std::shared_ptr<spdlog::logger>
        { return gx::detail::Logger::logger(); });
}

}  // namespace bind
