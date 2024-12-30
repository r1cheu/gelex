#pragma once
#include <memory>
#include <string>
#include <vector>

#include <spdlog/common.h>
#include <spdlog/logger.h>

namespace chenx
{
extern const std::string logger_name{"default"};
std::shared_ptr<spdlog::logger> setup_logger(
    std::vector<spdlog::sink_ptr> sinks);
}  // namespace chenx
