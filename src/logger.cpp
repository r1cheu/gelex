#include <memory>
#include <string>
#include <vector>

#include <spdlog/logger.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

namespace chenx
{

static const std::string logger_name = "defalut";
std::shared_ptr<spdlog::logger> setup_logger(
    std::vector<spdlog::sink_ptr> sinks)
{
    auto logger = spdlog::get(logger_name);
    if (not logger)
    {
        if (sinks.size() > 0)
        {
            logger = std::make_shared<spdlog::logger>(
                logger_name, std::begin(sinks), std::end(sinks));
            spdlog::register_logger(logger);
        }
        else
        {
            logger = spdlog::stdout_color_st(logger_name);
        }
    }

    return logger;
}
}  // namespace chenx
