#include "chenx/logger.h"

#include <iostream>
#include <memory>
#include <string>

#include <spdlog/logger.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/rotating_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

namespace chenx
{

// 构造函数实现，初始化 spdlog
Logger::Logger()
{
    try
    {
        auto console_sink
            = std::make_shared<spdlog::sinks::stdout_color_sink_st>();
        console_sink->set_level(spdlog::level::debug);
        console_sink->set_pattern("[%^%l%$] %v");

        auto file_sink = std::make_shared<spdlog::sinks::rotating_file_sink_st>(
            "chenx.log", 1024 * 1024, 5, false);  // maximum 5 files, 1MB each
        file_sink->set_level(spdlog::level::trace);
        file_sink->set_pattern("[%Y-%m-%d %H:%M:%S] [%l] %v");

        logger_ = std::make_shared<spdlog::logger>(
            "default", spdlog::sinks_init_list{console_sink, file_sink});
        spdlog::register_logger(logger_);
        logger_->set_level(spdlog::level::debug);
    }
    catch (const spdlog::spdlog_ex& ex)
    {
        std::cout << "Log initialization failed: " << ex.what() << "\n";
    }
}

std::shared_ptr<spdlog::logger> Logger::GetSpdLogger()
{
    return logger_;
}

std::shared_ptr<spdlog::logger> Logger::logger()
{
    static Logger instance;
    return instance.GetSpdLogger();
}

}  // namespace chenx
