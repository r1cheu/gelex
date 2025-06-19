#include "gelex/utils/utils.h"

#include <chrono>
#include <iostream>
#include <memory>
#include <string>

#include <fmt/chrono.h>
#include <fmt/format.h>
#include <spdlog/common.h>
#include <spdlog/formatter.h>
#include <spdlog/logger.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include "spdlog/pattern_formatter.h"

#include <armadillo>

namespace gelex
{

arma::dvec centralize(arma::dmat& x)
{
    arma::dvec x_mean(x.n_cols, arma::fill::zeros);
#pragma omp parallel for default(none) shared(x, x_mean)
    for (size_t i = 0; i < x.n_cols; ++i)
    {
        arma::dvec col = x.unsafe_col(i);
        double mean_i = arma::mean(col);
        x_mean.at(i) = mean_i;
        col -= mean_i;
    }
    return x_mean;
}

std::pair<arma::dvec, arma::dvec> standradize(arma::dmat& x)
{
    arma::dvec x_mean(x.n_cols, arma::fill::zeros);
    arma::dvec x_stddev(x.n_cols, arma::fill::zeros);
#pragma omp parallel for default(none) shared(x, x_mean, x_stddev)
    for (size_t i = 0; i < x.n_cols; ++i)
    {
        arma::dvec col = x.unsafe_col(i);
        double mean = arma::mean(col);
        double stddev = arma::stddev(col);

        x_mean.at(i) = mean;
        x_stddev.at(i) = stddev;
        col -= mean;
        if (stddev != 0)
        {
            col /= stddev;
        }
    }
    return {x_mean, x_stddev};
}
LevelFormatter::LevelFormatter()
{
    info_formatter_ = std::make_unique<spdlog::pattern_formatter>("%v");
    default_formatter_
        = std::make_unique<spdlog::pattern_formatter>("[%^%l%$] %v");
}

void LevelFormatter::format(
    const spdlog::details::log_msg& msg,
    spdlog::memory_buf_t& dest)
{
    if (msg.level == spdlog::level::info)
    {
        info_formatter_->format(msg, dest);
    }
    else
    {
        default_formatter_->format(msg, dest);
    }
}

std::unique_ptr<spdlog::formatter> LevelFormatter::clone() const
{
    return std::make_unique<LevelFormatter>();
}

Timer::~Timer()
{
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start_;
    elapsed_time_ = duration.count();
}

Logger::Logger()
{
    try
    {
        auto console_sink
            = std::make_shared<spdlog::sinks::stdout_color_sink_st>();
        console_sink->set_level(spdlog::level::debug);
        console_sink->set_formatter(std::make_unique<LevelFormatter>());

        auto now = std::chrono::system_clock::now();
        auto file_sink
            = std::make_shared<spdlog::sinks::basic_file_sink_st>(fmt::format(
                "gelex_{:%F_%H-%M-%S}.log",
                fmt::localtime(std::chrono::system_clock::to_time_t(now))));

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

std::string compute_time_left(
    const std::chrono::high_resolution_clock::time_point& start,
    size_t iter,
    size_t total_iter)
{
    // Handle edge cases
    if (iter == 0)
    {
        return "--:--:--";
    }

    if (iter >= total_iter)
    {
        return "00:00:00";
    }

    auto elapsed = std::chrono::high_resolution_clock::now() - start;
    auto time_per_iter = elapsed / iter;
    auto time_left = time_per_iter * (total_iter - iter);

    auto seconds_left
        = std::chrono::duration_cast<std::chrono::seconds>(time_left);
    return fmt::format("{:%H:%M:%S}", seconds_left);
}
}  // namespace gelex
