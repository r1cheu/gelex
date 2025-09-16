#include "gelex/utils/utils.h"

#include <chrono>
#include <string>

#include <fmt/chrono.h>
#include <fmt/format.h>
#include <spdlog/common.h>
#include <spdlog/formatter.h>
#include <spdlog/logger.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include "Eigen/Core"

#include <omp.h>
#include <Eigen/Dense>
#include <utility>

#include <armadillo>

namespace gelex
{
using Eigen::Ref;

using Eigen::MatrixXd;
using Eigen::RowVectorXd;
using Eigen::VectorXd;

RowVectorXd centralize(Ref<MatrixXd> x)
{
    RowVectorXd x_mean = RowVectorXd::Zero(x.cols());
#pragma omp parallel for default(none) shared(x, x_mean)
    for (int i = 0; i < x.cols(); ++i)
    {
        auto col = x.col(i);
        double mean_i = col.mean();
        x_mean(i) = mean_i;
        col.array() -= mean_i;
    }
    return x_mean;
}

std::pair<RowVectorXd, RowVectorXd> standardize(Ref<MatrixXd> x)
{
    RowVectorXd x_mean = RowVectorXd::Zero(x.cols());
    RowVectorXd x_stddev = RowVectorXd::Zero(x.cols());
#pragma omp parallel for default(none) shared(x, x_mean, x_stddev)
    for (int i = 0; i < x.cols(); ++i)
    {
        auto col = x.col(i);
        double mean = col.mean();
        double stddev = std::sqrt(
            (col.array() - mean).square().sum()
            / static_cast<double>(col.size() - 1));
        x_mean(i) = mean;
        x_stddev(i) = stddev;
        col.array() -= mean;
        if (stddev != 0)
        {
            col.array() /= stddev;
        }
    }
    return {x_mean, x_stddev};
}

Timer::~Timer()
{
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start_;
    elapsed_time_ = duration.count();
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
