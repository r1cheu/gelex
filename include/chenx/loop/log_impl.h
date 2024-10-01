#pragma once
#include "log.h"

namespace chenx {
namespace optim {

using namespace arma;
template <typename eT>
void Logger<eT>::log(int iterations, std::string_view method, double logLikelihood, const Col<eT>& variances) {
    // Calculate elapsed time in seconds
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = endTime - startTime;
    double timeCost = elapsed.count();
    if (iterations == 1) {
        updateColumnWidths(iterations, method, logLikelihood, variances, timeCost);
        printHeader();
    }
    logs.push_back({iterations, method, logLikelihood, variances, timeCost});
    printEntry(logs.back());
}

template <typename eT>
void Logger<eT>::updateColumnWidths(
    // Update the column widths based on the longest element
    int iterations,
    std::string_view method,
    double logLikelihood,
    const Col<double>& variances,
    double timeCost) {
    iterWidth = std::max(iterWidth, static_cast<int>(std::to_string(iterations).length()));
    logLikelihoodWidth = std::max(logLikelihoodWidth, static_cast<int>(fmt::format("{:.4f}", logLikelihood).length()));

    // Adjust variance width for each variance value

    int sum_variances = 0;
    sum_variances += variances.size() - 1;  // Add space for the separator
    for (const eT& var : variances) {
        sum_variances += static_cast<int>(fmt::format("{:.4f}", var).length());
    }
    varianceWidth = std::max(varianceWidth, sum_variances);
    timeCostWidth = std::max(timeCostWidth, static_cast<int>(fmt::format("{:.4f}", timeCost).length()));
}

// Print the header using fmt
template <typename eT>
void Logger<eT>::printHeader() const {
    fmt::print(
        "| {:^{}} | {:^{}} | {:^{}} | {:^{}} | {:^{}} |\n",
        "Iter.",
        iterWidth,
        "Method",
        methodWidth,
        "LogL.",
        logLikelihoodWidth,
        "Variance",
        varianceWidth,
        "Time(s)",
        timeCostWidth);
    int totalWidth = iterWidth + methodWidth + logLikelihoodWidth + varianceWidth + timeCostWidth
                     + 16;  // 16 is the number of spaces and separators
    fmt::print("{:-<{}}\n", "", totalWidth);
}
template <typename eT>
void Logger<eT>::printEntry(const LogEntry<eT>& entry) const {
    fmt::print(
        "| {:^{}} | {:^{}} | {:^{}} |",
        entry.iterations,
        iterWidth,
        entry.method,
        methodWidth,
        fmt::format("{:.4f}", entry.logLikelihood),
        logLikelihoodWidth);

    // Print multiple variances
    for (const eT& var : entry.variances) {
        fmt::print(" {}", fmt::format("{:.4f}", var));
    }
    fmt::print(" | {:^{}} |\n", fmt::format("{:.3f}", entry.timeCost), timeCostWidth);
};
}  // namespace optim
}  // namespace chenx
