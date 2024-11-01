#pragma once
#include <fmt/format.h>
#include "log.h"

namespace chenx {

using namespace arma;

template <typename eT>
void Logger<eT>::log(int iterations, std::string_view method, double logLikelihood, const Col<eT>& variances) {
    // Calculate elapsed time in seconds
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = endTime - startTime;
    double timeCost = elapsed.count();
    if (iterations == 1) {
        updateColumnWidths(logLikelihood, variances, timeCost);
        printHeader();
    }
    logs.push_back({iterations, method, logLikelihood, variances, timeCost});
    printEntry(logs.back());
}

template <typename eT>
void Logger<eT>::end() {
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = endTime - initTime;
    if (verbose) {
        std::cout << fmt::format("{:-<{}}", "", totalWidth) << std::endl;
    }
    std::cout << fmt::format("Total time cost: {:.3f} seconds", elapsed.count()) << std::endl;
}

template <typename eT>
void Logger<eT>::updateColumnWidths(
    // Update the column widths based on the longest element
    double logLikelihood,
    const Col<double>& variances,
    double timeCost) {
    logLikelihoodWidth = std::max(logLikelihoodWidth, static_cast<int>(fmt::format("{:.4f}", logLikelihood).length()));

    // Adjust variance width for each variance value
    for (const eT& var : variances) {
        varianceWidth = std::max(varianceWidth, static_cast<int>(fmt::format("{:.4f}", var).length()));
    }
    for (std::string_view var_name : _var_names) {
        varianceWidth = std::max(varianceWidth, static_cast<int>(var_name.length()));
    }
    timeCostWidth = std::max(timeCostWidth, static_cast<int>(fmt::format("{:.3f}", timeCost).length()));

    varianceWidth += 2;
    logLikelihoodWidth += 2;
    timeCostWidth += 2;
}

// Print the header using fmt
template <typename eT>
void Logger<eT>::printHeader() {
    std::string header = fmt::format("{:^{}}{:^{}}{:^{}}", "Iter.", 7, "Method", 8, " LogL.", logLikelihoodWidth);
    for (std::string_view var_name : _var_names) {
        header += fmt::format("{:>{}}", var_name, varianceWidth);
    }
    header += fmt::format("{:>{}}", "V(e)", varianceWidth);
    header += fmt::format("{:>{}}", "Time(s)", timeCostWidth);
    totalWidth = header.length();
    if (verbose) {
        std::cout << fmt::format("{:=<{}}", "", totalWidth) << std::endl;
        std::cout << header << std::endl;
        std::cout << fmt::format("{:-<{}}", "", totalWidth) << std::endl;
    }
}
template <typename eT>
void Logger<eT>::printEntry(const LogEntry<eT>& entry) const {
    std::string line
        = fmt::format("{:^7}{:^8}{:^{}.2f}", entry.iterations, entry.method, entry.logLikelihood, logLikelihoodWidth);
    for (const eT& var : entry.variances) {
        line += fmt::format("{:>{}.4f}", var, varianceWidth);
    }
    line += fmt::format("{:>{}.3f}", entry.timeCost, timeCostWidth);
    if (verbose) {
        std::cout << line << std::endl;
    }
};
}  // namespace chenx
