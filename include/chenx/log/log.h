#pragma once
#include <armadillo>
#include <chrono>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <string>
#include <vector>

namespace chenx
{

using namespace arma;
template <typename eT>
struct LogEntry
{
    int iterations;
    std::string_view method;
    double logLikelihood;
    Col<eT> variances;
    double timeCost;
};

template <typename eT>
class Logger
{
  public:
    // Call this to start the timer
    Logger(std::vector<std::string> var_names, bool verbose)
        : _var_names{var_names}, verbose{verbose}
    {
        for (auto& var_name : _var_names)
        {
            var_name = fmt::format("V({})", var_name);
        }
        initTime = std::chrono::high_resolution_clock::now();
    }
    void start()
    {
        startTime = std::chrono::high_resolution_clock::now();
    }

    // Add a log entry, and automatically calculate the time cost since start()
    // was called
    void log(
        int iterations,
        std::string_view method,
        double logLikelihood,
        const Col<eT>& variances);
    void end();

  private:
    std::vector<LogEntry<eT>> logs;
    bool verbose;
    std::chrono::time_point<std::chrono::high_resolution_clock> startTime;
    std::vector<std::string> _var_names;
    std::chrono::time_point<std::chrono::high_resolution_clock> initTime;

    // Column width settings
    int logLikelihoodWidth = 6;
    int varianceWidth = 4; // Width for each variance
    int timeCostWidth = 7;
    int totalWidth = 0;

    // Update the column widths based on the longest element
    void updateColumnWidths(
        double logLikelihood,
        const Col<double>& variances,
        double timeCost);
    // Print the header using fmt
    void printHeader();
    // Print a log entry using fmt
    void printEntry(const LogEntry<eT>& entry) const;
};

} // namespace chenx
#include "log_impl.h"
