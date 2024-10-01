#pragma once
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>  // To enable printing via ostream
#include <armadillo>
#include <chrono>
#include <string_view>
#include <vector>

namespace chenx {
namespace optim {

using namespace arma;
template <typename eT>
struct LogEntry {
    int iterations;
    std::string_view method;
    double logLikelihood;
    Col<eT> variances;
    double timeCost;
};

template <typename eT>
class Logger {
   public:
    // Call this to start the timer
    void start() { startTime = std::chrono::high_resolution_clock::now(); }

    // Add a log entry, and automatically calculate the time cost since start() was called
    void log(int iterations, std::string_view method, double logLikelihood, const Col<eT>& variances);

   private:
    std::vector<LogEntry<eT>> logs;
    std::chrono::time_point<std::chrono::high_resolution_clock> startTime;

    // Column width settings
    int iterWidth = 7;
    int methodWidth = 10;
    int logLikelihoodWidth = 11;
    int varianceWidth = 10;  // Width for each variance
    int timeCostWidth = 9;

    // Update the column widths based on the longest element
    void updateColumnWidths(
        int iterations,
        std::string_view method,
        double logLikelihood,
        const Col<double>& variances,
        double timeCost);
    // Print the header using fmt
    void printHeader() const;
    // Print a log entry using fmt
    void printEntry(const LogEntry<eT>& entry) const;
};

}  // namespace optim
}  // namespace chenx
#include "log_impl.h"
