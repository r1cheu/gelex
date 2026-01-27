#ifndef GELEX_LOGGER_LOCO_REML_LOGGER_H_
#define GELEX_LOGGER_LOCO_REML_LOGGER_H_

#include "reml_logger_base.h"

#include <string>
#include <vector>

#include "../src/types/freq_effect.h"

namespace gelex::detail
{

struct VarianceComponent
{
    freq::GrmType type{freq::GrmType::Unknown};
    double variance{};
    double heritability{};
};

struct LocoRemlResult
{
    std::string chr_name;
    double loglike{};
    std::vector<VarianceComponent> genetic;
    double residual_variance{};
    bool converged{true};
    double elapsed{};

    auto total_h2() const -> double
    {
        double sum = 0.0;
        for (const auto& g : genetic)
        {
            sum += g.heritability;
        }
        return sum;
    }
};

class LocoRemlLogger : public RemlLoggerBase
{
   public:
    explicit LocoRemlLogger(std::string chr_name);

    void set_verbose(bool verbose) override;

    void log_iteration(
        size_t iter,
        double loglike,
        const FreqState& state,
        double time_cost) override;

    void log_results(
        const FreqModel& model,
        const FreqState& state,
        double loglike,
        bool converged,
        size_t iter_count,
        size_t max_iter,
        double elapsed) override;

    auto result() const -> const LocoRemlResult& { return result_; }

   private:
    std::string chr_name_;
    LocoRemlResult result_;
};

void print_loco_reml_summary(const std::vector<LocoRemlResult>& results);

}  // namespace gelex::detail

#endif  // GELEX_LOGGER_LOCO_REML_LOGGER_H_
