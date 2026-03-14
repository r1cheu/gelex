#include "predict_reporter.h"

#include <fmt/base.h>

#include "gelex/infra/logger.h"

namespace gelex::cli
{

PredictReporter::PredictReporter() : logger_(gelex::logging::get()) {}

auto PredictReporter::on_event(const PredictParamsLoadedEvent& event) const
    -> void
{
    if (logger_)
    {
        logger_->info(
            "Parameters loaded: {} SNP effects, {} covariate terms, "
            "geno-method={}",
            event.num_snp_effects,
            event.num_covar_terms,
            fmt::format("{}", event.geno_method));
    }
}

auto PredictReporter::on_event(const PredictSnpSelectionEvent& event) const
    -> void
{
    if (logger_)
    {
        logger_->info(
            "SNP selection: {}/{} matched, {} missing",
            event.num_matched,
            event.num_total,
            event.num_missing);
    }
}

auto PredictReporter::on_event(const PredictDataLoadedEvent& event) const
    -> void
{
    if (logger_)
    {
        logger_->info(
            "Data loaded: {} samples, {} SNPs, {} covariate terms "
            "(design matrix {}x{})",
            event.num_samples,
            event.num_snps,
            event.num_covar_terms,
            event.design_rows,
            event.design_cols);
    }
}

auto PredictReporter::on_event(const PredictResultsWrittenEvent& event) const
    -> void
{
    if (logger_)
    {
        logger_->info(
            "Results written to '{}' ({} samples)",
            event.output_path,
            event.num_samples);
    }
}

}  // namespace gelex::cli
