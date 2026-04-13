#include "predict_reporter.h"

#include <fmt/format.h>

#include <fmt/base.h>

#include "config.h"
#include "gelex/infra/logger.h"
#include "gelex/infra/logging/formatter.h"

namespace gelex::cli
{

PredictReporter::PredictReporter() : logger_(gelex::logging::get()) {}

auto PredictReporter::on_event(const PredictBannerEvent& /*event*/) const
    -> void
{
    logger_->info(gelex::command_banner(PROJECT_VERSION, "Genomic Prediction"));
    logger_->info("");
}

auto PredictReporter::on_event(const PredictParamsLoadedEvent& event) const
    -> void
{
    logger_->info(gelex::section("[Config]"));
    logger_->info("  {:<12}: {}", "bfile", event.bfile_prefix);
    logger_->info("  {:<12}: {}", "gfile", event.gfile_prefix);
    logger_->info(
        "  {:<12}: {}", "Geno method", fmt::format("{}", event.geno_method));
    logger_->info("");
}

auto PredictReporter::on_event(const PredictSnpSelectionEvent& event) const
    -> void
{
    logger_->info(gelex::section("[SNP Alignment]"));
    logger_->info(
        "   {:<13}: {}/{}", "Matched", event.num_matched, event.num_total);
    logger_->info("   {:<13}: {}", "Missing", event.num_missing);
    logger_->info("   {:<13}: {}", "Mismatched", event.num_mismatched);

    if (event.num_mismatched > 0)
    {
        std::string plink_hint = fmt::format(
            "plink2 --bfile {} --alt1-allele {} 4 1 1 --make-bed --out "
            "<output>",
            event.bfile_path,
            event.snp_effect_path);
        logger_->warn(
            "Allele mismatch detected for {} SNPs. "
            "To fix, run:\n  {}",
            event.num_mismatched,
            plink_hint);
    }
    logger_->info("");
}

auto PredictReporter::on_event(const PredictDataLoadedEvent& event) const
    -> void
{
    logger_->info(gelex::section("[Dataset Summary]"));
    logger_->info("   {:<13}: {} samples", "Samples", event.num_samples);
    logger_->info("   {:<13}: {} markers", "SNPs", event.num_snps);
    logger_->info("   {:<13}: {}", "Covariates", event.num_covar_terms);
    logger_->info("");
}

auto PredictReporter::on_event(const PredictResultsWrittenEvent& event) const
    -> void
{
    logger_->info(
        gelex::success(
            "Results saved to '{}' ({} samples)",
            event.output_path,
            event.num_samples));
}

}  // namespace gelex::cli
