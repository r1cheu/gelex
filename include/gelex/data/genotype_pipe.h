#pragma once

#include <concepts>
#include <cstdint>
#include <expected>
#include <filesystem>
#include <memory>
#include <vector>

#include <gelex/barkeep.h>
#include <Eigen/Core>

#include "../src/data/binary_matrix_writer.h"
#include "../src/data/snp_stats_writer.h"
#include "../src/estimator/bayes/indicator.h"
#include "gelex/data/bed_pipe.h"
#include "gelex/data/genotype_mmap.h"
#include "gelex/data/sample_manager.h"
#include "gelex/error.h"
#include "gelex/logger.h"

namespace gelex
{

namespace bk = barkeep;

struct VariantStats
{
    double mean{0.0};
    double variance{0.0};
    bool is_monomorphic{false};
};

template <typename T>
concept VariantProcessor = requires(T processor, Eigen::VectorXd& variant) {
    { processor.process_variant(variant) } -> std::same_as<VariantStats>;
};

struct StandardizingProcessor
{
   public:
    static VariantStats process_variant(Eigen::Ref<Eigen::VectorXd> variant);
};

struct NonStandardizingProcessor
{
   public:
    static VariantStats process_variant(Eigen::Ref<Eigen::VectorXd> variant);
};

struct HardWenbergProcessor
{
   public:
    static VariantStats process_variant(Eigen::Ref<Eigen::VectorXd> variant);
};

class GenotypePipe
{
   public:
    static auto create(
        const std::filesystem::path& bed_path,
        std::shared_ptr<SampleManager> sample_manager,
        const std::filesystem::path& output_prefix,
        bool dominant = false) -> std::expected<GenotypePipe, Error>;

    GenotypePipe(const GenotypePipe&) = delete;
    GenotypePipe(GenotypePipe&&) noexcept = default;
    GenotypePipe& operator=(const GenotypePipe&) = delete;
    GenotypePipe& operator=(GenotypePipe&&) noexcept = default;
    ~GenotypePipe() = default;

    // Template method for processing with different strategies
    template <VariantProcessor Processor = StandardizingProcessor>
    auto process(size_t chunk_size = 10000) -> std::expected<GenotypeMap, Error>
    {
        auto logger = gelex::logging::get();
        global_snp_idx_ = 0;

        logger->info("");
        logger->info("Starting SNP processing for memory-map...");

        auto pbar = bk::ProgressBar(
            &global_snp_idx_,
            {.total = static_cast<uint64_t>(num_variants_),
             .format = "{bar} {value}/{total} ",
             .style = detail::BAR_STYLE});

        Processor processor;

        for (int64_t start_variant = 0; start_variant < num_variants_;)
        {
            int64_t end_variant = std::min(
                static_cast<int64_t>(start_variant + chunk_size),
                num_variants_);

            auto chunk
                = bed_pipe_.load_chunk(start_variant, end_variant, dominant_);
            if (!chunk)
            {
                return std::unexpected(chunk.error());
            }

            if (auto result
                = process_chunk(std::move(*chunk), start_variant, processor);
                !result)
            {
                return std::unexpected(result.error());
            }

            start_variant = end_variant;
        }
        pbar->done();

        return finalize();
    }

    Eigen::Index num_samples() const noexcept { return sample_size_; }
    Eigen::Index num_variants() const noexcept { return num_variants_; }
    size_t num_processed_variants() const noexcept { return means_.size(); }

   private:
    GenotypePipe(
        BedPipe&& bed_pipe,
        detail::BinaryMatrixWriter&& matrix_writer,
        detail::SnpStatsWriter&& stats_writer,
        bool dominant);

    template <VariantProcessor Processor>
    auto process_chunk(
        Eigen::MatrixXd&& chunk,
        size_t global_start,
        Processor& processor) -> std::expected<void, Error>
    {
        Eigen::MatrixXd matrix = std::move(chunk);
        const int64_t num_variants_in_chunk = matrix.cols();

        for (int64_t variant_idx = 0; variant_idx < num_variants_in_chunk;
             ++variant_idx)
        {
            global_snp_idx_++;
            const size_t global_idx = global_start + variant_idx;
            auto variant = matrix.col(variant_idx);

            VariantStats stats = processor.process_variant(variant);

            means_.push_back(stats.mean);
            variances_.push_back(stats.variance);

            if (stats.is_monomorphic)
            {
                monomorphic_indices_.push_back(
                    static_cast<int64_t>(global_idx));
            }
        }

        if (auto result = matrix_writer_.write(matrix); !result)
        {
            return std::unexpected(result.error());
        }

        return {};
    }

    auto finalize() -> std::expected<GenotypeMap, Error>;

    BedPipe bed_pipe_;
    bool dominant_{false};

    int64_t sample_size_{};
    int64_t num_variants_{};

    size_t global_snp_idx_{};

    std::vector<double> means_;
    std::vector<double> variances_;
    std::vector<int64_t> monomorphic_indices_;

    detail::BinaryMatrixWriter matrix_writer_;
    detail::SnpStatsWriter stats_writer_;
};

}  // namespace gelex
