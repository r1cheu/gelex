#pragma once

#include <cstdint>
#include <expected>
#include <filesystem>
#include <vector>

#include <Eigen/Core>

#include "../src/data/binary_matrix_writer.h"
#include "../src/data/snp_stats_writer.h"
#include "gelex/data/bed_pipe.h"
#include "gelex/data/sample_manager.h"
#include "gelex/error.h"

namespace gelex
{

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
    ~GenotypePipe();

    auto process(size_t chunk_size = 10000) -> std::expected<void, Error>;

    const std::vector<double>& means() const noexcept { return means_; }
    const std::vector<double>& stddevs() const noexcept { return stddevs_; }
    const std::vector<int64_t>& monomorphic_indices() const noexcept
    {
        return monomorphic_indices_;
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

    auto process_chunk(Eigen::MatrixXd&& chunk, size_t global_start)
        -> std::expected<void, Error>;

    auto finalize() -> std::expected<void, Error>;

    BedPipe bed_pipe_;
    double monomorphic_threshold_;
    bool dominant_{false};

    int64_t sample_size_{};
    int64_t num_variants_{};

    size_t global_snp_idx_{};

    std::vector<double> means_;
    std::vector<double> stddevs_;
    std::vector<int64_t> monomorphic_indices_;

    detail::BinaryMatrixWriter matrix_writer_;
    detail::SnpStatsWriter stats_writer_;
};

}  // namespace gelex
