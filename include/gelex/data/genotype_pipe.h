#pragma once

#include <expected>
#include <filesystem>
#include <memory>
#include <vector>

#include <Eigen/Core>

#include "../src/data/binary_matrix_writer.h"
#include "../src/data/snp_stats_writer.h"
#include "gelex/data/bed_pipe.h"
#include "gelex/error.h"

namespace gelex
{

class GenotypePipe
{
   public:
    static auto create(
        std::unique_ptr<BedPipe> bed_pipe,
        std::filesystem::path matrix_path,
        std::filesystem::path stats_path,
        std::optional<std::unordered_map<std::string, Eigen::Index>> id_map
        = std::nullopt,
        double monomorphic_threshold = 1e-6)
        -> std::expected<GenotypePipe, Error>;

    GenotypePipe(const GenotypePipe&) = delete;
    GenotypePipe(GenotypePipe&&) noexcept = default;
    GenotypePipe& operator=(const GenotypePipe&) = delete;
    GenotypePipe& operator=(GenotypePipe&&) noexcept = default;
    ~GenotypePipe();

    auto process(
        size_t chunk_size = 1000,
        std::optional<std::unordered_map<std::string, Eigen::Index>> id_map
        = std::nullopt) -> std::expected<void, Error>;

    const std::vector<double>& means() const noexcept { return means_; }
    const std::vector<double>& stddevs() const noexcept { return stddevs_; }
    const std::vector<size_t>& monomorphic_indices() const noexcept
    {
        return monomorphic_indices_;
    }

    Eigen::Index num_samples() const noexcept { return num_samples_; }
    Eigen::Index num_variants() const noexcept { return num_variants_; }
    size_t num_processed_variants() const noexcept { return means_.size(); }

   private:
    GenotypePipe(
        std::unique_ptr<BedPipe> bed_pipe,
        std::filesystem::path matrix_path,
        std::filesystem::path stats_path,
        std::optional<std::unordered_map<std::string, Eigen::Index>> id_map,
        double monomorphic_threshold);

    auto process_chunk(Eigen::MatrixXd&& chunk, size_t global_start)
        -> std::expected<void, Error>;

    auto finalize() -> std::expected<void, Error>;

    std::unique_ptr<BedPipe> bed_pipe_;
    std::filesystem::path matrix_path_;
    std::filesystem::path stats_path_;
    double monomorphic_threshold_;

    Eigen::Index num_samples_;
    Eigen::Index num_variants_;
    std::optional<std::unordered_map<std::string, Eigen::Index>> id_map_;

    std::vector<double> means_;
    std::vector<double> stddevs_;
    std::vector<size_t> monomorphic_indices_;

    std::unique_ptr<detail::BinaryMatrixWriter> matrix_writer_;
    std::unique_ptr<detail::SnpStatsWriter> stats_writer_;
};

}  // namespace gelex
