#ifndef GELEX_DATA_GRM_H_
#define GELEX_DATA_GRM_H_

#include <filesystem>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/bed_pipe.h"
#include "gelex/data/sample_manager.h"

namespace gelex
{

struct GrmResult
{
    Eigen::MatrixXd grm;
    double denominator;
};

class GRM
{
   public:
    explicit GRM(std::filesystem::path bed_path);
    GRM(const GRM&) = delete;
    GRM(GRM&&) noexcept = default;
    GRM& operator=(const GRM&) = delete;
    GRM& operator=(GRM&&) noexcept = default;

    ~GRM() = default;

    template <typename CodePolicy>
    auto compute(
        Eigen::Index chunk_size,
        bool additive,
        std::function<void(Eigen::Index, Eigen::Index)> progress_callback
        = nullptr) -> GrmResult;

    template <typename CodePolicy>
    auto compute(
        const std::vector<std::pair<Eigen::Index, Eigen::Index>>& ranges,
        Eigen::Index chunk_size,
        bool additive,
        const std::function<void(Eigen::Index, Eigen::Index)>& progress_callback
        = nullptr) -> GrmResult;

    [[nodiscard]] auto sample_ids() const -> const std::vector<std::string>&
    {
        return sample_manager_->common_ids();
    }

    [[nodiscard]] auto num_snps() const -> Eigen::Index
    {
        return bed_.num_snps();
    }

   private:
    std::shared_ptr<SampleManager> sample_manager_;
    BedPipe bed_;

    static auto update_grm(
        Eigen::Ref<Eigen::MatrixXd> grm,
        const Eigen::Ref<const Eigen::MatrixXd>& genotype) -> void;
};

template <typename CodePolicy>
auto GRM::compute(
    Eigen::Index chunk_size,
    bool add,
    std::function<void(Eigen::Index, Eigen::Index)> progress_callback)
    -> GrmResult
{
    return compute<CodePolicy>(
        {{0, bed_.num_snps()}}, chunk_size, add, progress_callback);
}

template <typename CodePolicy>
auto GRM::compute(
    const std::vector<std::pair<Eigen::Index, Eigen::Index>>& ranges,
    Eigen::Index chunk_size,
    bool add,
    const std::function<void(Eigen::Index, Eigen::Index)>& progress_callback)
    -> GrmResult
{
    const Eigen::Index n = bed_.num_samples();
    Eigen::MatrixXd grm = Eigen::MatrixXd::Zero(n, n);

    Eigen::Index total_snps_to_process = 0;
    for (const auto& [start, end] : ranges)
    {
        total_snps_to_process += (end - start);
    }

    Eigen::Index processed_snps = 0;
    for (const auto& [range_start, range_end] : ranges)
    {
        for (Eigen::Index start_col = range_start; start_col < range_end;
             start_col += chunk_size)
        {
            const Eigen::Index end_col
                = std::min(start_col + chunk_size, range_end);
            Eigen::MatrixXd genotype_chunk
                = bed_.load_chunk(start_col, end_col);

            CodePolicy{}(genotype_chunk, add);
            update_grm(grm, genotype_chunk);

            processed_snps += (end_col - start_col);
            if (progress_callback)
            {
                progress_callback(processed_snps, total_snps_to_process);
            }
        }
    }

    double denominator = grm.trace() / static_cast<double>(n);

    return {grm, denominator};
}
}  // namespace gelex

#endif  // GELEX_DATA_GRM_H_
