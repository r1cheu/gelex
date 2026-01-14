#ifndef GELEX_DATA_GRM_H_
#define GELEX_DATA_GRM_H_

#include <filesystem>

#include <Eigen/Core>

#include "gelex/data/bed_pipe.h"

namespace gelex
{

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
    auto compute(Eigen::Index chunk_size, bool additive) -> Eigen::MatrixXd;

   private:
    BedPipe bed_;
    static auto update_grm(
        Eigen::Ref<Eigen::MatrixXd> grm,
        const Eigen::Ref<const Eigen::MatrixXd>& genotype) -> void;
};

template <typename CodePolicy>
auto GRM::compute(Eigen::Index chunk_size, bool add) -> Eigen::MatrixXd
{
    const Eigen::Index n = bed_.num_samples();
    Eigen::MatrixXd grm = Eigen::MatrixXd::Zero(n, n);

    const Eigen::Index num_snps = bed_.num_snps();
    for (Eigen::Index start_col = 0; start_col < num_snps;
         start_col += chunk_size)
    {
        const Eigen::Index end_col = std::min(start_col + chunk_size, num_snps);
        Eigen::MatrixXd genotype_chunk = bed_.load_chunk(start_col, end_col);

        CodePolicy{}(genotype_chunk, add);
        update_grm(grm, genotype_chunk);
    }

    grm /= grm.trace() / static_cast<double>(n);

    return grm;
}
}  // namespace gelex

#endif  // GELEX_DATA_GRM_H_
