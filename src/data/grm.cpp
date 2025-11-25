#include "gelex/data/grm.h"

#include <algorithm>
#include <memory>

#include <Eigen/Core>

#include "gelex/data/bed_pipe.h"
#include "gelex/data/sample_manager.h"

namespace gelex
{
using Eigen::Index;

GRM::GRM(std::string_view bed_file, Index chunk_size) : chunk_size_(chunk_size)
{
    auto sample_manager = std::make_shared<SampleManager>(
        std::filesystem::path(bed_file).replace_extension(".fam"), false);
    bed_ = std::make_unique<BedPipe>(bed_file, sample_manager);

    p_major_ = Eigen::VectorXd::Zero(bed_->num_snps());
}

Eigen::MatrixXd GRM::compute(bool add)
{
    const auto n = static_cast<Eigen::Index>(
        id_map_.empty() ? bed_->num_samples() : id_map_.size());
    Eigen::MatrixXd grm = Eigen::MatrixXd::Zero(n, n);

    // Process in chunks
    const Index num_variants = bed_->num_snps();
    const Index chunk_size = 10000;  // Fixed chunk size for BedPipe

    for (Index start = 0; start < num_variants; start += chunk_size)
    {
        Index end = std::min(start + chunk_size, num_variants);
        Eigen::MatrixXd genotype = bed_->load_chunk(start, end);
        Eigen::VectorXd p_major_i = compute_p_major(genotype);
        code_method_varden(p_major_i, genotype, add);

        // Update p_major_ vector
        p_major_.segment(start, end - start) = p_major_i;

        grm.noalias() += genotype * genotype.transpose();
    }

    scale_factor_ = grm.trace() / static_cast<double>(grm.rows());
    return grm / scale_factor_;
}

CrossGRM::CrossGRM(
    std::string_view train_bed,
    Eigen::VectorXd p_major,
    double scale_factor,
    Index chunk_size)
    : scale_factor_{scale_factor}, chunk_size_{chunk_size}
{
    auto sample_manager = std::make_shared<SampleManager>(
        std::filesystem::path(train_bed).replace_extension(".fam"), false);

    bed_ = std::make_unique<BedPipe>(train_bed, sample_manager);

    if (p_major.size() != bed_->num_snps())
    {
        throw std::runtime_error(
            "p_major size does not match number of SNPs in training set.");
    }
    p_major_ = std::move(p_major);
}

Eigen::MatrixXd CrossGRM::compute(std::string_view test_bed, bool add)
{
    auto test_sample_manager = std::make_shared<SampleManager>(
        std::filesystem::path(test_bed).replace_extension(".fam"), false);

    BedPipe test_bed_pipe(test_bed, test_sample_manager);
    check_snp_consistency(test_bed_pipe);

    const Index test_n = test_bed_pipe.num_samples();
    const auto train_n = static_cast<Index>(
        id_map_.empty() ? bed_->num_samples() : id_map_.size());
    Eigen::MatrixXd grm = Eigen::MatrixXd::Zero(test_n, train_n);

    // Process in chunks
    const Index num_variants = bed_->num_snps();

    for (Index start = 0; start < num_variants; start += chunk_size_)
    {
        Index end = std::min(start + chunk_size_, num_variants);

        // Load training chunk
        Eigen::MatrixXd train_genotype = bed_->load_chunk(start, end);
        Eigen::MatrixXd test_genotype = test_bed_pipe.load_chunk(start, end);

        // Apply coding
        Eigen::VectorXd p_major_chunk = p_major_.segment(start, end - start);
        code_method_varden(p_major_chunk, train_genotype, add);
        code_method_varden(p_major_chunk, test_genotype, add);

        grm.noalias() += test_genotype * train_genotype.transpose();
    }

    return grm / scale_factor_;
}

void CrossGRM::check_snp_consistency(const BedPipe& test_bed) const {}

Eigen::VectorXd compute_p_major(
    const Eigen::Ref<const Eigen::MatrixXd>& genotype)
{
    Eigen::VectorXd p_major(genotype.cols());

#pragma omp parallel for default(none) shared(genotype, p_major)
    for (Eigen::Index i = 0; i < genotype.cols(); ++i)
    {
        p_major(i) = genotype.col(i).mean() / 2;
    }
    return p_major;
}

void code_method_varden(
    const Eigen::Ref<const Eigen::VectorXd>& p_major,
    Eigen::Ref<Eigen::MatrixXd> genotype,
    bool add)
{
    if (p_major.size() != genotype.cols())
    {
        throw std::runtime_error(
            "allele freq does not match genotype snp number");
    }

    if (add)
    {
        // For additive coding: subtract 2*p_major
#pragma omp parallel for default(none) shared(genotype, p_major)
        for (Eigen::Index i = 0; i < genotype.cols(); ++i)
        {
            genotype.col(i).array() -= 2 * p_major(i);
        }
    }
    else
    {
        // For dominant coding: use different formula
        // p_major = 2 * p_major % (1 - p_major) becomes:
        Eigen::VectorXd adjusted_p
            = 2 * p_major.array() * (1 - p_major.array());

        // Replace 2 with 0 (assuming 2 represents homozygous alternate)
        genotype = (genotype.array() == 2).select(0, genotype);

#pragma omp parallel for default(none) shared(genotype, adjusted_p)
        for (Eigen::Index i = 0; i < genotype.cols(); ++i)
        {
            genotype.col(i).array() -= adjusted_p(i);
        }
    }
}

}  // namespace gelex
