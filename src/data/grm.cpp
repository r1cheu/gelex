#include "gelex/data/grm.h"

#include <expected>
#include <memory>

#include <Eigen/Core>

#include "gelex/data/bed_pipe.h"

namespace gelex
{
using Eigen::Index;

GRM::GRM(
    std::string_view bed_file,
    Index chunk_size,
    const std::unordered_map<std::string, Eigen::Index>& target_order)
{
    auto bed_pipe = BedPipe::create(bed_file);
    if (!bed_pipe)
    {
        throw std::runtime_error(
            "Failed to create BedPipe: "
            + std::string(bed_pipe.error().message));
    }
    bed_ = std::make_unique<BedPipe>(std::move(*bed_pipe));

    // Use the provided id_map directly
    if (!target_order.empty())
    {
        id_map_ = target_order;

        // Validate that all IDs in id_map exist in the bed file
        for (const auto& [id, index] : id_map_)
        {
            if (!bed_->sample_map().contains(id))
            {
                throw std::runtime_error(
                    "Sample ID '" + id + "' not found in BED file");
            }
        }
    }
    chunk_size_ = chunk_size;

    p_major_ = Eigen::VectorXd::Zero(bed_->num_variants());
}

Eigen::MatrixXd GRM::compute(bool add)
{
    const auto n = static_cast<Eigen::Index>(
        id_map_.empty() ? bed_->sample_map().size() : id_map_.size());
    Eigen::MatrixXd grm = Eigen::MatrixXd::Zero(n, n);

    // Process in chunks
    const Index num_variants = bed_->num_variants();
    const Index chunk_size = 10000;  // Fixed chunk size for BedPipe

    for (Index start = 0; start < num_variants; start += chunk_size)
    {
        Index end = std::min(start + chunk_size, num_variants);
        auto genotype_result = bed_->load_chunk(start, end);
        if (!genotype_result)
        {
            throw std::runtime_error(
                "Failed to load genotype chunk: "
                + std::string(genotype_result.error().message));
        }

        Eigen::MatrixXd genotype = std::move(*genotype_result);
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
    Index chunk_size,
    const std::unordered_map<std::string, Eigen::Index>& target_order)
    : scale_factor_{scale_factor}, chunk_size_{chunk_size}
{
    auto bed_pipe = BedPipe::create(train_bed);
    if (!bed_pipe)
    {
        throw std::runtime_error(
            "Failed to create BedPipe: "
            + std::string(bed_pipe.error().message));
    }
    bed_ = std::make_unique<BedPipe>(std::move(*bed_pipe));

    // Use the provided id_map directly
    if (!target_order.empty())
    {
        id_map_ = target_order;

        // Validate that all IDs in id_map exist in the bed file
        for (const auto& [id, index] : id_map_)
        {
            if (!bed_->sample_map().contains(id))
            {
                throw std::runtime_error(
                    "Sample ID '" + id + "' not found in BED file");
            }
        }
    }

    if (p_major.size() != bed_->num_variants())
    {
        throw std::runtime_error(
            "p_major size does not match number of SNPs in training set.");
    }
    p_major_ = std::move(p_major);
}

Eigen::MatrixXd CrossGRM::compute(std::string_view test_bed, bool add)
{
    // Create test BedPipe
    auto test_bed_pipe = BedPipe::create(test_bed);
    if (!test_bed_pipe)
    {
        throw std::runtime_error(
            "Failed to create test BedPipe: "
            + std::string(test_bed_pipe.error().message));
    }

    check_snp_consistency(*test_bed_pipe);

    const Index test_n = test_bed_pipe->sample_size();
    const auto train_n = static_cast<Index>(
        id_map_.empty() ? bed_->sample_map().size() : id_map_.size());
    Eigen::MatrixXd grm = Eigen::MatrixXd::Zero(test_n, train_n);

    // Process in chunks
    const Index num_variants = bed_->num_variants();
    const Index chunk_size = 10000;

    for (Index start = 0; start < num_variants; start += chunk_size)
    {
        Index end = std::min(start + chunk_size, num_variants);

        // Load training chunk
        auto train_genotype_result = bed_->load_chunk(start, end);
        if (!train_genotype_result)
        {
            throw std::runtime_error(
                "Failed to load training genotype chunk: "
                + std::string(train_genotype_result.error().message));
        }
        Eigen::MatrixXd train_genotype = std::move(*train_genotype_result);

        // Load test chunk
        auto test_genotype_result = test_bed_pipe->load_chunk(start, end);
        if (!test_genotype_result)
        {
            throw std::runtime_error(
                "Failed to load test genotype chunk: "
                + std::string(test_genotype_result.error().message));
        }
        Eigen::MatrixXd test_genotype = std::move(*test_genotype_result);

        // Apply coding
        Eigen::VectorXd p_major_chunk = p_major_.segment(start, end - start);
        code_method_varden(p_major_chunk, train_genotype, add);
        code_method_varden(p_major_chunk, test_genotype, add);

        grm.noalias() += test_genotype * train_genotype.transpose();
    }

    return grm / scale_factor_;
}

void CrossGRM::check_snp_consistency(const BedPipe& test_bed) const
{
    if (bed_->num_variants() != test_bed.num_variants())
    {
        throw std::runtime_error{
            "Number of SNPs in training and test sets do not match."};
    }

    const auto& train_snps = bed_->snp_ids();
    const auto& test_snps = test_bed.snp_ids();

    for (Index i{0}; i < bed_->num_variants(); ++i)
    {
        if (train_snps[i] != test_snps[i])
        {
            throw std::runtime_error{
                "SNPs in training and test sets do not match."};
        }
    }
}

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
