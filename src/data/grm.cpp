#include "gelex/data/grm.h"

#include <Eigen/Core>

namespace gelex
{
using Eigen::Index;

GRM::GRM(
    std::string_view bed_file,
    Index chunk_size,
    const std::vector<std::string>& target_order)
    : bed_(bed_file, chunk_size, target_order),
      p_major_(Eigen::VectorXd::Zero(bed_.num_snps())) {};

Eigen::MatrixXd GRM::compute(bool add)
{
    const Index n{bed_.num_individuals()};
    Eigen::MatrixXd grm = Eigen::MatrixXd::Zero(n, n);
    while (bed_.has_next())
    {
        Eigen::MatrixXd genotype = bed_.read_chunk();
        Eigen::VectorXd p_major_i = compute_p_major(genotype);
        code_method_varden(p_major_i, genotype, add);

        // Update p_major_ vector
        Index start_idx = bed_.current_chunk_index() - p_major_i.size();
        for (Index i = 0; i < p_major_i.size(); ++i)
        {
            p_major_(start_idx + i) = p_major_i(i);
        }

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
    const std::vector<std::string>& target_order)
    : bed_(train_bed, chunk_size, target_order),
      p_major_(std::move(p_major)),
      scale_factor_{scale_factor},
      chunk_size_{chunk_size}
{
    if (p_major_.size() != bed_.num_snps())
    {
        throw std::runtime_error(
            "p_major size does not match number of SNPs in training set.");
    }
}

Eigen::MatrixXd CrossGRM::compute(std::string_view test_bed, bool add)
{
    BedReader test_bed_reader(test_bed, chunk_size_);
    individuals_ = test_bed_reader.individuals();
    check_snp_consistency(test_bed_reader);

    Eigen::MatrixXd grm = Eigen::MatrixXd::Zero(
        test_bed_reader.num_individuals(), bed_.num_individuals());
    while (bed_.has_next())
    {
        auto start = bed_.current_chunk_index();
        Eigen::MatrixXd train_genotype = bed_.read_chunk();
        Eigen::MatrixXd test_genotype = test_bed_reader.read_chunk();
        auto n = train_genotype.cols();
        code_method_varden(p_major_.segment(start, n), train_genotype, add);
        code_method_varden(p_major_.segment(start, n), test_genotype, add);
        grm.noalias() += test_genotype * train_genotype.transpose();
    }
    return grm / scale_factor_;
}

void CrossGRM::check_snp_consistency(const BedReader& test_bed) const
{
    if (bed_.num_snps() != test_bed.num_snps())
    {
        throw std::runtime_error{
            "Number of SNPs in training and test sets do not match."};
    }

    for (Index i{0}; i < bed_.num_snps(); ++i)
    {
        if (bed_.snps()[i] != test_bed.snps()[i])
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
