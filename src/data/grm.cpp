#include "gelex/data/grm.h"

#include <armadillo>

namespace gelex
{

GRM::GRM(
    std::string_view bed_file,
    size_t chunk_size,
    const std::vector<std::string>& target_order)
    : bed_(bed_file, chunk_size, target_order), p_major_(bed_.num_snps()) {};

arma::dmat GRM::compute(bool add)
{
    const size_t n{bed_.num_individuals()};
    arma::dmat grm(n, n, arma::fill::zeros);
    while (bed_.has_next())
    {
        arma::dmat genotype{bed_.read_chunk()};
        arma::dvec p_major_i{compute_p_major(genotype)};
        code_method_varden(p_major_i, genotype, add);
        p_major_.subvec(
            bed_.current_chunk_index() - p_major_i.n_elem,
            arma::size(p_major_i))
            = p_major_i;
        grm += genotype * genotype.t();
    }
    scale_factor_ = arma::trace(grm) / static_cast<double>(grm.n_rows);
    return grm / scale_factor_;
}

CrossGRM::CrossGRM(
    std::string_view train_bed,
    arma::dvec p_major,
    double scale_factor,
    size_t chunk_size,
    const std::vector<std::string>& target_order)
    : bed_(train_bed, chunk_size, target_order),
      p_major_(std::move(p_major)),
      scale_factor_{scale_factor},
      chunk_size_{chunk_size}
{
    if (p_major_.n_elem != bed_.num_snps())
    {
        throw std::runtime_error(
            "p_major size does not match number of SNPs in training set.");
    }
}

arma::dmat CrossGRM::compute(std::string_view test_bed, bool add)
{
    BedReader test_bed_reader(test_bed, chunk_size_);
    individuals_ = test_bed_reader.individuals();
    check_snp_consistency(test_bed_reader);

    arma::dmat grm(
        test_bed_reader.num_individuals(),
        bed_.num_individuals(),
        arma::fill::zeros);
    while (bed_.has_next())
    {
        auto start = bed_.current_chunk_index();
        arma::dmat train_genotype{bed_.read_chunk()};
        arma::dmat test_genotype{test_bed_reader.read_chunk()};
        auto end = start + train_genotype.n_cols - 1;
        code_method_varden(p_major_.subvec(start, end), train_genotype, add);
        code_method_varden(p_major_.subvec(start, end), test_genotype, add);
        grm += test_genotype * train_genotype.t();
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

    for (size_t i{0}; i < bed_.num_snps(); ++i)
    {
        if (bed_.snps()[i] != test_bed.snps()[i])
        {
            throw std::runtime_error{
                "SNPs in training and test sets do not match."};
        }
    }
}

arma::dvec compute_p_major(const arma::dmat& genotype)
{
    arma::dvec p_major(genotype.n_cols);

#pragma omp parallel for default(none) shared(genotype, p_major)
    for (size_t i = 0; i < genotype.n_cols; ++i)
    {
        p_major.at(i) = arma::mean(genotype.unsafe_col(i)) / 2;
    }
    return p_major;
}

void code_method_varden(arma::dvec p_major, arma::dmat& genotype, bool add)
{
    if (p_major.n_elem != genotype.n_cols)
    {
        throw std::runtime_error(
            "allele freq does not match genotype snp number");
    }

    if (add)
    {
        p_major *= 2;

#pragma omp parallel for default(none) shared(genotype, p_major)
        for (size_t i = 0; i < genotype.n_cols; ++i)
        {
            genotype.col(i) -= p_major(i);
        }
    }

    else
    {
        p_major = 2 * p_major % (1 - p_major);
        genotype.replace(2, 0);
#pragma omp parallel for default(none) shared(genotype, p_major)
        for (size_t i = 0; i < genotype.n_cols; ++i)
        {
            genotype.col(i) -= p_major(i);
        }
    }
}

}  // namespace gelex
