#include "chenx/data/grm.h"

#include <armadillo>
#include <cmath>

#include "chenx/data/bed_reader.h"

namespace chenx
{

CrossGrm::CrossGrm(
    const std::string& train_bed,
    const std::string& test_bed,
    rowvec&& center,
    double scale_factor,
    uint64_t chunk_size)
    : train_bed_{train_bed, chunk_size},
      test_bed_{test_bed, chunk_size},
      center_{std::move(center)},
      scale_factor_{scale_factor} {};

dmat CrossGrm::Compute(bool dom)
{
    CheckSnpConsistency();
    dmat grm{
        test_bed_.n_individuals(),
        train_bed_.n_individuals(),
        arma::fill::zeros};
    while (train_bed_.HasNext())
    {
        dmat train_genotype{train_bed_.ReadChunk()};
        dmat test_genotype{test_bed_.ReadChunk()};
        if (dom)
        {
            train_genotype.replace(2.0, 0.0);
            test_genotype.replace(2.0, 0.0);
        }
        train_genotype.each_row() -= center_;
        test_genotype.each_row() -= center_;
        grm += test_genotype * train_genotype.t();
    }
    grm /= scale_factor_;
    return grm;
}

void CrossGrm::CheckSnpConsistency()
{
    for (uint64_t i{0}; i < train_bed_.n_snps(); ++i)
    {
        if (train_bed_.snps()[i] != test_bed_.snps()[i])
        {
            throw std::runtime_error{
                "SNPs in training and test sets do not match."};
        }
    }
}

Grm::Grm(const std::string& bed_file, uint64_t chunk_size)
    : bed_{bed_file, chunk_size}
{
}

double Grm::scale_factor() const noexcept
{
    return scale_factor_;
}

const BedReader& Grm::bed() const noexcept
{
    return bed_;
}

dmat Grm::Compute()
{
    const uint64_t n{bed_.n_individuals()};
    dmat grm{n, n, arma::fill::zeros};

    while (bed_.HasNext())
    {
        dmat genotype{bed_.ReadChunk()};
        Centerlize(genotype);
        grm += genotype * genotype.t();
    };
    scale_factor_ = Scale(grm);
    return grm;
}

double Grm::Scale(dmat& grm)
{
    double scale_factor = arma::trace(grm) / static_cast<double>(grm.n_rows);
    grm /= scale_factor;
    return scale_factor;
}

AddGrm::AddGrm(const std::string& bed_file, uint64_t chunk_size)
    : Grm{bed_file, chunk_size}
{
    center_.zeros(bed().n_snps());
}

const arma::rowvec& AddGrm::center() const
{
    return center_;
}

void AddGrm::Centerlize(dmat& genotype)
{
    rowvec center{arma::mean(genotype, 0)};
    genotype.each_row() -= center;
    auto start_index = bed().current_chunk_index() - center.n_cols;
    center_.subvec(start_index, arma::size(center)) = center;
}

DomGrm::DomGrm(const std::string& bed_file, uint64_t chunk_size)
    : Grm{bed_file, chunk_size}
{
    center_.zeros(bed().n_snps());
}

const arma::rowvec& DomGrm::center() const
{
    return center_;
}

void DomGrm::Centerlize(dmat& genotype)
{
    rowvec pA = arma::mean(genotype, 0) / 2;
    rowvec center = 2 * (pA % (1 - pA));
    auto start_index = bed().current_chunk_index() - center.n_cols;
    genotype.each_row() -= center;
    center_.subvec(start_index, arma::size(center)) = center;
}

}  // namespace chenx
