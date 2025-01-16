#include "chenx/data/grm.h"

#include <armadillo>
#include <cmath>

#include "chenx/data/bed_reader.h"

namespace chenx
{
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
        dmat genotype{bed_.GetNextChunk()};
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
    center_.zeros(bed().n_individuals());
}

void AddGrm::Centerlize(dmat& genotype)
{
    rowvec center{arma::mean(genotype, 0)};
    genotype.each_row() -= center;
    center_.cols(
        bed().current_chunk_index(),
        bed().current_chunk_index() + genotype.n_cols - 1)
        = center;
}

DomGrm::DomGrm(const std::string& bed_file, uint64_t chunk_size)
    : Grm{bed_file, chunk_size}
{
    center_.zeros(bed().n_individuals());
}

void DomGrm::Centerlize(dmat& genotype)
{
    rowvec pA = arma::mean(genotype, 0) / 2;
    rowvec center = 2 * (pA % (1 - pA));
    genotype.each_row() -= center;
    center_.cols(
        bed().current_chunk_index(),
        bed().current_chunk_index() + genotype.n_cols - 1)
        = center;
}

}  // namespace chenx
