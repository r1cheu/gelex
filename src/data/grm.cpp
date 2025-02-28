#include "chenx/data/grm.h"

#include <armadillo>
#include <cmath>
#include <vector>

#include "chenx/data/bed_reader.h"

namespace chenx
{

void dom_encode(dmat& genotype)
{
    genotype.replace(2, 0);
}

IGrm::IGrm(
    std::string_view bed_file,
    std::vector<std::string>&& exclude_individuals,
    uint64_t chunk_size)
    : bed_{bed_file, std::move(exclude_individuals), chunk_size} {};

Grm::Grm(
    std::string_view bed_file,
    std::vector<std::string>&& exclude_individuals,
    uint64_t chunk_size)
    : IGrm{bed_file, std::move(exclude_individuals), chunk_size}
{
    set_center(rowvec{bed().num_snps(), arma::fill::zeros});
}

void Grm::Centerlize(dmat& genotype)
{
    rowvec center = ComputeCenter(genotype);
    auto start_index = bed().current_chunk_index() - center.n_cols;
    Encode(genotype);
    genotype.each_row() -= center;
    set_center(start_index, center);
}

dmat Grm::Compute()
{
    const uint64_t num_ind{bed().num_individuals()};
    dmat grm{num_ind, num_ind, arma::fill::zeros};

    while (bed().HasNext())
    {
        dmat genotype{bed().ReadChunk()};
        Centerlize(genotype);
        grm += genotype * genotype.t();
    };
    set_scale_factor(Scale(grm));
    return grm;
}

double Grm::Scale(dmat& grm)
{
    double scale_factor = arma::trace(grm) / static_cast<double>(grm.n_rows);
    grm /= scale_factor;
    return scale_factor;
}

void AddGrm::Encode(dmat& genotype) {
};  // add is the default encoding, so do nothing

rowvec AddGrm::ComputeCenter(const dmat& genotype)
{
    return arma::mean(genotype, 0);
}

void DomGrm::Encode(dmat& genotype)
{
    dom_encode(genotype);
}

rowvec DomGrm::ComputeCenter(const dmat& genotype)
{
    rowvec pA = arma::mean(genotype, 0) / 2;
    rowvec center = 2 * (pA % (1 - pA));
    return center;
}

}  // namespace chenx
