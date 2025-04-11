#include "gelex/data/grm.h"

#include <armadillo>
#include <cmath>
#include <vector>

#include "gelex/data/bed_reader.h"

namespace gelex
{

void dom_encode(dmat& genotype)
{
    genotype.replace(2, 0);
}

IGrm::IGrm(
    std::string_view bed_file,
    uint64_t chunk_size,
    const std::vector<std::string>& exclude_individuals)
    : bed_{bed_file, chunk_size, exclude_individuals} {};

Grm::Grm(
    std::string_view bed_file,
    uint64_t chunk_size,
    const std::vector<std::string>& exclude_individuals)
    : IGrm{bed_file, chunk_size, exclude_individuals}
{
    set_center(rowvec{bed().num_snps(), arma::fill::zeros});
}

void Grm::centerlize(dmat& genotype)
{
    rowvec center = compute_center(genotype);
    auto start_index = bed().current_chunk_index() - center.n_cols;
    encode(genotype);
    genotype.each_row() -= center;
    set_center(start_index, center);
}

dmat Grm::compute()
{
    reset();
    const uint64_t num_ind{bed().num_individuals()};
    dmat grm{num_ind, num_ind, arma::fill::zeros};

    while (bed().has_next())
    {
        dmat genotype{bed().read_chunk()};
        centerlize(genotype);
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

void AddGrm::encode(dmat& genotype) {
};  // add is the default encoding, so do nothing

rowvec AddGrm::compute_center(const dmat& genotype)
{
    return arma::mean(genotype, 0);
}

void DomGrm::encode(dmat& genotype)
{
    dom_encode(genotype);
}

rowvec DomGrm::compute_center(const dmat& genotype)
{
    rowvec pA = arma::mean(genotype, 0) / 2;
    rowvec center = 2 * (pA % (1 - pA));
    return center;
}

}  // namespace gelex
