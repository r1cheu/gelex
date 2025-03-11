#include "gelex/data/cross_grm.h"

#include <armadillo>
#include <string_view>

#include "gelex/data/bed_reader.h"
#include "gelex/data/grm.h"

namespace gelex
{

CrossGrm::CrossGrm(
    std::string_view train_bed_file,
    rowvec&& center,
    double scale_factor,
    uint64_t chunk_size,
    const std::vector<std::string>& exclude_individuals)
    : IGrm{train_bed_file, chunk_size, exclude_individuals}
{
    set_center(std::move(center));
    set_scale_factor(scale_factor);
}

void CrossGrm::CheckSnpConsistency(const BedReader& test_bed) const
{
    for (uint64_t i{0}; i < bed().num_snps(); ++i)
    {
        if (bed().snps()[i] != test_bed.snps()[i])
        {
            throw std::runtime_error{
                "SNPs in training and test sets do not match."};
        }
    }
}

void CrossGrm::Reset()
{
    bed().Reset();
}

dmat CrossGrm::Compute(std::string_view test_bed_path)
{
    Reset();
    BedReader test_bed{test_bed_path, bed().chunk_size()};
    test_individuals_ = test_bed.individuals();
    CheckSnpConsistency(test_bed);
    dmat grm{
        test_bed.num_individuals(), bed().num_individuals(), arma::fill::zeros};
    while (bed().HasNext())
    {
        auto start = bed().current_chunk_index();

        dmat train_genotype{bed().ReadChunk()};
        dmat test_genotype{test_bed.ReadChunk()};
        auto end = start + train_genotype.n_cols - 1;

        Encode(train_genotype);
        Encode(test_genotype);

        train_genotype.each_row() -= center().subvec(start, end);
        test_genotype.each_row() -= center().subvec(start, end);
        grm += test_genotype * train_genotype.t();
    }
    grm /= scale_factor();
    return grm;
}

void AddCrossGrm::Encode(dmat& genotype) {}
void DomCrossGrm::Encode(dmat& genotype)
{
    dom_encode(genotype);
}
}  // namespace gelex
