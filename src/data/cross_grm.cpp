#include "chenx/data/cross_grm.h"

#include <armadillo>
#include <string_view>

#include "chenx/data/bed_reader.h"
#include "chenx/data/grm.h"

namespace chenx
{

CrossGrm::CrossGrm(
    std::string_view train_bed_file,
    rowvec&& center,
    double scale_factor,
    std::vector<std::string>&& exclude_individuals,
    uint64_t chunk_size)
    : IGrm{train_bed_file, std::move(exclude_individuals), chunk_size}
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

dmat CrossGrm::Compute(std::string_view test_bed_path)
{
    BedReader test_bed{
        test_bed_path, bed().dropped_individuals(), bed().chunk_size()};
    CheckSnpConsistency(test_bed);
    dmat grm{
        test_bed.num_individuals(), bed().num_individuals(), arma::fill::zeros};
    while (bed().HasNext())
    {
        dmat train_genotype{bed().ReadChunk()};
        dmat test_genotype{test_bed.ReadChunk()};
        Encode(train_genotype);
        Encode(test_genotype);
        train_genotype.each_row() -= center();
        test_genotype.each_row() -= center();
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

CrossOnceGrm::CrossOnceGrm(
    std::string_view train_bed_file,
    rowvec&& center,
    double scale_factor,
    std::vector<std::string>&& exclude_individuals)
    : CrossGrm{
          train_bed_file,
          std::move(center),
          scale_factor,
          std::move(exclude_individuals),
          std::numeric_limits<uint64_t>::max()}

{
    train_genotype_ = bed().ReadChunk();
    train_genotype_.each_row() -= IGrm::center();
}

dmat CrossOnceGrm::Compute(std::string_view test_bed_path)
{
    BedReader test_bed{
        test_bed_path,
        bed().dropped_individuals(),
        std::numeric_limits<uint64_t>::max()};

    CheckSnpConsistency(test_bed);

    dmat grm{
        test_bed.num_individuals(), train_genotype_.n_cols, arma::fill::zeros};

    dmat test_genotype{test_bed.ReadChunk()};
    Encode(test_genotype);
    test_genotype.each_row() -= center();
    grm += test_genotype * train_genotype_.t();
    grm /= scale_factor();
    return grm;
}

void AddCrossOnceGrm::Encode(dmat& genotype) {}
void DomCrossOnceGrm::Encode(dmat& genotype)
{
    dom_encode(genotype);
}
}  // namespace chenx
