#include "chenx/data/cross_grm.h"

#include <armadillo>
#include <string_view>

#include "chenx/data/bed_reader.h"
#include "chenx/data/grm.h"

namespace chenx
{
void BaseCrossGrm::CheckSnpConsistency(const BedReader& test_bed)
{
    for (uint64_t i{0}; i < train_bed_.num_snps(); ++i)
    {
        if (train_bed_.snps()[i] != test_bed.snps()[i])
        {
            throw std::runtime_error{
                "SNPs in training and test sets do not match."};
        }
    }
}

CrossChunkGrm::CrossChunkGrm(
    std::string_view train_bed_path,
    rowvec&& center,
    double scale_factor,
    uint64_t chunk_size)
    : chunk_size_{chunk_size},
      BaseCrossGrm{
          BedReader{train_bed_path, chunk_size},
          std::move(center),
          scale_factor} {};

dmat CrossChunkGrm::Compute(std::string_view test_bed_path)
{
    BedReader test_bed{test_bed_path, chunk_size_};
    CheckSnpConsistency(test_bed);
    dmat grm{
        test_bed.num_individuals(),
        train_bed().num_individuals(),
        arma::fill::zeros};
    while (train_bed().HasNext())
    {
        dmat train_genotype{train_bed().ReadChunk()};
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

void AddCrossChunkGrm::Encode(dmat& genotype) {}
void DomCrossChunkGrm::Encode(dmat& genotype)
{
    dom_encode(genotype);
}

CrossGrm::CrossGrm(
    std::string_view train_bed_path,
    rowvec&& center,
    double scale_factor)
    : BaseCrossGrm{
          BedReader{train_bed_path, std::numeric_limits<uint64_t>::max()},
          std::move(center),
          scale_factor}
{
    train_genotype = train_bed().ReadChunk();
    train_genotype.each_row() -= BaseCrossGrm::center();
}

dmat CrossGrm::Compute(std::string_view test_bed_path)
{
    BedReader test_bed{test_bed_path, std::numeric_limits<uint64_t>::max()};
    CheckSnpConsistency(test_bed);
    dmat grm{
        test_bed.num_individuals(), train_genotype.n_cols, arma::fill::zeros};

    dmat test_genotype{test_bed.ReadChunk()};
    Encode(test_genotype);
    test_genotype.each_row() -= center();
    grm += test_genotype * train_genotype.t();
    grm /= scale_factor();
    return grm;
}

void AddCrossGrm::Encode(dmat& genotype) {}
void DomCrossGrm::Encode(dmat& genotype)
{
    dom_encode(genotype);
}
}  // namespace chenx
