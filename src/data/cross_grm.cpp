#include "gelex/data/cross_grm.h"

#include <armadillo>
#include <string_view>

#include "gelex/data/bed_reader.h"
#include "gelex/data/grm.h"

namespace gelex
{

CrossGrm::CrossGrm(
    std::string_view train_bed_file,
    rowvec center,
    double scale_factor,
    size_t chunk_size,
    const std::vector<std::string>& exclude_individuals)
    : IGrm{train_bed_file, chunk_size, exclude_individuals}
{
    set_center(std::move(center));
    set_scale_factor(scale_factor);
}

void CrossGrm::check_snp_consistency(const BedReader& test_bed) const
{
    for (size_t i{0}; i < bed().num_snps(); ++i)
    {
        if (bed().snps()[i] != test_bed.snps()[i])
        {
            throw std::runtime_error{
                "SNPs in training and test sets do not match."};
        }
    }
}

void CrossGrm::reset()
{
    bed().reset();
}

dmat CrossGrm::compute(std::string_view test_bed_path)
{
    reset();
    BedReader test_bed{test_bed_path, bed().chunk_size()};
    test_individuals_ = test_bed.individuals();
    check_snp_consistency(test_bed);
    dmat grm{
        test_bed.num_individuals(), bed().num_individuals(), arma::fill::zeros};
    while (bed().has_next())
    {
        auto start = bed().current_chunk_index();

        dmat train_genotype{bed().read_chunk()};
        dmat test_genotype{test_bed.read_chunk()};
        auto end = start + train_genotype.n_cols - 1;

        encode(train_genotype);
        encode(test_genotype);

        train_genotype.each_row() -= center().subvec(start, end);
        test_genotype.each_row() -= center().subvec(start, end);
        grm += test_genotype * train_genotype.t();
    }
    grm /= scale_factor();
    return grm;
}

void AddCrossGrm::encode(dmat& genotype) {}
void DomCrossGrm::encode(dmat& genotype)
{
    dom_encode(genotype);
}
}  // namespace gelex
