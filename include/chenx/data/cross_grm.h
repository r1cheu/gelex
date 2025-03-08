#pragma once

#include <cstdint>
#include <cstdio>
#include <string_view>

#include <armadillo>
#include <vector>

#include "chenx/data/bed_reader.h"
#include "chenx/data/grm.h"

namespace chenx
{
using arma::dmat;
using arma::rowvec;
class CrossGrm : public IGrm
{
   public:
    CrossGrm(
        std::string_view train_bed_file,
        rowvec&& center,
        double scale_factor,
        uint64_t chunk_size = DEFAULT_CHUNK_SIZE,
        const std::vector<std::string>& exclude_individuals = {});

    CrossGrm(CrossGrm&&) noexcept = default;
    CrossGrm(const CrossGrm&) = delete;
    CrossGrm& operator=(CrossGrm&&) noexcept = default;
    CrossGrm& operator=(const CrossGrm&) = delete;

    ~CrossGrm() override = default;
    dmat Compute(std::string_view test_bed_path);

    const std::vector<std::string>& test_individuals() const noexcept
    {
        return test_individuals_;
    }

   protected:
    void CheckSnpConsistency(const BedReader& test_bed) const;

   private:
    void Reset();
    std::vector<std::string> test_individuals_;
};

class AddCrossGrm : public CrossGrm
{
    using CrossGrm::CrossGrm;

   protected:
    void Encode(dmat& genotype) override;
};

class DomCrossGrm : public CrossGrm
{
    using CrossGrm::CrossGrm;

   protected:
    void Encode(dmat& genotype) override;
};
}  // namespace chenx
