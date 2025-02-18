#pragma once

#include <armadillo>
#include <cstdint>
#include <string_view>

#include "chenx/data/bed_reader.h"

namespace chenx
{
using arma::dmat;
using arma::rowvec;
class BaseCrossGrm
{
   public:
    BaseCrossGrm(
        BedReader&& train_bed,
        rowvec&& center,
        double scale_factor) noexcept
        : center_{std::move(center)},
          scale_factor_{scale_factor},
          train_bed_{std::move(train_bed)}
    {
    }
    BaseCrossGrm(BaseCrossGrm&&) noexcept = default;
    BaseCrossGrm(const BaseCrossGrm&) = delete;
    BaseCrossGrm& operator=(BaseCrossGrm&&) noexcept = default;
    BaseCrossGrm& operator=(const BaseCrossGrm&) = delete;

    const rowvec& center() const noexcept { return center_; }
    double scale_factor() const noexcept { return scale_factor_; }
    const BedReader& train_bed() const noexcept { return train_bed_; }
    BedReader& train_bed() noexcept { return train_bed_; }

    virtual ~BaseCrossGrm() = default;
    virtual dmat Compute(std::string_view test_bed_path) = 0;

    void CheckSnpConsistency(const BedReader& test_bed);

   private:
    BedReader train_bed_;
    rowvec center_;
    double scale_factor_{};
};

class CrossChunkGrm : public BaseCrossGrm
{
   public:
    CrossChunkGrm(
        std::string_view train_bed_path,
        rowvec&& center,
        double scale_factor,
        uint64_t chunk_size = 10000);

    dmat Compute(std::string_view test_bed_path) override;

   private:
    uint64_t chunk_size_;
    virtual void Encode(dmat& genotype) = 0;
};

class AddCrossChunkGrm : public CrossChunkGrm
{
    using CrossChunkGrm::CrossChunkGrm;

   private:
    void Encode(dmat& genotype) override;
};

class DomCrossChunkGrm : public CrossChunkGrm
{
    using CrossChunkGrm::CrossChunkGrm;

   private:
    void Encode(dmat& genotype) override;
};

class CrossGrm : public BaseCrossGrm
{
   public:
    CrossGrm(
        std::string_view train_bed_path,
        rowvec&& center,
        double scale_factor);
    dmat Compute(std::string_view test_bed_path) override;

   private:
    dmat train_genotype;
    virtual void Encode(dmat& genotype) = 0;
};

class AddCrossGrm : public CrossGrm
{
    using CrossGrm::CrossGrm;

   private:
    void Encode(dmat& genotype) override;
};

class DomCrossGrm : public CrossGrm
{
    using CrossGrm::CrossGrm;

   private:
    void Encode(dmat& genotype) override;
};
}  // namespace chenx
