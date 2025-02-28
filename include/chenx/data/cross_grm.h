#pragma once

#include <armadillo>
#include <cstdint>
#include <cstdio>
#include <string_view>

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
        std::vector<std::string>&& exclude_individuals = {},
        uint64_t chunk_size = DEFAULT_CHUNK_SIZE);

    CrossGrm(CrossGrm&&) noexcept = default;
    CrossGrm(const CrossGrm&) = delete;
    CrossGrm& operator=(CrossGrm&&) noexcept = default;
    CrossGrm& operator=(const CrossGrm&) = delete;

    ~CrossGrm() override = default;
    virtual dmat Compute(std::string_view test_bed_path);

   protected:
    void CheckSnpConsistency(const BedReader& test_bed) const;
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

class CrossOnceGrm : public CrossGrm
{
   public:
    CrossOnceGrm(
        std::string_view train_bed_file,
        rowvec&& center,
        double scale_factor,
        std::vector<std::string>&& exclude_individuals = {});

    CrossOnceGrm(CrossOnceGrm&&) noexcept = default;
    CrossOnceGrm(const CrossOnceGrm&) = delete;
    CrossOnceGrm& operator=(CrossOnceGrm&&) noexcept = default;
    CrossOnceGrm& operator=(const CrossOnceGrm&) = delete;
    ~CrossOnceGrm() override = default;

    dmat Compute(std::string_view test_bed_path) override;

   private:
    dmat train_genotype_;
};

class AddCrossOnceGrm : public CrossOnceGrm
{
    using CrossOnceGrm::CrossOnceGrm;

   protected:
    void Encode(dmat& genotype) override;
};

class DomCrossOnceGrm : public CrossOnceGrm
{
    using CrossOnceGrm::CrossOnceGrm;

   protected:
    void Encode(dmat& genotype) override;
};

}  // namespace chenx
