#pragma once

#include <armadillo>
#include <cstdint>
#include <string>
#include "chenx/data/bed_reader.h"

namespace chenx
{
using arma::dmat;
using arma::rowvec;
class CrossGrm
{
   public:
    CrossGrm(
        const std::string& train_bed,
        const std::string& test_bed,
        rowvec&& center,
        double scale_factor,
        uint64_t chunk_size = 10000);

    dmat Compute(bool dom);

   private:
    BedReader train_bed_;
    BedReader test_bed_;
    rowvec center_;
    double scale_factor_{};

    void CheckSnpConsistency();
};

class Grm
{
   public:
    Grm(const std::string& bed_file, uint64_t chunk_size = 10000);
    Grm& operator=(Grm&&) = delete;
    Grm(const Grm&) = delete;
    Grm(Grm&&) = delete;
    Grm& operator=(const Grm&) = delete;

    virtual ~Grm() = default;

    double scale_factor() const noexcept;
    const BedReader& bed() const noexcept;

    dmat Compute();

   private:
    BedReader bed_;
    double scale_factor_{};

    virtual void Centerlize(dmat& genotype) = 0;
    static double Scale(dmat& grm);
};

class AddGrm : public Grm
{
   public:
    AddGrm(const std::string& bed_file, uint64_t chunk_size);
    const rowvec& center() const;

   private:
    rowvec center_;
    void Centerlize(dmat& genotype) override;
};

class DomGrm : public Grm
{
   public:
    DomGrm(const std::string& bed_file, uint64_t chunk_size);
    const rowvec& center() const;

   private:
    rowvec center_;
    void Centerlize(dmat& genotype) override;
};

}  // namespace chenx
