#pragma once

#include <armadillo>
#include <cstdint>
#include "chenx/data/bed_reader.h"

namespace chenx
{
using arma::dmat;
using arma::rowvec;
class Grm
{
   public:
    Grm(const std::string& bed_file, uint64_t chunk_size);
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

   private:
    rowvec center_;
    void Centerlize(dmat& genotype) override;
};

class DomGrm : public Grm
{
   public:
    DomGrm(const std::string& bed_file, uint64_t chunk_size);

   private:
    rowvec center_;
    void Centerlize(dmat& genotype) override;
};

}  // namespace chenx
