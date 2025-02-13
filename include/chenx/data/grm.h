#pragma once

#include <armadillo>
#include <cstdint>
#include <string_view>
#include "chenx/data/bed_reader.h"

namespace chenx
{
using arma::dmat;
using arma::rowvec;
void dom_encode(dmat& genotype);

class Grm
{
   public:
    explicit Grm(std::string_view bed_file, uint64_t chunk_size = 10000);
    Grm& operator=(Grm&&) = delete;
    Grm(const Grm&) = delete;
    Grm(Grm&&) = delete;
    Grm& operator=(const Grm&) = delete;

    virtual ~Grm() = default;

    double scale_factor() const noexcept { return scale_factor_; }
    const BedReader& bed() const noexcept { return bed_; }
    const rowvec& center() const noexcept { return center_; }
    virtual dmat Compute();

   private:
    BedReader bed_;
    double scale_factor_{};
    rowvec center_;

    virtual rowvec ComputeCenter(const dmat& genotype) = 0;
    virtual void Encode(dmat& genotype) = 0;
    void Centerlize(dmat& genotype);
    static double Scale(dmat& grm);
};

class AddGrm : public Grm
{
    using Grm::Grm;

   private:
    void Encode(dmat& genotype) override;
    rowvec ComputeCenter(const dmat& genotype) override;
};

class DomGrm : public Grm
{
    using Grm::Grm;

   private:
    void Encode(dmat& genotype) override;
    rowvec ComputeCenter(const dmat& genotype) override;
};

}  // namespace chenx
