#include "bayes_effects.h"

#include <utility>

#include <Eigen/Core>

namespace gelex
{
using Eigen::Ref;

using Eigen::Index;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

namespace bayes
{

FixedEffect::FixedEffect(
    std::optional<std::vector<std::string>> levels,
    MatrixXd&& design_matrix)
    : design_matrix(std::move(design_matrix)), levels(std::move(levels))
{
    cols_norm = this->design_matrix.colwise().squaredNorm();
}

FixedState::FixedState(const FixedEffect& effect)
    : coeffs(VectorXd::Zero(effect.design_matrix.cols())) {};

RandomEffect::RandomEffect(
    std::optional<std::vector<std::string>> levels,
    MatrixXd&& design_matrix)
    : design_matrix(std::move(design_matrix)), levels(std::move(levels))
{
    cols_norm = this->design_matrix.colwise().squaredNorm();
}

RandomState::RandomState(const RandomEffect& effect)
    : coeffs(VectorXd::Zero(effect.design_matrix.cols())),
      variance{effect.init_variance}
{
}

AdditiveEffect::AdditiveEffect(GenotypeMap&& design_matrix)
    : design_matrix(std::move(design_matrix))
{
    cols_norm = get_matrix_ref(this->design_matrix).colwise().squaredNorm();
}

AdditiveEffect::AdditiveEffect(GenotypeMatrix&& design_matrix)
    : design_matrix(std::move(design_matrix))
{
    cols_norm = get_matrix_ref(this->design_matrix).colwise().squaredNorm();
}

AdditiveEffect::AdditiveEffect(GenotypeStorage&& design_matrix)
    : design_matrix(std::move(design_matrix))
{
    cols_norm = get_matrix_ref(this->design_matrix).colwise().squaredNorm();
}

AdditiveState::AdditiveState(const AdditiveEffect& effect)
    : coeffs(VectorXd::Zero(get_cols(effect.design_matrix))),
      u(VectorXd::Zero(get_rows(effect.design_matrix))),
      tracker(VectorXi::Zero(get_cols(effect.design_matrix))),
      pi{effect.pi, Eigen::VectorXi::Zero(effect.pi.size())},
      marker_variance(
          Eigen::VectorXd::Constant(
              effect.marker_variance_size,
              effect.init_marker_variance)) {};

DominantEffect::DominantEffect(GenotypeMap&& design_matrix)
    : design_matrix(std::move(design_matrix))
{
    cols_norm = get_matrix_ref(this->design_matrix).colwise().squaredNorm();
}

DominantEffect::DominantEffect(GenotypeMatrix&& design_matrix)
    : design_matrix(std::move(design_matrix))
{
    cols_norm = get_matrix_ref(this->design_matrix).colwise().squaredNorm();
}

DominantEffect::DominantEffect(GenotypeStorage&& design_matrix)
    : design_matrix(std::move(design_matrix))
{
    cols_norm = get_matrix_ref(this->design_matrix).colwise().squaredNorm();
}

DominantState::DominantState(const DominantEffect& effect)
    : coeffs(VectorXd::Zero(get_cols(effect.design_matrix))),
      ratios(VectorXd::Zero(get_cols(effect.design_matrix))),
      u(VectorXd::Zero(get_rows(effect.design_matrix))),
      ratio_mean(effect.ratio_mean),
      ratio_variance(effect.ratio_variance)
{
}

bool AdditiveEffect::is_monomorphic(Eigen::Index snp_index) const
{
    return is_monomorphic_variant(design_matrix, snp_index);
}

Index AdditiveEffect::num_mono() const
{
    return num_mono_variant(design_matrix);
}

bool DominantEffect::is_monomorphic(Eigen::Index snp_index) const
{
    return is_monomorphic_variant(design_matrix, snp_index);
}

Index DominantEffect::num_mono() const
{
    return num_mono_variant(design_matrix);
}
}  // namespace bayes
}  // namespace gelex
