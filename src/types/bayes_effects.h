#pragma once
#include <string>
#include <variant>
#include <vector>

#include <Eigen/Core>

#include "../src/model/bayes/distribution.h"
#include "gelex/data/genotype_matrix.h"
#include "gelex/data/genotype_mmap.h"

namespace gelex
{
namespace bayes
{

// Variant type for genotype storage
using GenotypeStorage = std::variant<GenotypeMap, GenotypeMatrix>;

// Helper to get matrix reference from either storage type
// Returns Eigen::Ref which can wrap both Map and MatrixXd
inline Eigen::Ref<const Eigen::MatrixXd> get_matrix_ref(
    const GenotypeStorage& storage)
{
    return std::visit(
        [](const auto& s) -> Eigen::Ref<const Eigen::MatrixXd>
        { return s.matrix(); },
        storage);
}

// Helper to get number of rows
inline Eigen::Index get_rows(const GenotypeStorage& storage)
{
    return std::visit([](const auto& s) { return s.rows(); }, storage);
}

// Helper to get number of columns
inline Eigen::Index get_cols(const GenotypeStorage& storage)
{
    return std::visit([](const auto& s) { return s.cols(); }, storage);
}

inline const Eigen::VectorXd& get_means(const GenotypeStorage& storage)
{
    return std::visit(
        [](const auto& s) -> const Eigen::VectorXd& { return s.mean(); },
        storage);
}

inline const Eigen::VectorXd& get_variances(const GenotypeStorage& storage)
{
    return std::visit(
        [](const auto& s) -> const Eigen::VectorXd& { return s.variance(); },
        storage);
}

// Helper to check monomorphic status
inline bool is_monomorphic_variant(
    const GenotypeStorage& storage,
    Eigen::Index idx)
{
    return std::visit(
        [idx](const auto& s) { return s.is_monomorphic(idx); }, storage);
}

// Helper to get number of monomorphic markers
inline Eigen::Index num_mono_variant(const GenotypeStorage& storage)
{
    return std::visit([](const auto& s) { return s.num_mono(); }, storage);
}

struct Pi
{
    Eigen::VectorXd prop;
    Eigen::VectorXi count;
};

struct FixedEffect
{
    FixedEffect(
        std::optional<std::vector<std::string>> levels,
        Eigen::MatrixXd&& design_matrix);

    Eigen::MatrixXd design_matrix;
    Eigen::VectorXd cols_norm;

    std::optional<std::vector<std::string>> levels;
};

struct FixedState
{
    explicit FixedState(const FixedEffect& effect);
    Eigen::VectorXd coeffs;
};

struct RandomEffect
{
    RandomEffect(
        std::optional<std::vector<std::string>> levels,
        Eigen::MatrixXd&& design_matrix);

    Eigen::MatrixXd design_matrix;
    Eigen::VectorXd cols_norm;

    std::optional<std::vector<std::string>> levels;

    detail::ScaledInvChiSqParams prior{4, 0};
    double init_variance{0.0};
};

struct RandomState
{
    explicit RandomState(const RandomEffect& effect);

    Eigen::VectorXd coeffs;
    double variance{0.0};
};

struct GeneticEffect
{
    explicit GeneticEffect(GenotypeMap&& design_matrix)
        : design_matrix(std::move(design_matrix))
    {
        cols_norm = get_matrix_ref(this->design_matrix).colwise().squaredNorm();
    }

    explicit GeneticEffect(GenotypeMatrix&& design_matrix)
        : design_matrix(std::move(design_matrix))
    {
        cols_norm = get_matrix_ref(this->design_matrix).colwise().squaredNorm();
    }

    explicit GeneticEffect(GenotypeStorage&& design_matrix)
        : design_matrix(std::move(design_matrix))
    {
        cols_norm = get_matrix_ref(this->design_matrix).colwise().squaredNorm();
    }

    GenotypeStorage design_matrix;
    Eigen::VectorXd cols_norm;

    detail::ScaledInvChiSqParams marker_variance_prior{4, 0};
    double init_marker_variance{0.0};
    Eigen::Index marker_variance_size{0};

    std::optional<Eigen::VectorXd> init_pi;
    std::optional<Eigen::VectorXd> scale;
    bool estimate_pi{false};

    bool is_monomorphic(Eigen::Index snp_index) const
    {
        return is_monomorphic_variant(design_matrix, snp_index);
    }

    Eigen::Index num_mono() const { return num_mono_variant(design_matrix); }
};

struct AdditiveEffect : GeneticEffect
{
    using GeneticEffect::GeneticEffect;
};

struct GeneticState
{
    explicit GeneticState(const GeneticEffect& effect)
        : coeffs(Eigen::VectorXd::Zero(bayes::get_cols(effect.design_matrix))),
          u(Eigen::VectorXd::Zero(bayes::get_rows(effect.design_matrix))),
          marker_variance(
              Eigen::VectorXd::Constant(
                  effect.marker_variance_size,
                  effect.init_marker_variance))

    {
        if (effect.init_pi)
        {
            tracker
                = Eigen::VectorXi::Zero(bayes::get_cols(effect.design_matrix));
            pi
                = {effect.init_pi.value(),
                   Eigen::VectorXi::Zero(effect.init_pi->size())};
        };
    }
    Eigen::VectorXd coeffs;
    Eigen::VectorXd u;

    Eigen::VectorXi tracker;
    Pi pi;

    double variance{};
    double heritability{};
    Eigen::VectorXd marker_variance;
};

struct AdditiveState : GeneticState
{
    using GeneticState::GeneticState;
};

struct DominantEffect : GeneticEffect
{
    using GeneticEffect::GeneticEffect;

    Eigen::VectorXd w;  // freq_q - freq_p
    double ratio_mean{};
    double ratio_variance{};

    detail::NormalParams mean_prior{0.2, 1};
    detail::ScaledInvChiSqParams var_prior{4, 0};
};

struct DominantState : GeneticState
{
    explicit DominantState(const DominantEffect& effect);

    Eigen::VectorXd ratios;

    double ratio_mean{};
    double ratio_variance{};
};

struct Residual
{
    detail::ScaledInvChiSqParams prior{-2, 0};
    double init_variance{0.0};
};

struct ResidualState
{
    Eigen::VectorXd y_adj;
    double variance{0.0};
};
}  // namespace bayes
}  // namespace gelex
