#pragma once
#include <string>
#include <variant>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "../src/model/effects_manager.h"
#include "distribution.h"
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

struct AdditiveEffect
{
    explicit AdditiveEffect(GenotypeMap&& design_matrix);
    explicit AdditiveEffect(GenotypeMatrix&& design_matrix);
    explicit AdditiveEffect(GenotypeStorage&& design_matrix);

    GenotypeStorage design_matrix;
    Eigen::VectorXd cols_norm;

    detail::ScaledInvChiSqParams prior{4, 0};
    double init_marker_variance{0.0};
    Eigen::Index marker_variance_size{0};

    Eigen::VectorXd pi;

    bool is_monomorphic(Eigen::Index snp_index) const;
    Eigen::Index num_mono() const;
};

struct AdditiveState
{
    explicit AdditiveState(const AdditiveEffect& effect, bool is_mixture_model);

    Eigen::VectorXd coeffs;
    Eigen::VectorXd u;
    Eigen::VectorXi tracker;

    Pi pi;
    double variance{};
    double heritability{};
    Eigen::VectorXd marker_variance;

    bool is_mixture_model_;
};

struct DominantEffect
{
    explicit DominantEffect(GenotypeMap&& design_matrix);
    explicit DominantEffect(GenotypeMatrix&& design_matrix);
    explicit DominantEffect(GenotypeStorage&& design_matrix);

    GenotypeStorage design_matrix;
    Eigen::VectorXd cols_norm;
    Eigen::VectorXd wj;  // freq_q - freq_p

    double ratio_mean{};
    double ratio_variance{};

    bool is_monomorphic(Eigen::Index snp_index) const;
    Eigen::Index num_mono() const;
};

struct DominantState
{
    explicit DominantState(const DominantEffect& effect);

    Eigen::VectorXd coeffs;
    Eigen::VectorXd ratios;
    Eigen::VectorXd u;

    double ratio_mean{};
    double ratio_variance{};

    double variance{};
    double heritability{};
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
