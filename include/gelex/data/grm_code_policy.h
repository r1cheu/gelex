#ifndef GELEX_DATA_GRM_CODE_POLICY_H
#define GELEX_DATA_GRM_CODE_POLICY_H

#include <Eigen/Core>

namespace gelex::grm
{

constexpr double EPSILON = 1e-10;

namespace detail
{
auto additive_mean_center(
    Eigen::Ref<Eigen::MatrixXd> genotype,
    Eigen::VectorXd* freqs = nullptr) -> void;
}  // namespace detail

struct Su
{
    auto operator()(
        Eigen::Ref<Eigen::MatrixXd> genotype,
        bool use_additive,
        Eigen::VectorXd* freqs = nullptr) const -> void;
};

struct Zeng
{
    auto operator()(
        Eigen::Ref<Eigen::MatrixXd> genotype,
        bool use_additive,
        Eigen::VectorXd* freqs = nullptr) const -> void;
};

struct Yang
{
    auto operator()(
        Eigen::Ref<Eigen::MatrixXd> genotype,
        bool use_additive,
        Eigen::VectorXd* freqs = nullptr) const -> void;
};

struct Vitezica
{
    auto operator()(
        Eigen::Ref<Eigen::MatrixXd> genotype,
        bool use_additive,
        Eigen::VectorXd* freqs = nullptr) const -> void;
};

}  // namespace gelex::grm

#endif  // GELEX_DATA_GRM_CODE_POLICY_H
