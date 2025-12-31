#pragma once

#include <concepts>

#include <Eigen/Dense>

namespace gelex
{

struct VariantStats
{
    double mean{0.0};
    double stddev{0.0};
    bool is_monomorphic{false};
};

template <typename T>
concept VariantProcessor
    = requires(T processor, Eigen::Ref<Eigen::VectorXd> variant) {
          { T::process_variant(variant) } -> std::same_as<VariantStats>;
      };

namespace detail
{
struct GenotypeCounts
{
    double p_AA{0.0};
    double p_Aa{0.0};
    double p_aa{0.0};
};

VariantStats compute_and_standardize(Eigen::Ref<Eigen::VectorXd> variant);
GenotypeCounts count_frequencies(
    const Eigen::Ref<const Eigen::VectorXd>& variant);

}  // namespace detail
}  // namespace gelex

namespace gelex
{
struct StandardizingProcessor
{
    static VariantStats process_variant(Eigen::Ref<Eigen::VectorXd> variant);
};

struct RawProcessor
{
    static VariantStats process_variant(Eigen::Ref<Eigen::VectorXd> variant);
};

struct HardWenbergProcessor
{
    static VariantStats process_variant(Eigen::Ref<Eigen::VectorXd> variant);
};

struct NOIAProcessor
{
    static VariantStats process_variant(Eigen::Ref<Eigen::VectorXd> variant);
};

struct DominantStandardizingProcessor
{
    static VariantStats process_variant(Eigen::Ref<Eigen::VectorXd> variant);
};

struct DominantRawProcessor
{
    static VariantStats process_variant(Eigen::Ref<Eigen::VectorXd> variant);
};

struct DominantOrthogonalHWEProcessor
{
    static VariantStats process_variant(Eigen::Ref<Eigen::VectorXd> variant);
};

struct DominantNOIAProcessor
{
    static VariantStats process_variant(Eigen::Ref<Eigen::VectorXd> variant);
};

}  // namespace gelex
