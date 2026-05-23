/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "gelex/algo/infer/mcmc/state.h"

#include <variant>

#include "gelex/exception.h"
#include "gelex/model/bayes/legacy_method.h"
#include "gelex/model/bayes/legacy_prior.h"

namespace gelex::bayes
{

namespace
{

// 从 GeneticSpec 中提取匹配 mode 的 GeneticSpec（JointSpec 按 mode 分支）
auto resolve_spec(const OldGeneticPrior& prior, GeneticMode mode)
    -> const GeneticSpec&
{
    if (const auto* gs = std::get_if<GeneticSpec>(&prior.spec))
    {
        return *gs;
    }
    const auto& js = std::get<JointSpec>(prior.spec);
    return (mode == GeneticMode::A) ? js.additive : js.dominance;
}

auto make_marker_variance(const GeneticSpec& spec, Eigen::Index num_markers)
    -> Eigen::VectorXd
{
    const Eigen::Index size
        = (spec.variance.scope == MarkerVarianceScope::per_marker) ? num_markers
                                                                   : 1;
    return Eigen::VectorXd::Constant(size, spec.variance.init);
}

auto make_group(
    const std::optional<Mixture>& mixture,
    Eigen::Index num_markers,
    Eigen::Index num_samples) -> std::optional<MarkerAllocation>
{
    if (!mixture)
    {
        return std::nullopt;
    }
    if (std::holds_alternative<SpikeSlab>(mixture->strategy))
    {
        return Assignment(num_markers, mixture->proportions.init);
    }
    if (std::holds_alternative<ScaledMixture>(mixture->strategy))
    {
        return ComponentAllocation(
            num_markers, num_samples, mixture->proportions.init);
    }
    throw GelexException("JointMixture unsupported in Step E");
}

auto make_sign(
    const std::optional<CategoricalSpec>& sign_spec,
    Eigen::Index num_markers) -> std::optional<Assignment>
{
    if (!sign_spec)
    {
        return std::nullopt;
    }
    // 内部保留 3-bin 表示：{0.0, positive_prob, negative_prob}
    Eigen::Vector3d sign_3bin{{0.0, sign_spec->init[0], sign_spec->init[1]}};
    return Assignment(num_markers, sign_3bin);
}

}  // namespace

LegacyGeneticState::LegacyGeneticState(
    const GeneticEffect& effect,
    const OldGeneticPrior& prior,
    GeneticMode mode)
    : type(mode),
      coeffs(Eigen::VectorXd::Zero(effect.X.cols())),
      u(Eigen::VectorXd::Zero(effect.X.rows()))
{
    const auto num_markers = effect.X.cols();
    const auto num_samples = effect.X.rows();
    const auto& spec = resolve_spec(prior, mode);

    marker_variance = make_marker_variance(spec, num_markers);
    group = make_group(prior.mixture, num_markers, num_samples);
    sign = make_sign(spec.sign, num_markers);
}

}  // namespace gelex::bayes
