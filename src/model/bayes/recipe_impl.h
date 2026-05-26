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

#ifndef GELEX_MODEL_BAYES_RECIPE_IMPL_H_
#define GELEX_MODEL_BAYES_RECIPE_IMPL_H_

#include <memory>
#include <string_view>
#include <vector>

#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/model/bayes/recipe_options.h"

namespace gelex
{
class BayesModel;
}

namespace gelex::bayes
{

class BayesRecipeImpl
{
   public:
    BayesRecipeImpl(const BayesRecipeImpl&) = delete;
    auto operator=(const BayesRecipeImpl&) -> BayesRecipeImpl& = delete;
    BayesRecipeImpl(BayesRecipeImpl&&) noexcept = delete;
    auto operator=(BayesRecipeImpl&&) noexcept -> BayesRecipeImpl& = delete;
    virtual ~BayesRecipeImpl() = default;

    virtual auto make_genetic_priors(const BayesModel& model) const
        -> std::vector<std::unique_ptr<GeneticPrior>> = 0;
    virtual auto make_genetic_prior_blocks(const BayesModel& model) const
        -> std::vector<GeneticPriorBlockV2> = 0;

   protected:
    BayesRecipeImpl(std::string_view name, const BayesRecipeConfig& options);

    auto name() const -> std::string_view { return name_; }
    auto options() const -> const BayesRecipeConfig& { return options_; }

    static auto default_heritability(GeneticMode mode)
        -> OpenUnitInterval<double>;
    static auto marker_variance_from_heritability(
        const BayesModel& model,
        GeneticMode mode,
        double heritability,
        double active_marker_weight) -> double;
    static auto make_marker_variance(
        MarkerVarianceLayout layout,
        double target_marker_variance) -> MarkerVariance;
    static auto make_mixture_proportion(
        const Simplex<double>& proportion,
        UpdatePolicy update) -> MixtureProportion;

    auto reject_dominance_positive_probability_override() const -> void;

   private:
    std::string_view name_;
    const BayesRecipeConfig& options_;
};

class IndependentMethod : public BayesRecipeImpl
{
   protected:
    IndependentMethod(std::string_view name, const BayesRecipeConfig& options);

    template <typename Fn>
    auto for_each_effect(Fn fn) const -> void
    {
        for (const auto mode : options().modes)
        {
            switch (mode)
            {
                case GeneticMode::A:
                    fn(GeneticMode::A, options().additive);
                    break;
                case GeneticMode::D:
                    fn(GeneticMode::D, options().dominance);
                    break;
            }
        }
    }

    auto reject_joint_overrides() const -> void;
    auto reject_per_effect_proportion() const -> void;
    auto reject_per_effect_multiplier() const -> void;
    auto require_paired_proportion_and_multiplier() const -> void;

   private:
    auto make_genetic_priors(const BayesModel& model) const
        -> std::vector<std::unique_ptr<GeneticPrior>> final;
    auto make_genetic_prior_blocks(const BayesModel& model) const
        -> std::vector<GeneticPriorBlockV2> final;

    virtual auto make_genetic_prior(
        GeneticMode mode,
        const EffectConfig& effect,
        const BayesModel& model) const -> std::unique_ptr<GeneticPrior> = 0;
    virtual auto make_single_genetic_prior(
        GeneticMode mode,
        const EffectConfig& effect,
        const BayesModel& model) const
        -> std::unique_ptr<SingleGeneticPrior> = 0;
};

class JointMethod : public BayesRecipeImpl
{
   protected:
    JointMethod(std::string_view name, const BayesRecipeConfig& options);

    auto require_both_modes() const -> void;
    auto reject_per_effect_proportion() const -> void;
    auto reject_per_effect_multiplier() const -> void;

   private:
    auto make_genetic_priors(const BayesModel& model) const
        -> std::vector<std::unique_ptr<GeneticPrior>> final;
    auto make_genetic_prior_blocks(const BayesModel& model) const
        -> std::vector<GeneticPriorBlockV2> final;

    virtual auto make_joint_prior(
        const BayesRecipeConfig& config,
        const BayesModel& model) const -> std::unique_ptr<GeneticPrior> = 0;
    virtual auto make_joint_prior_v2(
        const BayesRecipeConfig& config,
        const BayesModel& model) const
        -> std::unique_ptr<JointGeneticPrior> = 0;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_RECIPE_IMPL_H_
