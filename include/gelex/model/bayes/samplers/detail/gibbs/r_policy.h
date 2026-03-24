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

#ifndef GELEX_MODEL_BAYES_SAMPLERS_DETAIL_GIBBS_R_POLICY_H_
#define GELEX_MODEL_BAYES_SAMPLERS_DETAIL_GIBBS_R_POLICY_H_

#include <algorithm>
#include <cmath>
#include <numbers>
#include <random>
#include <span>
#include <utility>

#include "gelex/infra/stats/normal.h"
#include "gelex/infra/stats/truncated_normal.h"
#include "gelex/model/bayes/samplers/detail/common_op.h"
#include "gelex/model/bayes/states.h"

namespace gelex::detail::Gibbs::policy
{

struct Symmetric
{
    void init() {}

    static auto num_options(int num_scale_classes) -> int
    {
        return num_scale_classes;
    }

    static auto compute_log_probs(
        std::span<double> ll_buf,
        std::span<const double> logpi_buf,
        const std::array<PosteriorParams, kMaxMixtureComponents>& scale_pp,
        int num_scale_classes) -> double
    {
        ll_buf[0] = logpi_buf[0];
        double max_ll = ll_buf[0];
        for (int c = 1; c < num_scale_classes; ++c)
        {
            ll_buf[c] = scale_pp[c].log_likelihood_kernel + logpi_buf[c];
            max_ll = std::max(ll_buf[c], max_ll);
        }
        return max_ll;
    }

    static auto sample_effect(
        int dist_index,
        const std::array<PosteriorParams, kMaxMixtureComponents>& scale_pp,
        std::mt19937_64& rng) -> std::pair<int, double>
    {
        std::normal_distribution<double> normal{0, 1};
        const auto& pp = scale_pp[dist_index];
        double new_i = (normal(rng) * pp.stddev) + pp.mean;
        return {dist_index, new_i};
    }

    void finalize(std::mt19937_64& /*rng*/) {}
};

struct AT
{
    double log_p_{};
    double log_1mp_{};
    int n_positive_{};
    int n_negative_{};

    void init(const bayes::Assignment& sa)
    {
        log_p_ = std::log(sa.proportion(1));
        log_1mp_ = std::log(sa.proportion(2));
        n_positive_ = 0;
        n_negative_ = 0;
    }

    static auto num_options(int num_scale_classes) -> int
    {
        return (2 * (num_scale_classes - 1)) + 1;
    }

    // Maps 9 AT dist indices to 5 scale class indices:
    //   dist_index:  0    1   2   3   4   5   6   7   8
    //   class:       0    1   1   2   2   3   3   4   4
    static auto at_scale_index(int dist_index) -> int
    {
        return (dist_index == 0) ? 0 : (dist_index + 1) / 2;
    }

    auto compute_log_probs(
        std::span<double> ll_buf,
        std::span<const double> logpi_buf,
        const std::array<PosteriorParams, kMaxMixtureComponents>& scale_pp,
        int num_scale_classes) const -> double
    {
        ll_buf[0] = logpi_buf[0];
        double max_ll = ll_buf[0];

        for (int c = 1; c < num_scale_classes; ++c)
        {
            const double log_common = std::numbers::ln2
                                      + scale_pp[c].log_likelihood_kernel
                                      + logpi_buf[c];
            const double z = scale_pp[c].mean / scale_pp[c].stddev;

            const int idx_pos = (2 * c) - 1;
            ll_buf[idx_pos]
                = log_common + std::log(std::max(norm_cdf(z), 1e-300)) + log_p_;
            max_ll = std::max(ll_buf[idx_pos], max_ll);

            const int idx_neg = 2 * c;
            ll_buf[idx_neg] = log_common
                              + std::log(std::max(norm_cdf(-z), 1e-300))
                              + log_1mp_;
            max_ll = std::max(ll_buf[idx_neg], max_ll);
        }
        return max_ll;
    }

    auto sample_effect(
        int dist_index,
        const std::array<PosteriorParams, kMaxMixtureComponents>& scale_pp,
        std::mt19937_64& rng) -> std::pair<int, double>
    {
        const int class_index = at_scale_index(dist_index);
        const auto& pp = scale_pp[class_index];
        double new_i{};
        if (dist_index % 2 == 1)
        {
            new_i = sample_left_truncated_normal(pp.mean, pp.stddev, 0.0, rng);
            ++n_positive_;
        }
        else
        {
            new_i = sample_right_truncated_normal(pp.mean, pp.stddev, 0.0, rng);
            ++n_negative_;
        }
        return {class_index, new_i};
    }

    void finalize(bayes::Assignment& sa, std::mt19937_64& rng) const
    {
        std::gamma_distribution<double> gamma_a(1.0 + n_positive_, 1.0);
        std::gamma_distribution<double> gamma_b(1.0 + n_negative_, 1.0);
        const double xa = gamma_a(rng);
        const double xb = gamma_b(rng);
        const double p = xa / (xa + xb);
        sa.proportion(1) = p;
        sa.proportion(2) = 1.0 - p;
    }
};

}  // namespace gelex::detail::Gibbs::policy

#endif  // GELEX_MODEL_BAYES_SAMPLERS_DETAIL_GIBBS_R_POLICY_H_
