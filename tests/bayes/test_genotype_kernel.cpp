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

#include <Eigen/Core>
#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <random>
#include <span>
#include <vector>

#include "gelex/bayes/genotype_kernel.h"

#include "bayes/detail/genotype_kernel.h"

namespace
{

using RawDotKernel = gelex::bayes::detail::DotKernel;
using RawAxpyKernel = gelex::bayes::detail::AxpyKernel;
using RawMultiTargetAxpyKernel = gelex::bayes::detail::MultiTargetAxpyKernel;

constexpr std::array<std::size_t, 13>
    TEST_SIZES{0, 1, 3, 4, 7, 8, 15, 16, 17, 31, 32, 33, 517};
constexpr std::array<double, 3> test_scales{0.0, 0.125, -1.75};

auto adapt_raw_kernel(RawDotKernel kernel)
{
    return [kernel](
               std::span<const std::uint8_t> genotype_column,
               const Eigen::Ref<const Eigen::Array4d>& lut,
               std::span<const double> residual) noexcept -> double
    {
        return kernel(
            genotype_column.data(),
            lut.data(),
            residual.data(),
            residual.size());
    };
}

auto adapt_raw_kernel(RawAxpyKernel kernel)
{
    return [kernel](
               std::span<const std::uint8_t> genotype_column,
               const Eigen::Ref<const Eigen::Array4d>& lut,
               double scale,
               std::span<double> residual) noexcept -> void
    {
        kernel(
            genotype_column.data(),
            lut.data(),
            scale,
            residual.data(),
            residual.size());
    };
}

auto adapt_raw_kernel(RawMultiTargetAxpyKernel kernel)
{
    return
        [kernel](
            std::span<const std::uint8_t> genotype_column,
            const Eigen::Ref<const Eigen::Array4d>& lut,
            std::span<const gelex::bayes::AxpyTarget> targets) noexcept -> void
    {
        kernel(
            genotype_column.data(),
            lut.data(),
            targets,
            genotype_column.size());
    };
}

template <typename Dot>
auto check_dot_kernel(Dot kernel) -> void
{
    const Eigen::Array4d lut{{1.5, -0.75, 0.25, 2.0}};
    std::mt19937_64 random_engine{42};
    std::normal_distribution<double> normal_distribution{0.0, 1.0};

    for (const std::size_t size : TEST_SIZES)
    {
        INFO("size = " << size);
        std::vector<std::uint8_t> genotype_column(size);
        Eigen::VectorXd residual(static_cast<Eigen::Index>(size));
        Eigen::VectorXd decoded(static_cast<Eigen::Index>(size));
        for (std::size_t index = 0; index < size; ++index)
        {
            const auto code = static_cast<std::uint8_t>(index % 4);
            genotype_column[index] = code;
            residual[static_cast<Eigen::Index>(index)]
                = normal_distribution(random_engine);
            decoded[static_cast<Eigen::Index>(index)] = lut[code];
        }

        const double expected = decoded.dot(residual);
        const double actual = kernel(
            genotype_column,
            lut,
            std::span<const double>{
                residual.data(), static_cast<std::size_t>(residual.size())});

        REQUIRE(actual == Catch::Approx(expected).epsilon(1e-12).margin(1e-12));
    }
}

template <typename Axpy>
auto check_axpy_kernel(Axpy kernel) -> void
{
    const std::array<Eigen::Array4d, 3> luts{
        Eigen::Array4d{{2.0, 1.0, 1.0, 0.0}},
        Eigen::Array4d{{-0.5, 0.0, 1.5, -0.5}},
        Eigen::Array4d{{0.75, -0.125, -1.25, 2.5}}};
    std::mt19937_64 random_engine{42};
    std::normal_distribution<double> normal_distribution{0.0, 1.0};

    for (const auto& lut : luts)
    {
        for (const std::size_t size : TEST_SIZES)
        {
            INFO("size = " << size);
            std::vector<std::uint8_t> genotype_column(size);
            Eigen::VectorXd initial_residual(static_cast<Eigen::Index>(size));
            Eigen::VectorXd decoded(static_cast<Eigen::Index>(size));
            for (std::size_t index = 0; index < size; ++index)
            {
                const auto code = static_cast<std::uint8_t>(index % 4);
                genotype_column[index] = code;
                initial_residual[static_cast<Eigen::Index>(index)]
                    = normal_distribution(random_engine);
                decoded[static_cast<Eigen::Index>(index)] = lut[code];
            }

            for (const double scale : test_scales)
            {
                INFO("scale = " << scale);
                Eigen::VectorXd expected = initial_residual;
                expected.array() += scale * decoded.array();
                Eigen::VectorXd actual = initial_residual;
                kernel(
                    genotype_column,
                    lut,
                    scale,
                    std::span<double>{
                        actual.data(),
                        static_cast<std::size_t>(actual.size())});

                REQUIRE(actual.isApprox(expected, 1e-13));
            }
        }
    }
}

template <typename Axpy>
auto check_multi_target_axpy_kernel(Axpy kernel) -> void
{
    const Eigen::Array4d lut{{0.75, -0.125, -1.25, 2.5}};
    std::mt19937_64 random_engine{42};
    std::normal_distribution<double> normal_distribution{0.0, 1.0};

    for (const std::size_t size : TEST_SIZES)
    {
        INFO("size = " << size);
        std::vector<std::uint8_t> genotype_column(size);
        Eigen::VectorXd decoded(static_cast<Eigen::Index>(size));
        for (std::size_t index = 0; index < size; ++index)
        {
            const auto code = static_cast<std::uint8_t>(index % 4);
            genotype_column[index] = code;
            decoded[static_cast<Eigen::Index>(index)] = lut[code];
        }

        std::array<Eigen::VectorXd, test_scales.size()> expected;
        std::array<Eigen::VectorXd, test_scales.size()> actual;
        std::array<gelex::bayes::AxpyTarget, test_scales.size()> targets;
        for (std::size_t target_index = 0; target_index < test_scales.size();
             ++target_index)
        {
            actual.at(target_index).resize(static_cast<Eigen::Index>(size));
            for (Eigen::Index index = 0; index < actual.at(target_index).size();
                 ++index)
            {
                actual.at(target_index)(index)
                    = normal_distribution(random_engine);
            }
            expected.at(target_index) = actual.at(target_index);
            expected.at(target_index).array()
                += test_scales.at(target_index) * decoded.array();
            targets.at(target_index) = gelex::bayes::AxpyTarget{
                test_scales.at(target_index), actual.at(target_index)};
        }

        kernel(
            genotype_column, lut, std::span<const gelex::bayes::AxpyTarget>{});
        kernel(genotype_column, lut, targets);

        for (std::size_t target_index = 0; target_index < test_scales.size();
             ++target_index)
        {
            REQUIRE(actual.at(target_index)
                        .isApprox(expected.at(target_index), 1e-13));
        }
    }
}

}  // namespace

TEST_CASE(
    "genotype dot dispatch agrees with Eigen",
    "[bayes][genotype_kernel][dot]")
{
    check_dot_kernel(
        [](std::span<const std::uint8_t> genotype_column,
           const Eigen::Ref<const Eigen::Array4d>& lut,
           std::span<const double> residual) noexcept -> double
        { return gelex::bayes::dot(genotype_column, lut, residual); });
}

TEST_CASE(
    "genotype dot accepts a const LUT column",
    "[bayes][genotype_kernel][dot]")
{
    Eigen::Array<double, 4, Eigen::Dynamic> luts(4, 1);
    luts.col(0) = Eigen::Array4d{{1.5, -0.75, 0.25, 2.0}};
    const auto& const_luts = luts;
    const std::array<std::uint8_t, 4> genotype_column{0, 1, 2, 3};
    const Eigen::Vector4d residual{{0.5, -1.0, 2.0, -0.25}};
    const Eigen::Vector4d decoded = const_luts.col(0).matrix();

    const double expected = decoded.dot(residual);
    const double actual = gelex::bayes::dot(
        genotype_column,
        const_luts.col(0),
        std::span<const double>{residual.data(), 4});

    REQUIRE(actual == Catch::Approx(expected).epsilon(1e-12).margin(1e-12));
}

TEST_CASE(
    "genotype axpy dispatch agrees with Eigen",
    "[bayes][genotype_kernel][axpy]")
{
    check_axpy_kernel(
        [](std::span<const std::uint8_t> genotype_column,
           const Eigen::Ref<const Eigen::Array4d>& lut,
           double scale,
           std::span<double> residual) noexcept -> void
        { gelex::bayes::axpy(genotype_column, lut, scale, residual); });
}

TEST_CASE(
    "genotype axpy accepts a const LUT column",
    "[bayes][genotype_kernel][axpy]")
{
    Eigen::Array<double, 4, Eigen::Dynamic> luts(4, 1);
    luts.col(0) = Eigen::Array4d{{1.5, -0.75, 0.25, 2.0}};
    const auto& const_luts = luts;
    const std::array<std::uint8_t, 4> genotype_column{0, 1, 2, 3};
    const Eigen::Vector4d initial_residual{{0.5, -1.0, 2.0, -0.25}};
    const Eigen::Vector4d decoded = const_luts.col(0).matrix();
    Eigen::Vector4d expected = initial_residual;
    expected.array() += 0.25 * decoded.array();
    Eigen::Vector4d actual = initial_residual;

    gelex::bayes::axpy(
        genotype_column,
        const_luts.col(0),
        0.25,
        std::span<double>{actual.data(), 4});

    REQUIRE(actual.isApprox(expected));
}

TEST_CASE(
    "genotype multi-target axpy dispatch agrees with Eigen",
    "[bayes][genotype_kernel][axpy]")
{
    check_multi_target_axpy_kernel(
        [](std::span<const std::uint8_t> genotype_column,
           const Eigen::Ref<const Eigen::Array4d>& lut,
           std::span<const gelex::bayes::AxpyTarget> targets) noexcept -> void
        { gelex::bayes::axpy(genotype_column, lut, targets); });
}

TEST_CASE(
    "scalar genotype dot agrees with Eigen",
    "[bayes][genotype_kernel][dot]")
{
    check_dot_kernel(adapt_raw_kernel(gelex::bayes::detail::dot_scalar));
}

TEST_CASE(
    "scalar genotype axpy agrees with Eigen",
    "[bayes][genotype_kernel][axpy]")
{
    check_axpy_kernel(adapt_raw_kernel(gelex::bayes::detail::axpy_scalar));
}

TEST_CASE(
    "scalar genotype multi-target axpy agrees with Eigen",
    "[bayes][genotype_kernel][axpy]")
{
    check_multi_target_axpy_kernel(
        adapt_raw_kernel(gelex::bayes::detail::axpy_multi_target_scalar));
}

TEST_CASE(
    "AVX2 genotype dot agrees with Eigen",
    "[bayes][genotype_kernel][dot]")
{
    if (!gelex::bayes::detail::supports_avx2())
    {
        SKIP("AVX2 and FMA are required");
    }
    check_dot_kernel(adapt_raw_kernel(gelex::bayes::detail::dot_avx2));
}

TEST_CASE(
    "AVX2 genotype axpy agrees with Eigen",
    "[bayes][genotype_kernel][axpy]")
{
    if (!gelex::bayes::detail::supports_avx2())
    {
        SKIP("AVX2 and FMA are required");
    }
    check_axpy_kernel(adapt_raw_kernel(gelex::bayes::detail::axpy_avx2));
}

TEST_CASE(
    "AVX2 genotype multi-target axpy agrees with Eigen",
    "[bayes][genotype_kernel][axpy]")
{
    if (!gelex::bayes::detail::supports_avx2())
    {
        SKIP("AVX2 and FMA are required");
    }
    check_multi_target_axpy_kernel(
        adapt_raw_kernel(gelex::bayes::detail::axpy_multi_target_avx2));
}

TEST_CASE(
    "AVX-512 genotype dot agrees with Eigen",
    "[bayes][genotype_kernel][dot]")
{
    if (!gelex::bayes::detail::supports_avx512())
    {
        SKIP("AVX-512F and AVX-512BW are required");
    }
    check_dot_kernel(adapt_raw_kernel(gelex::bayes::detail::dot_avx512));
}

TEST_CASE(
    "AVX-512 genotype axpy agrees with Eigen",
    "[bayes][genotype_kernel][axpy]")
{
    if (!gelex::bayes::detail::supports_avx512())
    {
        SKIP("AVX-512F and AVX-512BW are required");
    }
    check_axpy_kernel(adapt_raw_kernel(gelex::bayes::detail::axpy_avx512));
}

TEST_CASE(
    "AVX-512 genotype multi-target axpy agrees with Eigen",
    "[bayes][genotype_kernel][axpy]")
{
    if (!gelex::bayes::detail::supports_avx512())
    {
        SKIP("AVX-512F and AVX-512BW are required");
    }
    check_multi_target_axpy_kernel(
        adapt_raw_kernel(gelex::bayes::detail::axpy_multi_target_avx512));
}
