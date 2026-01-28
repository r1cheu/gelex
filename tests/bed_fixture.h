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

#ifndef GELEX_TEST_BED_FIXTURE_H
#define GELEX_TEST_BED_FIXTURE_H

#include <cstddef>
#include <filesystem>
#include <random>
#include <span>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <utility>

#include "file_fixture.h"

namespace gelex::test
{
bool are_matrices_equal(
    const Eigen::Ref<Eigen::MatrixXd>& mat1,
    const Eigen::Ref<Eigen::MatrixXd>& mat2,
    double tol = 1e-8);

class BedFixture
{
   public:
    BedFixture();
    ~BedFixture() = default;

    BedFixture(const BedFixture&) = delete;
    BedFixture& operator=(const BedFixture&) = delete;
    BedFixture(BedFixture&&) noexcept = default;
    BedFixture& operator=(BedFixture&&) noexcept = default;

    std::pair<std::filesystem::path, Eigen::MatrixXd> create_bed_files(
        Eigen::Index num_samples,
        Eigen::Index num_snps,
        double missing_rate = 0.0,
        double maf_min = 0.05,
        double maf_max = 0.5,
        uint64_t seed = std::random_device{}());

    std::pair<std::filesystem::path, Eigen::MatrixXd> create_deterministic_bed_files(
        const Eigen::MatrixXd& genotypes,
        const std::vector<std::string>& sample_ids = {},
        const std::vector<std::string>& snp_ids = {},
        const std::vector<std::string>& chromosomes = {},
        const std::vector<std::pair<char, char>>& alleles = {});

    [[nodiscard]] std::filesystem::path get_prefix() const noexcept
    {
        return current_prefix_;
    }

    [[nodiscard]] FileFixture& get_file_fixture() noexcept
    {
        return file_fixture_;
    }

   private:
    FileFixture file_fixture_;
    std::filesystem::path current_prefix_;
    std::mt19937_64 rng_;

    static std::vector<std::byte> encode_variant(
        const Eigen::VectorXd& variant);

    void write_bed_file(const Eigen::MatrixXd& genotypes);

    void write_bim_file(
        Eigen::Index num_snps,
        std::span<const std::string> snp_ids,
        std::span<const std::string> chromosomes,
        std::span<const std::pair<char, char>> alleles = {});

    void write_fam_file(
        Eigen::Index num_samples,
        std::span<const std::string> sample_ids);

    static std::pair<char, char> generate_random_alleles(std::mt19937_64& rng);

    static std::vector<std::string> generate_random_sample_ids(
        Eigen::Index num_samples,
        std::mt19937_64& rng);

    static std::vector<std::string> generate_random_snp_ids(
        Eigen::Index num_snps,
        std::mt19937_64& rng);

    static std::vector<std::string> generate_random_chromosomes(
        Eigen::Index num_snps,
        std::mt19937_64& rng);
};

}  // namespace gelex::test

#endif  // GELEX_TEST_BED_FIXTURE_H
