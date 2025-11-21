#include "gelex/data/simulation.h"

#include <cmath>
#include <expected>
#include <format>
#include <fstream>
#include <memory>
#include <random>
#include <sstream>
#include <string_view>
#include <vector>

#include <Eigen/Core>

#include "../src/data/loader.h"
#include "../src/utils/math_utils.h"
#include "gelex/data/bed_pipe.h"
#include "gelex/data/sample_manager.h"

namespace gelex
{
using Eigen::MatrixXd;
using Eigen::VectorXd;

auto PhenotypeSimulator::create(const Config& config)
    -> std::expected<PhenotypeSimulator, Error>
{
    if (config.heritability <= 0.0 || config.heritability >= 1.0)
    {
        return std::unexpected(
            Error{
                ErrorCode::InvalidArgument, "Heritability must be in (0, 1)"});
    }

    if (!std::filesystem::is_regular_file(config.bed_path))
    {
        return std::unexpected(
            Error{
                ErrorCode::FileNotFound,
                std::format(
                    "BED file not found: {}", config.bed_path.string())});
    }

    if (!std::filesystem::is_regular_file(config.causal_variants_path))
    {
        return std::unexpected(
            Error{
                ErrorCode::FileNotFound,
                std::format(
                    "Causal variants file not found: {}",
                    config.causal_variants_path.string())});
    }

    return PhenotypeSimulator(config);
}

PhenotypeSimulator::PhenotypeSimulator(Config config)
    : config_(std::move(config))
{
}

auto PhenotypeSimulator::simulate() -> std::expected<void, Error>
{
    // Stage 1: Initialization
    if (auto result = initialize_rng(); !result)
    {
        return std::unexpected(result.error());
    }

    auto causal_effects_exp = load_or_generate_causal_effects();
    if (!causal_effects_exp)
    {
        return std::unexpected(causal_effects_exp.error());
    }

    // Stage 2: Data Processing Pipeline Setup
    auto fam_path = config_.bed_path;
    fam_path.replace_extension(".fam");

    auto sample_manager_result = SampleManager::create(fam_path);
    if (!sample_manager_result)
    {
        return std::unexpected(sample_manager_result.error());
    }
    sample_manager_result->finalize();
    auto sample_ptr = std::make_shared<SampleManager>(
        std::move(sample_manager_result.value()));

    auto bed_pipe_exp = BedPipe::create(config_.bed_path, sample_ptr);

    if (!bed_pipe_exp)
    {
        return std::unexpected(bed_pipe_exp.error());
    }
    auto& bed_pipe = *bed_pipe_exp;

    auto genetic_values_exp
        = calculate_genetic_values(bed_pipe, *causal_effects_exp);
    if (!genetic_values_exp)
    {
        return std::unexpected(genetic_values_exp.error());
    }

    VectorXd phenotypes = generate_phenotypes(*genetic_values_exp);

    return write_results(phenotypes, sample_ptr);
}

// --- Private Helper Methods ---

auto PhenotypeSimulator::initialize_rng() -> std::expected<void, Error>
{
    try
    {
        if (config_.seed == -1)
        {
            std::random_device rd;
            rng_.seed(rd());
        }
        else
        {
            rng_.seed(config_.seed);
        }
        return {};
    }
    catch (const std::exception& e)
    {
        return std::unexpected(
            Error{
                ErrorCode::InvalidArgument,
                std::format("Failed to initialize RNG: {}", e.what())});
    }
}

auto PhenotypeSimulator::load_or_generate_causal_effects()
    -> std::expected<std::unordered_map<std::string, double>, Error>
{
    auto file_result = detail::open_file<std::ifstream>(
        config_.causal_variants_path, std::ios::in);
    if (!file_result)
    {
        return std::unexpected(file_result.error());
    }
    auto& file = *file_result;

    std::unordered_map<std::string, double> variants;
    std::normal_distribution<double> dist(0.0, 1.0);
    std::string line;
    bool effects_were_generated = false;

    while (std::getline(file, line))
    {
        if (line.empty())
        {
            continue;
        }

        std::stringstream ss(line);
        std::string id;
        double effect{};
        ss >> id;

        if (ss >> effect)
        {
            variants.emplace(std::move(id), effect);
        }
        else
        {
            effects_were_generated = true;
            variants.emplace(std::move(id), dist(rng_));
        }
    }

    // If we generated effects, write them to a new file as a side-effect.
    if (effects_were_generated)
    {
        auto effects_path = config_.causal_variants_path;
        effects_path.concat(".effects");

        auto outfile_result
            = detail::open_file<std::ofstream>(effects_path, std::ios::out);
        if (!outfile_result)
        {
            return std::unexpected(
                Error{
                    ErrorCode::FileIOError,
                    std::format(
                        "Failed to create effects output file: {}",
                        effects_path.string())});
        }
        auto& outfile = *outfile_result;
        for (const auto& [id, eff] : variants)
        {
            outfile << std::format("{}\t{}\n", id, eff);
        }
    }

    return variants;
}

auto PhenotypeSimulator::calculate_genetic_values(
    BedPipe& bed_pipe,
    const std::unordered_map<std::string, double>& causal_effects)
    -> std::expected<VectorXd, Error>
{
    const auto n_individuals = bed_pipe.num_samples();
    const auto n_snps = bed_pipe.num_variants();
    const auto& snp_ids = bed_pipe.snp_ids();

    VectorXd genetic_values = VectorXd::Zero(n_individuals);

    // Helper lambda for efficient, in-place standardization
    auto standardize = [n_individuals](Eigen::Ref<VectorXd> vec)
    {
        if (vec.size() <= 1)
        {
            return;
        }
        const double mean = vec.mean();
        vec.array() -= mean;
        const double stddev
            = std::sqrt(vec.squaredNorm() / (n_individuals - 1));
        if (stddev > 1e-10)
        {  // Avoid division by zero for monomorphic variants
            vec.array() /= stddev;
        }
        else
        {
            vec.setZero();
        }
    };

    for (Eigen::Index start = 0; start < n_snps; start += SNP_CHUNK_SIZE)
    {
        const Eigen::Index end = std::min(start + SNP_CHUNK_SIZE, n_snps);

        auto genotype_chunk_exp = bed_pipe.load_chunk(start, end);
        if (!genotype_chunk_exp)
        {
            return std::unexpected(genotype_chunk_exp.error());
        }
        MatrixXd genotype_chunk = std::move(*genotype_chunk_exp);

        for (Eigen::Index i = 0; i < genotype_chunk.cols(); ++i)
        {
            const auto& current_snp_id = snp_ids[start + i];
            auto it = causal_effects.find(current_snp_id);
            if (it == causal_effects.end())
            {
                continue;
            }
            standardize(genotype_chunk.col(i));
            genetic_values += genotype_chunk.col(i) * it->second;
        }
    }
    return genetic_values;
}

auto PhenotypeSimulator::generate_phenotypes(const VectorXd& genetic_values)
    -> VectorXd
{
    // Calculate variance of genetic values
    const double genetic_variance = detail::var(genetic_values)(0);

    // Derive residual variance from heritability: h^2 = Vg / (Vg + Ve)
    // => Ve = Vg * (1/h^2 - 1)
    const double residual_variance
        = genetic_variance * (1.0 / config_.heritability - 1.0);

    std::normal_distribution<double> residual_dist(
        0.0, std::sqrt(residual_variance));

    VectorXd residuals(genetic_values.size());
    for (Eigen::Index i = 0; i < residuals.size(); ++i)
    {
        residuals(i) = residual_dist(rng_);
    }

    return genetic_values + residuals;
}

auto PhenotypeSimulator::write_results(
    const VectorXd& phenotypes,
    const std::shared_ptr<SampleManager>& sample_manager) const
    -> std::expected<void, Error>
{
    const auto output_path = config_.output_path.empty()
                                 ? std::filesystem::path(config_.bed_path)
                                       .replace_extension(".phen")
                                 : config_.output_path;

    auto output_result
        = detail::open_file<std::ofstream>(output_path, std::ios::out);
    if (!output_result)
    {
        return std::unexpected(
            Error{
                ErrorCode::FileIOError,
                std::format(
                    "Failed to create phenotype output file: {}",
                    output_path.string())});
    }
    auto& output = *output_result;

    output << "FID\tIID\tphenotype\n";

    const auto& sample_ids = sample_manager->common_ids();
    for (size_t i = 0; i < sample_ids.size(); ++i)
    {
        const std::string_view full_id(sample_ids[i]);
        std::string_view fid = full_id;
        std::string_view iid = full_id;

        if (const auto pos = full_id.find('_'); pos != std::string_view::npos)
        {
            fid = full_id.substr(0, pos);
            iid = full_id.substr(pos + 1);
        }

        output << std::format("{}\t{}\t{}\n", fid, iid, phenotypes[i]);
    }

    return {};
}

}  // namespace gelex
