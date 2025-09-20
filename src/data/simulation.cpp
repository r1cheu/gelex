#include "gelex/data/simulation.h"

#include <cmath>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <Eigen/Core>

#include "../src/data/loader.h"
#include "data/math_utils.h"
#include "gelex/data/bedpipe.h"

namespace gelex
{
using Eigen::Ref;

using Eigen::MatrixXd;
using Eigen::VectorXd;

PhenotypeSimulator::PhenotypeSimulator() : rng_(std::random_device{}()) {}

void PhenotypeSimulator::initialize_rng(int seed)
{
    if (seed == -1)
    {
        std::random_device rd;
        rng_.seed(rd());
    }
    else
    {
        rng_.seed(seed);
    }
}

std::unordered_map<std::string, double>
PhenotypeSimulator::load_causal_variants(
    const std::string& causal_variants_list)
{
    auto file = *detail::openfile<std::ifstream>(causal_variants_list);
    std::normal_distribution<double> dist(0, 1);

    std::unordered_map<std::string, double> variants;
    std::string line;

    bool no_effects{false};

    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        std::string id;
        double effect{};

        ss >> id;

        if (ss >> effect)
        {
            variants.emplace(id, effect);
        }
        else
        {
            no_effects = true;
            variants.emplace(id, dist(rng_));
        }
    }

    if (no_effects)
    {
        auto outfile = *detail::openfile<std::ofstream>(
            causal_variants_list + ".effects");
        for (const auto& [id, eff] : variants)
        {
            outfile << id << "\t" << eff << "\n";
        }
    }
    return variants;
}

void PhenotypeSimulator::simulate_qt_from_bed(
    const std::string& bfile,
    const std::string& causal_variants_list,
    double heritability,
    int seed)
{
    if (heritability <= 0.0 || heritability >= 1.0)
    {
        throw std::invalid_argument("Heritability must be between 0 and 1");
    }

    auto causal_variants = load_causal_variants(causal_variants_list);
    initialize_rng(seed);

    // Create BedPipe for genotype data processing
    auto bed_pipe_result = BedPipe::create(bfile);
    if (!bed_pipe_result)
    {
        throw std::runtime_error(
            "Failed to create BedPipe: "
            + std::string(bed_pipe_result.error().message));
    }
    auto bed_pipe = std::move(*bed_pipe_result);

    const auto& snp_ids = bed_pipe.snp_ids();
    const auto& sample_ids = bed_pipe.raw_sample_ids();

    const size_t n_individuals = bed_pipe.raw_sample_size();
    const size_t n_snps = bed_pipe.num_variants();

    Eigen::VectorXd genetic_values = Eigen::VectorXd::Zero(n_individuals);

    // Process SNPs in chunks for efficiency
    const Eigen::Index chunk_size = 10000;
    for (Eigen::Index start = 0; start < static_cast<Eigen::Index>(n_snps);
         start += chunk_size)
    {
        Eigen::Index end
            = std::min(start + chunk_size, static_cast<Eigen::Index>(n_snps));

        // Load genotype chunk
        auto genotype_result = bed_pipe.load_chunk(start, end);
        if (!genotype_result)
        {
            throw std::runtime_error(
                "Failed to load genotype chunk: "
                + std::string(genotype_result.error().message));
        }
        Eigen::MatrixXd genotype_chunk = std::move(*genotype_result);

        // Process each SNP in the chunk
        for (Eigen::Index chunk_idx = 0; chunk_idx < genotype_chunk.cols();
             ++chunk_idx)
        {
            const auto& current_snp = snp_ids[start + chunk_idx];
            auto it = causal_variants.find(current_snp);
            if (it == causal_variants.end())
            {
                continue;
            }

            // Extract raw genotype data for this SNP
            Eigen::VectorXd raw_genotype = genotype_chunk.col(chunk_idx);

            double mean = raw_genotype.mean();
            double sum_sq = (raw_genotype.array() - mean).square().sum();
            double stddev
                = std::sqrt(sum_sq / static_cast<double>(n_individuals - 1));

            Eigen::VectorXd standardized = raw_genotype;
            standardized.array() -= mean;
            if (stddev > 1e-10)
            {
                standardized.array() /= stddev;
            }
            genetic_values += standardized * it->second;
        }
    }

    double genetic_variance = detail::var(genetic_values)(0);
    double residual_variance = genetic_variance * (1.0 / heritability - 1.0);

    std::normal_distribution<double> normal_dist(
        0.0, std::sqrt(residual_variance));
    Eigen::VectorXd residuals(n_individuals);
    for (int i = 0; i < residuals.size(); ++i)
    {
        residuals(i) = normal_dist(rng_);
    }
    Eigen::VectorXd phenotype = genetic_values + residuals;

    auto output = *detail::openfile<std::ofstream>(bfile + ".phen");
    output << "FID\t" << "IID\t" << "phenotype\n";
    // Convert sample IDs to individual format
    std::vector<std::string> fam_ids;
    for (const auto& sample_id : sample_ids)
    {
        fam_ids.push_back(sample_id);
    }

    for (int i = 0; i < n_individuals; ++i)
    {
        std::string fid;
        std::string iid;
        size_t pos = fam_ids[i].find('_');
        if (pos != std::string::npos)
        {
            fid = fam_ids[i].substr(0, pos);
            iid = fam_ids[i].substr(pos + 1);
        }
        else
        {
            fid = fam_ids[i];
            iid = fam_ids[i];
        }
        output << fid << "\t" << iid << "\t" << phenotype[i] << "\n";
    }
}

}  // namespace gelex
