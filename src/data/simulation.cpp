#include "gelex/data/simulation.h"

#include <cmath>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <Eigen/Core>

#include "data/math_utils.h"
#include "gelex/data/bed_io.h"
#include "gelex/data/io.h"

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
    auto file = detail::open_or_throw<std::ifstream>(causal_variants_list);
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
        auto outfile = detail::open_or_throw<std::ofstream>(
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

    BedIO bed_io(bfile);
    const auto& snp_ids = bed_io.snp_ids();
    const auto& fam_ids = bed_io.fam_ids();

    auto bed_stream = bed_io.create_bed();

    const size_t n_individuals = bed_io.n_individuals();
    const size_t n_snps = bed_io.n_snp();

    const size_t buffer_size = (n_individuals + 3) / 4;
    std::vector<std::byte> buffer(buffer_size);
    Eigen::VectorXd raw_genotype(n_individuals);

    Eigen::VectorXd genetic_values = Eigen::VectorXd::Zero(n_individuals);

    for (size_t snp_idx = 0; snp_idx < n_snps; ++snp_idx)
    {
        const auto& current_snp = snp_ids[snp_idx];

        auto it = causal_variants.find(current_snp);
        bed_stream.read(reinterpret_cast<char*>(buffer.data()), buffer_size);
        if (it == causal_variants.end())
        {
            continue;
        }
        BedIO::read_locus(buffer, raw_genotype);

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
        break;
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

    auto output = detail::open_or_throw<std::ofstream>(bfile + ".phen");
    output << "FID\t" << "IID\t" << "phenotype\n";
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
