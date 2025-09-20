#include "../src/data/genotype_mmap.h"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include <gelex/barkeep.h>
#include <spdlog/spdlog.h>
#include <Eigen/Core>
#include <mio/mmap.hpp>

#include "../src/data/loader.h"
#include "../src/estimator/bayes/indicator.h"
#include "gelex/data/bedpipe.h"
#include "gelex/logger.h"

namespace bk = barkeep;

namespace gelex
{
namespace detail
{

using index_map = std::unordered_map<std::string, Eigen::Index>;
using Eigen::Index;

using Eigen::MatrixXd;
using Eigen::Ref;
using Eigen::VectorXd;

struct GenotypeStats
{
    double mean;
    double stddev;
    bool is_monomorphic;
};

struct ProcessedGenotype
{
    Eigen::VectorXd data;
    GenotypeStats stats;
};

void write_metadata(GenotypeFile& file, int64_t rows, int64_t cols)
{
    file.metadata.write(reinterpret_cast<const char*>(&rows), sizeof(int64_t));
    file.metadata.write(reinterpret_cast<const char*>(&cols), sizeof(int64_t));

    // Write mono indices count and data
    const auto mono_count = static_cast<int64_t>(file.mono_indices.size());
    file.metadata.write(
        reinterpret_cast<const char*>(&mono_count), sizeof(int64_t));
    if (mono_count > 0)
    {
        file.metadata.write(
            reinterpret_cast<const char*>(file.mono_indices.data()),
            mono_count * sizeof(int64_t));
    }

    // Write statistics arrays
    file.metadata.write(
        reinterpret_cast<const char*>(file.means.data()),
        cols * sizeof(double));
    file.metadata.write(
        reinterpret_cast<const char*>(file.stddevs.data()),
        cols * sizeof(double));
}

GenotypeStats calculate_genotype_stats(const Ref<const VectorXd>& genotype)
{
    constexpr double MONOMORPHIC_THRESHOLD = 1e-9;

    const double mean = genotype.mean();
    const double variance = (genotype.array() - mean).square().sum()
                            / static_cast<double>(genotype.size() - 1);
    const double stddev = std::sqrt(variance);
    const bool is_monomorphic = stddev < MONOMORPHIC_THRESHOLD;

    return {mean, stddev, is_monomorphic};
}

VectorXd standardize_genotype(
    std::span<const double> raw_genotype,
    std::span<const Index> indices,
    const GenotypeStats& stats)
{
    VectorXd result(indices.size());
    for (Index i = 0; i < static_cast<Index>(indices.size()); ++i)
    {
        result(i) = raw_genotype[indices[i]];
    }

    result.array() -= stats.mean;
    if (!stats.is_monomorphic)
    {
        result.array() /= stats.stddev;
    }

    return result;
}

std::pair<ProcessedGenotype, VectorXd> process_additive_genotype(
    std::span<const double> raw_genotype,
    std::span<const Index> g_indices,
    std::span<const Index> p_indices)
{
    // Reorder genotypes using BedPipe equivalent functionality
    VectorXd ordered_genotype(g_indices.size());
    for (Index i = 0; i < static_cast<Index>(g_indices.size()); ++i)
    {
        ordered_genotype(i) = raw_genotype[g_indices[i]];
    }

    const auto stats = calculate_genotype_stats(ordered_genotype);
    const auto standardized
        = standardize_genotype(raw_genotype, p_indices, stats);

    return {{standardized, stats}, ordered_genotype};
}

ProcessedGenotype process_dominance_genotype(
    const Ref<const VectorXd>& additive_ordered_genotype,
    std::span<const double> raw_genotype,
    std::span<const Index> p_indices)
{
    // Transform additive genotype to dominance (2 -> 0)
    Eigen::VectorXd ordered_dom_genotype = additive_ordered_genotype;
    ordered_dom_genotype = ordered_dom_genotype.unaryExpr(
        [](double val) { return (val == 2.0) ? 0.0 : val; });

    const auto stats = calculate_genotype_stats(ordered_dom_genotype);

    // Transform and standardize genotype for p_list
    Eigen::VectorXd final_genotype(p_indices.size());
    for (Index i = 0; i < static_cast<Index>(p_indices.size()); ++i)
    {
        const double val = raw_genotype[p_indices[i]];
        final_genotype(i) = (val == 2.0) ? 0.0 : val;
    }

    final_genotype.array() -= stats.mean;
    if (!stats.is_monomorphic)
    {
        final_genotype.array() /= stats.stddev;
    }
    return {final_genotype, stats};
}

void write_genotype_data(
    GenotypeFile& file,
    const ProcessedGenotype& processed,
    Index snp_index)
{
    // Write genotype data to .bin file
    file.genotype.write(
        reinterpret_cast<const char*>(processed.data.data()),
        static_cast<int64_t>(processed.data.size() * sizeof(double)));

    // Store statistics for later metadata writing
    file.means.push_back(processed.stats.mean);
    file.stddevs.push_back(processed.stats.stddev);

    // Track monomorphic SNPs
    if (processed.stats.is_monomorphic)
    {
        file.mono_indices.push_back(snp_index);
    }
}

std::pair<double, double> get_rows_and_cols(const std::string& bin_file)
{
    // Read dimensions from metadata file instead of binary file
    const std::string meta_file
        = bin_file.substr(0, bin_file.find_last_of('.')) + ".meta";
    auto file = *detail::openfile<std::ifstream>(
        meta_file, detail::file_type::binary);

    int64_t rows = 0;
    int64_t cols = 0;
    file.read(reinterpret_cast<char*>(&rows), sizeof(int64_t));
    file.read(reinterpret_cast<char*>(&cols), sizeof(int64_t));

    if (rows <= 0 || cols <= 0)
    {
        throw std::runtime_error("Matrix dimensions in metadata are invalid.");
    }

    return {static_cast<double>(rows), static_cast<double>(cols)};
}

double get_genetic_variance(const std::string& bin_file)
{
    const std::string meta_file
        = bin_file.substr(0, bin_file.find_last_of('.')) + ".meta";
    auto file = *detail::openfile<std::ifstream>(
        meta_file, detail::file_type::binary);
    constexpr int64_t size = sizeof(int64_t);

    int64_t rows = 0;
    int64_t cols = 0;
    int64_t n_mono = 0;

    file.read(reinterpret_cast<char*>(&rows), size);
    file.read(reinterpret_cast<char*>(&cols), size);
    file.read(reinterpret_cast<char*>(&n_mono), size);

    if (rows <= 0 || cols <= 0 || n_mono < 0)
    {
        throw std::runtime_error("Matrix dimensions in metadata are invalid.");
    }
    // skip mono indices and mean;
    file.seekg((n_mono * size) + (cols * size), std::ios::cur);

    VectorXd stddev(cols);
    file.read(reinterpret_cast<char*>(stddev.data()), size * cols);

    return stddev.squaredNorm();
}

void create_genotype_binary(
    const std::string& bfile,
    bool dom,
    std::span<const std::string> p_list,
    std::span<const std::string> g_list,
    bool iid_only)
{
    // Create BedPipe for genotype data processing
    auto bed_pipe_result = BedPipe::create(bfile, iid_only);
    if (!bed_pipe_result)
    {
        throw std::runtime_error(
            "Failed to create BedPipe: "
            + std::string(bed_pipe_result.error().message));
    }
    auto bed_pipe = std::move(*bed_pipe_result);

    const size_t n_snp = bed_pipe.num_variants();
    const size_t n_individuals = bed_pipe.sample_size();

    auto logger = gelex::logging::get();

    logger->info("Creating genotype binary files from: [{}].", bfile + ".bed");
    logger->info("Found {} SNPs in [{}].", n_snp, bfile + ".bim");

    GenotypeFile add_file(bfile, false, static_cast<Index>(n_snp));
    std::optional<GenotypeFile> dom_file;

    if (dom)
    {
        dom_file = std::make_optional<GenotypeFile>(bfile, true, n_snp);
        logger->info(
            "Generating both additive and dominance genotype files...");
    }
    else
    {
        logger->info("Generating additive genotype files...");
    }

    // Create index vectors using BedPipe's sample maps
    const auto& raw_sample_map = bed_pipe.sample_map();
    std::vector<Index> g_index;
    g_index.reserve(g_list.size());
    for (const auto& id : g_list)
    {
        auto it = raw_sample_map.find(id);
        if (it == raw_sample_map.end())
        {
            throw std::runtime_error("individual ID not found: " + id);
        }
        g_index.push_back(it->second);
    }

    std::vector<Index> p_index;
    p_index.reserve(p_list.size());
    for (const auto& id : p_list)
    {
        auto it = raw_sample_map.find(id);
        if (it == raw_sample_map.end())
        {
            throw std::runtime_error("individual ID not found: " + id);
        }
        p_index.push_back(it->second);
    }

    // Store dimensions for later header writing
    const auto rows = static_cast<int64_t>(p_list.size());
    const auto cols = static_cast<int64_t>(n_snp);

    std::atomic<int64_t> snp_index{0};

    // Create progress bar using shared style from Indicator
    auto progress_bar = bk::ProgressBar(
        &snp_index,
        bk::ProgressBarConfig<int64_t>{
            .total = cols,
            .format = "{bar} {value}/{total} ({speed:.1f}/s)",
            .speed = 0.1,
            .style = gelex::detail::Indicator::PROGRESS_BAR_STYLE,
            .show = false});
    progress_bar->show();

    // Process each SNP using BedPipe's chunk loading
    const Index chunk_size = 10000;
    for (Index start = 0; start < static_cast<Index>(n_snp);
         start += chunk_size)
    {
        Index end = std::min(start + chunk_size, static_cast<Index>(n_snp));

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
        for (Index chunk_idx = 0; chunk_idx < genotype_chunk.cols();
             ++chunk_idx)
        {
            // Extract raw genotype data for this SNP
            std::vector<double> raw_genotype(n_individuals);
            for (Index i = 0; i < static_cast<Index>(n_individuals); ++i)
            {
                raw_genotype[i] = genotype_chunk(i, chunk_idx);
            }

            // Process additive genotype
            const auto [add_processed, ordered_genotype]
                = process_additive_genotype(raw_genotype, g_index, p_index);
            write_genotype_data(add_file, add_processed, snp_index);

            // Process dominance genotype if requested
            if (dom_file)
            {
                const auto dom_processed = process_dominance_genotype(
                    ordered_genotype, raw_genotype, p_index);
                write_genotype_data(*dom_file, dom_processed, snp_index);
            }

            ++snp_index;
        }
    }

    // Close progress bar
    progress_bar->done();

    // Write metadata files with all collected statistics
    write_metadata(add_file, rows, cols);
    if (dom_file)
    {
        write_metadata(*dom_file, rows, cols);
    }

    // Log completion information
    const std::string add_bin_file = bfile + ".add.bin";
    const std::string add_meta_file = bfile + ".add.meta";

    logger->info("Additive genotype processing completed");
    logger->info(" - Monomorphic SNPs: {}", add_file.mono_indices.size());
    logger->info(" - Binary file: [{}]", add_bin_file);
    logger->info(" - Metadata file: [{}]", add_meta_file);

    if (dom_file)
    {
        const std::string dom_bin_file = bfile + ".dom.bin";
        const std::string dom_meta_file = bfile + ".dom.meta";

        logger->info("Dominance genotype processing completed");
        logger->info(" - Monomorphic SNPs: {}", dom_file->mono_indices.size());
        logger->info(" - Binary file: [{}]", dom_bin_file);
        logger->info(" - Metadata file: [{}]", dom_meta_file);
    }

    logger->info("Genotype binary creation completed successfully.");
}

GenotypeMap::GenotypeMap(const std::string& bin_file)
    : mmap(bin_file), mat(nullptr, 0, 0)
{
    // Read dimensions from metadata file
    const std::string meta_file
        = bin_file.substr(0, bin_file.find_last_of('.')) + ".meta";
    auto meta_stream = *detail::openfile<std::ifstream>(
        meta_file, detail::file_type::binary);

    int64_t rows = 0;
    int64_t cols = 0;
    meta_stream.read(reinterpret_cast<char*>(&rows), sizeof(int64_t));
    meta_stream.read(reinterpret_cast<char*>(&cols), sizeof(int64_t));

    if (rows <= 0 || cols <= 0)
    {
        throw std::runtime_error("Matrix dimensions in metadata are invalid.");
    }

    // Verify binary file size matches expected dimensions
    const size_t expected_size = static_cast<size_t>(rows)
                                 * static_cast<size_t>(cols) * sizeof(double);
    if (mmap.size() != expected_size)
    {
        throw std::runtime_error(
            "Binary file size does not match dimensions in metadata.");
    }

    // Binary file contains only genotype data, no header
    const double* data_ptr = reinterpret_cast<const double*>(mmap.data());
    new (&mat) decltype(mat)(data_ptr, rows, cols);
}

}  // namespace detail
}  // namespace gelex
