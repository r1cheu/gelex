#pragma once

#include <fstream>
#include <span>
#include <string>

#include <Eigen/Core>
#include <mio/mmap.hpp>

namespace gelex
{
namespace detail
{

#ifdef USE_AVX512
static constexpr size_t ALIGNMENT = 64;
#else
static constexpr size_t ALIGNMENT = 32;
#endif

struct GenotypeFile
{
    std::ofstream genotype;
    std::ofstream metadata;
    std::vector<int64_t> mono_indices;
    std::vector<double> means;
    std::vector<double> stddevs;

    GenotypeFile(const std::string& bfile, bool dom, Eigen::Index n_snp)
        : genotype(bfile + (dom ? ".dom.bin" : ".add.bin")),
          metadata(bfile + (dom ? ".dom.meta" : ".add.meta"))
    {
        means.reserve(n_snp);
        stddevs.reserve(n_snp);
        mono_indices.reserve(n_snp / 100);  // Reserve space for mono SNPs
    }

    GenotypeFile(const GenotypeFile&) = delete;
    GenotypeFile& operator=(const GenotypeFile&) = delete;
    GenotypeFile(GenotypeFile&&) noexcept = default;
    GenotypeFile& operator=(GenotypeFile&&) noexcept = default;
    ~GenotypeFile() = default;
};

std::pair<double, double> get_rows_and_cols(const std::string& bin_file);

double get_genetic_variance(const std::string& bin_file);

void create_genotype_binary(
    const std::string& bfile,
    bool dom,
    std::span<const std::string> p_list,
    std::span<const std::string> g_list,
    bool iid_only);

struct GenotypeMap
{
    // the order of these members is NOT ALLOWED TO CHANGE
    mio::mmap_source mmap;

#ifdef USE_AVX512
    Eigen::Map<
        const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,
        Eigen::Aligned64>
        mat;
#else
    Eigen::Map<
        const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,
        Eigen::Aligned32>
        mat;

#endif
    explicit GenotypeMap(const std::string& bin_file);
};

}  // namespace detail

}  // namespace gelex
