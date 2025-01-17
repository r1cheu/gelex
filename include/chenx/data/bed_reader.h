#pragma once

#include <omp.h>
#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

#include <armadillo>

// Structure to hold SNP information from .bim file
namespace chenx
{
class BedReader
{
   public:
    BedReader(const BedReader&) = delete;
    BedReader(BedReader&&) = delete;
    BedReader& operator=(const BedReader&) = delete;
    BedReader& operator=(BedReader&&) = delete;
    /**
     * @brief Construct a new BedReader object
     *
     * @param bed_file Path to the .bed file
     * @param chunk_size Number of SNPs to process per chunk
     * @param dosage If true, use "dosage" encoding
     * @param threads Number of threads for parallel processing (0 = let OpenMP
     * decide)
     */
    explicit BedReader(const std::string& bed_file, size_t chunk_size = 10000);

    /**
     * @brief Destroy the BedReader object
     *        Closes the BED file if open.
     */
    ~BedReader();

    /**
     * @brief Check if there are more SNPs to read
     * @return true if there are remaining SNPs in the BED file
     * @return false if no more SNPs are left
     */
    bool HasNext() const;

    /**
     * @brief Read the next chunk of SNP data
     * @return A matrix of size [chunk_size or less x n_individuals]
     *         Each element is an integer genotype code.
     */
    arma::dmat GetNextChunk();

    /**
     * @brief Get the total number of SNPs (from .bim file)
     */
    uint64_t n_snps() const { return snps_.size(); }
    const std::vector<std::string>& snps() const { return snps_; }
    /**
     * @brief Get the number of individuals (from .fam file)
     */
    uint64_t n_individuals() const { return n_individuals_; }

    /**
     * @brief Get the current chunk index
     * @return The index (0-based) of the current chunk
     */
    uint64_t current_chunk_index() const { return current_chunk_index_; }

   private:
    std::ifstream fin_;  ///< File stream for .bed
    std::string bed_file_;
    std::string bim_file_;
    std::string fam_file_;

    std::vector<std::string> snps_;
    uint64_t n_individuals_{};        ///< Number of individuals from .fam
    uint64_t chunk_size_;             ///< Number of SNPs to process per chunk
    uint64_t current_chunk_index_{};  ///< Which chunk are we on?
    uint64_t bytes_per_snp_{};        ///< Computed from n_individuals

    static constexpr double genotype_map[4] = {2.0, 1.0, 1.0, 0.0};

    static uint64_t parseFam(const std::string& fam_file);
    static std::vector<std::string> parseBim(const std::string& bim_file);
    void OpenBed();
};
}  // namespace chenx
