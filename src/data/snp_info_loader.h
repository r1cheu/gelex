#pragma once

#include <expected>
#include <filesystem>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "gelex/error.h"

namespace gelex
{

struct SnpInfo
{
    std::string id;
    std::string chrom;
    int position;
    std::string a1;
    std::string a2;
};

class SnpInfoLoader
{
   public:
    /**
     * @brief Create SnpInfoLoader from BIM file path
     *
     * @param bim_file_path Path to BIM file
     * @return std::expected<SnpInfoLoader, Error>
     */
    static auto create(const std::filesystem::path& bim_file_path)
        -> std::expected<SnpInfoLoader, Error>;

    /**
     * @brief Get all SNP information
     *
     * @return const std::vector<SnpInfo>&
     */
    const std::vector<SnpInfo>& snp_info() const { return snp_info_; }

    /**
     * @brief Get SNP information by index
     *
     * @param index SNP index
     * @return const SnpInfo&
     */
    const SnpInfo& operator[](size_t index) const { return snp_info_[index]; }

    /**
     * @brief Get number of SNPs
     *
     * @return size_t
     */
    size_t size() const { return snp_info_.size(); }

   private:
    explicit SnpInfoLoader(std::vector<SnpInfo>&& snp_info)
        : snp_info_(std::move(snp_info))
    {
    }

    static auto read_bim_file(const std::filesystem::path& path)
        -> std::expected<std::vector<SnpInfo>, Error>;

    std::vector<SnpInfo> snp_info_;
};

}  // namespace gelex
