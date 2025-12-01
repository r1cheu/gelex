#ifndef GELEX_DATA_LOADER_SNP_LOADER_H
#define GELEX_DATA_LOADER_SNP_LOADER_H

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <Eigen/Dense>

namespace gelex::detail
{

struct SnpInfo
{
    char chrom;
    std::string id;
    int base_coordinate;
    char A1;
    char A2;
};

class BimLoader
{
   public:
    explicit BimLoader(const std::filesystem::path& path);

    const std::vector<SnpInfo>& info() const { return snp_info_; }
    std::vector<SnpInfo>&& take_info() && { return std::move(snp_info_); }

    std::vector<std::string> get_ids() const;

    const SnpInfo& operator[](size_t index) const { return snp_info_[index]; }
    size_t size() const { return snp_info_.size(); }

   private:
    std::vector<SnpInfo> snp_info_;
    void set_snp_info(char delimiter, std::ifstream& file);
    static char detect_delimiter(std::ifstream& file);
};

}  // namespace gelex::detail

#endif  // GELEX_DATA_LOADER_SNP_LOADER_H
