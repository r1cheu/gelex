#ifndef GELEX_DATA_LOADER_SNP_LOADER_H
#define GELEX_DATA_LOADER_SNP_LOADER_H

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "gelex/types/snp_info.h"
namespace gelex::detail
{

class BimLoader
{
   public:
    explicit BimLoader(const std::filesystem::path& path);

    const SnpEffects& info() const { return snp_effects_; }
    SnpEffects& info() { return snp_effects_; }
    SnpEffects&& take_info() && { return std::move(snp_effects_); }

    std::vector<std::string> get_ids() const;

    size_t size() const { return snp_effects_.size(); }

   private:
    SnpEffects snp_effects_;
    void set_snp_info(char delimiter, std::ifstream& file);
};

}  // namespace gelex::detail

#endif  // GELEX_DATA_LOADER_SNP_LOADER_H
