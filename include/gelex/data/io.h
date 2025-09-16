#pragma once

#include <concepts>
#include <cstddef>
#include <expected>
#include <filesystem>
#include <format>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <Eigen/Core>

namespace gelex
{

class DataReader
{
    using encode_vector = std::vector<double>;
    using level_encode_map = std::unordered_map<std::string, encode_vector>;
    using covariate_encoding_map
        = std::unordered_map<std::string, level_encode_map>;
    using id_set = std::unordered_set<std::string>;

   public:
    [[nodiscard]] static DataReader Create(
        const std::filesystem::path& pheno_path,
        const std::filesystem::path& fam_path,
        const std::filesystem::path& qcovar_path,
        const std::filesystem::path& covar_path,
        size_t pheno_col,
        bool iid_only);

    Eigen::VectorXd&& take_phenotype() && { return std::move(phenotype_); }
    Eigen::MatrixXd&& take_fixed() && { return std::move(fixed_); }

    const Eigen::VectorXd& phenotype() const { return phenotype_; }
    const Eigen::MatrixXd& fixed() const { return fixed_; }

    const std::vector<std::string>& final_ids() const { return final_ids_; }

   private:
    DataReader(
        Eigen::VectorXd&& phenotype,
        Eigen::MatrixXd&& fixed,
        std::vector<std::string> final_ids,
        std::unordered_map<std::string, Eigen::Index> id_map)
        : phenotype_(std::move(phenotype)),
          fixed_(std::move(fixed)),
          final_ids_(std::move(final_ids)),
          id_to_idx_map_(std::move(id_map))
    {
    }

    static Eigen::VectorXd load_phenotype(
        const std::filesystem::path& pheno_path,
        const std::unordered_map<std::string, Eigen::Index>& id_map,
        bool iid_only,
        size_t pheno_col);

    static Eigen::MatrixXd load_qcovar(
        const std::filesystem::path& qcovar_path,
        const std::unordered_map<std::string, Eigen::Index>& id_map,
        bool iid_only);

    static Eigen::MatrixXd load_covar(
        const std::filesystem::path& covar_path,
        const std::unordered_map<std::string, Eigen::Index>& id_map,
        bool iid_only);

    static id_set parse_phenotype(
        const std::filesystem::path& path,
        size_t pheno_col,
        bool iid_only);

    static id_set parse_fam(const std::filesystem::path& path, bool iid_only);

    static id_set parse_covar(const std::filesystem::path& path, bool iid_only);

    static covariate_encoding_map encode_covar(
        const std::filesystem::path& path,
        const std::unordered_map<std::string, Eigen::Index>& id_map,
        bool iid_only);

    static std::vector<std::string_view> parse_header(
        std::string_view line,
        const std::filesystem::path& path,
        std::string_view delimiters = "\t");

    static auto parse_id(
        std::string_view line,
        bool iid_only,
        std::string_view delimiters = "\t")
        -> std::expected<std::string, std::errc>;

    static void intersect_in_place(
        std::unordered_set<std::string>& main_set,
        const std::unordered_set<std::string>& other_set);

    std::vector<std::string> final_ids_;
    std::unordered_map<std::string, Eigen::Index> id_to_idx_map_;

    Eigen::VectorXd phenotype_;
    Eigen::MatrixXd fixed_;

    std::filesystem::path pheno_path_;
    std::filesystem::path fam_path_;
    std::filesystem::path qcovar_path_;
    std::filesystem::path covar_path_;
};
}  // namespace gelex
