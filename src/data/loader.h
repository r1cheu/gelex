#pragma once

#include <expected>
#include <filesystem>
#include <format>
#include <fstream>
#include <istream>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <Eigen/Dense>

#include "gelex/error.h"

namespace gelex::detail
{

template <typename T>
concept CppStream = std::derived_from<T, std::ios_base>;

enum class file_type : uint8_t
{
    text,
    binary
};

template <CppStream StreamType>
[[nodiscard]] std::expected<StreamType, Error> openfile(
    const std::filesystem::path& path,
    file_type type = file_type::text)
{
    std::ios_base::openmode mode{};
    if constexpr (std::is_base_of_v<std::istream, StreamType>)
    {
        if (!std::filesystem::exists(path))
        {
            return std::unexpected(
                Error{
                    ErrorCode::FileNotFound,
                    std::format("file not found [{}].", path.string())});
        }
        mode |= std::ios_base::in;
    }
    if constexpr (std::is_base_of_v<std::ostream, StreamType>)
    {
        mode |= std::ios_base::out;
        if constexpr (
            std::is_same_v<StreamType, std::ofstream>
            || std::is_same_v<StreamType, std::wofstream>)
        {
            mode |= std::ios_base::trunc;
        }
    }

    if (type == file_type::binary)
    {
        mode |= std::ios_base::binary;
    }

    StreamType stream(path, mode);

    if (!stream.is_open())
    {
        return std::unexpected(
            Error{
                ErrorCode::FileIOError,
                std::format("failed to open file [{}].", path.string())});
    }

    return stream;
}

class PhenotypeLoader
{
   public:
    static auto create(std::string_view path, int pheno_column, bool iid_only)
        -> std::expected<PhenotypeLoader, Error>;

    void intersect(std::unordered_set<std::string>& id_set) const;

    Eigen::VectorXd load(
        const std::unordered_map<std::string, Eigen::Index>& id_map) const;
    const std::string& phenotype_name() const { return phenotype_name_; }
    const std::unordered_map<std::string, double>& phenotype_data() const
    {
        return phenotype_data_;
    }

   private:
    explicit PhenotypeLoader(
        std::string&& phenotype_name,
        std::unordered_map<std::string, double>&& phenotype_data)
        : phenotype_name_(std::move(phenotype_name)),
          phenotype_data_(std::move(phenotype_data))
    {
    }

    static auto get_header(std::ifstream& file, size_t target_column)
        -> std::expected<std::vector<std::string>, Error>;

    static auto read(
        std::ifstream& file,
        size_t target_column,
        size_t expected_columns,
        bool iid_only)
        -> std::expected<std::unordered_map<std::string, double>, Error>;

    std::string phenotype_name_;
    std::unordered_map<std::string, double> phenotype_data_;
};

class QcovarLoader
{
   public:
    static auto create(std::string_view path, bool iid_only)
        -> std::expected<QcovarLoader, Error>;

    void intersect(std::unordered_set<std::string>& id_set) const;

    Eigen::MatrixXd load(
        const std::unordered_map<std::string, Eigen::Index>& id_map) const;

    const std::vector<std::string>& covariate_names() const
    {
        return covariate_names_;
    }
    const std::unordered_map<std::string, std::vector<double>>& covariate_data()
        const
    {
        return covariate_data_;
    }

   private:
    explicit QcovarLoader(
        std::vector<std::string>&& covariate_names,
        std::unordered_map<std::string, std::vector<double>>&& covariate_data)
        : covariate_names_(std::move(covariate_names)),
          covariate_data_(std::move(covariate_data))
    {
    }

    static auto get_header(std::ifstream& file)
        -> std::expected<std::vector<std::string>, Error>;

    static auto
    read(std::ifstream& file, size_t expected_columns, bool iid_only) -> std::
        expected<std::unordered_map<std::string, std::vector<double>>, Error>;

    std::vector<std::string> covariate_names_;
    std::unordered_map<std::string, std::vector<double>> covariate_data_;
};

class CovarLoader
{
   public:
    static auto create(std::string_view path, bool iid_only)
        -> std::expected<CovarLoader, Error>;

    void intersect(std::unordered_set<std::string>& id_set);

    Eigen::MatrixXd load(
        const std::unordered_map<std::string, Eigen::Index>& id_map) const;

    const std::vector<std::string>& covariate_names() const
    {
        return covariate_names_;
    }
    const std::unordered_map<std::string, std::vector<std::string>>&
    covariate_data() const
    {
        return covariate_data_;
    }
    const std::vector<std::unordered_map<std::string, std::vector<int>>>&
    encode_maps() const
    {
        return encode_maps_;
    }

   private:
    explicit CovarLoader(
        std::vector<std::string>&& covariate_names,
        std::unordered_map<std::string, std::vector<std::string>>&&
            covariate_data,
        std::vector<std::unordered_map<std::string, std::vector<int>>>&&
            encode_maps)
        : covariate_names_(std::move(covariate_names)),
          covariate_data_(std::move(covariate_data)),
          encode_maps_(std::move(encode_maps))
    {
    }

    static auto get_header(std::ifstream& file)
        -> std::expected<std::vector<std::string>, Error>;

    static auto
    read(std::ifstream& file, size_t expected_columns, bool iid_only)
        -> std::expected<
            std::unordered_map<std::string, std::vector<std::string>>,
            Error>;

    static auto build_encode_maps(
        const std::unordered_map<std::string, std::vector<std::string>>&
            covariate_data,
        std::span<const std::string> covariate_names)
        -> std::vector<std::unordered_map<std::string, std::vector<int>>>;

    void rebuild_encode_maps();

    static auto create_encoding_for_one_covariate(
        const std::unordered_set<std::string_view>& unique_levels)
        -> std::unordered_map<std::string, std::vector<int>>;

    static auto collect_unique_levels(
        const std::unordered_map<std::string, std::vector<std::string>>& data,
        size_t num_covariates)
        -> std::vector<std::unordered_set<std::string_view>>;

    std::vector<std::string> covariate_names_;
    std::unordered_map<std::string, std::vector<std::string>> covariate_data_;
    std::vector<std::unordered_map<std::string, std::vector<int>>> encode_maps_;
};

class BimLoader
{
   public:
    static auto create(std::string_view path)
        -> std::expected<BimLoader, Error>;

    const std::vector<std::string>& snp_ids() const { return snp_ids_; }

   private:
    explicit BimLoader(std::vector<std::string>&& snp_ids)
        : snp_ids_(std::move(snp_ids))
    {
    }

    static auto read(std::ifstream& file)
        -> std::expected<std::vector<std::string>, Error>;

    std::vector<std::string> snp_ids_;
};

class FamLoader
{
   public:
    static auto create(std::string_view path, bool iid_only)
        -> std::expected<FamLoader, Error>;

    void intersect(std::unordered_set<std::string>& id_set) const;

    const std::unordered_set<std::string>& sample_ids() const
    {
        return sample_ids_;
    }

   private:
    explicit FamLoader(std::unordered_set<std::string>&& individual_ids)
        : sample_ids_(std::move(individual_ids))
    {
    }

    static auto read(std::ifstream& file, bool iid_only)
        -> std::expected<std::unordered_set<std::string>, Error>;

    std::unordered_set<std::string> sample_ids_;
};
}  // namespace gelex::detail
