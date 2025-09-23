#pragma once

#include <expected>
#include <filesystem>
#include <format>
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
concept FileStream = std::derived_from<T, std::ios_base>
                     && std::is_constructible_v<
                         T,
                         const std::filesystem::path&,
                         std::ios_base::openmode>;

template <FileStream StreamType>
[[nodiscard]] auto open_file(
    const std::filesystem::path& path,
    std::ios_base::openmode mode) -> std::expected<StreamType, Error>
{
    StreamType stream(path, mode);

    if (!stream.is_open())
    {
        if ((mode & std::ios_base::in))
        {
            if (!std::filesystem::exists(path))
            {
                return std::unexpected(
                    Error{
                        ErrorCode::FileNotFound,
                        std::format("File not found: '{}'", path.string())});
            }
        }

        return std::unexpected(
            Error{
                ErrorCode::FileIOError,
                std::format("Failed to open file: '{}'", path.string())});
    }
    if ((mode & std::ios_base::in) && std::filesystem::file_size(path) == 0)
    {
        return std::unexpected(
            Error{
                ErrorCode::InvalidFile,
                std::format("File is empty: '{}'", path.string())}

        );
    }
    return stream;
}

class PhenotypeLoader
{
   public:
    static auto create(
        const std::filesystem::path& path,
        int pheno_column,
        bool iid_only) -> std::expected<PhenotypeLoader, Error>;

    Eigen::VectorXd load(
        const std::unordered_map<std::string, Eigen::Index>& id_map) const;
    const std::string& name() const { return name_; }
    const std::unordered_map<std::string, double>& data() const
    {
        return data_;
    }

   private:
    explicit PhenotypeLoader(
        std::string&& name,
        std::unordered_map<std::string, double>&& data)
        : name_(std::move(name)), data_(std::move(data))
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

    std::string name_;
    std::unordered_map<std::string, double> data_;
};

class QcovarLoader
{
   public:
    static auto create(const std::filesystem::path& path, bool iid_only)
        -> std::expected<QcovarLoader, Error>;

    Eigen::MatrixXd load(
        const std::unordered_map<std::string, Eigen::Index>& id_map) const;

    const std::vector<std::string>& names() const { return names_; }
    const std::unordered_map<std::string, std::vector<double>>& data() const
    {
        return data_;
    }

   private:
    explicit QcovarLoader(
        std::vector<std::string>&& names,
        std::unordered_map<std::string, std::vector<double>>&& data)
        : names_(std::move(names)), data_(std::move(data))
    {
    }

    static auto get_header(std::ifstream& file)
        -> std::expected<std::vector<std::string>, Error>;

    static auto
    read(std::ifstream& file, size_t expected_columns, bool iid_only) -> std::
        expected<std::unordered_map<std::string, std::vector<double>>, Error>;

    std::vector<std::string> names_;
    std::unordered_map<std::string, std::vector<double>> data_;
};

class CovarLoader
{
   public:
    static auto create(const std::filesystem::path& path, bool iid_only)
        -> std::expected<CovarLoader, Error>;

    Eigen::MatrixXd load(
        const std::unordered_map<std::string, Eigen::Index>& id_map) const;

    const std::vector<std::string>& names() const { return names_; }
    const std::unordered_map<std::string, std::vector<std::string>>& data()
        const
    {
        return data_;
    }

   private:
    explicit CovarLoader(
        std::vector<std::string>&& names,
        std::unordered_map<std::string, std::vector<std::string>>&& data)
        : names_(std::move(names)), data_(std::move(data))
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

    static auto create_encoding_for_one_covariate(
        const std::unordered_set<std::string_view>& unique_levels)
        -> std::unordered_map<std::string, std::vector<int>>;

    static auto collect_unique_levels(
        const std::unordered_map<std::string, std::vector<std::string>>& data,
        size_t num_covariates)
        -> std::vector<std::unordered_set<std::string_view>>;

    std::vector<std::string> names_;
    std::unordered_map<std::string, std::vector<std::string>> data_;
};

class BimLoader
{
   public:
    static auto create(const std::filesystem::path& path)
        -> std::expected<BimLoader, Error>;

    const std::vector<std::string>& ids() const { return snp_ids_; }

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
    static auto create(const std::filesystem::path& path, bool iid_only)
        -> std::expected<FamLoader, Error>;

    const std::vector<std::string>& ids() const { return ids_; }

    const std::unordered_map<std::string, Eigen::Index>& data() const
    {
        return data_;
    }

   private:
    explicit FamLoader(
        std::vector<std::string>&& ids,
        std::unordered_map<std::string, Eigen::Index>&& data)
        : ids_(std::move(ids)), data_(std::move(data))
    {
    }

    static auto read(std::ifstream& file, bool iid_only)
        -> std::expected<std::vector<std::string>, Error>;

    std::vector<std::string> ids_;
    std::unordered_map<std::string, Eigen::Index> data_;
};
}  // namespace gelex::detail
