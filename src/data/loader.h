#pragma once

#include <filesystem>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

namespace gelex::detail
{

// -----------------------------------------------------------------------------
// PhenotypeLoader
// -----------------------------------------------------------------------------
class PhenotypeLoader
{
   public:
    PhenotypeLoader(
        const std::filesystem::path& path,
        int pheno_column,
        bool iid_only);

    [[nodiscard]] Eigen::VectorXd load(
        const std::unordered_map<std::string, Eigen::Index>& id_map) const;

    const std::string& name() const { return name_; }
    const std::unordered_map<std::string, double>& data() const
    {
        return data_;
    }

   private:
    std::string name_;
    std::unordered_map<std::string, double> data_;
};

// -----------------------------------------------------------------------------
// QcovarLoader (Quantitative Covariates)
// -----------------------------------------------------------------------------
class QcovarLoader
{
   public:
    QcovarLoader(const std::filesystem::path& path, bool iid_only);

    [[nodiscard]] Eigen::MatrixXd load(
        const std::unordered_map<std::string, Eigen::Index>& id_map) const;

    const std::vector<std::string>& names() const { return names_; }
    const std::unordered_map<std::string, std::vector<double>>& data() const
    {
        return data_;
    }

   private:
    std::vector<std::string> names_;
    std::unordered_map<std::string, std::vector<double>> data_;
};

// -----------------------------------------------------------------------------
// CovarLoader (Categorical Covariates)
// -----------------------------------------------------------------------------
class CovarLoader
{
   public:
    using EncodingMap = std::unordered_map<std::string, std::vector<int>>;

    CovarLoader(const std::filesystem::path& path, bool iid_only);

    [[nodiscard]] Eigen::MatrixXd load(
        const std::unordered_map<std::string, Eigen::Index>& id_map) const;

    const std::vector<std::string>& names() const { return names_; }

    const std::unordered_map<std::string, std::vector<std::string>>& data()
        const
    {
        return raw_data_;
    }

   private:
    struct EncodedCovariate
    {
        std::vector<int> encoding;  // The One-Hot vector
        int start_col_offset;       // Offset in the final matrix
    };

    // 预计算的编码表： CovarIndex -> (ValueString -> EncodedVector)
    std::vector<std::unordered_map<std::string, EncodedCovariate>>
        optimized_encodings_;

    std::vector<std::string> names_;
    std::unordered_map<std::string, std::vector<std::string>> raw_data_;

    // 缓存总的 dummy 变量数量，避免重复计算
    Eigen::Index total_dummy_vars_ = 0;

    void build_optimized_encodings();
};

// -----------------------------------------------------------------------------
// BimLoader (SNP Metadata)
// -----------------------------------------------------------------------------
struct SnpMeta
{
    std::string id;
    std::string chrom;
    int position;
    char A1;
    char A2;
};

class BimLoader
{
   public:
    explicit BimLoader(const std::filesystem::path& path);

    const std::vector<SnpMeta>& meta() const { return snp_meta_; }

    std::vector<std::string> get_ids() const;

    const SnpMeta& operator[](size_t index) const { return snp_meta_[index]; }
    size_t size() const { return snp_meta_.size(); }

    std::vector<SnpMeta>&& take_meta() && { return std::move(snp_meta_); }

   private:
    std::vector<SnpMeta> snp_meta_;
    static char detect_delimiter(std::ifstream& file);
};

// -----------------------------------------------------------------------------
// FamLoader (Sample Info)
// -----------------------------------------------------------------------------
class FamLoader
{
   public:
    explicit FamLoader(const std::filesystem::path& path, bool iid_only);

    const std::vector<std::string>& ids() const { return ids_; }
    const std::unordered_map<std::string, Eigen::Index>& data() const
    {
        return data_;
    }

    std::vector<std::string>&& take_ids() && { return std::move(ids_); }

   private:
    std::vector<std::string> ids_;
    std::unordered_map<std::string, Eigen::Index> data_;
};

}  // namespace gelex::detail
