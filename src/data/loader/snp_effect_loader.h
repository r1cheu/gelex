#ifndef GELEX_DATA_LOADER_SNP_EFFECT_LOADER_H
#define GELEX_DATA_LOADER_SNP_EFFECT_LOADER_H

#include <filesystem>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>

#include <Eigen/Core>

namespace gelex
{

struct SnpEffect
{
    Eigen::Index index;
    double A1freq;
    char A1;
    char A2;
    double add;
    double dom;
};

using SnpEffects = std::unordered_map<std::string, SnpEffect>;

}  // namespace gelex

namespace gelex::detail
{

/// Structure to store column indices for .snp.eff file parsing
struct ColumnIndices
{
    int id = -1;
    int a1 = -1;
    int a2 = -1;
    int a1frq = -1;
    int add = -1;
    int dom = -1;

    /// Check if all required columns are present
    [[nodiscard]] bool has_required_columns() const
    {
        return id != -1 && a1 != -1 && a2 != -1 && a1frq != -1 && add != -1;
    }

    /// Get the maximum index required to safely access the row vector
    [[nodiscard]] int max_required_index() const
    {
        int m = std::max({id, a1, a2, a1frq, add});
        if (dom != -1)
        {
            m = std::max(m, dom);
        }
        return m;
    }
};

class SnpEffectLoader
{
   public:
    explicit SnpEffectLoader(const std::filesystem::path& snp_effect_path);

    const SnpEffects& effects() const { return snp_effects_; }
    SnpEffects&& take_effects() && { return std::move(snp_effects_); }
    bool has_dom_effects() const { return has_dom_; }

   private:
    static ColumnIndices assign_column_indices(
        std::span<const std::string_view> header_columns);

    void load(const std::filesystem::path& snp_effect_path);
    static void parse_header(std::string_view line, ColumnIndices& indices);
    void parse_line(
        std::string_view line,
        int line_number,
        const ColumnIndices& indices);

    SnpEffects snp_effects_;
    bool has_dom_ = false;
    Eigen::Index current_index_ = 0;
};

/// Check if a .snp.eff file contains a dominance effect column
/// @param snp_effect_path Path to the .snp.eff file
/// @return true if the file contains a "Dom" column, false otherwise
/// @throws FileFormatException if the file cannot be read or has invalid header
bool has_dom_effect_column(const std::filesystem::path& snp_effect_path);

}  // namespace gelex::detail
#endif  // GELEX_DATA_LOADER_SNP_EFFECT_LOADER_H
