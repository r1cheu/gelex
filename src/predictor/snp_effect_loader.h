#pragma once

#include <filesystem>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>

#include <Eigen/Core>

namespace gelex
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

class SnpEffectLoader
{
   public:
    static SnpEffects load(const std::filesystem::path& snp_effect_path);

    static bool has_dom_effects(const std::filesystem::path& snp_effect_path);

   private:
    SnpEffectLoader() = default;

    static ColumnIndices assign_column_indices(
        std::span<const std::string_view> header_columns);
};

}  // namespace gelex
