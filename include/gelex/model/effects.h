#pragma once
#include <cstdint>

#include <fmt/base.h>
#include <armadillo>
#include <vector>

namespace gelex
{

enum class effect_type : uint8_t
{
    random,
    genetic,
    gxe,
    residual,
};

enum class BayesAlphabet : uint8_t
{
    A,
    RR,
    B,
    Bpi,
    C,
    Cpi,
    R,
    None,
    Count
};

template <typename Effect>
class Effects
{
   public:
    template <typename... Args>
    void add(Args&&... args)
    {
        effects_.emplace_back(std::forward<Args>(args)...);
        index_map_[effects_.back().name] = effects_.size() - 1;
    }

    const Effect* get(const std::string& name) const
    {
        if (auto it = index_map_.find(name); it != index_map_.end())
        {
            return &effects_[it->second];
        }
        return nullptr;
    }
    Effect* get(const std::string& name)
    {
        if (auto it = index_map_.find(name); it != index_map_.end())
        {
            return &effects_[it->second];
        }
        return nullptr;
    }

    size_t size() const { return effects_.size(); }

    const std::vector<Effect>& effects() const { return effects_; }
    std::vector<Effect>& effects() { return effects_; }

    std::vector<std::string> keys() const
    {
        std::vector<std::string> names;
        names.reserve(effects_.size());
        for (const auto& effect : effects_)
        {
            names.push_back(effect.name);
        }
        return names;
    }

    std::vector<double> values() const
    {
        std::vector<double> vals;
        vals.reserve(effects_.size());
        for (const auto& effect : effects_)
        {
            vals.push_back(effect.sigma);
        }
        return vals;
    }

    void clear()
    {
        effects_.clear();
        index_map_.clear();
    }

    const Effect& operator[](size_t index) const { return effects_[index]; }
    Effect& operator[](size_t index) { return effects_[index]; }
    const Effect& back() const { return effects_.back(); }
    Effect& back() { return effects_.back(); }

    auto begin() { return effects_.begin(); }
    auto end() { return effects_.end(); }
    auto begin() const { return effects_.begin(); }
    auto end() const { return effects_.end(); }

    explicit operator bool() const { return !effects_.empty(); }

   private:
    std::vector<Effect> effects_;
    std::unordered_map<std::string, size_t> index_map_;
};
}  // namespace gelex

namespace fmt
{
template <>
struct formatter<gelex::BayesAlphabet> : formatter<string_view>
{
    auto format(gelex::BayesAlphabet t, format_context& ctx) const
        -> format_context::iterator
    {
        string_view name = "unknown";
        switch (t)
        {
            case gelex::BayesAlphabet::A:
                name = "BayesA";
                break;
            case gelex::BayesAlphabet::RR:
                name = "BayesRR";
                break;
            case gelex::BayesAlphabet::B:
                name = "BayesB";
                break;
            case gelex::BayesAlphabet::Bpi:
                name = "BayesBpi";
                break;
            case gelex::BayesAlphabet::C:
                name = "BayesC";
                break;
            case gelex::BayesAlphabet::Cpi:
                name = "BayesCpi";
                break;
            case gelex::BayesAlphabet::R:
                name = "BayesR";
                break;
            default:
                break;
        }
        return formatter<string_view>::format(name, ctx);
    }
};

}  // namespace fmt
