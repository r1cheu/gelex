#pragma once

#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace gelex
{

namespace detail
{
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
    auto cbegin() const { return effects_.cbegin(); }
    auto cend() const { return effects_.cend(); }

    explicit operator bool() const { return !effects_.empty(); }
    bool empty() const { return effects_.empty(); }

   private:
    std::vector<Effect> effects_;
    std::unordered_map<std::string, size_t> index_map_;
};
}  // namespace detail
}  // namespace gelex
