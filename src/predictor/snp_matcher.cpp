#include "snp_matcher.h"

#include <ranges>

#include "gelex/exception.h"

#include "../src/data/loader/bim_loader.h"

namespace gelex
{

SnpMatcher::SnpMatcher(const std::filesystem::path& snp_effect_path)
    : snp_effects_(SnpEffectLoader::load(snp_effect_path))
{
}

MatchPlan SnpMatcher::match(const std::filesystem::path& predict_bed_path) const
{
    auto bim_path = predict_bed_path;
    bim_path.replace_extension(".bim");

    if (!std::filesystem::exists(bim_path))
    {
        throw FileNotFoundException(bim_path);
    }

    detail::BimLoader bim_loader(bim_path);

    std::vector<detail::SnpInfo> predict_snp_meta
        = std::move(bim_loader).take_info();

    MatchPlan match_info(predict_snp_meta.size());

    for (auto [meta, plan] : std::views::zip(predict_snp_meta, match_info))
    {
        auto it = snp_effects_.find(meta.id);

        if (it != snp_effects_.end())
        {
            const auto& model_meta = it->second;
            plan.type = determine_match_type(model_meta, meta);

            if (plan.type != MatchType::skip)
            {
                plan.target_col = model_meta.index;
            }
        }
    }

    return match_info;
}

constexpr char SnpMatcher::normalize_allele(char allele) noexcept
{
    return allele & 0xDF;
}

MatchType SnpMatcher::determine_match_type(
    const SnpEffect& model,
    const detail::SnpInfo& predict) noexcept
{
    const char model_a1 = normalize_allele(model.A1);
    const char model_a2 = normalize_allele(model.A2);
    const char predict_a1 = normalize_allele(predict.A1);
    const char predict_a2 = normalize_allele(predict.A2);

    if (model_a1 == predict_a1 && model_a2 == predict_a2)
    {
        return MatchType::keep;
    }

    if (model_a1 == predict_a2 && model_a2 == predict_a1)
    {
        return MatchType::reverse;
    }

    return MatchType::skip;
}

}  // namespace gelex
