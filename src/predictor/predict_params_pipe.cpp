#include "predict_params_pipe.h"

#include <filesystem>
#include <memory>

#include "gelex/logger.h"

namespace gelex
{

PredictParamsPipe::PredictParamsPipe(const Config& config)
{
    if (!config.snp_effect_path.empty())
    {
        load_snp_effects(config.snp_effect_path);
    }

    if (!config.covar_effect_path.empty())
    {
        load_covar_effects(config.covar_effect_path);
    }

    auto logger = logging::get();
    if (logger)
    {
        logger->info(
            "Loaded parameters: SNP effects={}, covariate effects={}",
            has_snp_effects_,
            has_covar_effects());
    }
}

void PredictParamsPipe::load_snp_effects(const std::filesystem::path& path)
{
    detail::SnpEffectLoader loader(path);
    snp_effects_ = std::move(loader).take_effects();
    has_snp_effects_ = true;
}

void PredictParamsPipe::load_covar_effects(const std::filesystem::path& path)
{
    covar_effect_loader_ = std::make_unique<detail::CovarEffectLoader>(path);
}

}  // namespace gelex
