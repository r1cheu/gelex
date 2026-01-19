#include "predict_params_pipe.h"

#include <filesystem>

#include "gelex/exception.h"

namespace gelex
{

PredictParamsPipe::PredictParamsPipe(const Config& config)
{
    if (config.snp_effect_path.empty())
    {
        throw InvalidInputException("SNP effect path must be provided");
    }
    load_snp_effects(config.snp_effect_path);

    if (config.covar_effect_path.empty())
    {
        throw InvalidInputException("params effect path must be provided");
    }
    load_covar_effects(config.covar_effect_path);
}

void PredictParamsPipe::load_snp_effects(const std::filesystem::path& path)
{
    detail::SnpEffectLoader loader(path);
    snp_effects_ = std::move(loader).take_effects();
}

void PredictParamsPipe::load_covar_effects(const std::filesystem::path& path)
{
    detail::CovarEffectLoader loader(path);
    covar_effects_ = std::move(loader).take_effects();
}

}  // namespace gelex
