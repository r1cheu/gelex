#ifndef GELEX_PREDICTOR_PREDICT_PARAMS_PIPE_H
#define GELEX_PREDICTOR_PREDICT_PARAMS_PIPE_H

#include <filesystem>
#include <memory>

#include "../data/loader/snp_effect_loader.h"
#include "covar_effect_loader.h"
#include "gelex/exception.h"

namespace gelex
{

class PredictParamsPipe
{
   public:
    struct Config
    {
        std::filesystem::path snp_effect_path;
        std::filesystem::path covar_effect_path;
    };

    explicit PredictParamsPipe(const Config& config);

    const SnpEffects& snp_effects() const
    {
        if (!has_snp_effects_)
        {
            throw InvalidOperationException("SNP effects not available");
        }
        return snp_effects_;
    }
    const detail::CovarEffectLoader& covar_effect_loader() const
    {
        if (!covar_effect_loader_)
        {
            throw InvalidOperationException(
                "Covariate effect loader not available");
        }
        return *covar_effect_loader_;
    }
    bool has_snp_effects() const { return has_snp_effects_; }
    bool has_covar_effects() const { return covar_effect_loader_ != nullptr; }

    SnpEffects take_snp_effects() &&
    {
        if (!has_snp_effects_)
        {
            throw InvalidOperationException("SNP effects not available");
        }
        return std::move(snp_effects_);
    }
    std::unique_ptr<detail::CovarEffectLoader> take_covar_effect_loader() &&
    {
        if (!covar_effect_loader_)
        {
            throw InvalidOperationException(
                "Covariate effect loader not available");
        }
        return std::move(covar_effect_loader_);
    }

   private:
    void load_snp_effects(const std::filesystem::path& path);
    void load_covar_effects(const std::filesystem::path& path);

    SnpEffects snp_effects_;
    std::unique_ptr<detail::CovarEffectLoader> covar_effect_loader_;
    bool has_snp_effects_ = false;
};

}  // namespace gelex

#endif  // GELEX_PREDICTOR_PREDICT_PARAMS_PIPE_H
