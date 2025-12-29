#ifndef GELEX_PREDICT_SNP_PREDICTOR_H
#define GELEX_PREDICT_SNP_PREDICTOR_H

#include <Eigen/Core>

#include "gelex/types/snp_info.h"

namespace gelex
{

struct SnpComputeResult
{
    Eigen::VectorXd add;
    Eigen::VectorXd dom;
    [[nodiscard]] Eigen::VectorXd total() const
    {
        if (dom.size() > 0)
        {
            return add + dom;
        }
        return add;
    }
};

class SnpPredictor
{
   public:
    explicit SnpPredictor(const SnpEffects& effects);

    SnpComputeResult compute(
        const Eigen::Ref<const Eigen::MatrixXd>& genotype) const;

   private:
    SnpEffects effects_;
    void validate_dimensions(
        const Eigen::Ref<const Eigen::MatrixXd>& genotype) const;
};

}  // namespace gelex

#endif  // GELEX_PREDICT_SNP_PREDICTOR_H
