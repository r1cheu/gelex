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
};

class SnpPredictor
{
   public:
    explicit SnpPredictor(const SnpEffects& effects);

    Eigen::VectorXd compute(const Eigen::MatrixXd& genotype);
    SnpComputeResult compute_split(const Eigen::MatrixXd& genotype);

   private:
    SnpEffects effects_;
};

}  // namespace gelex

#endif  // GELEX_PREDICT_SNP_PREDICTOR_H
