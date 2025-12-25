#include "predict_engine.h"

#include <cmath>
#include <filesystem>
#include <format>

#include <omp.h>

#include "../src/data/parser.h"
#include "predict/genotype_aligner.h"
#include "predict_params_pipe.h"
#include "predict_pipe.h"

namespace gelex
{

void PredictEngine::Config::validate() const
{
    if (bed_path.empty())
    {
        throw InvalidInputException("BED path must be provided");
    }
    if (snp_effect_path.empty())
    {
        throw InvalidInputException("SNP effect path must be provided");
    }
    if (covar_effect_path.empty())
    {
        throw InvalidInputException("Covariate effect path must be provided");
    }
    if (output_path.empty())
    {
        throw InvalidInputException("Output path must be provided");
    }
}

PredictEngine::PredictEngine(const Config& config) : config_(config)
{
    config_.validate();

    load_parameters();
    load_data();
    validate_dimensions();
}

void PredictEngine::load_parameters()
{
    PredictParamsPipe::Config params_config{
        .snp_effect_path = config_.snp_effect_path,
        .covar_effect_path = config_.covar_effect_path};

    PredictParamsPipe params_pipe(params_config);

    snp_effects_ = std::move(params_pipe).take_snp_effects();
    covar_effects_ = std::move(params_pipe).take_covar_effects();
}

void PredictEngine::load_data()
{
    PredictDataPipe::Config data_config{
        .bed_path = config_.bed_path,
        .qcovar_path = config_.qcovar_path,
        .dcovar_path = config_.dcovar_path,
        .iid_only = config_.iid_only};

    PredictDataPipe data_pipe(data_config);
    data_ = std::move(data_pipe).take_data();
    sample_ids_ = data_.sample_ids;
    GenotypeAligner genotype_filter(config_.bed_path, snp_effects_);
    data_.genotype = genotype_filter.load(std::move(data_.genotype));
}

void PredictEngine::validate_dimensions()
{
    const Eigen::Index n_snps = data_.genotype.cols();
    const auto n_snp_effects = static_cast<Eigen::Index>(snp_effects_.size());

    if (n_snps != n_snp_effects)
    {
        throw InvalidInputException(
            std::format(
                "Dimension mismatch: genotype matrix has {} SNPs, but SNP "
                "effects has {}",
                n_snps,
                n_snp_effects));
    }
}

void PredictEngine::run()
{
    compute_predictions();

    if (sample_ids_.size() != static_cast<size_t>(predictions_.size()))
    {
        throw InvalidInputException(
            std::format(
                "Dimension mismatch: {} sample IDs but {} predictions",
                sample_ids_.size(),
                predictions_.size()));
    }

    write_output();
}

void PredictEngine::compute_predictions()
{
    compute_covar_predictions();

    snp_predictions_ = compute_snp_predictions();

    predictions_ = snp_predictions_ + covar_predictions_.rowwise().sum();
}

void PredictEngine::compute_covar_predictions()
{
    const Eigen::Index n_samples = data_.genotype.rows();
    const size_t n_cont = data_.qcovariate_names.size();
    const size_t n_cat = data_.dcovariate_names.size();

    covar_predictions_ = Eigen::MatrixXd::Zero(n_samples, 1 + n_cont + n_cat);
    covar_prediction_names_.clear();
    covar_prediction_names_.reserve(1 + n_cont + n_cat);

    if (std::isnan(covar_effects_.intercept))
    {
        throw InvalidInputException("Intercept coefficient is missing or NaN");
    }
    covar_predictions_.col(0).setConstant(covar_effects_.intercept);
    covar_prediction_names_.emplace_back("Intercept");

    if (data_.qcovariates.cols() != static_cast<Eigen::Index>(n_cont + 1))
    {
        throw InvalidInputException(
            std::format(
                "qcovariates matrix has {} columns, expected {} ({} continuous "
                "+ intercept)",
                data_.qcovariates.cols(),
                n_cont + 1,
                n_cont));
    }

    for (size_t i = 0; i < n_cont; ++i)
    {
        const std::string& var_name = data_.qcovariate_names[i];
        auto it = covar_effects_.continuous_coeffs.find(var_name);
        if (it == covar_effects_.continuous_coeffs.end())
        {
            throw InvalidInputException(
                std::format(
                    "Missing coefficient for continuous variable '{}'",
                    var_name));
        }
        covar_predictions_.col(1 + i)
            = data_.qcovariates.col(i + 1) * it->second;
        covar_prediction_names_.push_back(var_name);
    }

    for (size_t i = 0; i < n_cat; ++i)
    {
        const std::string& var_name = data_.dcovariate_names[i];
        auto cat_it = covar_effects_.categorical_coeffs.find(var_name);
        if (cat_it == covar_effects_.categorical_coeffs.end())
        {
            throw InvalidInputException(
                std::format(
                    "Missing coefficient for categorical variable '{}'",
                    var_name));
        }

        const auto& level_coeffs = cat_it->second;
        const auto& levels_vec = data_.dcovariates.at(var_name);

        for (Eigen::Index j = 0; j < n_samples; ++j)
        {
            const std::string& level = levels_vec[j];
            auto level_it = level_coeffs.find(level);
            if (level_it == level_coeffs.end())
            {
                throw InvalidInputException(
                    std::format(
                        "Missing coefficient for level '{}' of variable '{}'",
                        level,
                        var_name));
            }
            covar_predictions_(j, 1 + n_cont + i) = level_it->second;
        }
        covar_prediction_names_.push_back(var_name);
    }
}

Eigen::VectorXd PredictEngine::compute_snp_predictions()
{
    const Eigen::Index n_samples = data_.genotype.rows();
    const Eigen::Index n_snps = data_.genotype.cols();

    Eigen::VectorXd frequencies = snp_effects_.frequencies();
    Eigen::VectorXd additive_effects = snp_effects_.additive_effects();
    const bool has_dominant = snp_effects_.dominance_effects().size() > 0;
    Eigen::VectorXd dominant_effects
        = has_dominant ? snp_effects_.dominance_effects() : Eigen::VectorXd();

    if (frequencies.size() != n_snps || additive_effects.size() != n_snps)
    {
        throw InvalidInputException(
            std::format(
                "Dimension mismatch: genotype has {} SNPs, but frequencies has "
                "{} and additive effects has {}",
                n_snps,
                frequencies.size(),
                additive_effects.size()));
    }
    if (has_dominant && dominant_effects.size() != n_snps)
    {
        throw InvalidInputException(
            std::format(
                "Dimension mismatch: genotype has {} SNPs, but dominant "
                "effects has {}",
                n_snps,
                dominant_effects.size()));
    }

    const double eps = 1e-10;

    Eigen::VectorXd predictions = Eigen::VectorXd::Zero(n_samples);

#pragma omp parallel
    {
        Eigen::VectorXd thread_predictions = Eigen::VectorXd::Zero(n_samples);

#pragma omp for schedule(dynamic) nowait
        for (Eigen::Index j = 0; j < n_snps; ++j)
        {
            const double p = frequencies[j];
            const double q = 1.0 - p;

            const double scale_add = std::sqrt(std::max(2.0 * p * q, eps));
            const double mean_add = 2.0 * p;
            Eigen::VectorXd std_add
                = (data_.genotype.col(j).array() - mean_add) / scale_add;
            thread_predictions += std_add * additive_effects[j];

            if (has_dominant)
            {
                const double mean_dom = 2.0 * p * p;
                const double scale_dom = std::max(2.0 * p * q, eps);

                Eigen::VectorXd dom_transformed
                    = data_.genotype.col(j).unaryExpr(
                        [p](double x) -> double
                        {
                            if (x == 1.0)
                            {
                                return 2.0 * p;
                            }
                            if (x == 2.0)
                            {
                                return (4.0 * p) - 2.0;
                            }
                            return 0.0;
                        });

                Eigen::VectorXd std_dom
                    = (dom_transformed.array() - mean_dom) / scale_dom;
                thread_predictions += std_dom * dominant_effects[j];
            }
        }

#pragma omp critical
        predictions += thread_predictions;
    }

    return predictions;
}

void PredictEngine::write_output()
{
    auto stream
        = detail::open_file<std::ofstream>(config_.output_path, std::ios::out);

    // 写入表头
    stream << "sample_id\ttotal_prediction";
    for (const auto& name : covar_prediction_names_)
    {
        stream << "\t" << name;
    }
    stream << "\tsnp_contribution\n";

    // 写入每行数据
    const Eigen::Index n_samples = predictions_.size();
    for (Eigen::Index i = 0; i < n_samples; ++i)
    {
        stream << sample_ids_[i];
        stream << std::format("\t{:.6f}", predictions_[i]);

        // 协变量贡献
        for (Eigen::Index j = 0; j < covar_predictions_.cols(); ++j)
        {
            stream << std::format("\t{:.6f}", covar_predictions_(i, j));
        }

        // SNP贡献
        stream << std::format("\t{:.6f}", snp_predictions_[i]);
        stream << "\n";
    }
}

}  // namespace gelex
