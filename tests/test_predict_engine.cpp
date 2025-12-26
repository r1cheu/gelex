#include <cmath>
#include <string>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../src/predict/predict_engine.h"
#include "gelex/exception.h"
#include "predict_engine_fixture.h"

using namespace gelex;
using gelex::test::PredictEngineTestFixture;

namespace
{

double compute_std_additive(double geno, double p)
{
    const double q = 1.0 - p;
    const double eps = 1e-10;
    const double scale = std::sqrt(std::max(2.0 * p * q, eps));
    const double mean = 2.0 * p;
    return (geno - mean) / scale;
}

double compute_std_dominance(double geno, double p)
{
    const double q = 1.0 - p;
    const double eps = 1e-10;
    const double scale_dom = std::max(2.0 * p * q, eps);

    double dom = 0.0;
    if (geno == 0.0)
    {
        dom = 0.0;
    }
    else if (geno == 1.0)
    {
        dom = 2.0 * p;
    }
    else
    {
        dom = (4.0 * p) - 2.0;
    }

    const double mean_dom = 2.0 * p * p;
    return (dom - mean_dom) / scale_dom;
}

}  // namespace

TEST_CASE(
    "PredictEngine - SNP only prediction",
    "[predict][predict_engine][numerical]")
{
    PredictEngineTestFixture fixture;

    Eigen::MatrixXd genotypes(2, 2);
    genotypes << 0.0, 2.0, 1.0, 1.0;

    std::vector<std::string> snp_ids = {"rs1", "rs2"};
    std::vector<std::pair<char, char>> alleles = {{'A', 'C'}, {'T', 'G'}};
    std::vector<std::vector<std::string>> snp_rows
        = {{"1", "1000", "rs1", "A", "C", "0.30", "0.10"},
           {"1", "2000", "rs2", "T", "G", "0.40", "-0.05"}};

    auto [bed_prefix, _] = fixture.create_deterministic_bed_files(
        genotypes,
        {"sample1", "sample2"},
        snp_ids,
        std::vector<std::string>(snp_ids.size(), "1"),
        alleles);

    auto snp_path = fixture.create_snp_effects_file(snp_rows, false);
    auto param_path = fixture.create_param_intercept_only(1.0);

    PredictEngine::Config config{
        .bed_path = bed_prefix,
        .snp_effect_path = snp_path,
        .covar_effect_path = param_path,
        .qcovar_path = "",
        .dcovar_path = "",
        .output_path = "test.predictions",
        .iid_only = false};

    PredictEngine engine(config);
    REQUIRE_NOTHROW(engine.run());

    const double intercept = 1.0;
    const std::vector<double> freqs = {0.3, 0.4};
    const std::vector<double> effects = {0.1, -0.05};
    const auto n_samples = static_cast<size_t>(genotypes.rows());
    const auto n_snps = static_cast<size_t>(genotypes.cols());
    std::vector<double> expected_snp(n_samples, 0.0);
    for (size_t snp = 0; snp < n_snps; ++snp)
    {
        for (size_t sample = 0; sample < n_samples; ++sample)
        {
            const double geno = genotypes(
                static_cast<Eigen::Index>(sample),
                static_cast<Eigen::Index>(snp));
            expected_snp[sample]
                += compute_std_additive(geno, freqs[snp]) * effects[snp];
        }
    }

    const auto& predictions = engine.predictions();
    const auto& snp_preds = engine.snp_predictions();
    const auto& covar_preds = engine.covar_predictions();

    REQUIRE(predictions.size() == 2);
    REQUIRE(snp_preds.size() == 2);
    REQUIRE(covar_preds.rows() == 2);
    REQUIRE(covar_preds.cols() == 1);

    const double eps = 1e-8;
    for (size_t i = 0; i < 2; ++i)
    {
        REQUIRE_THAT(
            predictions[i],
            Catch::Matchers::WithinAbs(intercept + expected_snp[i], eps));
        REQUIRE_THAT(
            snp_preds[i], Catch::Matchers::WithinAbs(expected_snp[i], eps));
        REQUIRE_THAT(
            covar_preds(i, 0), Catch::Matchers::WithinAbs(intercept, eps));
    }
}

TEST_CASE(
    "PredictEngine - Quantitative covariate",
    "[predict][predict_engine][numerical][covariate]")
{
    PredictEngineTestFixture fixture;

    Eigen::MatrixXd genotypes(2, 2);
    genotypes << 0.0, 2.0, 1.0, 1.0;

    std::vector<std::string> snp_ids = {"rs1", "rs2"};
    std::vector<std::pair<char, char>> alleles = {{'A', 'C'}, {'T', 'G'}};
    std::vector<std::vector<std::string>> snp_rows
        = {{"1", "1000", "rs1", "A", "C", "0.30", "0.10"},
           {"1", "2000", "rs2", "T", "G", "0.40", "-0.05"}};

    std::vector<std::string> iids = {"sample1", "sample2"};

    auto [bed_prefix, _] = fixture.create_deterministic_bed_files(
        genotypes,
        iids,
        snp_ids,
        std::vector<std::string>(snp_ids.size(), "1"),
        alleles);

    auto fam_path = bed_prefix;
    fam_path.replace_extension(".fam");
    auto [fids, loaded_iids] = PredictEngineTestFixture::read_fam(fam_path);

    auto snp_path = fixture.create_snp_effects_file(snp_rows, false);
    auto qcovar_path = fixture.create_qcovar_file(
        fids, loaded_iids, {{"Age", {25.0, 30.0}}});
    auto param_path = fixture.create_param_with_qcovar(1.0, {{"Age", 0.2}});

    PredictEngine::Config config{
        .bed_path = bed_prefix,
        .snp_effect_path = snp_path,
        .covar_effect_path = param_path,
        .qcovar_path = qcovar_path.string(),
        .dcovar_path = "",
        .output_path = "test.predictions",
        .iid_only = false};

    PredictEngine engine(config);
    REQUIRE_NOTHROW(engine.run());

    const double intercept = 1.0;
    const double age_coef = 0.2;
    const std::vector<double> ages = {25.0, 30.0};
    const std::vector<double> freqs = {0.3, 0.4};
    const std::vector<double> effects = {0.1, -0.05};
    const auto n_samples = static_cast<size_t>(genotypes.rows());
    const auto n_snps = static_cast<size_t>(genotypes.cols());
    std::vector<double> expected_snp(n_samples, 0.0);
    for (size_t snp = 0; snp < n_snps; ++snp)
    {
        for (size_t sample = 0; sample < n_samples; ++sample)
        {
            const double geno = genotypes(
                static_cast<Eigen::Index>(sample),
                static_cast<Eigen::Index>(snp));
            expected_snp[sample]
                += compute_std_additive(geno, freqs[snp]) * effects[snp];
        }
    }

    const auto& predictions = engine.predictions();
    const auto& covar_preds = engine.covar_predictions();

    REQUIRE(predictions.size() == 2);
    REQUIRE(covar_preds.cols() == 2);

    const double eps = 1e-8;
    for (size_t i = 0; i < 2; ++i)
    {
        double expected_covar = intercept + (age_coef * ages[i]);
        double expected_total = expected_covar + expected_snp[i];

        REQUIRE_THAT(
            predictions[i], Catch::Matchers::WithinAbs(expected_total, eps));
        REQUIRE_THAT(
            covar_preds(i, 0), Catch::Matchers::WithinAbs(intercept, eps));
        REQUIRE_THAT(
            covar_preds(i, 1),
            Catch::Matchers::WithinAbs(age_coef * ages[i], eps));
    }

    const auto& covar_names = engine.covar_prediction_names();
    REQUIRE(covar_names[0] == "Intercept");
    REQUIRE(covar_names[1] == "Age");
}

TEST_CASE(
    "PredictEngine - Categorical covariate",
    "[predict][predict_engine][numerical][covariate]")
{
    PredictEngineTestFixture fixture;

    Eigen::MatrixXd genotypes(3, 2);
    genotypes << 0.0, 1.0, 1.0, 2.0, 2.0, 0.0;

    std::vector<std::string> snp_ids = {"rs1", "rs2"};
    std::vector<std::pair<char, char>> alleles = {{'A', 'C'}, {'T', 'G'}};
    std::vector<std::vector<std::string>> snp_rows
        = {{"1", "1000", "rs1", "A", "C", "0.30", "0.10"},
           {"1", "2000", "rs2", "T", "G", "0.40", "0.05"}};

    std::vector<std::string> iids = {"s1", "s2", "s3"};

    auto [bed_prefix, _] = fixture.create_deterministic_bed_files(
        genotypes,
        iids,
        snp_ids,
        std::vector<std::string>(snp_ids.size(), "1"),
        alleles);

    auto fam_path = bed_prefix;
    fam_path.replace_extension(".fam");
    auto [fids, loaded_iids] = PredictEngineTestFixture::read_fam(fam_path);

    auto snp_path = fixture.create_snp_effects_file(snp_rows, false);
    auto dcovar_path = fixture.create_dcovar_file(
        fids, loaded_iids, {{"Sex", {"M", "F", "M"}}});
    auto param_path = fixture.create_param_with_dcovar(
        1.5, {{"Sex_M", -0.3}, {"Sex_F", 0.2}});

    PredictEngine::Config config{
        .bed_path = bed_prefix,
        .snp_effect_path = snp_path,
        .covar_effect_path = param_path,
        .qcovar_path = "",
        .dcovar_path = dcovar_path.string(),
        .output_path = "test.predictions",
        .iid_only = false};

    PredictEngine engine(config);
    REQUIRE_NOTHROW(engine.run());

    const double intercept = 1.5;
    const double sex_m_coef = -0.3;
    const double sex_f_coef = 0.2;
    const std::vector<double> expected_sex_effects
        = {sex_m_coef, sex_f_coef, sex_m_coef};

    const auto& predictions = engine.predictions();
    const auto& covar_preds = engine.covar_predictions();

    REQUIRE(predictions.size() == 3);
    REQUIRE(covar_preds.cols() == 2);

    const double eps = 1e-8;
    for (size_t i = 0; i < 3; ++i)
    {
        REQUIRE_THAT(
            covar_preds(i, 0), Catch::Matchers::WithinAbs(intercept, eps));
        REQUIRE_THAT(
            covar_preds(i, 1),
            Catch::Matchers::WithinAbs(expected_sex_effects[i], eps));
    }

    const auto& covar_names = engine.covar_prediction_names();
    REQUIRE(covar_names[0] == "Intercept");
    REQUIRE(covar_names[1] == "Sex");
}

TEST_CASE(
    "PredictEngine - Full model",
    "[predict][predict_engine][numerical][covariate]")
{
    PredictEngineTestFixture fixture;

    Eigen::MatrixXd genotypes(3, 2);
    genotypes << 0.0, 2.0, 1.0, 1.0, 2.0, 0.0;

    std::vector<std::string> snp_ids = {"rs1", "rs2"};
    std::vector<std::pair<char, char>> alleles = {{'A', 'C'}, {'T', 'G'}};
    std::vector<std::vector<std::string>> snp_rows
        = {{"1", "1000", "rs1", "A", "C", "0.30", "0.10"},
           {"1", "2000", "rs2", "T", "G", "0.40", "-0.05"}};

    std::vector<std::string> iids = {"s1", "s2", "s3"};

    auto [bed_prefix, _] = fixture.create_deterministic_bed_files(
        genotypes,
        iids,
        snp_ids,
        std::vector<std::string>(snp_ids.size(), "1"),
        alleles);

    auto fam_path = bed_prefix;
    fam_path.replace_extension(".fam");
    auto [fids, loaded_iids] = PredictEngineTestFixture::read_fam(fam_path);

    auto snp_path = fixture.create_snp_effects_file(snp_rows, false);
    auto qcovar_path = fixture.create_qcovar_file(
        fids, loaded_iids, {{"Age", {25.0, 30.0, 35.0}}});
    auto dcovar_path = fixture.create_dcovar_file(
        fids, loaded_iids, {{"Sex", {"M", "F", "M"}}});
    auto param_path = fixture.create_param_full(
        1.0, {{"Age", 0.2}}, {{"Sex_M", -0.3}, {"Sex_F", 0.1}});

    PredictEngine::Config config{
        .bed_path = bed_prefix,
        .snp_effect_path = snp_path,
        .covar_effect_path = param_path,
        .qcovar_path = qcovar_path.string(),
        .dcovar_path = dcovar_path.string(),
        .output_path = "test.predictions",
        .iid_only = false};

    PredictEngine engine(config);
    REQUIRE_NOTHROW(engine.run());

    const double intercept = 1.0;
    const double age_coef = 0.2;
    const double sex_m_coef = -0.3;
    const double sex_f_coef = 0.1;
    const std::vector<double> ages = {25.0, 30.0, 35.0};
    const std::vector<double> sex_effects
        = {sex_m_coef, sex_f_coef, sex_m_coef};

    const auto& predictions = engine.predictions();
    const auto& covar_preds = engine.covar_predictions();

    REQUIRE(predictions.size() == 3);
    REQUIRE(covar_preds.cols() == 3);

    const double eps = 1e-8;
    for (size_t i = 0; i < 3; ++i)
    {
        REQUIRE_THAT(
            covar_preds(i, 0), Catch::Matchers::WithinAbs(intercept, eps));
        REQUIRE_THAT(
            covar_preds(i, 1),
            Catch::Matchers::WithinAbs(age_coef * ages[i], eps));
        REQUIRE_THAT(
            covar_preds(i, 2), Catch::Matchers::WithinAbs(sex_effects[i], eps));
    }

    const auto& covar_names = engine.covar_prediction_names();
    REQUIRE(covar_names[0] == "Intercept");
    REQUIRE(covar_names[1] == "Age");
    REQUIRE(covar_names[2] == "Sex");
}

TEST_CASE(
    "PredictEngine - Dominance effect",
    "[predict][predict_engine][numerical][dominance]")
{
    PredictEngineTestFixture fixture;

    Eigen::MatrixXd genotypes(3, 2);
    genotypes << 0.0, 1.0, 1.0, 2.0, 2.0, 0.0;

    std::vector<std::string> snp_ids = {"rs1", "rs2"};
    std::vector<std::pair<char, char>> alleles = {{'A', 'C'}, {'T', 'G'}};
    std::vector<std::vector<std::string>> snp_rows
        = {{"1", "1000", "rs1", "A", "C", "0.30", "0.10", "0.02"},
           {"1", "2000", "rs2", "T", "G", "0.40", "-0.05", "0.01"}};

    auto [bed_prefix, _] = fixture.create_deterministic_bed_files(
        genotypes,
        {"s1", "s2", "s3"},
        snp_ids,
        std::vector<std::string>(snp_ids.size(), "1"),
        alleles);

    auto snp_path = fixture.create_snp_effects_file(snp_rows, true);
    auto param_path = fixture.create_param_intercept_only(1.0);

    PredictEngine::Config config{
        .bed_path = bed_prefix,
        .snp_effect_path = snp_path,
        .covar_effect_path = param_path,
        .qcovar_path = "",
        .dcovar_path = "",
        .output_path = "test.predictions",
        .iid_only = false};

    PredictEngine engine(config);
    REQUIRE_NOTHROW(engine.run());

    const std::vector<double> p_values = {0.3, 0.4};
    const std::vector<double> add_effects = {0.1, -0.05};
    const std::vector<double> dom_effects = {0.02, 0.01};

    const auto n_samples = static_cast<size_t>(genotypes.rows());
    const auto n_snps = static_cast<size_t>(genotypes.cols());
    std::vector<double> expected_snp(n_samples, 0.0);

    for (size_t snp = 0; snp < n_snps; ++snp)
    {
        const double p = p_values[snp];
        const double add = add_effects[snp];
        const double dom = dom_effects[snp];

        for (size_t sample = 0; sample < n_samples; ++sample)
        {
            const double geno = genotypes(
                static_cast<Eigen::Index>(sample),
                static_cast<Eigen::Index>(snp));

            const double std_add = compute_std_additive(geno, p);
            const double std_dom = compute_std_dominance(geno, p);

            expected_snp[sample] += (std_add * add) + (std_dom * dom);
        }
    }

    const auto& predictions = engine.predictions();
    const auto& snp_preds = engine.snp_predictions();

    REQUIRE(predictions.size() == 3);
    REQUIRE(snp_preds.size() == 3);

    const double eps = 1e-8;
    for (size_t i = 0; i < 3; ++i)
    {
        REQUIRE_THAT(
            snp_preds[i], Catch::Matchers::WithinAbs(expected_snp[i], eps));
    }
}

TEST_CASE("PredictEngine - iid_only mode", "[predict][predict_engine][mode]")
{
    PredictEngineTestFixture fixture;

    Eigen::MatrixXd genotypes(2, 2);
    genotypes << 0.0, 1.0, 1.0, 2.0;

    std::vector<std::string> snp_ids = {"rs1", "rs2"};
    std::vector<std::pair<char, char>> alleles = {{'A', 'C'}, {'T', 'G'}};
    std::vector<std::vector<std::string>> snp_rows
        = {{"1", "1000", "rs1", "A", "C", "0.30", "0.10"},
           {"1", "2000", "rs2", "T", "G", "0.40", "-0.05"}};

    auto [bed_prefix, _] = fixture.create_deterministic_bed_files(
        genotypes,
        {"sample1", "sample2"},
        snp_ids,
        std::vector<std::string>(snp_ids.size(), "1"),
        alleles);

    auto snp_path = fixture.create_snp_effects_file(snp_rows, false);
    auto param_path = fixture.create_param_intercept_only(1.0);

    PredictEngine::Config config{
        .bed_path = bed_prefix,
        .snp_effect_path = snp_path,
        .covar_effect_path = param_path,
        .qcovar_path = "",
        .dcovar_path = "",
        .output_path = "test.predictions",
        .iid_only = true};

    PredictEngine engine(config);
    REQUIRE_NOTHROW(engine.run());

    const auto& sample_ids = engine.sample_ids();
    REQUIRE(sample_ids.size() == 2);

    for (size_t i = 0; i < sample_ids.size(); ++i)
    {
        REQUIRE(sample_ids[i].find('_') == std::string::npos);
    }
}

TEST_CASE("PredictEngine - Error handling", "[predict][predict_engine][error]")
{
    PredictEngineTestFixture fixture;

    Eigen::MatrixXd genotypes(2, 2);
    genotypes << 0.0, 1.0, 1.0, 2.0;

    std::vector<std::string> snp_ids = {"rs1", "rs2"};
    std::vector<std::pair<char, char>> alleles = {{'A', 'C'}, {'T', 'G'}};
    std::vector<std::vector<std::string>> snp_rows
        = {{"1", "1000", "rs1", "A", "C", "0.30", "0.10"},
           {"1", "2000", "rs2", "T", "G", "0.40", "-0.05"}};

    auto [bed_prefix, _] = fixture.create_deterministic_bed_files(
        genotypes,
        {"s1", "s2"},
        snp_ids,
        std::vector<std::string>(snp_ids.size(), "1"),
        alleles);

    auto fam_path = bed_prefix;
    fam_path.replace_extension(".fam");
    auto [fids, iids] = PredictEngineTestFixture::read_fam(fam_path);

    auto snp_path = fixture.create_snp_effects_file(snp_rows, false);
    auto param_path = fixture.create_param_intercept_only(1.0);

    SECTION("Missing output path")
    {
        PredictEngine::Config config{
            .bed_path = bed_prefix,
            .snp_effect_path = snp_path,
            .covar_effect_path = param_path,
            .qcovar_path = "",
            .dcovar_path = "",
            .output_path = "",
            .iid_only = false};

        REQUIRE_THROWS_AS(PredictEngine(config), InvalidInputException);
    }

    SECTION("Missing covariate coefficient")
    {
        auto qcovar_path
            = fixture.create_qcovar_file(fids, iids, {{"Age", {25.0, 30.0}}});

        PredictEngine::Config config{
            .bed_path = bed_prefix,
            .snp_effect_path = snp_path,
            .covar_effect_path = param_path,
            .qcovar_path = qcovar_path.string(),
            .dcovar_path = "",
            .output_path = "test.predictions",
            .iid_only = false};

        PredictEngine engine(config);
        REQUIRE_THROWS_AS(engine.run(), InvalidInputException);
    }

    SECTION("Missing categorical level coefficient")
    {
        auto dcovar_path
            = fixture.create_dcovar_file(fids, iids, {{"Sex", {"M", "F"}}});
        auto param_no_f
            = fixture.create_param_with_dcovar(1.0, {{"Sex_M", -0.3}});

        PredictEngine::Config config{
            .bed_path = bed_prefix,
            .snp_effect_path = snp_path,
            .covar_effect_path = param_no_f,
            .qcovar_path = "",
            .dcovar_path = dcovar_path.string(),
            .output_path = "test.predictions",
            .iid_only = false};

        PredictEngine engine(config);
        REQUIRE_THROWS_AS(engine.run(), InvalidInputException);
    }
}
