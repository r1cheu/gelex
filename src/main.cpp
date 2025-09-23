#include <argparse/argparse.hpp>

#include "gelex/data/data_pipe.h"
#include "gelex/data/genotype_pipe.h"
#include "gelex/data/sample_manager.h"
#include "gelex/logger.h"

int main(int argc, char* argv[])
{
    argparse::ArgumentParser program("gelex", "0.4.0");
    argparse::ArgumentParser fit("fit");
    fit.add_description(
        "Fit a genomic prediction model using Bayesian or GBLUP methods");
    fit.add_argument("--bfile")
        .help("prefix for PLINK1 binary files")
        .required();
    fit.add_argument("--pheno")
        .help("phenotype file with columns [FID, IID, A, B...], sep by tab")
        .required();
    fit.add_argument("--pheno-col")
        .help(
            "specify which phenotype column to use, default is the 3rd column")
        .default_value(3)
        .scan<'i', int>();
    fit.add_argument("--chunk-size")
        .help("chunk size for processing snps, default is 10000")
        .default_value(10000)
        .scan<'i', int>();
    fit.add_argument("--dom")
        .help("enable estimation of dominance effects")
        .flag();
    fit.add_argument("--qcovar")
        .default_value("")
        .help(
            "quantitative covariate file with columns [FID, IID, C1, C2...], "
            "sep by tab");
    fit.add_argument("--covar").default_value("").help(
        "categorical covariate file with columns [FID, IID, C1, C2...], sep by "
        "tab");
    fit.add_argument("-o", "--out").help("output prefix").required();
    fit.add_argument("-m", "--method")
        .help("genomic prediction method: A, B(pi), C(pi), R, RR, GBLUP")
        .default_value("RR")
        .choices("A", "B", "Bpi", "C", "Cpi", "R", "RR", "GBLUP")
        .required();
    fit.add_argument("--iid_only")
        .help("use IID for sample ID, default is false, which will use FID_IID")
        .flag();

    program.add_subparser(fit);

    try
    {
        program.parse_args(argc, argv);
    }
    catch (const std::exception& err)
    {
        std::cerr << err.what() << "\n";
        if (program.is_subcommand_used("fit"))
        {
            std::cerr << fit;
        }
        else
        {
            std::cerr << program;
        }
        std::exit(1);
    }

    if (program.is_subcommand_used("fit"))
    {
        gelex::logging::initialize(fit.get("--out"));
        auto logger = gelex::logging::get();

        auto bed = gelex::valid_bed(fit.get("--bfile"));

        if (!bed)
        {
            logger->error(bed.error().message);
            std::exit(1);
        }

        // Create shared SampleManager
        auto sample_manager_result = gelex::SampleManager::create(
            bed->replace_extension(".fam"), fit.get<bool>("--iid_only"));

        if (!sample_manager_result)
        {
            logger->error(sample_manager_result.error().message);
            std::exit(1);
        }

        auto sample_manager = std::make_shared<gelex::SampleManager>(
            std::move(*sample_manager_result));

        gelex::DataPipe::Config config{
            .phenotype_path = fit.get("--pheno"),
            .phenotype_column = fit.get<int>("--pheno-col"),
            .qcovar_path = fit.get("--qcovar"),
            .covar_path = fit.get("--covar"),
            .iid_only = fit.get<bool>("--iid_only"),
            .output_prefix = fit.get("--out")};

        auto data_pipe = gelex::DataPipe::create(config, sample_manager);

        if (!data_pipe)
        {
            logger->error(data_pipe.error().message);
            std::exit(1);
        }
        {
            auto genotype_pipe = gelex::GenotypePipe::create(
                bed->replace_extension(".bed"),
                sample_manager,
                fit.get("--out"));

            if (!genotype_pipe)
            {
                logger->error(genotype_pipe.error().message);
                std::exit(1);
            }

            if (auto process_result
                = genotype_pipe->process(fit.get<int>("--chunk-size"));
                !process_result)

            {
                logger->error(process_result.error().message);
                std::exit(1);
            }
        }

        if (fit.get<bool>("--dom"))
        {
            auto genotype_pipe = gelex::GenotypePipe::create(
                bed->replace_extension(".bed"),
                sample_manager,
                fit.get("--out"),
                true);

            if (!genotype_pipe)
            {
                logger->error(genotype_pipe.error().message);
                std::exit(1);
            }

            if (auto process_result
                = genotype_pipe->process(fit.get<int>("--chunk-size"));
                !process_result)

            {
                logger->error(process_result.error().message);
                std::exit(1);
            }
        }
    }
    else
    {
        if (argc <= 1)
        {
            std::cerr << program;
            std::exit(1);
        }
    }

    return 0;
}
