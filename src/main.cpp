#include <argparse/argparse.hpp>

#include "Eigen/Core"

#include "gelex/data/io.h"
#include "gelex/model/bayes/model.h"

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
        auto data = gelex::DataReader::Create(
            fit.get("--pheno"),
            fit.get("--bfile") + ".fam",
            fit.get("--qcovar"),
            fit.get("--covar"),
            fit.get<int>("--pheno-col"),
            fit.get<bool>("--iid_only"));
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
