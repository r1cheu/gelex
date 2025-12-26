#include <argparse.h>

#include "config.h"

#include "gelex/cli/fit_command.h"
#include "gelex/cli/predict_command.h"
#include "gelex/cli/simulation_command.h"
#include "gelex/cli/utils.h"
#include "gelex/logger.h"

int main(int argc, char* argv[])
{
    argparse::ArgumentParser program(PROJECT_NAME, PROJECT_VERSION);
    argparse::ArgumentParser fit("fit");
    fit_command(fit);
    argparse::ArgumentParser simulate("simulate");
    simulate_command(simulate);
    argparse::ArgumentParser predict("predict");
    predict_command(predict);

    program.add_subparser(fit);
    program.add_subparser(simulate);
    program.add_subparser(predict);
    constexpr std::string_view error_marker = "[\033[31merror\033[0m] ";

    try
    {
        program.parse_args(argc, argv);
    }
    catch (const std::exception& err)
    {
        std::cerr << error_marker << err.what() << "\n";

        if (program.is_subcommand_used("fit"))
        {
            std::cerr << fit;
        }
        else if (program.is_subcommand_used("simulate"))
        {
            std::cerr << simulate;
        }
        else if (program.is_subcommand_used("predict"))
        {
            std::cerr << predict;
        }
        else
        {
            std::cerr << program;
        }
        return 1;
    }

    if (program.is_subcommand_used("fit"))
    {
        try
        {
            gelex::cli::log_command(fit, gelex::cli::parse_command(argc, argv));
            int code = fit_execute(fit);
            return code;
        }
        catch (const std::exception& e)
        {
            auto logger = gelex::logging::get();
            if (logger)
            {
                logger->error("{}", e.what());
            }
            else
            {
                std::cerr << error_marker << e.what() << "\n";
            }
            return 1;
        }
    }

    if (program.is_subcommand_used("simulate"))
    {
        try
        {
            gelex::cli::log_command(
                simulate, gelex::cli::parse_command(argc, argv));
            int code = simulate_execute(simulate);
            return code;
        }
        catch (const std::exception& e)
        {
            auto logger = gelex::logging::get();
            if (logger)
            {
                logger->error("{}", e.what());
            }
            else
            {
                std::cerr << error_marker << e.what() << "\n";
            }
            return 1;
        }
    }

    if (program.is_subcommand_used("predict"))
    {
        try
        {
            gelex::cli::log_command(
                predict, gelex::cli::parse_command(argc, argv));
            int code = predict_execute(predict);
            return code;
        }
        catch (const std::exception& e)
        {
            auto logger = gelex::logging::get();
            if (logger)
            {
                logger->error("{}", e.what());
            }
            else
            {
                std::cerr << error_marker << e.what() << "\n";
            }
            return 1;
        }
    }

    if (argc <= 1)
    {
        gelex::cli::print_banner_message(PROJECT_VERSION);
        std::cerr << "\n" << program;
        return 1;
    }

    return 0;
}
