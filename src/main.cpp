#include <gelex/argparse.h>

#include "config.h"

#include "gelex/cli/fit_command.h"
#include "gelex/cli/simulation_command.h"
#include "gelex/cli/utils.h"

int main(int argc, char* argv[])
{
    argparse::ArgumentParser program(PROJECT_NAME, PROJECT_VERSION);
    argparse::ArgumentParser fit("fit");
    fit_command(fit);
    argparse::ArgumentParser simulate("simulate");
    simulate_command(simulate);

    program.add_subparser(fit);
    program.add_subparser(simulate);
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
        else
        {
            std::cerr << program;
        }
        return 1;
    }

    if (program.is_subcommand_used("fit"))
    {
        gelex::cli::log_command(fit, gelex::cli::parse_command(argc, argv));
        int code = fit_execute(fit);
        return code;
    }

    if (program.is_subcommand_used("simulate"))
    {
        gelex::cli::log_command(
            simulate, gelex::cli::parse_command(argc, argv));
        int code = simulate_execute(simulate);
        return code;
    }

    if (argc <= 1)
    {
        gelex::cli::print_banner_message(PROJECT_VERSION);
        std::cerr << "\n" << program;
        return 1;
    }

    return 0;
}
