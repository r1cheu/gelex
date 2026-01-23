#include <array>
#include <memory>
#include <string_view>

#include <argparse.h>

#include "config.h"

#include "cli/assoc_args.h"
#include "cli/assoc_command.h"
#include "cli/cli_helper.h"
#include "cli/fit_args.h"
#include "cli/fit_command.h"
#include "cli/grm_args.h"
#include "cli/grm_command.h"
#include "cli/predict_args.h"
#include "cli/predict_command.h"
#include "cli/simulation_args.h"
#include "cli/simulation_command.h"
#include "gelex/logger.h"

struct CommandDescriptor
{
    std::string_view name;
    void (*setup_fn)(argparse::ArgumentParser&);
    int (*execute_fn)(argparse::ArgumentParser&);
};

auto execute_command(
    argparse::ArgumentParser& parser,
    int (*execute_fn)(argparse::ArgumentParser&)) -> int
{
    constexpr std::string_view error_marker = "[\033[31merror\033[0m] ";
    try
    {
        gelex::logging::initialize(parser.get("--out"));
        return execute_fn(parser);
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

int main(int argc, char* argv[])
{
    constexpr std::array commands = {
        CommandDescriptor{"fit", setup_fit_args, fit_execute},
        CommandDescriptor{"simulate", setup_simulation_args, simulate_execute},
        CommandDescriptor{"predict", setup_predict_args, predict_execute},
        CommandDescriptor{"grm", setup_grm_args, grm_execute},
        CommandDescriptor{"assoc", setup_assoc_args, assoc_execute}};

    constexpr std::string_view error_marker = "[\033[31merror\033[0m] ";

    argparse::ArgumentParser program(PROJECT_NAME, PROJECT_VERSION);
    argparse::ArgumentParser fit("fit");
    setup_fit_args(fit);
    argparse::ArgumentParser simulate("simulate");
    setup_simulation_args(simulate);
    argparse::ArgumentParser predict("predict");
    setup_predict_args(predict);
    argparse::ArgumentParser grm("grm");
    setup_grm_args(grm);
    argparse::ArgumentParser assoc("assoc");
    setup_assoc_args(assoc);

    program.add_subparser(fit);
    program.add_subparser(simulate);
    program.add_subparser(predict);
    program.add_subparser(grm);
    program.add_subparser(assoc);

    std::array<argparse::ArgumentParser*, 5> parsers
        = {&fit, &simulate, &predict, &grm, &assoc};

    try
    {
        program.parse_args(argc, argv);
    }
    catch (const std::exception& err)
    {
        std::cerr << error_marker << err.what() << "\n";

        bool found = false;
        for (size_t i = 0; i < commands.size(); ++i)
        {
            if (program.is_subcommand_used(commands[i].name))
            {
                std::cerr << *parsers[i];
                found = true;
                break;
            }
        }
        if (!found)
        {
            std::cerr << program;
        }
        return 1;
    }

    for (size_t i = 0; i < commands.size(); ++i)
    {
        if (program.is_subcommand_used(commands[i].name))
        {
            return execute_command(*parsers[i], commands[i].execute_fn);
        }
    }

    if (argc <= 1)
    {
        gelex::cli::print_gelex_banner_message(PROJECT_VERSION);
        std::cerr << "\n" << program;
        return 1;
    }

    return 0;
}
