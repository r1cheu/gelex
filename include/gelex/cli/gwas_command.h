#ifndef GELEX_CLI_GWAS_COMMAND_H
#define GELEX_CLI_GWAS_COMMAND_H

#include <argparse.h>

auto gwas_command(argparse::ArgumentParser& cmd) -> void;
auto gwas_execute(argparse::ArgumentParser& cmd) -> int;

#endif  // GELEX_CLI_GWAS_COMMAND_H
