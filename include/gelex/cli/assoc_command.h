#ifndef GELEX_CLI_ASSOC_COMMAND_H
#define GELEX_CLI_ASSOC_COMMAND_H

#include <argparse.h>

auto assoc_command(argparse::ArgumentParser& cmd) -> void;
auto assoc_execute(argparse::ArgumentParser& cmd) -> int;

#endif  // GELEX_CLI_ASSOC_COMMAND_H
