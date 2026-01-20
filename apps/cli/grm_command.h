#ifndef GELEX_CLI_GRM_COMMAND_H_
#define GELEX_CLI_GRM_COMMAND_H_

#include <argparse.h>

void grm_command(argparse::ArgumentParser& cmd);
auto grm_execute(argparse::ArgumentParser& cmd) -> int;

#endif  // GELEX_CLI_GRM_COMMAND_H_
