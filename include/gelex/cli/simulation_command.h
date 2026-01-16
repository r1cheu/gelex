#ifndef GELEX_CLI_SIMULATE_COMMAND_H_
#define GELEX_CLI_SIMULATE_COMMAND_H_
#include <argparse.h>

void simulate_command(argparse::ArgumentParser& cmd);
int simulate_execute(argparse::ArgumentParser& sim);

#endif  // GELEX_CLI_SIMULATE_COMMAND_H_
