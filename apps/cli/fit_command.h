#ifndef GELEX_CLI_FIT_COMMAND_H_
#define GELEX_CLI_FIT_COMMAND_H_

#include <argparse.h>

void fit_command(argparse::ArgumentParser& cmd);
int fit_execute(argparse::ArgumentParser& fit);

#endif  // GELEX_CLI_FIT_COMMAND_H_
