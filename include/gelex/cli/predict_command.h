#ifndef GELEX_CLI_PREDICT_COMMAND_H_
#define GELEX_CLI_PREDICT_COMMAND_H_

#include <argparse.h>

void predict_command(argparse::ArgumentParser& cmd);
int predict_execute(argparse::ArgumentParser& predict);

#endif  // GELEX_CLI_PREDICT_COMMAND_H_
