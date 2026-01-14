#pragma once

#include <argparse.h>

void grm_command(argparse::ArgumentParser& cmd);
auto grm_execute(argparse::ArgumentParser& cmd) -> int;
