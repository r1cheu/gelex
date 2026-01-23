#ifndef GELEX_CLI_GRM_ARGS_H_
#define GELEX_CLI_GRM_ARGS_H_

#include <argparse.h>

struct GrmConfig
{
    std::filesystem::path bed_path;
    std::string out_prefix;
    std::string method;
    int chunk_size;
    bool do_additive;
    bool do_dominant;
    bool do_loco;
    int threads;
};

void setup_grm_args(argparse::ArgumentParser& cmd);

#endif  // GELEX_CLI_GRM_ARGS_H_
