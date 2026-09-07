// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "cli/runtime.h"

#include <Eigen/Core>
#include <cstdio>
#include <omp.h>

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

extern "C"
{
    void openblas_set_num_threads(int) __attribute__((weak));
    void MKL_Set_Num_Threads(int) __attribute__((weak));
}

namespace cli
{

auto is_tty() -> bool
{
#ifdef _WIN32
    return _isatty(_fileno(stdout)) != 0;
#else
    return isatty(fileno(stdout)) != 0;
#endif
}

auto setup_parallelization(int num_threads) -> void
{
    if (num_threads <= 0)
    {
        return;
    }
    omp_set_num_threads(num_threads);
    Eigen::setNbThreads(num_threads);
    if (openblas_set_num_threads != nullptr)
    {
        openblas_set_num_threads(num_threads);
    }
    if (MKL_Set_Num_Threads != nullptr)
    {
        MKL_Set_Num_Threads(num_threads);
    }
}

}  // namespace cli
