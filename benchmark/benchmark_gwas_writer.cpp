/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <nanobench.h>

#include <array>
#include <cstdio>
#include <fstream>
#include <string>

#include "gelex/algo/gwas/assoc_tester.h"
#include "gelex/data/reader.h"
#include "gelex/io/gwas_writer.h"

#include <catch2/catch_test_macros.hpp>

using namespace gelex;

TEST_CASE("GWAS writer output", "[!benchmark][io][gwas]")
{
    std::string temp_filename = "/tmp/bench_test_ignore";
    std::string bim_file = "/tmp/bench_test_ignore.bim";
    {
        std::ofstream ofs(bim_file);
        ofs << "1\trs123456\t0\t100000\tA\tG\n";
    }
    auto bim = gelex::read_bim(bim_file);

    std::array freq{0.25};
    std::array beta{0.0123};
    std::array se{0.0045};
    std::array p{1.23e-8};
    std::array pve{0.01};
    gelex::TestResults results{
        .freq = freq,
        .additive = {.beta = beta, .se = se, .p = p, .pve = pve},
    };

    GwasWriter writer(temp_filename, bim);

    ankerl::nanobench::Bench().run(
        "GwasWriter WriteResult", [&]() { writer.write(0, results); });

    std::remove((temp_filename + ".gwas.tsv").c_str());
    std::remove(bim_file.c_str());
}
