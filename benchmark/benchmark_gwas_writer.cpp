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

#include <cstdio>
#include <fstream>
#include <string>

#include "gelex/data/reader.h"
#include "gelex/io/gwas_writer.h"

using namespace gelex::gwas;

int main()
{
    // 写往 /tmp 避免磁盘 I/O 波动影响格式化性能测试
    std::string temp_filename = "/tmp/bench_test_ignore";
    std::string bim_file = "/tmp/bench_test_ignore.bim";
    {
        std::ofstream ofs(bim_file);
        ofs << "1\trs123456\t0\t100000\tA\tG\n";
    }
    auto bim = gelex::read_bim(bim_file);

    GwasWriter::AssocResult result;
    result.beta = 0.0123;
    result.se = 0.0045;
    result.p_value = 1.23e-8;

    GwasWriter writer(temp_filename, bim);

    ankerl::nanobench::Bench().run(
        "GwasWriter WriteResult", [&]() { writer.write(0, result); });

    std::remove((temp_filename + ".gwas.tsv").c_str());
    std::remove(bim_file.c_str());
}
