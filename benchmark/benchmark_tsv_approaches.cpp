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

#include <benchmark/benchmark.h>

#include <format>
#include <fstream>

#include <fmt/compile.h>
#include <fmt/format.h>

// Approach A: std::format + stream << (current batch writers pattern)
static void BM_StdFormatStream(benchmark::State& state)
{
    std::ofstream ofs("/tmp/bench_tsv_a.tsv", std::ios::out | std::ios::trunc);
    for (auto _ : state)
    {
        ofs << std::format(
            "{}\t{}\t{:.6f}\t{:.6f}\t{:.6e}\n",
            "rs1234567",
            1,
            0.012345,
            0.006789,
            1.23e-8);
    }
}

// Approach B: fmt::format_to + fmt::memory_buffer + 64KB threshold flush
// (GwasWriter pattern)
static void BM_FmtMemoryBuffer(benchmark::State& state)
{
    std::ofstream ofs(
        "/tmp/bench_tsv_b.tsv",
        std::ios::out | std::ios::binary | std::ios::trunc);
    fmt::memory_buffer buf;
    buf.reserve(64 * 1024);
    for (auto _ : state)
    {
        fmt::format_to(
            std::back_inserter(buf),
            FMT_COMPILE("{}\t{}\t{:.6f}\t{:.6f}\t{:.6e}\n"),
            "rs1234567",
            1,
            0.012345,
            0.006789,
            1.23e-8);
        if (buf.size() >= 64 * 1024)
        {
            ofs.write(buf.data(), static_cast<std::streamsize>(buf.size()));
            buf.clear();
        }
    }
    if (buf.size() > 0)
        ofs.write(buf.data(), static_cast<std::streamsize>(buf.size()));
}

BENCHMARK(BM_StdFormatStream)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_FmtMemoryBuffer)->Unit(benchmark::kNanosecond);

BENCHMARK_MAIN();
