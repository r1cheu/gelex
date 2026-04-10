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

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <fmt/format.h>

#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/reader.h"

struct BimStruct
{
    std::vector<std::string> chrom;
    std::vector<std::int32_t> pos;
    std::vector<std::string> a1;
    std::vector<std::string> a2;
};

auto write_fake_bim(const std::filesystem::path& path, std::size_t n) -> void
{
    std::ofstream ofs(path);
    for (std::size_t i = 0; i < n; ++i)
    {
        // chr  snp_id  cm  bp  A1  A2
        fmt::format_to(
            std::ostreambuf_iterator(ofs),
            "{}\trs{}\t0\t{}\tA\tG\n",
            (i % 22) + 1,
            i,
            i * 1000);
    }
}

auto to_struct(const gelex::df::DataFrame<std::string>& df) -> BimStruct
{
    auto n = df.rows();
    BimStruct s;
    s.chrom.reserve(n);
    s.pos.reserve(n);
    s.a1.reserve(n);
    s.a2.reserve(n);

    auto chrom = df["chrom"].as<std::string>();
    auto pos = df["pos"].as<std::int32_t>();
    auto a1 = df["A1"].as<std::string>();
    auto a2 = df["A2"].as<std::string>();
    for (std::size_t i = 0; i < n; ++i)
    {
        s.chrom.push_back(std::string(chrom[i]));
        s.pos.push_back(pos[i]);
        s.a1.push_back(std::string(a1[i]));
        s.a2.push_back(std::string(a2[i]));
    }
    return s;
}

int main()
{
    constexpr std::size_t kNumSnps = 500'000;

    auto tmp = std::filesystem::temp_directory_path() / "bench_bim.bim";
    write_fake_bim(tmp, kNumSnps);
    auto df = gelex::read_bim(tmp);
    auto bim = to_struct(df);
    std::filesystem::remove(tmp);

    ankerl::nanobench::Bench b;
    b.warmup(3).minEpochIterations(3);

    // 1) current pattern: hash lookup + variant dispatch per row per column
    b.run(
        "df_lookup_per_row",
        [&]()
        {
            std::int64_t sink = 0;
            for (std::size_t i = 0; i < kNumSnps; ++i)
            {
                ankerl::nanobench::doNotOptimizeAway(
                    df["chrom"].as<std::string>()[i]);
                sink += df["pos"].as<std::int32_t>()[i];
                ankerl::nanobench::doNotOptimizeAway(
                    df["A1"].as<std::string>()[i]);
                ankerl::nanobench::doNotOptimizeAway(
                    df["A2"].as<std::string>()[i]);
            }
            ankerl::nanobench::doNotOptimizeAway(sink);
        });

    // 2) hoisted: lookup once, loop over spans
    b.run(
        "df_span_hoisted",
        [&]()
        {
            auto chrom = df["chrom"].as<std::string>();
            auto pos = df["pos"].as<std::int32_t>();
            auto a1 = df["A1"].as<std::string>();
            auto a2 = df["A2"].as<std::string>();
            std::int64_t sink = 0;
            for (std::size_t i = 0; i < kNumSnps; ++i)
            {
                ankerl::nanobench::doNotOptimizeAway(chrom[i]);
                sink += pos[i];
                ankerl::nanobench::doNotOptimizeAway(a1[i]);
                ankerl::nanobench::doNotOptimizeAway(a2[i]);
            }
            ankerl::nanobench::doNotOptimizeAway(sink);
        });

    // 3) typed struct: direct vector access
    b.run(
        "struct_direct",
        [&]()
        {
            std::int64_t sink = 0;
            for (std::size_t i = 0; i < kNumSnps; ++i)
            {
                ankerl::nanobench::doNotOptimizeAway(bim.chrom[i]);
                sink += bim.pos[i];
                ankerl::nanobench::doNotOptimizeAway(bim.a1[i]);
                ankerl::nanobench::doNotOptimizeAway(bim.a2[i]);
            }
            ankerl::nanobench::doNotOptimizeAway(sink);
        });

    return 0;
}
