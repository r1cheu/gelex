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

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <fstream>
#include <ios>
#include <sstream>
#include <stdexcept>
#include <string>

#include "gelex/exception.h"
#include "gelex/io/detail/atomic_output_stream.h"

#include "file_fixture.h"

namespace fs = std::filesystem;

using gelex::GelexException;
using gelex::detail::AtomicOutputStream;
using gelex::test::FileFixture;

namespace
{

auto read_file(const fs::path& p) -> std::string
{
    std::ifstream ifs(p);
    std::ostringstream oss;
    oss << ifs.rdbuf();
    return oss.str();
}

auto tmp_of(const fs::path& final_path) -> fs::path
{
    fs::path t = final_path;
    t += ".tmp";
    return t;
}

}  // namespace

TEST_CASE(
    "AtomicOutputStream - destructor discards uncommitted tmp file",
    "[atomic_output_stream][tmp]")
{
    FileFixture files;
    auto final_path = files.generate_random_file_path(".txt");

    {
        AtomicOutputStream ofs(final_path, std::ios::out);
        ofs.write("hello");
    }

    REQUIRE_FALSE(fs::exists(final_path));
    REQUIRE_FALSE(fs::exists(tmp_of(final_path)));
}

TEST_CASE(
    "AtomicOutputStream - tmp file exists before destruction",
    "[atomic_output_stream][tmp]")
{
    FileFixture files;
    auto final_path = files.generate_random_file_path(".txt");
    const auto tmp = tmp_of(final_path);

    AtomicOutputStream ofs(final_path, std::ios::out);
    ofs.write("partial");

    REQUIRE(fs::exists(tmp));
    REQUIRE_FALSE(fs::exists(final_path));
}

TEST_CASE(
    "AtomicOutputStream - explicit commit makes final file visible",
    "[atomic_output_stream][commit]")
{
    FileFixture files;
    auto final_path = files.generate_random_file_path(".txt");

    AtomicOutputStream ofs(final_path, std::ios::out);
    ofs.write("hello");
    ofs.commit();

    REQUIRE(fs::exists(final_path));
    REQUIRE_FALSE(fs::exists(tmp_of(final_path)));
    REQUIRE(read_file(final_path) == "hello");
}

TEST_CASE(
    "AtomicOutputStream - in-flight exception abandons tmp file",
    "[atomic_output_stream][exception]")
{
    FileFixture files;
    auto final_path = files.generate_random_file_path(".txt");
    const auto tmp = tmp_of(final_path);

    try
    {
        AtomicOutputStream ofs(final_path, std::ios::out);
        ofs.write("will be discarded");
        throw std::runtime_error("simulated failure");
    }
    catch (const std::runtime_error&)  // NOLINT(bugprone-empty-catch)
    {
    }

    REQUIRE_FALSE(fs::exists(final_path));
    REQUIRE_FALSE(fs::exists(tmp));
}

TEST_CASE(
    "AtomicOutputStream - commit atomically replaces existing file",
    "[atomic_output_stream][overwrite]")
{
    FileFixture files;
    auto final_path
        = files.create_text_file("old content", ".txt");  // pre-existing
    REQUIRE(read_file(final_path) == "old content");

    {
        AtomicOutputStream ofs(final_path, std::ios::out);
        ofs.write("new content");
        // pre-existing final file remains intact until rename.
        REQUIRE(read_file(final_path) == "old content");
        REQUIRE(fs::exists(tmp_of(final_path)));
        ofs.commit();
    }

    REQUIRE(fs::exists(final_path));
    REQUIRE_FALSE(fs::exists(tmp_of(final_path)));
    REQUIRE(read_file(final_path) == "new content");
}

TEST_CASE(
    "AtomicOutputStream - is_directory path throws",
    "[atomic_output_stream][error]")
{
    FileFixture files;
    const auto& dir = files.get_test_dir();
    REQUIRE_THROWS_AS(AtomicOutputStream(dir, std::ios::out), GelexException);
}

TEST_CASE(
    "AtomicOutputStream - opening in non-existent directory throws",
    "[atomic_output_stream][error]")
{
    FileFixture files;
    auto bad_path = files.get_test_dir() / "no_such_dir" / "file.txt";
    REQUIRE_THROWS_AS(
        AtomicOutputStream(bad_path, std::ios::out), GelexException);
    REQUIRE_FALSE(fs::exists(bad_path));
    REQUIRE_FALSE(fs::exists(tmp_of(bad_path)));
}

TEST_CASE(
    "AtomicOutputStream - binary mode roundtrip",
    "[atomic_output_stream][binary]")
{
    FileFixture files;
    auto final_path = files.generate_random_file_path(".bin");

    constexpr std::array<char, 4> data{'\x00', '\x01', '\x02', '\x03'};

    {
        AtomicOutputStream ofs(
            final_path, std::ios::out | std::ios::binary | std::ios::trunc);
        ofs.write(data.data(), static_cast<std::streamsize>(data.size()));
        ofs.commit();
    }

    REQUIRE(fs::exists(final_path));
    REQUIRE(fs::file_size(final_path) == data.size());

    std::ifstream ifs(final_path, std::ios::binary);
    std::array<char, 4> read_back{};
    ifs.read(read_back.data(), static_cast<std::streamsize>(read_back.size()));
    REQUIRE(read_back == data);
}
