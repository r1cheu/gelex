// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef TESTS_TEST_INIT_H_
#define TESTS_TEST_INIT_H_

#include <filesystem>
#include <fmt/format.h>
#include <unistd.h>

namespace test_utils
{

struct TestInitializer
{
    TestInitializer() { clean_test_files(); }

    ~TestInitializer() { clean_test_files(); }

    static auto log_prefix() -> std::string
    {
        return fmt::format("test_output_{}", ::getpid());
    }

    static void clean_test_files()
    {
        try
        {
            std::filesystem::remove(log_prefix() + ".log");
        }
        catch (const std::filesystem::filesystem_error&)
        {
        }
    }
};

// Global test initializer - will be constructed before any tests run
inline TestInitializer global_test_initializer;

}  // namespace test_utils

#endif  // TESTS_TEST_INIT_H_
