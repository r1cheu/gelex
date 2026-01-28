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

#ifndef TESTS_TEST_INIT_H_
#define TESTS_TEST_INIT_H_

#include "gelex/logger.h"
#include <filesystem>

namespace test_utils
{

struct TestInitializer
{
    TestInitializer()
    {
        // Initialize logger with test output
        gelex::logging::initialize("test_output");

        clean_test_files();
    }

    ~TestInitializer()
    {
        clean_test_files();
    }

    static void clean_test_files()
    {
        try
        {
            std::filesystem::remove("test_output.log");
        }
        catch (const std::filesystem::filesystem_error&)
        {
        }
    }
};

// Global test initializer - will be constructed before any tests run
inline TestInitializer global_test_initializer;

} // namespace test_utils

#endif  // TESTS_TEST_INIT_H_