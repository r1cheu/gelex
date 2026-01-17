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