#include "file_fixture.h"

#include <ctime>
#include <format>
#include <fstream>
#include <ios>
#include <utility>

#include "gelex/exception.h"

namespace fs = std::filesystem;

namespace gelex::test
{

FileFixture::FileFixture()
{
    auto dir_name = std::format("gelex_test_{}", std::time(nullptr));
    test_dir_ = fs::temp_directory_path() / dir_name;

    if (fs::exists(test_dir_))
    {
        fs::remove_all(test_dir_);
    }
    fs::create_directories(test_dir_);
}

FileFixture::~FileFixture()
{
    if (!test_dir_.empty())
    {
        std::error_code ec;
        fs::remove_all(test_dir_, ec);
    }
}

FileFixture::FileFixture(FileFixture&& other) noexcept
    : test_dir_(std::move(other.test_dir_)), file_counter_(other.file_counter_)
{
    other.test_dir_.clear();
    other.file_counter_ = 0;
}

FileFixture& FileFixture::operator=(FileFixture&& other) noexcept
{
    if (this != &other)
    {
        if (!test_dir_.empty())
        {
            std::error_code ec;
            fs::remove_all(test_dir_, ec);
        }
        test_dir_ = std::move(other.test_dir_);
        file_counter_ = other.file_counter_;
        other.test_dir_.clear();
        other.file_counter_ = 0;
    }
    return *this;
}

fs::path FileFixture::create_text_file(
    std::string_view content,
    std::string_view suffix)
{
    auto path = next_path(suffix);
    write_to_file(path, content);
    return path;
}

fs::path FileFixture::create_binary_file(
    std::span<const std::byte> content,
    std::string_view suffix)
{
    auto path = next_path(suffix);
    write_to_file(path, content);
    return path;
}

fs::path FileFixture::create_named_text_file(
    std::string_view filename,
    std::string_view content)
{
    auto path = test_dir_ / filename;
    write_to_file(path, content);
    return path;
}

fs::path FileFixture::create_named_binary_file(
    std::string_view filename,
    std::span<const std::byte> content)
{
    auto path = test_dir_ / filename;
    write_to_file(path, content);
    return path;
}

fs::path FileFixture::create_empty_file(std::string_view suffix)
{
    return create_text_file("", suffix);
}

const fs::path& FileFixture::get_test_dir() const noexcept
{
    return test_dir_;
}

fs::path FileFixture::next_path(std::string_view suffix)
{
    return test_dir_ / std::format("test_{}{}", file_counter_++, suffix);
}

void FileFixture::write_to_file(
    const fs::path& filepath,
    std::string_view content)
{
    if (auto parent = filepath.parent_path(); !parent.empty())
    {
        fs::create_directories(parent);
    }

    std::ofstream file(filepath);
    if (!file)
    {
        throw FileWriteException(
            std::format("{}:Failed to create file", filepath.string()));
    }
    file << content;
}

void FileFixture::write_to_file(
    const fs::path& filepath,
    std::span<const std::byte> content)
{
    if (auto parent = filepath.parent_path(); !parent.empty())
    {
        fs::create_directories(parent);
    }

    std::ofstream file(filepath, std::ios::binary);
    if (!file)
    {
        throw FileWriteException(
            std::format("{}:Failed to create binary file", filepath.string()));
    }
    file.write(
        reinterpret_cast<const char*>(content.data()),  // NOLINT
        static_cast<std::streamsize>(content.size()));
}

}  // namespace gelex::test
