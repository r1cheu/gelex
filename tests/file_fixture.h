#ifndef GELEX_TEST_FILE_FIXTURE_H
#define GELEX_TEST_FILE_FIXTURE_H

#include <cstddef>
#include <filesystem>
#include <span>
#include <string_view>

namespace gelex::test
{

class FileFixture
{
   public:
    FileFixture();
    ~FileFixture();

    FileFixture(const FileFixture&) = delete;
    FileFixture& operator=(const FileFixture&) = delete;

    FileFixture(FileFixture&& other) noexcept;
    FileFixture& operator=(FileFixture&& other) noexcept;

    [[nodiscard]] std::filesystem::path create_text_file(
        std::string_view content,
        std::string_view suffix = "");

    [[nodiscard]] std::filesystem::path create_binary_file(
        std::span<const std::byte> content,
        std::string_view suffix = "");

    [[nodiscard]] std::filesystem::path create_named_text_file(
        std::string_view filename,
        std::string_view content);

    [[nodiscard]] std::filesystem::path create_named_binary_file(
        std::string_view filename,
        std::span<const std::byte> content);

    [[nodiscard]] std::filesystem::path create_empty_file(
        std::string_view suffix = "");
    [[nodiscard]] std::filesystem::path generate_random_file_path(
        std::string_view suffix = "");

    [[nodiscard]] const std::filesystem::path& get_test_dir() const noexcept;

   private:
    std::filesystem::path test_dir_;
    int file_counter_ = 0;

    std::filesystem::path next_path(std::string_view suffix);

    static void write_to_file(
        const std::filesystem::path& filepath,
        std::string_view content);
    static void write_to_file(
        const std::filesystem::path& filepath,
        std::span<const std::byte> content);
};

}  // namespace gelex::test

#endif  // GELEX_TEST_FILE_FIXTURE_H
