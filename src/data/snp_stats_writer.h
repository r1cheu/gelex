#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <span>
#include <string_view>
#include <vector>

namespace gelex::detail
{

class SnpStatsWriter
{
   public:
    static constexpr auto kDefaultBufferSize
        = static_cast<const size_t>(64 * 1024);

    explicit SnpStatsWriter(const std::filesystem::path& file_path);

    SnpStatsWriter(const SnpStatsWriter&) = delete;
    SnpStatsWriter(SnpStatsWriter&&) noexcept = default;
    SnpStatsWriter& operator=(const SnpStatsWriter&) = delete;
    SnpStatsWriter& operator=(SnpStatsWriter&&) noexcept = default;
    ~SnpStatsWriter() = default;

    void write(
        int64_t num_samples,
        std::span<const int64_t> monomorphic_indices,
        std::span<const double> means,
        std::span<const double> stddevs);

    [[nodiscard]] auto path() const noexcept -> const std::filesystem::path&
    {
        return path_;
    }

   private:
    void write_data(
        const void* data,
        std::streamsize size,
        std::string_view error_msg);
    template <typename T>
    void write_data(std::span<const T> data, std::string_view error_msg)
    {
        write_data(
            data.data(),
            static_cast<std::streamsize>(data.size_bytes()),
            error_msg);
    }

    static void check_monomorphic_indices(
        std::span<const int64_t> monomorphic_indices,
        int64_t num_variants);

    std::filesystem::path path_;
    std::vector<char> io_buffer_;
    std::ofstream file_;
};

}  // namespace gelex::detail
