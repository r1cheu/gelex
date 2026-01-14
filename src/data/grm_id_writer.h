#ifndef GELEX_DATA_GRM_ID_WRITER_H
#define GELEX_DATA_GRM_ID_WRITER_H

#include <filesystem>
#include <fstream>
#include <span>
#include <string>
#include <string_view>
#include <utility>

namespace gelex::detail
{

class GrmIdWriter
{
   public:
    explicit GrmIdWriter(const std::filesystem::path& file_path);

    GrmIdWriter(const GrmIdWriter&) = delete;
    GrmIdWriter(GrmIdWriter&&) noexcept = default;
    GrmIdWriter& operator=(const GrmIdWriter&) = delete;
    GrmIdWriter& operator=(GrmIdWriter&&) noexcept = default;
    ~GrmIdWriter() = default;

    // Write sample IDs, each line contains: FID\tIID
    // Input ids are in "FID_IID" format, will be split by first '_'
    auto write(std::span<const std::string> ids) -> void;

    [[nodiscard]] auto path() const noexcept -> const std::filesystem::path&
    {
        return path_;
    }

   private:
    std::filesystem::path path_;
    std::ofstream file_;

    // Split "FID_IID" into (FID, IID) by the first '_'
    // If no '_' found, both FID and IID are set to the original id
    static auto split_id(std::string_view id)
        -> std::pair<std::string_view, std::string_view>;
};

}  // namespace gelex::detail

#endif  // GELEX_DATA_GRM_ID_WRITER_H
