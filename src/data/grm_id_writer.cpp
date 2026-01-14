#include "grm_id_writer.h"

#include <format>

#include "gelex/exception.h"
#include "parser.h"

namespace gelex::detail
{

GrmIdWriter::GrmIdWriter(const std::filesystem::path& file_path)
    : path_(file_path)
{
    file_ = detail::open_file<std::ofstream>(
        path_, std::ios::out | std::ios::trunc);
}

auto GrmIdWriter::split_id(std::string_view id)
    -> std::pair<std::string_view, std::string_view>
{
    auto pos = id.find('_');
    if (pos == std::string_view::npos)
    {
        // No '_' found, use the same value for both FID and IID
        return {id, id};
    }
    return {id.substr(0, pos), id.substr(pos + 1)};
}

auto GrmIdWriter::write(std::span<const std::string> ids) -> void
{
    for (const auto& id : ids)
    {
        auto [fid, iid] = split_id(id);
        file_ << fid << '\t' << iid << '\n';
    }

    if (!file_.good())
    {
        throw FileWriteException(
            std::format("{}: failed to write ID data", path_.string()));
    }
}

}  // namespace gelex::detail
