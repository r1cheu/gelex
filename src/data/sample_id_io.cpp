// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/data/sample_id_io.h"

#include <fmt/format.h>
#include <span>
#include <string>

#include "gelex/data/dataframe/index.h"
#include "gelex/data/dataframe/reader.h"
#include "gelex/data/sample_id.h"
#include "gelex/exception.h"
#include "gelex/io/detail/text_writer.h"

namespace gelex
{

auto write_sample_ids(
    const std::string& prefix,
    std::span<const std::string> ids) -> void
{
    detail::TextWriter writer(prefix + ".id");
    for (const auto& id : ids)
    {
        const auto [fid, iid] = split_sample_id(id);
        writer.write(fmt::format("{}\t{}", fid, iid));
    }
}

auto read_sample_ids(const std::string& prefix) -> DataFrameIndex<std::string>
{
    const std::string path = prefix + ".id";
    ReadOptions options;
    options.delimiter = '\t';
    options.header = false;
    options.index_cols = {0, 1};
    auto index = read_index<std::string>(path, options);
    if (index.size() == 0)
    {
        throw GelexException(fmt::format("{}: no sample IDs found", path));
    }
    return index;
}

}  // namespace gelex
