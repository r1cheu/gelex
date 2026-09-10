// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/bayes/diagnostics.h"

#include <filesystem>
#include <fmt/format.h>
#include <span>
#include <string_view>

#include "gelex/bayes/genetic/diagnostics_traits.h"
#include "gelex/io/detail/atomic_output_stream.h"

namespace gelex
{

auto write_diagnostics(
    std::string_view path,
    std::span<const DiagnosticEntry> entries) -> void
{
    detail::AtomicOutputStream file{std::filesystem::path{path}};
    file.write(
        "id\tindex\tmean\tsd\tmedian\thpdi_lower\thpdi_upper\tess\tmcse\t"
        "split_rhat\n");
    for (const auto& entry : entries)
    {
        const auto& s = entry.stats;
        file.write(
            fmt::format(
                "{}\t{}\t{:.10g}\t{:.10g}\t{:.10g}\t{:.10g}\t{:.10g}\t{:.10g}"
                "\t{:.10g}\t{:.10g}\n",
                entry.id,
                entry.index,
                s.mean,
                s.sd,
                s.median,
                s.hpdi_lower,
                s.hpdi_upper,
                s.ess,
                s.mcse,
                s.split_rhat));
    }
    file.commit();
}

}  // namespace gelex
