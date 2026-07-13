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

#include "color_formatter.h"

#include <CLI/CLI.hpp>
#include <algorithm>
#include <fmt/format.h>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "theme.h"

namespace cli
{

namespace
{

auto section_heading(std::string text) -> std::string
{
    return colorize(ColorRole::heading, std::move(text));
}

auto option_name_text(std::string text) -> std::string
{
    return colorize(ColorRole::option_name, std::move(text));
}

auto option_value_text(std::string text) -> std::string
{
    return colorize(ColorRole::option_value, std::move(text));
}

class ColorFormatter : public CLI::Formatter
{
   public:
    auto make_help(
        const CLI::App* app,
        std::string name,
        CLI::AppFormatMode mode) const -> std::string override
    {
        if (mode == CLI::AppFormatMode::Sub)
        {
            return CLI::Formatter::make_help(app, std::move(name), mode);
        }

        std::ostringstream out;
        out << '\n';
        auto description = make_description(app);
        if (!description.empty())
        {
            while (!description.empty() && description.back() == '\n')
            {
                description.pop_back();
            }
            out << description << "\n\n";
        }
        out << make_usage(app, std::move(name));
        out << make_positionals(app);
        out << make_subcommands(app, mode);
        out << make_groups(app, mode);
        out << make_footer(app);
        return out.str();
    }

    auto make_group(
        std::string group,
        bool is_positional,
        std::vector<const CLI::Option*> options) const -> std::string override
    {
        if (group == "OPTIONS")
        {
            group = "Options";
        }
        else if (group == "POSITIONALS")
        {
            group = "Positionals";
        }

        std::ostringstream out;
        out << "\n" << section_heading(std::move(group) + ":") << "\n";
        for (const CLI::Option* option : options)
        {
            out << make_option(option, is_positional);
        }
        return out.str();
    }

    auto make_groups(const CLI::App* app, CLI::AppFormatMode mode) const
        -> std::string
    {
        std::ostringstream out;
        const auto groups = app->get_groups();

        for (const auto& group : groups)
        {
            if (group == "OPTIONS")
            {
                continue;
            }
            out << make_group_options(app, mode, group);
        }

        for (const auto& group : groups)
        {
            if (group != "OPTIONS")
            {
                continue;
            }
            out << make_group_options(app, mode, group);
        }

        return out.str();
    }

    auto make_subcommands(const CLI::App* app, CLI::AppFormatMode mode) const
        -> std::string override
    {
        auto subcommands = app->get_subcommands(
            [](const CLI::App* subcommand)
            { return should_show_subcommand(subcommand); });
        if (subcommands.empty())
        {
            return {};
        }

        std::ostringstream out;
        out << "\n" << section_heading("Commands:") << "\n";
        for (const CLI::App* subcommand : subcommands)
        {
            if (mode == CLI::AppFormatMode::All)
            {
                out << subcommand->help(subcommand->get_name(), mode) << "\n";
            }
            else
            {
                out << make_subcommand(subcommand);
            }
        }
        return out.str();
    }

    auto make_usage(const CLI::App* app, std::string name) const
        -> std::string override
    {
        std::ostringstream out;
        out << section_heading("Usage:") << "\n  ";
        if (app->get_parent() != nullptr
            && (name.empty() || name == app->get_name()))
        {
            std::vector<std::string> command_path;
            for (const auto* current = app; current != nullptr;
                 current = current->get_parent())
            {
                if (!current->get_name().empty())
                {
                    command_path.push_back(current->get_name());
                }
            }
            std::ranges::reverse(command_path);
            out << CLI::detail::join(command_path, " ");
        }
        else
        {
            out << (name.empty() ? app->get_name() : name);
        }

        auto required_options = app->get_options(
            [](const CLI::Option* option)
            { return option->nonpositional() && option->get_required(); });
        for (const CLI::Option* option : required_options)
        {
            const auto names = CLI::detail::split(
                CLI::Formatter::make_option_name(option, false), ',');
            auto option_name = names.front();
            for (const auto& name : names)
            {
                if (name.find("--", 0) == std::string::npos)
                {
                    option_name = name;
                    break;
                }
            }
            out << " " << option_name_text(option_name)
                << make_option_opts(option);
        }

        auto optional_options = app->get_options(
            [](const CLI::Option* option)
            { return option->nonpositional() && !option->get_required(); });
        if (!optional_options.empty())
        {
            out << " " << option_name_text("[OPTIONS]");
        }

        auto positionals = app->get_options(
            [](const CLI::Option* option) { return option->get_positional(); });
        for (const CLI::Option* positional : positionals)
        {
            out << " " << make_option_usage(positional);
        }

        auto subcommands = app->get_subcommands(
            [](const CLI::App* subcommand)
            { return should_show_subcommand(subcommand); });
        if (!subcommands.empty())
        {
            out << " " << option_name_text("[COMMAND]");
        }

        out << "\n";
        return out.str();
    }

    auto make_option_name(const CLI::Option* option, bool is_positional) const
        -> std::string override
    {
        return option_name_text(
            CLI::Formatter::make_option_name(option, is_positional));
    }

    auto make_option(const CLI::Option* option, bool is_positional) const
        -> std::string override
    {
        if (is_positional)
        {
            return make_positional_option(option);
        }

        return make_nonpositional_option(option);
    }

    auto make_option_opts(const CLI::Option* option) const
        -> std::string override
    {
        return option_value_text(make_option_suffix(option));
    }

    auto make_option_desc(const CLI::Option* option) const
        -> std::string override
    {
        auto desc = CLI::Formatter::make_option_desc(option);
        if (option->get_default_str().empty())
        {
            return desc;
        }
        if (!desc.empty())
        {
            desc += " ";
        }
        desc += make_default_text(option);
        return desc;
    }

    auto make_subcommand(const CLI::App* subcommand) const
        -> std::string override
    {
        auto name = "  " + subcommand->get_display_name(true);
        if (subcommand->get_required())
        {
            name += " " + get_label("REQUIRED");
        }
        std::ostringstream out;
        out << option_name_text(name);
        if (name.length() < column_width_)
        {
            out << std::string(column_width_ - name.length(), ' ');
        }
        CLI::detail::streamOutAsParagraph(
            out,
            subcommand->get_description(),
            right_column_width_,
            std::string(column_width_, ' '),
            true);
        out << '\n';
        return out.str();
    }

    auto make_footer(const CLI::App* app) const -> std::string override
    {
        if (app->get_footer().empty())
        {
            return {};
        }

        std::istringstream lines(app->get_footer());
        std::ostringstream out;
        out << '\n';
        std::string line;
        while (std::getline(lines, line))
        {
            if (line.ends_with(':'))
            {
                out << section_heading(std::move(line));
            }
            else
            {
                out << line;
            }
            out << '\n';
        }
        return out.str();
    }

   private:
    auto make_group_options(
        const CLI::App* app,
        CLI::AppFormatMode mode,
        const std::string& group) const -> std::string
    {
        auto options = app->get_options(
            [app, mode, &group](const CLI::Option* option)
            { return should_show_group_option(app, mode, group, option); });
        if (group.empty() || options.empty())
        {
            return {};
        }
        return make_group(group, false, std::move(options));
    }

    static auto should_show_group_option(
        const CLI::App* app,
        CLI::AppFormatMode mode,
        const std::string& group,
        const CLI::Option* option) -> bool
    {
        return option->get_group() == group && option->nonpositional()
               && (mode != CLI::AppFormatMode::Sub
                   || (app->get_help_ptr() != option
                       && app->get_help_all_ptr() != option));
    }

    static auto should_show_subcommand(const CLI::App* subcommand) -> bool
    {
        // An empty group hides a subcommand, mirroring option-group semantics.
        return !subcommand->get_disabled() && !subcommand->get_name().empty()
               && !subcommand->get_group().empty();
    }

    auto make_positional_option(const CLI::Option* option) const -> std::string
    {
        std::ostringstream out;
        const auto name = CLI::Formatter::make_option_name(option, true);
        const auto opts = make_option_suffix(option);
        const auto left = "  " + name + opts;
        out << "  " << option_name_text(name) << option_value_text(opts);
        if (left.length() < column_width_)
        {
            out << std::string(column_width_ - left.length(), ' ');
        }
        append_option_description(out, option, left.length());
        out << '\n';
        return out.str();
    }

    auto make_nonpositional_option(const CLI::Option* option) const
        -> std::string
    {
        std::ostringstream out;
        const auto names = CLI::detail::split(
            CLI::Formatter::make_option_name(option, false), ',');
        std::vector<std::string> short_names;
        std::vector<std::string> long_names;
        for (const auto& name : names)
        {
            if (name.contains("--"))
            {
                long_names.push_back(name);
            }
            else
            {
                short_names.push_back(name);
            }
        }

        auto short_text = CLI::detail::join(short_names, ", ");
        auto long_text = CLI::detail::join(long_names, ", ");
        const auto opts = make_option_suffix(option);

        const auto short_width = static_cast<std::size_t>(
            static_cast<double>(column_width_) * long_option_alignment_ratio_);
        const auto long_width = column_width_ - short_width;

        if (!short_text.empty())
        {
            short_text = "  " + short_text;
            if (long_text.empty() && !opts.empty())
            {
                short_text += opts;
            }
            if (!long_text.empty())
            {
                short_text += ",";
            }
            out << option_name_text(short_text);
        }
        if (short_text.length() < short_width)
        {
            out << std::string(short_width - short_text.length(), ' ');
        }

        if (!long_text.empty())
        {
            out << option_name_text(long_text);
            if (!opts.empty())
            {
                out << option_value_text(opts);
                long_text += opts;
            }
        }
        if (long_text.length() < long_width)
        {
            out << std::string(long_width - long_text.length(), ' ');
        }

        const auto left_length
            = std::max(short_text.length(), short_width) + long_text.length();
        append_option_description(out, option, left_length);
        out << '\n';
        return out.str();
    }

    auto make_option_suffix(const CLI::Option* option) const -> std::string
    {
        std::ostringstream out;
        if (!option->get_option_text().empty())
        {
            out << " " << option->get_option_text();
        }
        else if (option->get_type_size() != 0)
        {
            if (enable_option_type_names_ && !option->get_type_name().empty())
            {
                out << " " << get_label(make_option_type_name(option));
            }
            if (option->get_expected_max()
                == CLI::detail::expected_max_vector_size)
            {
                out << " ...";
            }
            else if (option->get_expected_min() > 1)
            {
                out << " x " << option->get_expected();
            }
        }
        return out.str();
    }

    static auto make_option_type_name(const CLI::Option* option) -> std::string
    {
        auto type_name = option->get_type_name();
        // CLI11 appends validator descriptions as ":<desc>" after the "<...>"
        // placeholder; strip them while preserving any ':' inside the
        // placeholder itself (e.g. "<A:B>").
        const auto base_end = type_name.find('>');
        const auto from = base_end == std::string::npos ? 0 : base_end + 1;
        if (const auto pos = type_name.find(':', from);
            pos != std::string::npos)
        {
            type_name.erase(pos);
        }
        return type_name;
    }

    auto append_option_description(
        std::ostringstream& out,
        const CLI::Option* option,
        std::size_t left_length) const -> void
    {
        const auto desc = make_option_desc(option);
        if (desc.empty())
        {
            return;
        }
        bool skip_first_line_prefix = true;
        if (left_length >= column_width_)
        {
            out << '\n';
            skip_first_line_prefix = false;
        }
        std::ostringstream desc_out;
        CLI::detail::streamOutAsParagraph(
            desc_out,
            desc,
            right_column_width_,
            std::string(column_width_, ' '),
            skip_first_line_prefix);
        out << color_default_text(option, desc_out.str());
    }

    static auto color_default_text(
        const CLI::Option* option,
        std::string desc_text) -> std::string
    {
        if (option->get_default_str().empty())
        {
            return desc_text;
        }
        const auto default_text = make_default_text(option);
        if (const auto pos = desc_text.rfind(default_text);
            pos != std::string::npos)
        {
            desc_text.replace(
                pos, default_text.size(), option_name_text(default_text));
        }
        return desc_text;
    }

    static auto make_default_text(const CLI::Option* option) -> std::string
    {
        return "[default: " + option->get_default_str() + "]";
    }
};

}  // namespace

auto make_cli_formatter() -> std::shared_ptr<CLI::FormatterBase>
{
    auto formatter = std::make_shared<ColorFormatter>();
    formatter->column_width(20);
    formatter->right_column_width(80);
    formatter->enable_footer_formatting(false);
    return formatter;
}

auto format_parse_error(const CLI::App* app, const CLI::Error& err)
    -> std::string
{
    const CLI::App* active = app;
    if (const auto entered = app->get_subcommands(); !entered.empty())
    {
        active = entered.back();
    }
    const auto formatter
        = std::dynamic_pointer_cast<CLI::Formatter>(app->get_formatter());
    return fmt::format(
        "{}{}\n{}For more information, try '{}'\n",
        error_marker(),
        err.what(),
        formatter->make_usage(active, ""),
        option_name_text("--help"));
}

}  // namespace cli
