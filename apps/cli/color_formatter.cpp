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

#include <algorithm>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <fmt/color.h>
#include <CLI/CLI.hpp>

#include "cli_helper.h"

namespace cli
{

namespace
{

constexpr fmt::rgb CATPPUCCIN_SKY{0x89DCEB};
constexpr fmt::rgb CATPPUCCIN_GREEN{0xA6E3A1};
constexpr fmt::rgb CATPPUCCIN_YELLOW{0xF9E2AF};

auto colored(std::string text, fmt::text_style style) -> std::string
{
    if (text.empty() || !is_tty())
    {
        return text;
    }
    return fmt::format(style, "{}", text);
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
        if (app->get_parent() != nullptr && !description.empty())
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
        out << "\n"
            << colored(
                   std::move(group) + ":",
                   fmt::emphasis::bold | fmt::emphasis::underline
                       | fmt::fg(CATPPUCCIN_GREEN))
            << "\n";
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
            auto options = app->get_options(
                [app, mode, &group](const CLI::Option* option)
                {
                    return option->get_group() == group
                           && option->nonpositional()
                           && (mode != CLI::AppFormatMode::Sub
                               || (app->get_help_ptr() != option
                                   && app->get_help_all_ptr() != option));
                });
            if (!group.empty() && !options.empty())
            {
                out << make_group(group, false, std::move(options));
            }
        }

        for (const auto& group : groups)
        {
            if (group != "OPTIONS")
            {
                continue;
            }
            auto options = app->get_options(
                [app, mode, &group](const CLI::Option* option)
                {
                    return option->get_group() == group
                           && option->nonpositional()
                           && (mode != CLI::AppFormatMode::Sub
                               || (app->get_help_ptr() != option
                                   && app->get_help_all_ptr() != option));
                });
            if (!group.empty() && !options.empty())
            {
                out << make_group(group, false, std::move(options));
            }
        }

        return out.str();
    }

    auto make_subcommands(const CLI::App* app, CLI::AppFormatMode mode) const
        -> std::string override
    {
        auto subcommands = app->get_subcommands(
            [](const CLI::App* subcommand)
            {
                return !subcommand->get_disabled()
                       && !subcommand->get_name().empty();
            });
        if (subcommands.empty())
        {
            return {};
        }

        std::ostringstream out;
        out << "\n"
            << colored(
                   "Commands:",
                   fmt::emphasis::bold | fmt::emphasis::underline
                       | fmt::fg(CATPPUCCIN_GREEN))
            << "\n";
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
        out << colored(
            "Usage:",
            fmt::emphasis::bold | fmt::emphasis::underline
                | fmt::fg(CATPPUCCIN_GREEN))
            << "\n  ";
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
            out << " " << colored(option_name, fmt::fg(CATPPUCCIN_SKY))
                << make_option_opts(option);
        }

        auto optional_options = app->get_options(
            [](const CLI::Option* option)
            { return option->nonpositional() && !option->get_required(); });
        if (!optional_options.empty())
        {
            out << " " << colored("[OPTIONS]", fmt::fg(CATPPUCCIN_SKY));
        }

        auto positionals = app->get_options(
            [](const CLI::Option* option) { return option->get_positional(); });
        for (const CLI::Option* positional : positionals)
        {
            out << " " << make_option_usage(positional);
        }

        auto subcommands = app->get_subcommands(
            [](const CLI::App* subcommand)
            {
                return !subcommand->get_disabled()
                       && !subcommand->get_name().empty();
            });
        if (!subcommands.empty())
        {
            out << " " << colored("[COMMAND]", fmt::fg(CATPPUCCIN_SKY));
        }

        out << "\n";
        return out.str();
    }

    auto make_option_name(const CLI::Option* option, bool is_positional) const
        -> std::string override
    {
        return colored(
            CLI::Formatter::make_option_name(option, is_positional),
            fmt::fg(CATPPUCCIN_SKY));
    }

    auto make_option(const CLI::Option* option, bool is_positional) const
        -> std::string override
    {
        std::ostringstream out;
        if (is_positional)
        {
            const auto name = CLI::Formatter::make_option_name(option, true);
            const auto opts = [&option, this]()
            {
                std::ostringstream out;
                if (!option->get_option_text().empty())
                {
                    out << " " << option->get_option_text();
                }
                else if (option->get_type_size() != 0)
                {
                    if (enable_option_type_names_
                        && !option->get_type_name().empty())
                    {
                        auto type_name = option->get_type_name();
                        if (const auto pos = type_name.find(':');
                            pos != std::string::npos)
                        {
                            type_name.erase(pos);
                        }
                        out << " " << get_label(type_name);
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
            }();
            const auto left = "  " + name + opts;
            out << "  " << colored(name, fmt::fg(CATPPUCCIN_SKY))
                << colored(opts, fmt::fg(CATPPUCCIN_YELLOW));
            if (left.length() < column_width_)
            {
                out << std::string(column_width_ - left.length(), ' ');
            }

            const auto desc = make_option_desc(option);
            if (!desc.empty())
            {
                bool skip_first_line_prefix = true;
                if (left.length() >= column_width_)
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
                auto desc_text = desc_out.str();
                if (!option->get_default_str().empty())
                {
                    const auto default_text
                        = "[default: " + option->get_default_str() + "]";
                    if (const auto pos = desc_text.rfind(default_text);
                        pos != std::string::npos)
                    {
                        desc_text.replace(
                            pos,
                            default_text.size(),
                            colored(default_text, fmt::fg(CATPPUCCIN_SKY)));
                    }
                }
                out << desc_text;
            }
            out << '\n';
            return out.str();
        }

        const auto names = CLI::detail::split(
            CLI::Formatter::make_option_name(option, false), ',');
        std::vector<std::string> short_names;
        std::vector<std::string> long_names;
        for (const auto& name : names)
        {
            if (name.find("--", 0) != std::string::npos)
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
        const auto opts = [&option, this]()
        {
            std::ostringstream out;
            if (!option->get_option_text().empty())
            {
                out << " " << option->get_option_text();
            }
            else if (option->get_type_size() != 0)
            {
                if (enable_option_type_names_
                    && !option->get_type_name().empty())
                {
                    auto type_name = option->get_type_name();
                    if (const auto pos = type_name.find(':');
                        pos != std::string::npos)
                    {
                        type_name.erase(pos);
                    }
                    out << " " << get_label(type_name);
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
        }();

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
            out << colored(short_text, fmt::fg(CATPPUCCIN_SKY));
        }
        if (short_text.length() < short_width)
        {
            out << std::string(short_width - short_text.length(), ' ');
        }

        if (!long_text.empty())
        {
            out << colored(long_text, fmt::fg(CATPPUCCIN_SKY));
            if (!opts.empty())
            {
                out << colored(opts, fmt::fg(CATPPUCCIN_YELLOW));
                long_text += opts;
            }
        }
        if (long_text.length() < long_width)
        {
            out << std::string(long_width - long_text.length(), ' ');
        }

        const auto left_length
            = std::max(short_text.length(), short_width) + long_text.length();
        const auto desc = make_option_desc(option);
        if (!desc.empty())
        {
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
            auto desc_text = desc_out.str();
            if (!option->get_default_str().empty())
            {
                const auto default_text
                    = "[default: " + option->get_default_str() + "]";
                if (const auto pos = desc_text.rfind(default_text);
                    pos != std::string::npos)
                {
                    desc_text.replace(
                        pos,
                        default_text.size(),
                        colored(default_text, fmt::fg(CATPPUCCIN_SKY)));
                }
            }
            out << desc_text;
        }
        out << '\n';
        return out.str();
    }

    auto make_option_opts(const CLI::Option* option) const
        -> std::string override
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
                auto type_name = option->get_type_name();
                if (const auto pos = type_name.find(':');
                    pos != std::string::npos)
                {
                    type_name.erase(pos);
                }
                out << " " << get_label(type_name);
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
        return colored(out.str(), fmt::fg(CATPPUCCIN_YELLOW));
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
        desc += "[default: " + option->get_default_str() + "]";
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
        out << colored(name, fmt::fg(CATPPUCCIN_SKY));
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
                out << colored(
                    std::move(line),
                    fmt::emphasis::bold | fmt::emphasis::underline
                        | fmt::fg(CATPPUCCIN_GREEN));
            }
            else
            {
                out << line;
            }
            out << '\n';
        }
        return out.str();
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

}  // namespace cli
