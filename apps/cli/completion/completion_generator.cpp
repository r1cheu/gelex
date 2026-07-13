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

#include "cli/completion/completion_generator.h"

#include <CLI/CLI.hpp>
#include <algorithm>
#include <array>
#include <ranges>
#include <string>
#include <string_view>
#include <vector>

#include "cli/completion/choice_descriptions.h"

namespace cli
{
namespace
{
// Path-like option type-name tokens receive file completion; every other
// value option (numeric, column, plain text) is completed with nothing.
constexpr std::array<std::string_view, 11> PATH_TYPE_TOKENS{
    "BFILE",
    "OUT",
    "GRM",
    "GFILE",
    "PREFIX",
    "QRAND",
    "QCOVAR",
    "DRAND",
    "DCOVAR",
    "PHENOTYPE",
    "CKPT"};

struct Choice
{
    std::string value;
    std::string_view desc;
};

struct ArgInfo
{
    std::vector<std::string> longs;
    std::vector<std::string> shorts;
    std::string desc;
    bool takes_arg{};
    std::vector<Choice> choices;
    bool file{};
};

struct CmdInfo
{
    std::string name;
    std::string desc;
    std::vector<ArgInfo> args;
};

struct Model
{
    std::vector<ArgInfo> root_args;
    std::vector<CmdInfo> cmds;
};

auto base_type_token(std::string_view type_name) -> std::string_view
{
    const auto open = type_name.find('<');
    const auto close = type_name.find('>');
    if (open == std::string_view::npos || close == std::string_view::npos
        || close < open)
    {
        return {};
    }
    return type_name.substr(open + 1, close - open - 1);
}

auto parse_choices(std::string_view type_name) -> std::vector<std::string>
{
    const auto open = type_name.find('{');
    const auto close = type_name.find('}', open);
    if (open == std::string_view::npos || close == std::string_view::npos)
    {
        return {};
    }
    const auto inner = type_name.substr(open + 1, close - open - 1);
    std::vector<std::string> choices;
    for (const auto token : std::views::split(inner, ','))
    {
        if (!token.empty())
        {
            choices.emplace_back(token.begin(), token.end());
        }
    }
    return choices;
}

auto is_path_token(std::string_view token) -> bool
{
    return std::ranges::contains(PATH_TYPE_TOKENS, token);
}

auto to_arg_info(const CLI::Option& opt) -> ArgInfo
{
    ArgInfo info;
    info.longs = opt.get_lnames();
    info.shorts = opt.get_snames();
    info.desc = opt.get_description();
    info.takes_arg = opt.get_expected_max() != 0;
    if (info.takes_arg)
    {
        const auto type_name = opt.get_type_name();
        const auto type_token = base_type_token(type_name);
        const auto group = choice_group(type_token);
        for (const auto& token : parse_choices(type_name))
        {
            std::string_view desc;
            for (const auto& entry : group)
            {
                if (entry.first == token)
                {
                    desc = entry.second;
                    break;
                }
            }
            info.choices.push_back({token, desc});
        }
        info.file = info.choices.empty() && is_path_token(type_token);
    }
    return info;
}

auto collect_args(const CLI::App& app) -> std::vector<ArgInfo>
{
    std::vector<ArgInfo> args;
    for (const auto* opt : app.get_options())
    {
        if (!opt->nonpositional() || opt->get_group().empty())
        {
            continue;
        }
        args.push_back(to_arg_info(*opt));
    }
    return args;
}

auto build_model(const CLI::App& root) -> Model
{
    Model model;
    model.root_args = collect_args(root);
    for (const auto* sub :
         root.get_subcommands([](const CLI::App*) { return true; }))
    {
        if (sub->get_group().empty())
        {
            continue;
        }
        model.cmds.push_back(
            {sub->get_name(), sub->get_description(), collect_args(*sub)});
    }
    return model;
}

// Fish single-quoted strings only escape backslash and the quote itself.
auto escape_fish(std::string_view text) -> std::string
{
    std::string out;
    out.reserve(text.size());
    for (const char ch : text)
    {
        if (ch == '\\' || ch == '\'')
        {
            out += '\\';
        }
        out += ch;
    }
    return out;
}

auto fish_option_line(
    std::string_view program,
    std::string_view condition,
    const ArgInfo& arg) -> std::string
{
    std::string flags = "complete -c ";
    flags += program;
    flags += " -n '";
    flags += condition;
    flags += '\'';
    for (const auto& name : arg.longs)
    {
        flags += " -l " + name;
    }
    for (const auto& name : arg.shorts)
    {
        flags += " -s " + name;
    }

    // Enum options expand to one rule per candidate so fish shows each member's
    // own description instead of the option-level help. A -r rule falls back to
    // file completion unless it forbids files, so -f suppresses that here.
    if (arg.takes_arg && !arg.choices.empty())
    {
        std::string out;
        for (const auto& choice : arg.choices)
        {
            out += flags;
            out += " -r -f -a '";
            out += choice.value;
            out += '\'';
            if (!choice.desc.empty())
            {
                out += " -d '";
                out += escape_fish(choice.desc);
                out += '\'';
            }
            out += '\n';
        }
        return out;
    }

    std::string line = flags;
    if (!arg.desc.empty())
    {
        line += " -d '" + escape_fish(arg.desc) + '\'';
    }
    if (arg.takes_arg)
    {
        line += arg.file ? " -r -F" : " -r -f";
    }
    line += '\n';
    return line;
}

// Dashed CLI tokens for one option, e.g. {"--method", "-m"}.
auto arg_tokens(const ArgInfo& arg) -> std::vector<std::string>
{
    std::vector<std::string> tokens;
    tokens.reserve(arg.longs.size() + arg.shorts.size());
    for (const auto& name : arg.longs)
    {
        tokens.push_back("--" + name);
    }
    for (const auto& name : arg.shorts)
    {
        tokens.push_back("-" + name);
    }
    return tokens;
}

auto flag_tokens(const std::vector<ArgInfo>& args) -> std::vector<std::string>
{
    std::vector<std::string> tokens;
    for (const auto& arg : args)
    {
        for (auto& token : arg_tokens(arg))
        {
            tokens.push_back(std::move(token));
        }
    }
    return tokens;
}

// Alternation pattern of every flag token, e.g. "--method|-m".
auto prev_pattern(const ArgInfo& arg) -> std::string
{
    return CLI::detail::join(arg_tokens(arg), "|");
}

}  // namespace

auto generate_fish_completion(const CLI::App& root) -> std::string
{
    const auto model = build_model(root);
    const auto& program = root.get_name();

    std::string out = "# fish completion for " + program + "\n";
    out += "complete -c " + program + " -f\n";

    for (const auto& cmd : model.cmds)
    {
        out += "complete -c " + program + " -n '__fish_use_subcommand' -a '"
               + cmd.name + '\'';
        if (!cmd.desc.empty())
        {
            out += " -d '" + escape_fish(cmd.desc) + '\'';
        }
        out += '\n';
    }
    for (const auto& arg : model.root_args)
    {
        out += fish_option_line(program, "__fish_use_subcommand", arg);
    }
    for (const auto& cmd : model.cmds)
    {
        const std::string condition = "__fish_seen_subcommand_from " + cmd.name;
        for (const auto& arg : cmd.args)
        {
            out += fish_option_line(program, condition, arg);
        }
    }
    return out;
}

auto generate_bash_completion(const CLI::App& root) -> std::string
{
    const auto model = build_model(root);
    const auto& program = root.get_name();

    std::vector<std::string> cmd_names;
    cmd_names.reserve(model.cmds.size());
    for (const auto& cmd : model.cmds)
    {
        cmd_names.push_back(cmd.name);
    }

    std::string out = "# bash completion for " + program + "\n";
    out += "_" + program + "() {\n";
    out += "    local cur prev cmd i\n";
    out += "    cur=\"${COMP_WORDS[COMP_CWORD]}\"\n";
    out += "    prev=\"${COMP_WORDS[COMP_CWORD-1]}\"\n";
    out += "    cmd=\"\"\n";
    out += "    for ((i=1; i<COMP_CWORD; i++)); do\n";
    out += "        case \"${COMP_WORDS[i]}\" in\n";
    out += "            -*) ;;\n";
    out += "            *) cmd=\"${COMP_WORDS[i]}\"; break;;\n";
    out += "        esac\n";
    out += "    done\n";
    out += "    case \"$cmd\" in\n";

    for (const auto& cmd : model.cmds)
    {
        out += "        " + cmd.name + ")\n";

        std::vector<const ArgInfo*> value_args;
        for (const auto& arg : cmd.args)
        {
            if (arg.takes_arg)
            {
                value_args.push_back(&arg);
            }
        }
        if (!value_args.empty())
        {
            out += "            case \"$prev\" in\n";
            for (const auto* arg : value_args)
            {
                out += "                " + prev_pattern(*arg) + ")\n";
                if (!arg->choices.empty())
                {
                    out += "                    COMPREPLY=($(compgen -W \""
                           + CLI::detail::join(
                               arg->choices,
                               [](const Choice& choice)
                               { return choice.value; },
                               " ")
                           + "\" -- \"$cur\")); return;;\n";
                }
                else if (arg->file)
                {
                    out += "                    COMPREPLY=($(compgen -f -- "
                           "\"$cur\")); return;;\n";
                }
                else
                {
                    out += "                    COMPREPLY=(); return;;\n";
                }
            }
            out += "            esac\n";
        }
        out += "            COMPREPLY=($(compgen -W \""
               + CLI::detail::join(flag_tokens(cmd.args), " ")
               + "\" -- \"$cur\"))\n";
        out += "            return;;\n";
    }

    out += "        \"\")\n";
    out += "            if [[ \"$cur\" == -* ]]; then\n";
    out += "                COMPREPLY=($(compgen -W \""
           + CLI::detail::join(flag_tokens(model.root_args), " ")
           + "\" -- \"$cur\"))\n";
    out += "            else\n";
    out += "                COMPREPLY=($(compgen -W \""
           + CLI::detail::join(cmd_names, " ") + "\" -- \"$cur\"))\n";
    out += "            fi\n";
    out += "            return;;\n";
    out += "    esac\n";
    out += "}\n";
    out += "complete -F _" + program + ' ' + program + '\n';
    return out;
}

}  // namespace cli
