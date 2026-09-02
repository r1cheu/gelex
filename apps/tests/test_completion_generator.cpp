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

#include <CLI/CLI.hpp>
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <ranges>
#include <string>
#include <string_view>
#include <vector>

#include "gelex/bayes/builtin_method.h"
#include "gelex/data/genotype_method.h"
#include "gelex/data/rank_inverse_norm_transform.h"
#include "gelex/genetic_mode.h"

#include "cli/completion/choice_descriptions.h"
#include "cli/completion/completion_generator.h"

using Catch::Matchers::ContainsSubstring;

namespace
{
// CLI::App is neither copyable nor movable, so the app and the storage its
// options bind to must live together in one fixture.
struct AppFixture
{
    CLI::App app{"desc", "gelex"};
    std::string method;
    std::string bfile;
    std::string pheno;
    int threads{};
    bool loco{};

    AppFixture()
    {
        auto* mcmc = app.add_subcommand("mcmc", "Fit models with MCMC");
        mcmc->add_option("-m,--method", method)
            ->type_name("<MODEL>")
            ->check(CLI::IsMember(std::vector<std::string>{"RR", "A", "B"}));
        mcmc->add_option("-b,--bfile", bfile)->type_name("<BFILE>");
        mcmc->add_option("--pheno", pheno)
            ->type_name("<PHENOTYPE>")
            ->check(CLI::ExistingFile);
        mcmc->add_option("-t,--threads", threads)->type_name("<N>");
        mcmc->add_flag("--loco", loco);
    }
};
}  // namespace

TEST_CASE("fish completion introspects the app tree", "[cli][completion]")
{
    const AppFixture fx;
    const auto script = cli::generate_fish_completion(fx.app);

    SECTION("subcommand is offered")
    {
        REQUIRE_THAT(
            script, ContainsSubstring("-n '__fish_use_subcommand' -a 'mcmc'"));
    }
    SECTION("enum choices expand to one rule per candidate with a description")
    {
        REQUIRE_THAT(
            script,
            ContainsSubstring("-l method -s m -r -f -a 'RR' -d 'BayesRR'"));
        REQUIRE_THAT(
            script,
            ContainsSubstring("-l method -s m -r -f -a 'A' -d 'BayesA'"));
        REQUIRE_THAT(script, !ContainsSubstring("-a 'RR A B'"));
    }
    SECTION("path options get file completion")
    {
        REQUIRE_THAT(script, ContainsSubstring("-l bfile -s b -r -F"));
        REQUIRE_THAT(script, ContainsSubstring("-l pheno -r -F"));
    }
    SECTION("numeric option takes an arg but forbids files and candidates")
    {
        REQUIRE_THAT(script, ContainsSubstring("-l threads -s t -r -f"));
        REQUIRE_THAT(script, !ContainsSubstring("-l threads -s t -r -f -a"));
        REQUIRE_THAT(script, !ContainsSubstring("-l threads -s t -r -F"));
    }
    SECTION("flag takes no arg")
    {
        REQUIRE_THAT(script, ContainsSubstring("-l loco"));
        REQUIRE_THAT(script, !ContainsSubstring("-l loco -r"));
    }
}

TEST_CASE("bash completion introspects the app tree", "[cli][completion]")
{
    const AppFixture fx;
    const auto script = cli::generate_bash_completion(fx.app);

    REQUIRE_THAT(script, ContainsSubstring("complete -F _gelex gelex"));
    REQUIRE_THAT(script, ContainsSubstring("mcmc)"));
    REQUIRE_THAT(script, ContainsSubstring("--method|-m)"));
    REQUIRE_THAT(
        script, ContainsSubstring("compgen -W \"RR A B\" -- \"$cur\""));
    REQUIRE_THAT(script, ContainsSubstring("compgen -f -- \"$cur\""));
}

// The hand-written choice_descriptions table must cover exactly the tokens of
// the core name-map constant it mirrors: no missing member, no stale extra.
TEST_CASE(
    "choice descriptions stay in sync with the core name maps",
    "[cli][completion]")
{
    auto sorted = [](auto range)
    {
        std::vector<std::string_view> out(range.begin(), range.end());
        std::ranges::sort(out);
        return out;
    };

    // Each group's tokens must equal the mirrored name-map's names. The maps
    // store the name in different tuple positions, so each call names its own
    // projection.
    auto expect_sync
        = [&](std::string_view type_token, const auto& name_map, auto project)
    {
        REQUIRE(
            sorted(
                cli::choice_group(type_token)
                | std::views::transform([](const auto& entry)
                                        { return entry.first; }))
            == sorted(name_map | std::views::transform(project)));
    };

    expect_sync(
        "CODING",
        gelex::GENOTYPE_METHOD_CODES,
        [](const auto& entry) { return entry.first; });
    expect_sync(
        "MODEL",
        gelex::bayes_method_names,
        [](const auto& entry) { return entry.second; });
    expect_sync(
        "MODE",
        gelex::genetic_mode_set_names,
        [](const auto& entry) { return entry.second; });
    expect_sync(
        "TRANSFORM",
        gelex::rint_type_names,
        [](const auto& entry) { return entry.second; });
}
