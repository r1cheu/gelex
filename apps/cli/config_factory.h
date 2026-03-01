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

#ifndef GELEX_CLI_CONFIG_FACTORY_H_
#define GELEX_CLI_CONFIG_FACTORY_H_

#include <concepts>

namespace argparse
{
class ArgumentParser;
}

namespace gelex::cli
{

template <typename T>
concept CliConfig = requires(argparse::ArgumentParser& args, const T& config) {
    { T::make(args) } -> std::same_as<T>;
    { config.validate() } -> std::same_as<void>;
};

template <CliConfig T>
auto make_config(argparse::ArgumentParser& args) -> T
{
    T config = T::make(args);
    config.validate();
    return config;
}

}  // namespace gelex::cli

#endif  // GELEX_CLI_CONFIG_FACTORY_H_
