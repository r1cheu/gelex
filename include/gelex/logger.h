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

#ifndef GELEX_LOGGER_H_
#define GELEX_LOGGER_H_

#include <memory>
#include <string_view>

#include <spdlog/formatter.h>
#include <spdlog/logger.h>

namespace gelex::logging
{
void initialize(std::string_view output_prefix = "output");

std::shared_ptr<spdlog::logger>& get();

class LevelFormatter : public spdlog::formatter
{
   public:
    LevelFormatter();

    void format(const spdlog::details::log_msg& msg, spdlog::memory_buf_t& dest)
        override;

    [[nodiscard]] std::unique_ptr<formatter> clone() const override;

   private:
    std::unique_ptr<spdlog::formatter> info_formatter_;
    std::unique_ptr<spdlog::formatter> default_formatter_;
};
}  // namespace gelex::logging

#endif  // GELEX_LOGGER_H_
