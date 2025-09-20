#pragma once

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
