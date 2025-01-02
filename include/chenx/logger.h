#pragma once
#include <spdlog/spdlog.h>
#include <memory>

namespace chenx
{

class Logger
{
   public:
    static std::shared_ptr<spdlog::logger> logger();

   private:
    Logger();
    std::shared_ptr<spdlog::logger> GetSpdLogger();
    std::shared_ptr<spdlog::logger> logger_;
};

}  // namespace chenx
