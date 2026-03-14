#ifndef GELEX_CLI_PREDICT_REPORTER_H_
#define GELEX_CLI_PREDICT_REPORTER_H_

#include <memory>

#include "gelex/infra/logging/predict_event.h"

namespace spdlog
{
class logger;
}

namespace gelex::cli
{

class PredictReporter
{
   public:
    PredictReporter();

    auto on_event(const PredictParamsLoadedEvent& event) const -> void;
    auto on_event(const PredictSnpSelectionEvent& event) const -> void;
    auto on_event(const PredictDataLoadedEvent& event) const -> void;
    auto on_event(const PredictResultsWrittenEvent& event) const -> void;

    auto as_observer() -> PredictObserver
    {
        return [this](const PredictEvent& e)
        { std::visit([this](const auto& ev) { this->on_event(ev); }, e); };
    }

   private:
    std::shared_ptr<spdlog::logger> logger_;
};

}  // namespace gelex::cli

#endif  // GELEX_CLI_PREDICT_REPORTER_H_
