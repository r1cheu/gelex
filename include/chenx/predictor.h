#include <armadillo>
#include <cstdint>
#include <string_view>

#include "chenx/data/cross_grm.h"

namespace chenx
{
using arma::rowvec;
class Predictor
{
   public:
    Predictor(
        std::string_view train_bed_file,
        std::string_view method,
        rowvec&& center,
        double scale_factor,
        uint64_t chunk_size = 10000);
    Predictor(
        std::string_view train_bed_file,
        std::string_view method,
        rowvec&& center,
        double scale_factor);

   private:
    std::unique_ptr<BaseCrossGrm> cross_grm_;
};
}  // namespace chenx
