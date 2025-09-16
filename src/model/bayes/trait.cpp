#include "trait.h"
#include "model/bayes/traits/BayesA.h"

namespace gelex
{
std::unique_ptr<GeneticTrait> create_genetic_trait(BayesAlphabet type)
{
    switch (type)
    {
        case BayesAlphabet::A:
            return std::make_unique<BayesATrait>();
        case BayesAlphabet::B:
            return std::make_unique<BayesBTrait>();
        case BayesAlphabet::Bpi:
            return std::make_unique<BayesBpiTrait>();
        case BayesAlphabet::C:
            return std::make_unique<BayesCTrait>();
        case BayesAlphabet::Cpi:
            return std::make_unique<BayesCpiTrait>();
        case BayesAlphabet::RR:
            return std::make_unique<BayesRRTrait>();
        default:
            throw std::invalid_argument("Unknown Bayes trait type");
    }
}
}  // namespace gelex
