#pragma once

#include <memory>

#include "gelex/model/effects.h"
#include "traits/BayesA.h"
#include "traits/BayesB.h"
#include "traits/BayesBpi.h"
#include "traits/BayesC.h"
#include "traits/BayesCpi.h"
#include "traits/BayesRR.h"
#include "traits/base_trait.h"

namespace gelex
{
std::unique_ptr<GeneticTrait> create_genetic_trait(BayesAlphabet type);
}
