// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef APPS_CLI_MCMC_DIAGNOSTICS_REPORTER_H_
#define APPS_CLI_MCMC_DIAGNOSTICS_REPORTER_H_

#include <span>

#include "gelex/bayes/genetic/diagnostics_traits.h"

namespace cli
{

// Prints the model-level entries (everything except `*/coefficients`) as one
// table of mean, sd, ESS and R-hat; the remaining statistics live only in the
// summary file.
auto show_diagnostics(std::span<const gelex::DiagnosticEntry> entries) -> void;

}  // namespace cli

#endif  // APPS_CLI_MCMC_DIAGNOSTICS_REPORTER_H_
