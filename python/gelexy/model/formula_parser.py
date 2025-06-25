import re
from dataclasses import dataclass
from typing import Final

import pandas as pd


@dataclass
class GeneticTerm:
    """Represents a genetic term in the formula with its associated GRM."""

    name: str
    genetic: str
    env: str | None = None


class FormulaParser:
    """Parses mixed model formulas with genetic and GxE terms."""

    _GENETIC_PATTERN: Final = re.compile(r"g\[(.*?)\]")

    def __init__(self, formula: str, genetic_input: list[str]) -> None:
        """Initialize parser with formula and capture caller's locals.
        Args:
            formula: Mixed model formula (e.g., "y ~ x + g(GRM)")
        """
        self.genetic_input = genetic_input
        self.common = ""
        self.genetic_terms: list[GeneticTerm] = []
        self.gxe_terms: list[GeneticTerm] = []

        self.formula = ""
        self._parse(formula)

    def _parse(self, formula: str) -> None:
        """Parse the formula into components."""
        self._check_formula(formula)
        cleaned_formula = formula.replace(" ", "")
        lhs, rhs = cleaned_formula.split("~", 1)

        self.response = lhs
        self._build_terms(lhs, rhs)
        self.formula = self._format_formula(cleaned_formula)

    def _check_formula(self, formula: str) -> None:
        """Validate formula structure."""
        if "~" not in formula:
            msg = "Formula must contain '~' operator"
            raise ValueError(msg)

    def _check_term(self, term: str) -> tuple[str, pd.DataFrame]:
        """Parse genetic term into name and GRM DataFrame.
        Args:
            term: Either a variable name or file path
        Returns:
            Tuple of (name, GRM DataFrame)
        """
        if term not in self.genetic_input:
            msg = f"Genetic term '{term}' not found in provided GRMs."
            raise ValueError(msg)
        return term

    def _build_terms(self, lhs: str, rhs: str) -> str:
        """Build the common (non-genetic) terms portion."""
        terms = rhs.split("+")
        common_terms = []

        for term in terms:
            if genetic_match := self._GENETIC_PATTERN.match(term):
                self._process_genetic_term(genetic_match.group(1))
            else:
                common_terms.append(term)
        self.common = f"{lhs}~{'+'.join(common_terms)}"

    def _process_genetic_term(self, term: str) -> None:
        """Process genetic term (either simple or GxE)."""
        if ":" in term:
            genetic, env = term.split(":", 1)
            genetic = self._check_term(genetic)
            self.gxe_terms.append(GeneticTerm(term, genetic, env))
        else:
            genetic = self._check_term(term)
            self.genetic_terms.append(GeneticTerm(term, genetic))

    def _format_formula(self, common_term: str) -> str:
        """Format common terms with consistent spacing."""
        operators = r"([~+\-=*\/^])"
        return re.sub(r"\s*" + operators + r"\s*", r" \1 ", common_term).strip()
