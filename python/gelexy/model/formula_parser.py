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

    _GENETIC_PATTERN: Final = re.compile(r"g\(((?:[^()]|\([^()]*\))*)\)")
    _FORMAT_PATTERN: Final = re.compile(r"\s*([+~])\s*")
    _GROUP_PATTERNL: Final = re.compile(r"\(.*?\|.*?\)")
    _GXE_DELIMITER: Final = ":"

    def __init__(self, formula: str, grms: list[str]) -> None:
        """Initialize parser with formula and capture caller's locals.
        Args:
            formula: Mixed model formula (e.g., "y ~ x + g(GRM)")
        """
        self.grms = grms
        self.genetic_terms: list[GeneticTerm] = []
        self.gxe_terms: list[GeneticTerm] = []
        self._parse(formula)

    def _parse(self, formula: str) -> None:
        """Parse the formula into components."""
        self._validate_formula(formula)
        cleaned_formula = formula.replace(" ", "")
        lhs, rhs = cleaned_formula.split("~", 1)

        self.response = lhs
        self.common = self._build_common_term(lhs, rhs)
        self.format_common = self._format_common_term(self.common)

    def _validate_formula(self, formula: str) -> None:
        """Validate formula structure."""
        if "~" not in formula:
            msg = "Formula must contain '~' operator"
            raise ValueError(msg)

    def _build_common_term(self, lhs: str, rhs: str) -> str:
        """Build the common (non-genetic) terms portion."""
        terms = rhs.split("+")
        common_terms = []

        for term in terms:
            if genetic_match := self._GENETIC_PATTERN.search(term):
                self._process_genetic_term(genetic_match.group(1))
            else:
                common_terms.append(term)

        return f"{lhs}~{'+'.join(common_terms)}"

    def _process_genetic_term(self, term_content: str) -> None:
        """Process genetic term (either simple or GxE)."""
        if self._GXE_DELIMITER in term_content:
            genetic, env = term_content.split(self._GXE_DELIMITER, 1)
            genetic = self._parse_genetic(genetic)
            self.gxe_terms.append(GeneticTerm(term_content, genetic, env))
        else:
            genetic = self._parse_genetic(term_content)
            self.genetic_terms.append(GeneticTerm(term_content, genetic))

    def _parse_genetic(self, term: str) -> tuple[str, pd.DataFrame]:
        """Parse genetic term into name and GRM DataFrame.
        Args:
            term: Either a variable name or file path
        Returns:
            Tuple of (name, GRM DataFrame)
        """
        if term not in self.grms:
            msg = f"Genetic term '{term}' not found in provided GRMs."
            raise ValueError(msg)
        return term

    def _format_common_term(self, common_term: str) -> str:
        """Format common terms with consistent spacing."""
        # First add consistent spacing
        formatted_term = self._FORMAT_PATTERN.sub(r" \1 ", common_term)
        formatted_term = self._GROUP_PATTERNL.sub("", formatted_term)
        formatted_term = re.sub(r"\+\++", "+", formatted_term)
        if not formatted_term.endswith("+ "):
            formatted_term += " + "
        return formatted_term
