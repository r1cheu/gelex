import h5py
import pytest
from gelexy.model.formula_parser import FormulaParser


@pytest.fixture
def a_grm(tmp_path):
    path = tmp_path / "a.grm"
    with h5py.File(path, "w") as f:
        f.create_dataset("grm", data=[[1, 0], [0, 1]])
        f.create_dataset("individuals", data=["ind1", "ind2"])
    return path


@pytest.fixture
def d_grm(tmp_path):
    path = tmp_path / "d.grm"
    with h5py.File(path, "w") as f:
        f.create_dataset("grm", data=[[1, 0], [0, 1]])
        f.create_dataset("individuals", data=["ind1", "ind2"])
    return path


def test_genetic_term_parsing():
    """Test parsing with genetic terms"""

    formula = "y ~ x1 + g(a) + g(d)"
    parser = FormulaParser(formula, ["a", "d"])
    assert parser.response == "y"
    assert parser.common == "y~x1"
    assert len(parser.genetic_terms) == 2
    assert parser.genetic_terms[0].name == "a"
    assert parser.genetic_terms[1].name == "d"
    assert parser.genetic_terms[0].genetic == "a"
    assert parser.genetic_terms[1].genetic == "d"
    assert len(parser.gxe_terms) == 0


def test_gxe_terms_parsing():
    """Test parsing with mixed term types"""

    formula = "y ~ x1 + g(a:scale(e1)) + g(a:e1)"
    parser = FormulaParser(formula, ["a"])
    assert parser.response == "y"
    assert parser.common == "y~x1"
    assert len(parser.genetic_terms) == 0
    assert len(parser.gxe_terms) == 2

    assert parser.gxe_terms[0].name == "a:scale(e1)"
    assert parser.gxe_terms[0].genetic == "a"
    assert parser.gxe_terms[0].env == "scale(e1)"

    assert parser.gxe_terms[1].name == "a:e1"
    assert parser.gxe_terms[1].genetic == "a"
    assert parser.gxe_terms[1].env == "e1"


def test_all_terms_parsing():
    """Test parsing with mixed term types"""

    formula = "y ~ x1 + g(a) + g(a:e1) + x2"
    parser = FormulaParser(formula, ["a"])
    assert parser.response == "y"
    assert parser.common == "y~x1+x2"
    assert len(parser.genetic_terms) == 1
    assert parser.genetic_terms[0].name == "a"
    assert parser.genetic_terms[0].genetic == "a"

    assert len(parser.gxe_terms) == 1
    assert parser.gxe_terms[0].name == "a:e1"
    assert parser.gxe_terms[0].genetic == "a"
    assert parser.gxe_terms[0].env == "e1"
