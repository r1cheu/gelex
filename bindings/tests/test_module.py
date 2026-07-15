import gelex


def test_enums_exposed():
    assert gelex.GeneticMode.A != gelex.GeneticMode.D
    assert gelex.GenotypeMethod.Standardize is not None
