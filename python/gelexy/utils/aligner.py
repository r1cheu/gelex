import pandas as pd


class Matcher:
    def __init__(
        self,
        phenotype: pd.DataFrame,
        genotypes: dict[str, pd.DataFrame],
        axis: int | list[int] = 0,
    ):
        if "id" not in phenotype.columns:
            msg = "Phenotype DataFrame must have an 'id' column."
            raise ValueError(msg)
        self._axis = axis if isinstance(axis, list) else [axis]

        self.common_set = set(phenotype["id"].astype(str))
        self._intersect(genotypes)
        self.common_order = sorted(self.common_set)

        self.phenotype = self._match_phenotype(phenotype)
        self.genotypes = self._match_genotype(genotypes)

        self.phenotype_ids = self.phenotype["id"].astype(str).tolist()
        self.genotype_ids = self.common_order

    def _match_phenotype(self, phenotype: pd.DataFrame) -> pd.DataFrame:
        return phenotype[phenotype["id"].isin(self.common_set)]

    def _match_genotype(
        self, genotypes: dict[str, pd.DataFrame]
    ) -> dict[str, pd.DataFrame]:
        for k, _ in genotypes.items():
            for a in self._axis:
                if a == 0:
                    genotypes[k] = genotypes[k].loc[self.common_order, :]
                elif a == 1:
                    genotypes[k] = genotypes[k].loc[:, self.common_order]
                else:
                    msg = (
                        f"Axis {a} is not supported. Only 0 and 1 are allowed."
                    )
                    raise ValueError(msg)
        return genotypes

    def _intersect(self, genotypes: dict[str, pd.DataFrame]):
        for k, g in genotypes.items():
            if not isinstance(g, pd.DataFrame):
                msg = f"Genotype for {k} must be a DataFrame."
                raise ValueError(msg)
            if g.index.name != "id":
                msg = f"Genotype DataFrame for {k} must have an 'id' index."
                raise ValueError(msg)
            self.common_set.intersection_update(set(g.index.astype(str)))
