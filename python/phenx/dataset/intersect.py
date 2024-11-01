from .genotype import Genotypes
from .phenotype import Phenotypes


def intersect(obj1, obj2):
    if isinstance(obj1, Genotypes) and isinstance(obj2, Phenotypes):
        commom_sample = obj1.data.columns.intersection(obj2.data.index)
        obj1.data = obj1.data.loc[:, commom_sample]
        obj2.data = obj2.data.loc[commom_sample, :]
        return obj1, obj2

    if isinstance(obj1, Phenotypes) and isinstance(obj2, Genotypes):
        commom_sample = obj1.data.index.intersection(obj2.data.columns)
        obj1.data = obj1.data.loc[commom_sample, :]
        obj2.data = obj2.data.loc[:, commom_sample]
        return obj1, obj2
    msg = f"Intersection is only supported between Genotypes and Phenotypes objects. But got {type(obj1)} and {type(obj2)}"
    raise ValueError(msg)
