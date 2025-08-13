import pandas as pd

from .._logger import setup_logger


def align_gblup(
    phenotype: pd.DataFrame,
    genotypes: dict[str, pd.DataFrame],
):
    logger = setup_logger(__name__)
    if "id" not in phenotype.columns:
        msg = "Phenotype DataFrame must have an 'id' column."
        raise ValueError(msg)

    pid = set(phenotype["id"].astype(str))
    gid_list = next(iter(genotypes.values())).index.to_list()
    gid = set(gid_list)
    drop = pid - gid
    if drop:
        pid = pid - drop
        msg = f"The following IDs {', '.join(list(drop))} have phenotype data but no genotype data. They will be dropped."
        logger.warning(msg)

    return phenotype[phenotype["id"].isin(pid)], gid_list


def align_bayes(
    phenotype: pd.DataFrame,
    genotype: pd.DataFrame,
):
    if "id" not in phenotype.columns:
        msg = "Phenotype DataFrame must have an 'id' column."
        raise ValueError(msg)

    pid = set(phenotype["id"].astype(str))
    gid = set(genotype.index.to_list())
    common_ids = pid.intersection(gid)
    common_ids_list = list(common_ids)

    return phenotype[phenotype["id"].isin(common_ids)], genotype.loc[
        common_ids_list
    ]
