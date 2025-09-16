import pandas as pd

from gelexy._logger import logger


def align_gblup(
    phenotype: pd.DataFrame,
    genotypes: dict[str, pd.DataFrame],
):
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


def intersection(
    phenotype: pd.DataFrame,
    genotype: pd.Index,
    iid_only: bool = False,
) -> tuple[pd.DataFrame, list[str], list[str]]:
    id_column = "IID" if iid_only else "FID_IID"
    pid = pd.Index(phenotype[id_column])

    common_ids = genotype.intersection(pid)
    logger.info(
        f"Found {len(common_ids)} individuals with both phenotype and genotype data."
    )

    phenotype = phenotype[phenotype[id_column].isin(common_ids)]

    return (
        phenotype,
        phenotype[id_column].tolist(),
        genotype[genotype.isin(common_ids)].tolist(),
    )
