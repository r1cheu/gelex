import pandas as pd


class Aligner:
    def __init__(
        self,
        index_ids: list[str],
    ):
        str_ids = [str(i) for i in index_ids]
        self.index_ids = list(dict.fromkeys(str_ids))
        self.index_set = set(self.index_ids)

    def align(
        self, other: pd.DataFrame, axis: int | list[int] = 0
    ) -> pd.DataFrame:
        axis = [axis] if isinstance(axis, int) else axis
        for a in axis:
            if a == 0:
                other.index = other.index.astype(str)
                set_other = set(other.index)
                missing_index = self.index_set.difference(set_other)
                if missing_index:
                    msg = f"The following index IDs are missing in the other DataFrame: {missing_index}"
                    raise KeyError(msg)
                index_for_reindex = [
                    id for id in self.index_ids if id in set_other
                ]
                other = other.reindex(index=index_for_reindex)
            elif a == 1:
                other.columns = other.columns.astype(str)
                set_other = set(other.columns)
                missing_columns = self.index_set.difference(set_other)
                if missing_columns:
                    msg = f"The following column IDs are missing in the other DataFrame: {missing_columns}"
                    raise KeyError(msg)
                columns_for_reindex = [
                    id for id in self.index_ids if id in set_other
                ]
                other = other.reindex(columns=columns_for_reindex)
            else:
                msg = "only 0, 1 are supported for axis."
                raise ValueError(msg)
        return other

    def difference(self, other: list[str]):
        set_other = set(other)
        return self.index_set.difference(set_other)
