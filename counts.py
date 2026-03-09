from __future__ import annotations

from pathlib import Path

import pandas as pd


def _read_count_file(path: str | Path, colname: str) -> pd.DataFrame:
    if not Path(path).exists() or Path(path).stat().st_size == 0:
        return pd.DataFrame(columns=["GuideID", colname])
    df = pd.read_csv(path, sep=r"\s+", header=None, names=[colname, "GuideID"])
    return df


def merge_counts(candidate_tsv: str | Path, mm0: str | Path, mm1: str | Path, mm2: str | Path) -> pd.DataFrame:
    base = pd.read_csv(candidate_tsv, sep="\t")
    c0 = _read_count_file(mm0, "n0")
    c1 = _read_count_file(mm1, "n1")
    c2 = _read_count_file(mm2, "n2")
    out = base.merge(c0, on="GuideID", how="left")\
              .merge(c1, on="GuideID", how="left")\
              .merge(c2, on="GuideID", how="left")
    for col in ["n0", "n1", "n2"]:
        out[col] = out[col].fillna(0).astype(int)
    return out
