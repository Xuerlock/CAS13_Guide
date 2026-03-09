from __future__ import annotations

import pandas as pd


def _score_row(row: pd.Series) -> float:
    # higher is better
    gc_dev = abs(row["gc"] - 50.0)
    pos_dev = abs(row["start"] - 115)
    return -(gc_dev) - 0.02 * pos_dev


def _pick_spread(group: pd.DataFrame, per_gene: int, min_dist: int) -> pd.DataFrame:
    group = group.copy()
    group["score"] = group.apply(_score_row, axis=1)
    group = group.sort_values(["score", "start"], ascending=[False, True])
    picked = []
    picked_starts = []
    for _, row in group.iterrows():
        if len(picked) >= per_gene:
            break
        if not picked_starts or min(abs(row["start"] - s) for s in picked_starts) >= min_dist:
            picked.append(row)
            picked_starts.append(int(row["start"]))
    if len(picked) < per_gene:
        used = {r["GuideID"] for r in picked}
        for _, row in group.iterrows():
            if len(picked) >= per_gene:
                break
            if row["GuideID"] not in used:
                picked.append(row)
                used.add(row["GuideID"])
    return pd.DataFrame(picked)


def select_guides(
    merged_df: pd.DataFrame,
    per_gene: int = 5,
    strict_mm0: int = 1,
    strict_mm1: int = 1,
    strict_mm2: int = 1,
    min_dist: int = 25,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    filt = merged_df[
        (merged_df["n0"] == strict_mm0)
        & (merged_df["n1"] <= strict_mm1)
        & (merged_df["n2"] <= strict_mm2)
        & (merged_df["gc"] >= 40)
        & (merged_df["gc"] <= 60)
    ].copy()

    final = []
    for gene, sub in filt.groupby("GeneID", sort=True):
        picked = _pick_spread(sub, per_gene=per_gene, min_dist=min_dist)
        final.append(picked)
    final_df = pd.concat(final, ignore_index=True) if final else pd.DataFrame(columns=filt.columns)
    final_df["rank"] = final_df.groupby("GeneID").cumcount() + 1
    qc = final_df.groupby("GeneID", as_index=False).size().rename(columns={"size": "n"})
    lt5 = qc[qc["n"] < per_gene].copy()
    return final_df, qc, lt5
