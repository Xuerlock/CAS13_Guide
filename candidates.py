from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd

from fasta import read_fasta
from qc import gc_percent, has_dinuc_repeat, has_homopolymer


@dataclass
class Candidate:
    guide_id: str
    gene_id: str
    seq_id: str
    start: int
    spacer: str
    gc: float
    length: int


def generate_candidates(
    fasta_path: str | Path,
    guide_len: int,
    window_start: int,
    window_end: int,
    gc_min: float,
    gc_max: float,
    homopoly_k: int = 5,
    remove_dinuc: bool = True,
) -> pd.DataFrame:
    rows: list[dict] = []
    for rec in read_fasta(fasta_path):
        seq = rec.sequence
        left = max(0, window_start)
        right = min(len(seq), window_end)
        if right - left < guide_len:
            continue
        for pos in range(left, right - guide_len + 1):
            spacer = seq[pos:pos + guide_len]
            gc = gc_percent(spacer)
            if gc < gc_min or gc > gc_max:
                continue
            if has_homopolymer(spacer, k=homopoly_k):
                continue
            if remove_dinuc and has_dinuc_repeat(spacer):
                continue
            guide_id = f"{rec.header}__{pos}"
            rows.append(
                {
                    "GuideID": guide_id,
                    "GeneID": rec.gene_id,
                    "SeqID": rec.header,
                    "start": pos,
                    "spacer": spacer,
                    "gc": gc,
                    "len": guide_len,
                }
            )
    return pd.DataFrame(rows)


def write_candidate_outputs(df: pd.DataFrame, output_prefix: str | Path) -> tuple[str, str]:
    output_prefix = str(output_prefix)
    tsv = f"{output_prefix}.candidates.tsv"
    fa = f"{output_prefix}.candidates.fa"
    df.to_csv(tsv, sep="\t", index=False)
    with open(fa, "w") as handle:
        for row in df.itertuples(index=False):
            handle.write(f">{row.GuideID}\n{row.spacer}\n")
    return tsv, fa
