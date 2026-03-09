from __future__ import annotations

import subprocess
from pathlib import Path


def _run_count(fasta: str, index: str, mismatches: int, out_tsv: str) -> None:
    cmd = (
        f"bowtie -f -v {mismatches} -a {index} {fasta} /dev/stdout "
        f"| cut -f1 | sort | uniq -c > {out_tsv}"
    )
    subprocess.run(cmd, shell=True, check=True)


def count_mismatches(fasta: str | Path, index: str | Path, output_prefix: str | Path) -> tuple[str, str, str]:
    output_prefix = str(output_prefix)
    mm0 = f"{output_prefix}.mm0.counts.tsv"
    mm1 = f"{output_prefix}.mm1.counts.tsv"
    mm2 = f"{output_prefix}.mm2.counts.tsv"
    _run_count(str(fasta), str(index), 0, mm0)
    _run_count(str(fasta), str(index), 1, mm1)
    _run_count(str(fasta), str(index), 2, mm2)
    return mm0, mm1, mm2
