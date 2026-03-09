from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass
class FastaRecord:
    header: str
    sequence: str

    @property
    def gene_id(self) -> str:
        return self.header.split("|")[0].split()[0]


def read_fasta(path: str | Path) -> Iterable[FastaRecord]:
    header = None
    seq = []
    with open(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield FastaRecord(header=header, sequence="".join(seq).upper())
                header = line[1:]
                seq = []
            else:
                seq.append(line)
    if header is not None:
        yield FastaRecord(header=header, sequence="".join(seq).upper())
