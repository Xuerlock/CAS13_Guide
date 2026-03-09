from __future__ import annotations

import re


def gc_percent(seq: str) -> float:
    gc = seq.count("G") + seq.count("C")
    return round(100.0 * gc / len(seq), 2)


def has_homopolymer(seq: str, k: int = 5) -> bool:
    return bool(re.search(rf"A{{{k},}}|T{{{k},}}|C{{{k},}}|G{{{k},}}", seq))


def has_dinuc_repeat(seq: str, repeats: int = 4) -> bool:
    patterns = [
        "AT", "TA", "AC", "CA", "AG", "GA",
        "CT", "TC", "CG", "GC", "GT", "TG",
    ]
    return any((p * repeats) in seq for p in patterns)
