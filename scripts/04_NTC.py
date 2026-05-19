#!/usr/bin/env python3
"""
04_NTC.py (Enhanced Version for Publication)
Non-Target Control (NTC) gRNA Generator with Advanced Features
Core Upgrades:
1. GC content control (min-gc / max-gc, consistent with 03_filter_merge.py)
2. Genome-wide 0/1/2 mismatch hit count quantification
3. Optional standalone mismatch statistics table output
4. Strict NTC filtering (0 perfect matches to reference genome)

Outputs:
- NTC.fa: Validated NTC sequences (FASTA)
- NTC_statistics.tsv: Full NTC metadata (GC, mismatches, quality)
- NTC_mismatch_counts.tsv: Standalone mismatch hit table (--mismatch-info enabled)
"""
import argparse
import random
import subprocess
from pathlib import Path
from typing import List

def calculate_gc_content(seq: str) -> float:
    """Calculate GC content percentage of a sequence"""
    if len(seq) == 0:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return (gc / len(seq)) * 100

def run_bowtie_alignment(bowtie_index: str, temp_fa: Path, max_mismatch: int) -> int:
    """Run bowtie alignment and return total hit counts with up to max_mismatch errors"""
    cmd = [
        "bowtie",
        "-v", str(max_mismatch),
        "-a",
        "--norc",
        "--quiet",
        bowtie_index,
        str(temp_fa)
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    hits = len(result.stdout.strip().splitlines()) if result.stdout.strip() else 0
    return hits

def main():
    parser = argparse.ArgumentParser(
        description="Generate high-quality NTCs with matched GC distribution and mismatch quantification",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Basic Parameters
    parser.add_argument("--length", type=int, default=28, help="Length of NTC spacer sequences")
    parser.add_argument("--count", type=int, default=100, help="Number of NTC sequences to generate")
    parser.add_argument("--bowtie-index", required=True, help="Path to bowtie genome index prefix")
    parser.add_argument("--quiet", action="store_true", help="Run without progress messages")
    parser.add_argument("--min-gc", type=float, default=30, help="Minimum GC percentage of spacer")
    parser.add_argument("--max-gc", type=float, default=70, help="Maximum GC percentage of space")
    parser.add_argument("--out-fa", default="NTC.fa", help="Output NTC FASTA file")
    parser.add_argument("--out-stats", default="NTC_statistics.tsv", help="Output full statistics file")
    parser.add_argument("--mismatch-info", action="store_true", help="Output standalone NTC mismatch hit counts table")
    parser.add_argument("--out-mismatch", default="NTC_mismatch_counts.tsv", help="Mismatch counts output file")
    args = parser.parse_args()

    # Validate GC range
    if args.min_gc < 0 or args.max_gc > 100 or args.min_gc > args.max_gc:
        raise ValueError("Invalid GC range: min-gc must be ≤ max-gc (0-100)")

    # Path setup
    root = Path(__file__).parent.parent
    output_dir = root / "output"
    output_dir.mkdir(exist_ok=True)

    out_fa = output_dir / args.out_fa
    out_stats = output_dir / args.out_stats
    out_mismatch = output_dir / args.out_mismatch

    bases = ["A", "T", "C", "G"]
    ntc_list: List[dict] = []
    generated = 0

    print("=== Enhanced NTC Generator (Publication Ready) ===")
    print(f"Target NTCs: {args.count} | Length: {args.length} bp")
    print(f"GC Range: {args.min_gc}% - {args.max_gc}%")
    print(f"Bowtie Index: {args.bowtie_index}")
    print(f"Standalone Mismatch Table: {'ENABLED' if args.mismatch_info else 'DISABLED'}")
    print("==================================================\n")

    while generated < args.count:
        # Generate sequence with valid GC content
        while True:
            seq = "".join(random.choices(bases, k=args.length))
            gc = calculate_gc_content(seq)
            if args.min_gc <= gc <= args.max_gc:
                break

        seq_id = f"NTC_{generated + 1}"
        tmp = output_dir / f"tmp_{seq_id}.fa"

        # Write temp FASTA
        with open(tmp, "w") as f:
            f.write(f">{seq_id}\n{seq}\n")

        # Calculate mismatch hits
        total_0mm = run_bowtie_alignment(args.bowtie_index, tmp, 0)
        total_1mm = run_bowtie_alignment(args.bowtie_index, tmp, 1)
        total_2mm = run_bowtie_alignment(args.bowtie_index, tmp, 2)

        # Exact mismatch counts
        mm0 = total_0mm
        mm1 = total_1mm - total_0mm
        mm2 = total_2mm - total_1mm

        tmp.unlink(missing_ok=True)

        # Strict NTC rule: 0 perfect matches
        if mm0 == 0:
            ntc_list.append({
                "id": seq_id,
                "seq": seq,
                "gc": round(gc, 2),
                "mm0": mm0,
                "mm1": mm1,
                "mm2": mm2
            })
            generated += 1
            if not args.quiet:
                print(f"Generated {generated}/{args.count} | {seq_id} | GC: {gc:.1f}% | mm1: {mm1} | mm2: {mm2}")

    # ==================== Write Output Files ====================
    # 1. FASTA file
    with open(out_fa, "w") as f:
        for n in ntc_list:
            f.write(f">{n['id']} GC={n['gc']}%\n{n['seq']}\n")

    # 2. Full statistics table
    with open(out_stats, "w") as f:
        f.write("NTC_ID\tSequence\tGC_Pct\tmm0_hits\tmm1_hits\tmm2_hits\n")
        for n in ntc_list:
            f.write(f"{n['id']}\t{n['seq']}\t{n['gc']}\t{n['mm0']}\t{n['mm1']}\t{n['mm2']}\n")

    # 3. STANDALONE MISMATCH TABLE (only if --mismatch-info is used)
    if args.mismatch_info:
        with open(out_mismatch, "w") as f:
            f.write("NTC_ID\tmm0_hit_count\tmm1_hit_count\tmm2_hit_count\n")
            for n in ntc_list:
                f.write(f"{n['id']}\t{n['mm0']}\t{n['mm1']}\t{n['mm2']}\n")

    # ==================== Summary ====================
    avg_gc = sum(n["gc"] for n in ntc_list) / len(ntc_list)
    print("\n=== NTC Generation Complete ===")
    print(f"Total valid NTCs: {len(ntc_list)}")
    print(f"Average GC: {avg_gc:.2f}%")
    print(f"FASTA: {out_fa}")
    print(f"Full Stats: {out_stats}")
    if args.mismatch_info:
        print(f"Mismatch Count Table: {out_mismatch}")
    print("===============================")

if __name__ == "__main__":
    main()
