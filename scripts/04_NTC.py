#!/usr/bin/env python3
"""
04_NTC.py
Generate non-target control (NTC) gRNAs with random sequences and genome-wide alignment check.
Output: NTC.fa (FASTA format)
"""
import argparse
import random
import subprocess
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Generate NTC gRNAs with no genome-wide matches.")
    parser.add_argument("--length", type=int, default=28, help="Length of gRNA sequences (default: 28)")
    parser.add_argument("--count", type=int, default=100, help="Number of NTC sequences to generate (default: 100)")
    parser.add_argument("--bowtie-index", required=True, help="Path to bowtie index prefix")
    parser.add_argument("--output", default="NTC.fa", help="Output FASTA file path (default: ./NTC.fa)")
    parser.add_argument("--quiet", action="store_true", help="Quiet mode: suppress per-sequence progress messages")
    args = parser.parse_args()

    # Dynamic paths
    script_dir = Path(__file__).parent
    root_dir = script_dir.parent
    output_path = Path(args.output)
    if not output_path.is_absolute():
        output_path = root_dir / "output" / output_path
    output_path.parent.mkdir(exist_ok=True, parents=True)

    bases = ["A", "T", "C", "G"]
    generated = 0
    ntc_seqs = []

    print(f"Generating {args.count} NTC sequences (length={args.length}) with no genome matches...")
    while generated < args.count:
        # Generate random sequence
        seq = "".join(random.choices(bases, k=args.length))
        seq_id = f"NTC_{generated+1}"
        temp_fa = root_dir / "tmp_ntc_check.fa"
        with open(temp_fa, "w") as f:
            f.write(f">{seq_id}\n{seq}\n")

        # Bowtie check: 0 mismatches, no alignments allowed
        cmd = [
            "bowtie",
            "-v", "0", "-a", "--norc", "--quiet",
            args.bowtie_index,
            str(temp_fa)
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if not result.stdout.strip():  # No alignment output → valid NTC
            ntc_seqs.append((seq_id, seq))
            generated += 1
            if not args.quiet:
                print(f"Valid NTC found: {generated}/{args.count}")

        temp_fa.unlink(missing_ok=True)

    # Write final output
    with open(output_path, "w") as f:
        for seq_id, seq in ntc_seqs:
            f.write(f">{seq_id}\n{seq}\n")

    print(f"\nDone! NTC sequences saved to: {output_path}")

if __name__ == "__main__":
    main()
