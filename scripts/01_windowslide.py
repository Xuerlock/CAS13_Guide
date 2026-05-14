# -*- coding: utf-8 -*-
"""
Sliding window tool for CRISPR spacer candidate generation
Supports:
1. Filter spacers by GC content, homopolymer, N content
2. Generate TSV with spacer metadata
3. Generate CRISPOR-compatible batch FASTA
4. Generate bowtie-compatible individual FASTA (candidates.fa)
"""
from __future__ import print_function
import sys
import argparse
import os
from pathlib import Path
from typing import Generator, List, Tuple


def parse_fasta(path: str) -> Generator[Tuple[str, str], None, None]:
    """
    Parse FASTA file and yield (sequence_id, uppercase_sequence) pairs
    
    Args:
        path: Path to FASTA file
    
    Yields:
        Tuple of (sequence ID, uppercase sequence string)
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"FASTA file not found: {path}")
    
    name = None
    seq_chunks = []
    with open(path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_chunks).upper()
                # Extract first field as sequence ID (ignore comments)
                name = line[1:].split()[0]
                seq_chunks = []
            else:
                # Validate sequence characters (basic check)
                if not all(c in 'ATCGNatcgn' for c in line):
                    sys.stderr.write(f"Warning: Invalid characters in line {line_num} of {path}\n")
                seq_chunks.append(line)
        # Yield the last sequence
        if name is not None:
            yield name, "".join(seq_chunks).upper()


def gc_percent(seq: str) -> float:
    """
    Calculate GC percentage of a DNA sequence
    
    Args:
        seq: DNA sequence string (uppercase)
    
    Returns:
        GC percentage (0-100)
    """
    if not seq:
        return 0.0
    gc_count = sum(1 for c in seq if c in ('G', 'C'))
    return 100.0 * gc_count / len(seq)


def has_homopolymer(seq: str, k: int = 5) -> bool:
    """
    Check if sequence contains homopolymer runs of length >= k
    
    Args:
        seq: DNA sequence string (uppercase)
        k: Minimum homopolymer run length to detect
    
    Returns:
        True if homopolymer run found, False otherwise
    """
    if k < 2 or len(seq) < k:
        return False
    
    run_length = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i-1]:
            run_length += 1
            if run_length >= k:
                return True
        else:
            run_length = 1
    return False


def filter_spacers(
    fasta_path: str,
    guide_len: int = 28,
    step: int = 1,
    gc_min: float = 30.0,
    gc_max: float = 70.0,
    homopolymer_k: int = 5,
    max_per_gene: int = 0
) -> List[Tuple[str, str, str, int, float, int]]:
    """
    Filter spacers from FASTA sequences using sliding window approach
    
    Args:
        fasta_path: Path to input FASTA file
        guide_len: Length of spacer to extract
        step: Sliding window step size
        gc_min: Minimum GC percentage (0-100)
        gc_max: Maximum GC percentage (0-100)
        homopolymer_k: Filter homopolymers >= k (0 to disable)
        max_per_gene: Max spacers to keep per gene (0 to disable)
    
    Returns:
        List of tuples: (guide_id, gene_id, start_pos, spacer_seq, gc_percent, spacer_len)
    """
    if guide_len < 1:
        raise ValueError("Guide length must be positive integer")
    if step < 1:
        raise ValueError("Step size must be positive integer")
    if not (0 <= gc_min <= gc_max <= 100):
        raise ValueError("Invalid GC range: 0 ≤ gc_min ≤ gc_max ≤ 100")
    
    filtered_spacers = []
    
    for gene_id, seq in parse_fasta(fasta_path):
        seq_len = len(seq)
        if seq_len < guide_len:
            sys.stderr.write(f"Warning: {gene_id} sequence length ({seq_len}) < guide length ({guide_len}) - skipped\n")
            continue
        
        kept_count = 0
        # Sliding window iteration
        for start in range(0, seq_len - guide_len + 1, step):
            spacer = seq[start:start+guide_len]
            
            # Filter out spacers with N
            if "N" in spacer:
                continue
            
            # Filter by GC content
            gc = gc_percent(spacer)
            if not (gc_min <= gc <= gc_max):
                continue
            
            # Filter by homopolymer
            if homopolymer_k > 0 and has_homopolymer(spacer, k=homopolymer_k):
                continue
            
            # Generate guide ID
            guide_id = f"{gene_id}__{start}"
            filtered_spacers.append((guide_id, gene_id, start, spacer, gc, len(spacer)))
            
            # Enforce max per gene
            kept_count += 1
            if max_per_gene > 0 and kept_count >= max_per_gene:
                break
    
    return filtered_spacers


def write_tsv(output_path: str, spacers: List[Tuple[str, str, int, str, float, int]]) -> None:
    """
    Write spacer metadata to TSV file
    
    Args:
        output_path: Path to output TSV file
        spacers: List of filtered spacer tuples
    """
    with open(output_path, 'w') as f:
        # Write header
        f.write("guideId\tgeneId\tstart\tspacer\tgc\tlen\n")
        # Write data rows
        for guide_id, gene_id, start, spacer, gc, length in spacers:
            f.write(f"{guide_id}\t{gene_id}\t{start}\t{spacer}\t{gc:.2f}\t{length}\n")
    
    print(f"Successfully wrote TSV: {output_path}")


def write_crispor_batch_fasta(
    output_path: str,
    spacers: List[Tuple[str, str, int, str, float, int]],
    batch_size: int = 80,
    prefix: str = "GGG",
    sep_n: str = "NNNN"
) -> None:
    """
    Write CRISPOR-compatible batch FASTA file
    
    Args:
        output_path: Path to output FASTA file
        spacers: List of filtered spacer tuples
        batch_size: Number of guides per batch record
        prefix: Prefix to add before each guide
        sep_n: Separator between guides in batch
    """
    if not spacers:
        sys.stderr.write("Warning: No spacers to write to CRISPOR batch FASTA\n")
        return
    
    with open(output_path, 'w') as f:
        batch_idx = 1
        idx = 0
        total_spacers = len(spacers)
        
        while idx < total_spacers:
            # Get current batch chunk
            chunk = spacers[idx:idx+batch_size]
            # Write batch header
            f.write(f">batch_{batch_idx}\n")
            # Build batch sequence
            seq_pieces = [prefix + spacer for _, _, _, spacer, _, _ in chunk]
            batch_seq = sep_n.join(seq_pieces)
            f.write(f"{batch_seq}\n")
            
            # Increment counters
            batch_idx += 1
            idx += batch_size
    
    print(f"Successfully wrote CRISPOR batch FASTA: {output_path}")


def write_bowtie_fasta(output_path: str, spacers: List[Tuple[str, str, int, str, float, int]]) -> None:
    """
    Write bowtie-compatible FASTA file (one record per spacer)
    
    Args:
        output_path: Path to output FASTA file
        spacers: List of filtered spacer tuples
    """
    if not spacers:
        sys.stderr.write("Warning: No spacers to write to bowtie FASTA\n")
        return
    
    with open(output_path, 'w') as f:
        for guide_id, _, _, spacer, _, _ in spacers:
            f.write(f">{guide_id}\n")
            f.write(f"{spacer}\n")
    
    print(f"Successfully wrote bowtie-compatible FASTA: {output_path}")


def main():
    """Main function for command line execution"""
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Sliding window tool for CRISPR spacer candidate generation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python 01windowslide.py -i input_genes.fa -o cas13
  python 01windowslide.py -i input_genes.fa -o cas13 --guide_len 30 --gc_min 25 --gc_max 75
  python 01windowslide.py -i input_genes.fa -o cas13 --max_per_gene 10 --batch_size 100
        """
    )
    parser.add_argument("-i", "--in_fa", required=True, help="Input FASTA file (one record per gene)")
    parser.add_argument("-o", "--out_prefix", default="cas13", help="Output file prefix (default: cas13)")
    parser.add_argument("--guide_len", type=int, default=28, help="Spacer length (default: 28)")
    parser.add_argument("--step", type=int, default=1, help="Sliding window step size (default: 1)")
    parser.add_argument("--gc_min", type=float, default=30.0, help="Minimum GC percentage (default: 30.0)")
    parser.add_argument("--gc_max", type=float, default=70.0, help="Maximum GC percentage (default: 70.0)")
    parser.add_argument("--no_homopolymer_k", type=int, default=5, help="Filter homopolymer runs >=k (0 to disable, default:5)")
    parser.add_argument("--max_per_gene", type=int, default=0, help="Max candidates per gene (0 to disable, default:0)")
    parser.add_argument("--batch_size", type=int, default=80, help="Guides per CRISPOR batch (default:80)")
    parser.add_argument("--sepN", default="NNNN", help="Separator between guides in CRISPOR batch (default:NNNN)")
    parser.add_argument("--prefix", default="GGG", help="Prefix before each guide in CRISPOR batch (default:GGG)")
    
    args = parser.parse_args()
    
    try:
        # Step 1: Filter spacers using sliding window
        print("Starting spacer filtering...")
        spacers = filter_spacers(
            fasta_path=args.in_fa,
            guide_len=args.guide_len,
            step=args.step,
            gc_min=args.gc_min,
            gc_max=args.gc_max,
            homopolymer_k=args.no_homopolymer_k,
            max_per_gene=args.max_per_gene
        )
        
        if not spacers:
            sys.stderr.write("Error: No spacers passed all filters! Check your input and parameters.\n")
            sys.exit(1)
        
        # Step 2: Define output paths
        script_dir = Path(__file__).parent
        output_dir = script_dir.parent / "output"
        output_dir.mkdir(exist_ok=True, parents=True)
                
        tsv_path = os.path.join(output_dir, f"{args.out_prefix}_candidates.tsv")
        crispor_fa_path = os.path.join(output_dir, f"{args.out_prefix}_spacer_candidate.fa")
        bowtie_fa_path = os.path.join(output_dir, f"{args.out_prefix}_candidates.fa")
        
        # Step 3: Write output files
        write_tsv(tsv_path, spacers)
        write_crispor_batch_fasta(crispor_fa_path, spacers, args.batch_size, args.prefix, args.sepN)
        write_bowtie_fasta(bowtie_fa_path, spacers)
        
        # Step 4: Print summary
        print("\n=== Summary ===")
        print(f"Total valid spacers: {len(spacers)}")
        print(f"TSV file: {tsv_path}")
        print(f"CRISPOR batch FASTA: {crispor_fa_path}")
        print(f"Bowtie-compatible FASTA: {bowtie_fa_path}")
        
    except Exception as e:
        sys.stderr.write(f"Error during execution: {str(e)}\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
