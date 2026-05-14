# -*- coding: utf-8 -*-
"""
03_filter_merge.py
Cas13 28mer gRNA filtering
"""
import os
import re
import argparse
import pandas as pd
import numpy as np
from pathlib import Path

# ===================== Parameter Setting=====================
parser = argparse.ArgumentParser(description="gRNA filtering script")
parser.add_argument("--keep-all-spacers", action="store_true", help="Retain all merged spacer information without deduplication")
parser.add_argument("--max-position", type=int, default=None, help="Only keep spacers located within the specified range from the mRNA start codon. For example, --max-position 1000 retains only spacers within the first 1000 nucleotides.")
parser.add_argument("--prioritize-position", action="store_true", help="Enable position-priority mode (higher weight given to positional criteria)")
parser.add_argument("--top-n", type=int, default=5, help="Number of top-ranked gRNAs to retain per gene (default: 5)")
args = parser.parse_args()

# ===================== 0. Path Configuration=====================
script_dir = Path(__file__).parent
input_dir = script_dir.parent / "output" / "bowtie"
candidate_file = script_dir.parent / "output" / "cas13_candidates.tsv"
output_dir = script_dir.parent / "output" / "filtered"
output_dir.mkdir(exist_ok=True, parents=True)

# ===================== 1. Read Bowtie Count Files =====================
print("=== Reading Bowtie alignment results ===")
def read_bowtie_counts(file_path):
    df = pd.read_csv(
        file_path, 
        sep=r"\s+", 
        header=None, 
        names=["count", "GuideID"],
        dtype={"count": int, "GuideID": str}
    )
    return df

mm0 = read_bowtie_counts(input_dir / "cas13_mm0_counts.txt")
mm1 = read_bowtie_counts(input_dir / "cas13_mm1_counts.txt")
mm2 = read_bowtie_counts(input_dir / "cas13_mm2_counts.txt")

# ===================== 2. Read Candidate Sequences =====================
print("=== Reading sliding window candidate sequences ===")
cand = pd.read_csv(
    candidate_file,
    sep=r"\s+", 
    header=0,
    names=["GuideID", "GeneID", "start", "spacer", "gc", "len"],
    dtype={
        "GuideID": str,
        "GeneID": str,
        "start": int,
        "spacer": str,
        "gc": float,
        "len": int
    },
    skip_blank_lines=True
)

# ===================== 3. Merge Bowtie Results =====================
print("=== Merging data ===")
hit = cand.merge(mm0, on="GuideID", how="left").rename(columns={"count": "n0"})
hit = hit.merge(mm1, on="GuideID", how="left").rename(columns={"count": "n1"})
hit = hit.merge(mm2, on="GuideID", how="left").rename(columns={"count": "n2"})
hit[["n0", "n1", "n2"]] = hit[["n0", "n1", "n2"]].fillna(0).astype(int)

# ===================== Optional: Output All Merged Spacers =====================
if args.keep_all_spacers:
    all_out = output_dir / "all_merged_spacers_with_counts.tsv"
    hit.to_csv(all_out, sep="\t", index=False)
    print(f"All merged spacer information saved to: {all_out}")

# ===================== 4. Core Filtering: No Off-Target=====================
print("=== Filtering off-target-free gRNAs ===")
good = hit[(hit["n0"] == 1) & (hit["n1"] == 1) & (hit["n2"] == 1)].copy()
print(f"Number of off-target-free gRNAs: {len(good)}")

# ===================== 5. Sequence Quality Filtering =====================
print("=== Applying sequence quality filters ===")
good["spacer"] = good["spacer"].str.upper()
good["pass_gc"] = (good["gc"] >= 30) & (good["gc"] <= 70)

def has_homopoly5(seq):
    return bool(re.search(r"A{5,}|T{5,}|C{5,}|G{5,}", seq)) 
good["has_homopoly5"] = good["spacer"].apply(has_homopoly5)

dinuc_pattern = r"(AT){4,}|(TA){4,}|(AC){4,}|(CA){4,}|(AG){4,}|(GA){4,}|(CT){4,}|(TC){4,}|(CG){4,}|(GC){4,}|(GT){4,}|(TG){4,}"
good["has_dinuc"] = good["spacer"].apply(lambda x: bool(re.search(dinuc_pattern, x)))
good["has_N"] = good["spacer"].str.contains("N")
good["pass_pos"] = good["start"] >= 30

final = good[
    (good["pass_gc"]) &
    (~good["has_homopoly5"]) &
    (~good["has_dinuc"]) &
    (~good["has_N"]) &
    (good["pass_pos"])
].copy()

# ===================== Optional: Restrict Maximum Position =====================
if args.max_position is not None:
    final = final[final["start"] <= args.max_position].copy()
    print(f"Filtered positions to ≤ {args.max_position} bp, remaining: {len(final)}")

print(f"Number of gRNAs passing all quality filters: {len(final)}")

# ===================== Three-Dimensional Scoring =====================
def score_position(pos):
    if pos <= 500:
        return 10
    elif pos <= 1000:
        return 8
    elif pos <= 1500:
        return 5
    else:
        return 2

def score_quality(gc):
    dev = abs(gc - 50)
    return max(0, 10 - dev * 0.2)

def score_off_target(n0, n1, n2):
    if n0 == 1 and n1 == 1 and n2 == 1:
        return 10
    else:
        return 0

# Calculate three scores
final["score_pos"] = final["start"].apply(score_position)
final["score_qual"] = final["gc"].apply(score_quality)
final["score_off"] = final.apply(lambda row: score_off_target(row["n0"], row["n1"], row["n2"]), axis=1)

# ===================== Weighted Total Score =====================
if args.prioritize_position:
    final["total_score"] = (
        0.6 * final["score_pos"] +
        0.2 * final["score_qual"] +
        0.2 * final["score_off"]
    )
    print("Using Position-Priority scoring mode")
else:
    final["total_score"] = (
        0.3 * final["score_pos"] +
        0.4 * final["score_qual"] +
        0.3 * final["score_off"]
    )
    print("Using default Balanced scoring mode")

# ===================== Select Top N Guides Per Gene =====================
print("=== Selecting top-ranked gRNAs per gene ===")

def pick_guides(sub_df, k=args.top_n, min_dist=20):
    sub_df = sub_df.sort_values(["total_score", "start"], ascending=[False, True]).reset_index(drop=True)
    picked = []
    for idx, row in sub_df.iterrows():
        if len(picked) >= k:
            break
        if not picked:
            picked.append(row)
        else:
            min_distance = min(abs(row["start"] - p["start"]) for p in picked)
            if min_distance >= min_dist:
                picked.append(row)
    return pd.DataFrame(picked)

res_list = []
for gene, group in final.groupby("GeneID"):
    picked = pick_guides(group)
    picked["rank"] = range(1, len(picked)+1)
    res_list.append(picked)

res = pd.concat(res_list, ignore_index=True)

# ===================== Output Final Results  =====================
print("=== Writing final results ===")
out_cols = ["GeneID", "GuideID", "start", "spacer", "gc", "n0", "n1", "n2",
            "score_pos", "score_qual", "score_off", "total_score", "rank"]

res[out_cols].to_csv(
    output_dir / "final_filtered_gRNA.tsv",
    sep="\t",
    index=False,
    encoding="utf-8"
)

qc = res.groupby("GeneID").size().reset_index(name="total_selected")
qc.to_csv(output_dir / "filter_qc_report.tsv", sep="\t", index=False)

# ===================== Summary =====================
print("\n====== Filtering Complete ======")
print(f"Total input genes: {len(cand['GeneID'].unique())}")
print(f"Genes with successful gRNAs: {len(res['GeneID'].unique())}")
print(f"Total filtered gRNAs:{len(res)}")
print(f"Top {args.top_n} gRNAs retained per gene")
print(f"Results directory: {output_dir}")
print("===================================================")
