# -*- coding: utf-8 -*-
"""
03_filter_merge.py
Cas13 28mer gRNA filtering
升级：三维评分 + 位置优先 + 可选输出全部spacer + 可选最大位置限制
"""
import os
import re
import argparse
import pandas as pd
import numpy as np
from pathlib import Path

# ===================== 0. 命令行参数（新增你要的功能） =====================
parser = argparse.ArgumentParser(description="gRNA筛选脚本")
parser.add_argument("--keep-all-spacers", action="store_true", help="保留所有合并后的spacer信息（含全部字段）")
parser.add_argument("--max-position", type=int, default=None, help="仅保留该位置以内的spacer，例如 --max-position 1000")
parser.add_argument("--prioritize-position", action="store_true", help="开启位置优先模式（权重更高）")
args = parser.parse_args()

# ===================== 0. 路径设置（完全不变） =====================
script_dir = Path(__file__).parent
input_dir = script_dir.parent / "output" / "bowtie"
candidate_file = script_dir.parent / "output" / "cas13_28mer_candidates.tsv"
output_dir = script_dir.parent / "output" / "filtered"
output_dir.mkdir(exist_ok=True, parents=True)

# ===================== 1. 读取bowtie counts文件 =====================
print("=== 读取bowtie比对结果 ===")
def read_bowtie_counts(file_path):
    df = pd.read_csv(
        file_path, 
        sep="\s+", 
        header=None, 
        names=["count", "GuideID"],
        dtype={"count": int, "GuideID": str}
    )
    return df

mm0 = read_bowtie_counts(input_dir / "cas13_mm0_counts.txt")
mm1 = read_bowtie_counts(input_dir / "cas13_mm1_counts.txt")
mm2 = read_bowtie_counts(input_dir / "cas13_mm2_counts.txt")

# ===================== 2. 读取candidates.tsv =====================
print("=== 读取滑窗候选序列 ===")
cand = pd.read_csv(
    candidate_file,
    sep="\s+",
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

# ===================== 3. 合并bowtie结果 =====================
print("=== 合并数据 ===")
hit = cand.merge(mm0, on="GuideID", how="left").rename(columns={"count": "n0"})
hit = hit.merge(mm1, on="GuideID", how="left").rename(columns={"count": "n1"})
hit = hit.merge(mm2, on="GuideID", how="left").rename(columns={"count": "n2"})
hit[["n0", "n1", "n2"]] = hit[["n0", "n1", "n2"]].fillna(0).astype(int)

# ===================== 【可选】输出全部合并的spacer =====================
if args.keep_all_spacers:
    all_out = output_dir / "all_merged_spacers_with_counts.tsv"
    hit.to_csv(all_out, sep="\t", index=False)
    print(f"✅ 全部合并spacer信息已保存：{all_out}")

# ===================== 4. 核心筛选：无脱靶（你原来的规则） =====================
print("=== 筛选无脱靶gRNA ===")
good = hit[(hit["n0"] == 1) & (hit["n1"] == 1) & (hit["n2"] == 1)].copy()
print(f"无脱靶gRNA数量：{len(good)}")

# ===================== 5. 序列质量筛选（你原来的规则） =====================
print("=== 序列质量筛选 ===")
good["spacer"] = good["spacer"].str.upper()
good["pass_gc"] = (good["gc"] >= 30) & (good["gc"] <= 70)

def has_homopoly5(seq):
    return bool(re.search(r"A{5,}|T{5,}|C{5,}|G{5,}", seq))
good["has_homopoly5"] = good["spacer"].apply(has_homopoly5)

dinuc_pattern = r"(AT){4,}|(TA){4,}|(AC){4,}|(CA){4,}|(AG){4,}|(GA){4,}|(CT){4,}|(TC){4,}|(CG){4,}|(GC){4,}|(GT){4,}|(TG){4,}"
good["has_dinuc"] = good["spacer"].apply(lambda x: bool(re.search(dinuc_pattern, x)))
good["has_N"] = good["spacer"].str.contains("N")
good["pass_pos"] = good["start"] >= 30

# 应用筛选
final = good[
    (good["pass_gc"]) &
    (~good["has_homopoly5"]) &
    (~good["has_dinuc"]) &
    (~good["has_N"]) &
    (good["pass_pos"])
].copy()

# ===================== 【新增】可选：限制最大位置 =====================
if args.max_position is not None:
    final = final[final["start"] <= args.max_position].copy()
    print(f"✅ 已过滤位置 > {args.max_position} bp，剩余：{len(final)}")

print(f"通过所有质量筛选的gRNA数量：{len(final)}")

# ===================== 【新增】三维度评分：位置 + 质量 + 脱靶 =====================
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
    # 30-70 满分10，偏离50越远分越低
    dev = abs(gc - 50)
    return max(0, 10 - dev * 0.2)

def score_off_target(n0, n1, n2):
    # 无脱靶满分10，有脱靶直接扣分
    if n0 == 1 and n1 == 1 and n2 == 1:
        return 10
    else:
        return 0

# 计算三个分数
final["score_pos"] = final["start"].apply(score_position)
final["score_qual"] = final["gc"].apply(score_quality)
final["score_off"] = final.apply(lambda row: score_off_target(row["n0"], row["n1"], row["n2"]), axis=1)

# ===================== 【新增】加权总分 =====================
if args.prioritize_position:
    final["total_score"] = (
        0.6 * final["score_pos"] +
        0.2 * final["score_qual"] +
        0.2 * final["score_off"]
    )
    print("✅ 已开启【位置优先】模式")
else:
    final["total_score"] = (
        0.3 * final["score_pos"] +
        0.4 * final["score_qual"] +
        0.3 * final["score_off"]
    )
    print("✅ 默认【平衡模式】")

# ===================== 6. 按基因挑选TOP5（改用总分排序） =====================
print("=== 挑选每条基因最优5条gRNA ===")

def pick_guides(sub_df, k=5, min_dist=20):
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

# ===================== 7. 输出结果 =====================
print("=== 输出最终结果 ===")
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

# ===================== 汇总 =====================
print("\n==================== 筛选完成 ====================")
print(f"输入基因总数：{len(cand['GeneID'].unique())}")
print(f"成功获得gRNA的基因数：{len(res['GeneID'].unique())}")
print(f"最终筛选gRNA总数：{len(res)}")
print(f"结果文件：{output_dir}")
print("===================================================")
