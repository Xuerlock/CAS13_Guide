# -*- coding: utf-8 -*-
"""
03_filter_merge.py
Cas13 28mer gRNA filtering
"""
import os
import re
import pandas as pd
import numpy as np
from pathlib import Path

# ===================== 0. 路径设置（适配你的实际路径） =====================
script_dir = Path(__file__).parent
input_dir = script_dir.parent / "output" / "bowtie"
candidate_file = script_dir.parent / "output" / "cas13_28mer_candidates.tsv"
output_dir = script_dir.parent / "output" / "filtered"
output_dir.mkdir(exist_ok=True, parents=True)

# ===================== 1. 读取bowtie counts文件（适配你的格式） =====================
print("=== 读取bowtie比对结果 ===")
def read_bowtie_counts(file_path):
    # 适配你的counts文件：第一列是数字，第二列是GuideID（可能含空格，用\s+分隔）
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

# ===================== 2. 读取candidates.tsv（关键修复：读表头+适配多空格分隔） =====================
print("=== 读取滑窗候选序列 ===")
cand = pd.read_csv(
    candidate_file,
    sep="\s+",  # 适配你的文件：多空格/制表符混合分隔
    header=0,   # 第一行是表头（guideId/geneId/start/spacer/gc/len）
    names=["GuideID", "GeneID", "start", "spacer", "gc", "len"],  # 强制列名统一
    dtype={
        "GuideID": str,
        "GeneID": str,
        "start": int,
        "spacer": str,
        "gc": float,
        "len": int
    },
    skip_blank_lines=True  # 跳过空行
)

# ===================== 3. 合并bowtie结果 =====================
print("=== 合并数据 ===")
# 合并3个错配counts
hit = cand.merge(mm0, on="GuideID", how="left").rename(columns={"count": "n0"})
hit = hit.merge(mm1, on="GuideID", how="left").rename(columns={"count": "n1"})
hit = hit.merge(mm2, on="GuideID", how="left").rename(columns={"count": "n2"})

# 填充NA为0
hit[["n0", "n1", "n2"]] = hit[["n0", "n1", "n2"]].fillna(0).astype(int)

# ===================== 4. 核心筛选：无脱靶 =====================
print("=== 筛选无脱靶gRNA ===")
# 0/1/2错配都必须=1（唯一靶向）
good = hit[(hit["n0"] == 1) & (hit["n1"] == 1) & (hit["n2"] == 1)].copy()
print(f"无脱靶gRNA数量：{len(good)}")

# ===================== 5. 序列质量筛选 =====================
print("=== 序列质量筛选 ===")
# 转大写
good["spacer"] = good["spacer"].str.upper()

# 5.1 GC 30-70%
good["pass_gc"] = (good["gc"] >= 30) & (good["gc"] <= 70)

# 5.2 无5个以上连续单碱基
def has_homopoly5(seq):
    return bool(re.search(r"A{5,}|T{5,}|C{5,}|G{5,}", seq))
good["has_homopoly5"] = good["spacer"].apply(has_homopoly5)

# 5.3 无二核苷酸重复≥4次
dinuc_pattern = r"(AT){4,}|(TA){4,}|(AC){4,}|(CA){4,}|(AG){4,}|(GA){4,}|(CT){4,}|(TC){4,}|(CG){4,}|(GC){4,}|(GT){4,}|(TG){4,}"
good["has_dinuc"] = good["spacer"].apply(lambda x: bool(re.search(dinuc_pattern, x)))

# 5.4 无N碱基
good["has_N"] = good["spacer"].str.contains("N")

# 5.5 位置≥30bp
good["pass_pos"] = good["start"] >= 30

# 应用所有筛选
final = good[
    (good["pass_gc"]) & 
    (~good["has_homopoly5"]) & 
    (~good["has_dinuc"]) & 
    (~good["has_N"]) & 
    (good["pass_pos"])
].copy()
print(f"通过所有质量筛选的gRNA数量：{len(final)}")

# ===================== 6. 打分 + 每条基因选TOP5 =====================
print("=== 挑选每条基因最优5条gRNA ===")
# 打分：GC越接近50分越高
final["gc_dev"] = abs(final["gc"] - 50)
final["score"] = -final["gc_dev"]  # 分数越高越好

# 按基因挑选，控制间距≥20bp
def pick_guides(sub_df, k=5, min_dist=20):
    # 按分数降序、起始位置升序排序
    sub_df = sub_df.sort_values(["score", "start"], ascending=[False, True]).reset_index(drop=True)
    picked = []
    for idx, row in sub_df.iterrows():
        if len(picked) >= k:
            break
        # 第一个直接选
        if not picked:
            picked.append(row)
        else:
            # 检查与已选gRNA的间距
            start_pos = row["start"]
            min_distance = min(abs(start_pos - p["start"]) for p in picked)
            if min_distance >= min_dist:
                picked.append(row)
    return pd.DataFrame(picked)

# 按GeneID分组挑选
res_list = []
for gene, group in final.groupby("GeneID"):
    picked = pick_guides(group)
    picked["rank"] = range(1, len(picked)+1)
    res_list.append(picked)

res = pd.concat(res_list, ignore_index=True)

# ===================== 7. 输出结果 =====================
print("=== 输出最终结果 ===")
# 选择输出列
out_cols = ["GeneID", "GuideID", "start", "spacer", "gc", "n0", "n1", "n2", "score", "rank"]
res[out_cols].to_csv(
    output_dir / "final_filtered_gRNA.tsv",
    sep="\t",
    index=False,
    encoding="utf-8"
)

# QC报告
qc = res.groupby("GeneID").size().reset_index(name="total_selected")
qc.to_csv(
    output_dir / "filter_qc_report.tsv",
    sep="\t",
    index=False,
    encoding="utf-8"
)

# 统计汇总
total_genes = len(cand["GeneID"].unique())
genes_with_gRNA = len(res["GeneID"].unique())
total_gRNAs = len(res)

print("\n==================== 筛选完成 ====================")
print(f"输入基因总数：{total_genes}")
print(f"成功获得gRNA的基因数：{genes_with_gRNA}")
print(f"最终筛选gRNA总数：{total_gRNAs}")
print("===================================================")
print(f"结果文件已保存至：{output_dir}")
