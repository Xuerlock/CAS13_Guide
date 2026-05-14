# CAS13_Guide: Genome-scale Cas13 gRNA Design Tookit
A robust, user-friendly pipeline for designing high-specificity CRISPR-Cas13 gRNAs, including sliding window candidate generation, off-target filtering, multi-criteria optimization, and non-target control (NTC) sequence generation.

## Overview
This pipeline automates the entire process of Cas13 guide RNA design:

1. Extract high-quality spacer candidates from input sequences using a sliding window approach.
2. Perform genome-wide off-target analysis using Bowtie.
3. Filter and rank spacers by sequence quality, position efficiency, and off-target specificity.
4. Generate validated non-target control (NTC) sequences with no matches in the reference genome.

## Highlight
* Sliding-window generation of Cas13 spacer candidates from transcript or CDS FASTA
* Sequence-level QC filtering (GC, homopolymer, low complexity)
* Bowtie-based mismatch counting against genome or transcriptome indexes
* Strict off-target filtering using mm0/mm1/mm2 uniqueness rules
* Position-aware selection of a designate number of guides per gene
* Outputs ready for downstream library synthesis

## Pipeline Workflow
1. **Candidate Generation**: Sliding window scans input FASTA sequences to generate potential gRNA spacers.
2. **Off-target Detection**: Bowtie aligns candidates to the reference genome to identify off-target matches.
3. **Filtering & Ranking**: Multi-criteria scoring (position, quality, off-target) filters and ranks optimal gRNAs.
4. **NTC Generation**: Creates validated control sequences with no perfect matches in the reference genome.

## File Structure
```
root/
├── scripts/
│   ├── 01_windowslide.py       # Spacer candidate generation
│   ├── 02_bowtie.sh            # Bowtie alignment for off-target checks
│   ├── 03_filter_merge.py      # Filtering and ranking of spacers
│   └── 04_NTC.py               # NTC sequence generation
├── output/                     # Auto-generated (all results stored here)
│   ├── bowtie/                 # Bowtie alignment counts
│   ├── filtered/               # Final filtered/ranked gRNAs
│   └── NTC.fa                  # NTC sequences (if generated)
└── README.md
```


## Prerequisites
### Software
```
Python 3.7 or higher
Bowtie 1.x (not Bowtie 2)
Bash shell (Linux / macOS / WSL)
```

### Python Dependencies
```{bash}
pip install pandas numpy
```

## Step-by-Step Usage

### 1. Generate gRNA Candidates

Extract candidate spacers from input FASTA using a sliding window approach.

**Script**: 01_windowslide.py

**Usage**：

```{bash}
cd scripts/
python 01_windowslide.py -i input_genes.fa -o cas13
```

**Parameters**:

-i/--input: Path to input FASTA file containing target genes.

-o/--output-prefix: Prefix for output candidate file (saved to ../output/).

**Output**: 


../output/cas13_candidates.fa (FASTA file of spacer candidates with position metadata).


### 2. Off-Target Analysis with Bowtie

Build a reference genome index and align candidates to detect off-targets.

**Script**: 02_bowtie.sh

**Usage**：

```{bash}
cd scripts/
bash 02_bowtie.sh reference_genome.fa ../output/cas13_candidates.fa cas13
```

**Parameters**:

Path to reference genome FASTA (for Bowtie index building).

Path to candidate spacer FASTA (output from 01_windowslide.py).

Prefix for Bowtie output files (saved to ../output/bowtie/).

**Output**: 

../output/bowtie/cas13_index.*: Bowtie reference index files.

../output/bowtie/cas13_alignments.txt: Off-target alignment results.

### 3. Filter and Select Top gRNAs

Merge alignment results and select optimal spacers using multi-criteria scoring.

**Script**: 03_filter_merge.py

**Usage**：

```{bash}
cd scripts/
python 03_filter_merge.py --top-n 5 --max-position 1000 --prioritize-position
```

**Parameters**:

--top-n: Number of top-ranked gRNAs to retain per gene (default: 5).

--max-position: Maximum allowable start position of spacer (default: 1000).

--prioritize-position: Enable position-prioritized scoring mode (default: standard mode).

--min-gc: Minimum GC percentage of spacer (default: 30).

--max-gc: Maximum GC percentage of spacer (default: 70).

**Scoring System**:

| Category          | Standard Mode Weight | Position-Prioritized Mode Weight | Scoring Rules                                                   |
|--------------------|----------------------|-----------------------------------|-----------------------------------------------------------------|
| Position Score     | 30%                  | 60%                               | ≤500 bp = 10; ≤1000 bp = 8; ≤1500 bp = 5; >1500 bp = 2          |
| Quality Score      | 40%                  | 20%                               | Based on GC proximity to 50% (0-10 scale)                      |
| Off-Target Score   | 30%                  | 20%                               | Perfect specificity = 10; Any off-targets = 0                  |

**Output**: 

../output/filtered/final_filtered_gRNA.tsv: Ranked gRNAs with full scoring metadata.

../output/filtered/filter_qc_report.tsv: QC report (valid gRNA count per gene).

### 4. Generate Non-Target Controls (NTCs)

Create random sequences with zero perfect matches in the reference genome.

**Script**: 04_NTC.py

**Usage**：

```{bash}
cd scripts/
python 04_NTC.py --bowtie-index ../output/bowtie/bowtie_index --count 100 --length 28 --quiet
```

**Parameters**:

--bowtie-index: Path to pre-built Bowtie index (from 02_bowtie.sh).

--count: Number of NTC sequences to generate (default: 100).

--length: Length of NTC sequences (matches Cas13 spacer length, default: 28).

--quiet: Suppress progress messages (optional).


**Output**：

../output/NTC.fa (FASTA file of validated NTC sequences).

## Output File Details

**Final gRNAs (output/filtered/final_filtered_gRNA.tsv)**

| Column                     | Description                                                                 |
|----------------------------|-----------------------------------------------------------------------------|
| GeneID                     | ID of target gene from input FASTA                                          |
| GuideID                    | Unique identifier for the gRNA (GeneID + position)                          |
| start position             | Start position of spacer on target gene                                     |
| spacer sequence            | 28 bp Cas13 gRNA spacer sequence                                            |
| GC%                        | GC percentage of spacer sequence                                            |
| n0/n1/n2 mismatch counts   | Number of perfect (n0), 1-mismatch (n1), 2-mismatch (n2) off-targets        |
| position score             | Score based on spacer position (0–10)                                        |
| quality score              | Score based on GC content (0–10)                                            |
| off-target score           | Score based on specificity (0 or 10)                                         |
| total score                | Weighted sum of position + quality + off-target scores (0–10)               |
| rank                       | Rank of gRNA within the target gene (1 = highest)                            |


**QC Report (output/filtered/filter_qc_report.tsv)**

| Column       | Description                                              |
|--------------|----------------------------------------------------------|
| GeneID       | ID of target gene from input FASTA                       |
| valid_gRNAs  | Number of gRNAs passing all filter criteria              |


**NTC Sequences (output/NTC.fa)**

FASTA-formatted file with validated non-target control sequences:

```
>NTC_001
GGTCTTTCTGGTTGGACGGTTGGCGG
>NTC_002
GGACGGTTGGCGGAGAAGAGGCCGGT
>NTC_003
GGCCGGTCAGACGGAGAAGCTGCGGG
......
```

## Troubleshooting

| Issue                      | Solution                                                                 |
|----------------------------|--------------------------------------------------------------------------|
| bowtie not found           | Install Bowtie 1.x (not Bowtie 2) and add to system PATH.                |
| No gRNAs pass filters      | Relax --min-gc/--max-gc thresholds or increase --max-position.           |
| NTC generation slow        | Reduce --count or use a transcriptome index (instead of full genome).    |
| Path errors (e.g., "file not found") | Always run scripts from the scripts/ directory (per file structure). |
| Low number of candidates   | Check input FASTA format (ensure valid gene sequences, no extra whitespace). |

## License
This pipeline is available for academic and research use only.

## Citation
If you use this pipeline in your work, please cite this repository.

