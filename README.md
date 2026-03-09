# CAS13-Guide

A lightweight command-line toolkit for genome-scale Cas13 guide RNA design in designated species.

## Features

- sliding-window generation of Cas13 spacer candidates from transcript or CDS FASTA
- sequence-level QC filtering (GC, homopolymer, low complexity)
- bowtie-based mismatch counting against genome or transcriptome indexes
- strict off-target filtering using mm0/mm1/mm2 uniqueness rules
- position-aware selection of a fixed number of guides per gene
- outputs ready for downstream library synthesis

## Installation

```bash
git clone https://github.com/Xuerlock/CAS13_Guide.git
cd CAS13_Guide
python -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -e .


#External dependency
#CAS13-Guide requires Bowtie 1 in PATH.
#Check installation:
bowtie --version
bowtie-build --version

