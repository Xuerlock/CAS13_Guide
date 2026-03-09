from __future__ import annotations
import click
from candidates import generate_candidates, write_candidate_outputs
from bowtie import count_mismatches
from counts import merge_counts
from guide_select import select_guides


@click.group()
def main() -> None:
    """CAS13-Guide: microbial Cas13 gRNA design toolkit."""


@main.command("generate")
@click.option("--input", "input_fasta", required=True, type=click.Path(exists=True))
@click.option("--output-prefix", required=True, type=click.Path())
@click.option("--guide-len", default=28, show_default=True, type=int)
@click.option("--window-start", default=30, show_default=True, type=int)
@click.option("--window-end", default=200, show_default=True, type=int)
@click.option("--gc-min", default=40.0, show_default=True, type=float)
@click.option("--gc-max", default=60.0, show_default=True, type=float)
def generate_cmd(input_fasta: str, output_prefix: str, guide_len: int, window_start: int, window_end: int, gc_min: float, gc_max: float) -> None:
    df = generate_candidates(
        fasta_path=input_fasta,
        guide_len=guide_len,
        window_start=window_start,
        window_end=window_end,
        gc_min=gc_min,
        gc_max=gc_max,
    )
    tsv, fa = write_candidate_outputs(df, output_prefix)
    click.echo(f"Wrote {tsv}")
    click.echo(f"Wrote {fa}")


@main.command("count")
@click.option("--fasta", required=True, type=click.Path(exists=True))
@click.option("--index", required=True, type=click.Path())
@click.option("--output-prefix", required=True, type=click.Path())
def count_cmd(fasta: str, index: str, output_prefix: str) -> None:
    mm0, mm1, mm2 = count_mismatches(fasta=fasta, index=index, output_prefix=output_prefix)
    click.echo(f"Wrote {mm0}")
    click.echo(f"Wrote {mm1}")
    click.echo(f"Wrote {mm2}")


@main.command("merge")
@click.option("--candidates", required=True, type=click.Path(exists=True))
@click.option("--mm0", required=True, type=click.Path(exists=True))
@click.option("--mm1", required=True, type=click.Path(exists=True))
@click.option("--mm2", required=True, type=click.Path(exists=True))
@click.option("--output", required=True, type=click.Path())
def merge_cmd(candidates: str, mm0: str, mm1: str, mm2: str, output: str) -> None:
    df = merge_counts(candidates, mm0, mm1, mm2)
    df.to_csv(output, sep="\t", index=False)
    click.echo(f"Wrote {output}")


@main.command("select")
@click.option("--merged", required=True, type=click.Path(exists=True))
@click.option("--output-prefix", required=True, type=click.Path())
@click.option("--per-gene", default=5, show_default=True, type=int)
@click.option("--strict-mm0", default=1, show_default=True, type=int)
@click.option("--strict-mm1", default=1, show_default=True, type=int)
@click.option("--strict-mm2", default=1, show_default=True, type=int)
@click.option("--min-dist", default=25, show_default=True, type=int)
def select_cmd(merged: str, output_prefix: str, per_gene: int, strict_mm0: int, strict_mm1: int, strict_mm2: int, min_dist: int) -> None:
    import pandas as pd
    df = pd.read_csv(merged, sep="\t")
    final_df, qc, lt5 = select_guides(
        merged_df=df,
        per_gene=per_gene,
        strict_mm0=strict_mm0,
        strict_mm1=strict_mm1,
        strict_mm2=strict_mm2,
        min_dist=min_dist,
    )
    final = f"{output_prefix}.final.tsv"
    qc_out = f"{output_prefix}.qc.tsv"
    lt5_out = f"{output_prefix}.genes_lt5.tsv"
    final_df.to_csv(final, sep="\t", index=False)
    qc.to_csv(qc_out, sep="\t", index=False)
    lt5.to_csv(lt5_out, sep="\t", index=False)
    click.echo(f"Wrote {final}")
    click.echo(f"Wrote {qc_out}")
    click.echo(f"Wrote {lt5_out}")
