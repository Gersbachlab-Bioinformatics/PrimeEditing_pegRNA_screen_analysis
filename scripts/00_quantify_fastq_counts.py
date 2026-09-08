#!/usr/bin/env python3
"""Quantify per-oligo read counts and reporter editing rates from demultiplexed
PE HCR-FlowFISH screen FASTQs, for every bin/replicate sample in a directory.

This is the step upstream of `01_generate_count_table.R`, which currently
assumes reads have already been resolved to a single best-matching oligo_id
and (for bulk samples) reporter editing has already been quantified. This
script does that resolution directly: read 1 carries the pegRNA spacer, read
2 carries a 7 bp barcode followed by the reporter; a read is kept only if
its independently-matched spacer (by substring search in read 1) and its
barcode-matched library row (by exact 7 bp prefix match) agree ("perfect
match"). For each perfect-matched read, the barcode is stripped from read 2
and the remainder is checked against the library's unedited vs. edited
reporter sequence to call per-read editing status.

Generalized from an interactive, single-sample, hardcoded-path notebook
(`Demultiplex_redo/claude_files_PE/030626_pegSTAG1_screen.ipynb`) that did
this for one STAG1 bulk replicate. Assumes FAM120A's and SV2A's full
PRIDICT2.0-output library CSVs have the same column schema as STAG1's
(this was not independently verified locally -- only a reduced 5-column
reference library is checked into this repo; the full-schema library CSVs
used here live on HPC).

Expected pegRNA library CSV columns: barcode, EditedAllele, OriginalAllele,
min_index, max_index, wide_mutated_target, Spacer-Sequence, class, pegRNA,
reporter, oligo_id.

Usage:
    python 00_quantify_fastq_counts.py \\
        --library <pegRNA_library.csv> \\
        --fastq-dir <demultiplexed_fastqs/> \\
        --output-dir <library_QC/>

Expects FASTQ files named `<sample_name>_read1.fastq` / `<sample_name>_read2.fastq`
in --fastq-dir, one pair per bin/replicate (e.g. `pegSTAG1_bulk_r1_read1.fastq`).
Use --strip-prefix if your demultiplexed files carry an extra prefix (e.g.
"i3N_") not part of the sample name used downstream.

Outputs, mirroring the original notebook's directory layout under
--output-dir, one row per file except recombination summaries (one row per
quality bucket):
    read_mapping_filtered/<sample_name>_mapped_reads_perfect_matches.csv
        Perfect-matched reads with their assigned oligo_id -- this is the
        direct input `01_generate_count_table.R` expects for every bin.
    recombination/<sample_name>_recombination.csv
        QC: counts of reads by (barcode matched?, spacer-vs-barcode agree?).
    diversity/<sample_name>_diversity.csv
        QC: per-oligo fraction of perfect-matched reads (library representation).
    editing_quantification/<sample_name>_reporter_editing.csv
        Bulk samples only (see --editing-sample-pattern): per-oligo
        match_edited/match_unedited/other read counts, percent_editing, and
        total_counts -- the direct input `01_generate_count_table.R` expects
        for reporter editing rate.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator

REQUIRED_LIBRARY_COLUMNS = [
    "barcode", "EditedAllele", "OriginalAllele", "min_index", "max_index",
    "wide_mutated_target", "Spacer-Sequence", "class", "pegRNA", "reporter",
    "oligo_id",
]


def load_library(library_path: str) -> pd.DataFrame:
    """Load and preprocess the pegRNA library, once, shared across samples."""
    lib = pd.read_csv(library_path)
    missing = [c for c in REQUIRED_LIBRARY_COLUMNS if c not in lib.columns]
    if missing:
        raise ValueError(
            f"Library file {library_path} is missing expected column(s): {missing}. "
            "This script assumes the full PRIDICT2.0-output schema, not the "
            "reduced reference schema checked into pegRNA_library/."
        )

    lib["barcode_complement"] = lib["barcode"].apply(lambda x: str(Seq(x).reverse_complement()))
    lib.replace("-", "", inplace=True)
    lib["edit_size"] = lib["EditedAllele"].str.len() - lib["OriginalAllele"].str.len()
    lib["min_index"] = lib["min_index"].astype("Int64")
    lib["max_index"] = lib["max_index"].astype("Int64")
    lib["edit_size"] = lib["edit_size"].astype("Int64")

    # Non-targeting controls have no target site, so no designed spacer column
    # value; fall back to the first 20 bases of the full pegRNA oligo.
    lib.loc[lib["class"] == "non-targeting", "Spacer-Sequence"] = (
        lib.loc[lib["class"] == "non-targeting", "pegRNA"].str[:20]
    )

    lib["edited_reporter"] = lib.apply(_extract_edited_reporter, axis=1)
    lib["barcode_prefix7"] = lib["barcode_complement"].astype(str).str[:7]
    return lib


def _extract_edited_reporter(row: pd.Series) -> str | pd.NA:
    seq = row["wide_mutated_target"]
    min_idx, max_idx, edit_size = row["min_index"], row["max_index"], row["edit_size"]
    if (
        isinstance(seq, str)
        and pd.notna(min_idx) and pd.notna(max_idx) and pd.notna(edit_size)
        and isinstance(min_idx, (int, np.integer))
        and isinstance(max_idx, (int, np.integer))
        and isinstance(edit_size, (int, np.integer))
    ):
        return seq[min_idx - 1 : max_idx + edit_size]
    return pd.NA


def _extract_seq_lines(fastq_path: str) -> list[str]:
    lines = []
    with open(fastq_path) as f:
        for i, line in enumerate(f):
            if i % 4 == 1:
                lines.append(line.strip())
    return lines


def _get_avg_quality(fastq_path: str) -> list[float]:
    avg_quals = []
    with open(fastq_path) as handle:
        for _, _seq, qual in FastqGeneralIterator(handle):
            phred_scores = [ord(ch) - 33 for ch in qual]
            avg_quals.append(np.mean(phred_scores))
    return avg_quals


def _match_spacer(seq: str, spacers: list[str], start: int, end: int) -> str | None:
    region = seq[start - 1 : end]
    for spacer in spacers:
        if spacer in region:
            return spacer
    return None


def quantify_sample(
    read1_path: str,
    read2_path: str,
    library: pd.DataFrame,
    quality_threshold: float,
    spacer_window: tuple[int, int],
    barcode_length: int,
) -> dict[str, pd.DataFrame]:
    """Run the full per-read quantification for one bin/replicate sample."""
    read1 = _extract_seq_lines(read1_path)
    read2 = _extract_seq_lines(read2_path)
    fastq_df = pd.DataFrame({"read_1": read1, "read_2": read2})

    fastq_df["average_read1_quality"] = _get_avg_quality(read1_path)
    fastq_df["average_read2_quality"] = _get_avg_quality(read2_path)
    fastq_df = fastq_df[
        (fastq_df["average_read1_quality"] > quality_threshold)
        & (fastq_df["average_read2_quality"] > quality_threshold)
    ].copy()

    spacer_seqs = [str(s) for s in library["Spacer-Sequence"]]
    start, end = spacer_window
    fastq_df["matched_id"] = fastq_df["read_1"].apply(lambda x: _match_spacer(x, spacer_seqs, start, end))
    fastq_df["read2_prefix7"] = fastq_df["read_2"].str[:barcode_length]

    merged_df = fastq_df.merge(library, left_on="read2_prefix7", right_on="barcode_prefix7", how="left")

    recombination = (
        merged_df.groupby([merged_df["read2_prefix7"].isna(), merged_df["matched_id"] == merged_df["Spacer-Sequence"]])
        .size()
        .reset_index(name="n")
    )
    recombination.columns = ["barcode_unmatched", "spacer_barcode_agree", "n"]

    perfect_match_df = merged_df[
        merged_df["read2_prefix7"].notna() & (merged_df["matched_id"] == merged_df["Spacer-Sequence"])
    ].copy()

    diversity = perfect_match_df["oligo_id"].value_counts(normalize=True).reset_index()
    diversity.columns = ["oligo_id", "fraction"]

    perfect_match_df["read_2_nobarcode"] = perfect_match_df["read_2"].str[barcode_length:]
    perfect_match_df["match_unedited"] = perfect_match_df.apply(
        lambda row: row["read_2_nobarcode"].lower().startswith(str(row["reporter"]).lower()), axis=1
    )
    perfect_match_df["match_edited"] = perfect_match_df.apply(
        lambda row: row["read_2_nobarcode"].lower().startswith(str(row["edited_reporter"]).lower()), axis=1
    )

    grouped = (
        perfect_match_df.groupby(["oligo_id", "match_unedited", "match_edited"]).size().reset_index(name="count")
    )
    pivoted = grouped.pivot_table(
        index="oligo_id", columns=["match_unedited", "match_edited"], values="count", fill_value=0
    )
    pivoted.columns = [f"{int(k1)}_{int(k2)}" for k1, k2 in pivoted.columns]
    pivoted.rename(columns={"1_0": "match_unedited", "0_1": "match_edited", "0_0": "other"}, inplace=True)
    for col in ["match_unedited", "match_edited", "other"]:
        if col not in pivoted.columns:
            pivoted[col] = 0
    pivoted["total_counts"] = pivoted["other"] + pivoted["match_edited"] + pivoted["match_unedited"]
    pivoted["percent_editing"] = pivoted["match_edited"] / pivoted["total_counts"] * 100
    pivoted = pivoted.reset_index()

    return {
        "mapped_reads_perfect_matches": perfect_match_df,
        "recombination": recombination,
        "diversity": diversity,
        "reporter_editing": pivoted,
    }


def find_sample_pairs(fastq_dir: str, strip_prefix: str | None) -> list[tuple[str, str, str]]:
    """Return (sample_name, read1_path, read2_path) for every matched pair."""
    fastq_dir_path = Path(fastq_dir)
    pairs = []
    for read1_path in sorted(fastq_dir_path.glob("*_read1.fastq")):
        stem = read1_path.name[: -len("_read1.fastq")]
        read2_path = fastq_dir_path / f"{stem}_read2.fastq"
        if not read2_path.exists():
            print(f"  WARNING: no matching read2 for {read1_path.name}, skipping")
            continue
        sample_name = stem
        if strip_prefix and sample_name.startswith(strip_prefix):
            sample_name = sample_name[len(strip_prefix):]
        pairs.append((sample_name, str(read1_path), str(read2_path)))
    return pairs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--library", required=True, help="Full PRIDICT2.0-output pegRNA library CSV")
    parser.add_argument("--fastq-dir", required=True, help="Directory of <sample>_read1.fastq / _read2.fastq pairs")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--quality-threshold", type=float, default=30.0)
    parser.add_argument("--spacer-window-start", type=int, default=18)
    parser.add_argument("--spacer-window-end", type=int, default=38)
    parser.add_argument("--barcode-length", type=int, default=7)
    parser.add_argument("--strip-prefix", default=None,
                         help='Prefix to remove from FASTQ filenames when deriving sample names, e.g. "i3N_"')
    parser.add_argument("--editing-sample-pattern", default=r"_bulk_r\d+$",
                         help="Regex (matched against sample name) selecting which samples get a "
                              "reporter_editing.csv output; 01_generate_count_table.R only needs this for bulk samples")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    for sub in ["read_mapping_filtered", "recombination", "diversity", "editing_quantification"]:
        (out_dir / sub).mkdir(parents=True, exist_ok=True)

    print(f"Loading library {args.library}...")
    library = load_library(args.library)
    print(f"  {len(library)} library rows")

    pairs = find_sample_pairs(args.fastq_dir, args.strip_prefix)
    print(f"Found {len(pairs)} sample(s) in {args.fastq_dir}")

    editing_pattern = re.compile(args.editing_sample_pattern)

    for sample_name, read1_path, read2_path in pairs:
        print(f"\n{sample_name}:")
        results = quantify_sample(
            read1_path, read2_path, library,
            quality_threshold=args.quality_threshold,
            spacer_window=(args.spacer_window_start, args.spacer_window_end),
            barcode_length=args.barcode_length,
        )

        n_perfect = len(results["mapped_reads_perfect_matches"])
        print(f"  {n_perfect} perfect-matched reads")
        results["mapped_reads_perfect_matches"].to_csv(
            out_dir / "read_mapping_filtered" / f"{sample_name}_mapped_reads_perfect_matches.csv", index=False
        )
        results["recombination"]["sample"] = sample_name
        results["recombination"].to_csv(out_dir / "recombination" / f"{sample_name}_recombination.csv", index=False)
        results["diversity"].to_csv(out_dir / "diversity" / f"{sample_name}_diversity.csv", index=False)

        if editing_pattern.search(sample_name):
            results["reporter_editing"].to_csv(
                out_dir / "editing_quantification" / f"{sample_name}_reporter_editing.csv", index=False
            )
            print(f"  wrote reporter_editing.csv ({len(results['reporter_editing'])} oligos)")

    print(f"\nDone. Outputs under {out_dir}")


if __name__ == "__main__":
    main()
