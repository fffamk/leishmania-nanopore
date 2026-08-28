#!/usr/bin/env python3
"""Recompute adaptive-sampling results from read evidence and run totals."""

import csv
import gzip
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "analysis_inputs"
OUTPUT = ROOT / "reproduced_source_data"
BIOPSIES = ["Biopsy L", "Biopsy A", "Biopsy O"]
LOCUS_WINDOW_BP = 10_000
FIRST_24H_SECONDS = 86_400


def read(name):
    opener = gzip.open if name.endswith(".gz") else open
    with opener(INPUT / name, "rt", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(name, fields, rows):
    OUTPUT.mkdir(exist_ok=True)
    with (OUTPUT / name).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def clean(value, digits=10):
    return f"{value:.{digits}f}".rstrip("0").rstrip(".")


def validate_assignment(row):
    checks = {
        "read length": int(row["read_length"]) >= 500,
        "mean Q score": float(row["qscore"]) >= 10,
        "mapping quality": int(row["best_mapq"]) >= 20,
        "alignment identity": float(row["best_identity"]) >= 0.85,
        "query coverage": float(row["best_query_coverage"]) >= 0.80,
        "genus score margin": float(row["genus_score_margin"]) > 0,
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise ValueError(f"Read {row['read_id']} failed: {', '.join(failed)}")


def detected(rows):
    loci = {
        (row["best_target"], int(row["best_target_start"]) // LOCUS_WINDOW_BP)
        for row in rows
    }
    chromosomes = {row["best_chromosome"] for row in rows}
    return len(rows) >= 3 and len(loci) >= 3 and len(chromosomes) >= 2


def summarize_target_reads(rows):
    nuclear = sorted(
        (
            row
            for row in rows
            if row["is_nuclear"] == "true"
            and float(row["elapsed_seconds"]) <= FIRST_24H_SECONDS
        ),
        key=lambda row: (float(row["elapsed_seconds"]), row["read_id"]),
    )
    detection_seconds = None
    for end in range(1, len(nuclear) + 1):
        if detected(nuclear[:end]):
            detection_seconds = float(nuclear[end - 1]["elapsed_seconds"])
            break
    if detection_seconds is None:
        raise ValueError("Matched sample did not meet the sequencing-positive rule")
    return {
        "qualifying_nuclear_reads": len(nuclear),
        "qualifying_nuclear_bases": sum(int(row["read_length"]) for row in nuclear),
        "time_to_detection_seconds": detection_seconds,
    }


def main():
    evidence = read("matched_target_reads.tsv.gz")
    for row in evidence:
        validate_assignment(row)
    grouped = defaultdict(list)
    for row in evidence:
        grouped[(row["biopsy"], row["mode"])].append(row)
    summaries = {key: summarize_target_reads(rows) for key, rows in grouped.items()}

    totals = {
        (row["biopsy"], row["mode"]): row for row in read("primary_24h_totals.tsv")
    }
    panel_a = []
    panel_b = []
    panel_c = []
    for biopsy in BIOPSIES:
        standard = summaries[(biopsy, "standard")]
        enrichment = summaries[(biopsy, "enrichment")]
        standard_total = totals[(biopsy, "standard")]
        enrichment_total = totals[(biopsy, "enrichment")]
        standard_bases = standard["qualifying_nuclear_bases"]
        enrichment_bases = enrichment["qualifying_nuclear_bases"]
        standard_fraction = float(standard_total["reported_target_base_fraction"])
        enrichment_fraction = float(enrichment_total["reported_target_base_fraction"])
        calculated_standard = standard_bases / int(
            float(standard_total["primary_bases"])
        )
        calculated_enrichment = enrichment_bases / int(
            float(enrichment_total["primary_bases"])
        )
        if abs(calculated_standard - standard_fraction) >= 5e-10:
            raise ValueError(f"Standard target fraction mismatch for {biopsy}")
        if abs(calculated_enrichment - enrichment_fraction) >= 5e-10:
            raise ValueError(f"Enrichment target fraction mismatch for {biopsy}")
        standard_percent = 100 * standard_fraction
        enrichment_percent = 100 * enrichment_fraction
        standard_minutes = standard["time_to_detection_seconds"] / 60
        enrichment_minutes = enrichment["time_to_detection_seconds"] / 60
        panel_a.append(
            {
                "biopsy": biopsy,
                "standard_target_base_percent": clean(standard_percent),
                "enrichment_target_base_percent": clean(enrichment_percent),
                "fold_change": clean(enrichment_percent / standard_percent),
            }
        )
        panel_b.append(
            {
                "biopsy": biopsy,
                "standard_target_bases": standard_bases,
                "enrichment_target_bases": enrichment_bases,
                "enrichment_over_standard": clean(enrichment_bases / standard_bases),
            }
        )
        panel_c.append(
            {
                "biopsy": biopsy,
                "standard_minutes": clean(standard_minutes),
                "enrichment_minutes": clean(enrichment_minutes),
                "enrichment_minus_standard_minutes": clean(
                    enrichment_minutes - standard_minutes
                ),
            }
        )
    panel_c.sort(key=lambda row: float(row["standard_minutes"]))
    write("panel_a.tsv", list(panel_a[0]), panel_a)
    write("panel_b.tsv", list(panel_b[0]), panel_b)
    write("panel_c.tsv", list(panel_c[0]), panel_c)

    checkpoints = read("run_yield_checkpoints.tsv")
    write("panel_d.tsv", list(checkpoints[0]), checkpoints)

    panel_e = []
    for row in read("adaptive_decision_counts.tsv"):
        ejected = int(row["unblock"])
        continued = int(row["stop_receiving"]) + int(row["no_decision"])
        no_log = int(row["absent_from_decision_log"])
        total = ejected + continued + no_log
        for outcome, count in (
            ("Ejected", ejected),
            ("Continued", continued),
            ("No log", no_log),
        ):
            panel_e.append(
                {
                    "mode": row["mode"],
                    "origin": row["origin"],
                    "outcome": outcome,
                    "read_count": count,
                    "percent_of_all_classified_reads": clean(100 * count / total),
                }
            )
    write("panel_e.tsv", list(panel_e[0]), panel_e)


if __name__ == "__main__":
    main()
