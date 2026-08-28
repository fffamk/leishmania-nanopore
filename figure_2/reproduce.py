#!/usr/bin/env python3
"""Recompute cohort calls, positive agreement, and detection times from read evidence."""

import csv
import gzip
from collections import defaultdict
from pathlib import Path

from scipy.stats import beta


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "analysis_inputs"
OUTPUT = ROOT / "reproduced_source_data"
LOCUS_WINDOW_BP = 10_000


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


def detection_state(rows):
    loci = {
        (row["best_target"], int(row["best_target_start"]) // LOCUS_WINDOW_BP)
        for row in rows
    }
    chromosomes = {row["best_chromosome"] for row in rows}
    return (
        len(rows),
        len(loci),
        len(chromosomes),
        (len(rows) >= 3 and len(loci) >= 3 and len(chromosomes) >= 2),
    )


def summarize(rows):
    nuclear = sorted(
        (row for row in rows if row["is_nuclear"] == "true"),
        key=lambda row: (float(row["elapsed_seconds"]), row["read_id"]),
    )
    read_count, locus_count, chromosome_count, positive = detection_state(nuclear)
    detection_seconds = None
    for end in range(1, len(nuclear) + 1):
        if detection_state(nuclear[:end])[3]:
            detection_seconds = float(nuclear[end - 1]["elapsed_seconds"])
            break
    return {
        "high_confidence_nuclear_reads": read_count,
        "distinct_nuclear_loci": locus_count,
        "distinct_nuclear_chromosomes": chromosome_count,
        "sequencing_positive": "true" if positive else "false",
        "time_to_detection_minutes": (
            str(detection_seconds / 60) if detection_seconds is not None else ""
        ),
    }


def main():
    evidence = read("cohort_target_reads.tsv.gz")
    for row in evidence:
        validate_assignment(row)
    by_biopsy = defaultdict(list)
    for row in evidence:
        by_biopsy[row["biopsy"]].append(row)

    cohort = []
    for sample in read("sample_metadata.tsv"):
        cohort.append(sample | summarize(by_biopsy[sample["biopsy"]]))

    positives = [row for row in cohort if row["pcr_status"] == "positive"]
    detected = sum(row["sequencing_positive"] == "true" for row in positives)
    total = len(positives)
    ppa = detected / total
    lower = beta.ppf(0.025, detected, total - detected + 1)
    upper = beta.ppf(0.975, detected + 1, total - detected)

    ranked = sorted(cohort, key=lambda row: int(row["high_confidence_nuclear_reads"]))
    fields = [
        "biopsy",
        "rank",
        "high_confidence_nuclear_reads",
        "distinct_nuclear_loci",
        "distinct_nuclear_chromosomes",
        "sequencing_positive",
        "pcr_status",
        "ppa",
        "exact_95ci_lower",
        "exact_95ci_upper",
    ]
    panel_a = []
    for rank, row in enumerate(ranked, 1):
        panel_a.append(
            {
                "biopsy": row["biopsy"],
                "rank": rank,
                "high_confidence_nuclear_reads": row["high_confidence_nuclear_reads"],
                "distinct_nuclear_loci": row["distinct_nuclear_loci"],
                "distinct_nuclear_chromosomes": row["distinct_nuclear_chromosomes"],
                "sequencing_positive": row["sequencing_positive"],
                "pcr_status": row["pcr_status"],
                "ppa": f"{ppa:.10f}",
                "exact_95ci_lower": f"{lower:.10f}",
                "exact_95ci_upper": f"{upper:.10f}",
            }
        )
    write("panel_a.tsv", fields, panel_a)

    run_end = float(read("run_metadata.tsv")[0]["run_duration_hours"]) * 60
    events = sorted(
        (row for row in positives if row["sequencing_positive"] == "true"),
        key=lambda row: float(row["time_to_detection_minutes"]),
    )
    event_order = {row["biopsy"]: index for index, row in enumerate(events, 1)}
    fields = [
        "biopsy",
        "pcr_status",
        "sequencing_positive",
        "time_to_detection_minutes",
        "censor_time_minutes",
        "event_order",
    ]
    panel_b = []
    for row in positives:
        called = row["sequencing_positive"] == "true"
        panel_b.append(
            {
                "biopsy": row["biopsy"],
                "pcr_status": row["pcr_status"],
                "sequencing_positive": row["sequencing_positive"],
                "time_to_detection_minutes": row["time_to_detection_minutes"]
                if called
                else "",
                "censor_time_minutes": ""
                if called
                else f"{run_end:.8f}".rstrip("0").rstrip("."),
                "event_order": event_order.get(row["biopsy"], ""),
            }
        )
    write("panel_b.tsv", fields, panel_b)


if __name__ == "__main__":
    main()
