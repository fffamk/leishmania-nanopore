#!/usr/bin/env python3
"""Recompute marker evidence from public reads and biopsy-sample BLAST hits."""

import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "analysis_inputs"
OUTPUT = ROOT / "reproduced_source_data"
MARKERS = [
    ("its1", "ITS1"),
    ("its2", "ITS2"),
    ("miniexon", "Miniexon"),
    ("hsp20", "HSP20"),
    ("hsp70", "HSP70 coding"),
    ("hsp70_3utr", "HSP70 3′ UTR"),
    ("cytb", "Cytochrome b"),
]
MIN_ALIGNED = {
    "its1": 180,
    "its2": 250,
    "miniexon": 180,
    "hsp20": 220,
    "hsp70": 700,
    "hsp70_3utr": 350,
    "cytb": 650,
}
BIOPSIES = [
    *(f"Biopsy {letter}" for letter in "ABCDEFGHIJKLMNO"),
    "PCR-negative biopsy",
]


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


def public_panel():
    public = read("public_marker_scores.tsv.gz")
    runs = list(dict.fromkeys(row["run_accession"] for row in public))
    grouped = defaultdict(list)
    for row in public:
        best_score = float(row["best_bitscore"])
        second_score = float(row["second_bitscore"] or 0)
        margin = best_score - second_score
        required_margin = max(15.0, 0.02 * best_score)
        row["called_species"] = (
            row["best_species"] if margin >= required_margin else "unresolved"
        )
        row["outcome"] = (
            "unresolved"
            if row["called_species"] == "unresolved"
            else "correct"
            if row["called_species"] == row["true_species"]
            else "wrong"
        )
        grouped[(row["run_accession"], row["marker"])].append(row)

    panel = []
    for run in runs:
        for marker, label in MARKERS:
            rows = grouped[(run, marker)]
            counts = Counter(
                row["called_species"]
                for row in rows
                if row["called_species"] != "unresolved"
            )
            ranked = counts.most_common()
            tied = len(ranked) > 1 and ranked[0][1] == ranked[1][1]
            plurality = "" if not ranked or tied else ranked[0][0]
            support = ranked[0][1] if plurality else 0
            resolved = sum(counts.values())
            panel.append(
                {
                    "run_accession": run,
                    "known_species": rows[0]["true_species"],
                    "marker_label": label,
                    "evidence_reads": len(rows),
                    "resolved_reads": resolved,
                    "unresolved_reads": len(rows) - resolved,
                    "plurality_species": plurality,
                    "plurality_support": support,
                    "cell_state": "resolved" if plurality else "unresolved",
                }
            )
    return panel


def marker_target(value):
    fields = value.split("|", 3)
    if len(fields) != 4:
        raise ValueError(f"Invalid marker target identifier: {value}")
    marker, species, accession, role = fields
    return marker, species, accession, role


def biopsy_calls():
    grouped = defaultdict(list)
    for row in read("biopsy_marker_blast.tsv.gz"):
        marker, species, accession, role = marker_target(row["sseqid"])
        if marker not in MIN_ALIGNED:
            raise ValueError(f"Unknown marker: {marker}")
        if int(row["length"]) < MIN_ALIGNED[marker] or float(row["pident"]) < 80:
            continue
        row |= {
            "marker": marker,
            "species": species,
            "accession": accession,
            "role": role,
        }
        grouped[(row["biopsy"], row["qseqid"], marker)].append(row)

    calls = []
    for (biopsy, read_id, marker), hits in sorted(grouped.items()):
        best_by_species = {}
        for hit in hits:
            old = best_by_species.get(hit["species"])
            key = (float(hit["bitscore"]), int(hit["length"]), float(hit["pident"]))
            old_key = (
                (float(old["bitscore"]), int(old["length"]), float(old["pident"]))
                if old
                else None
            )
            if old is None or key > old_key:
                best_by_species[hit["species"]] = hit
        ranked = sorted(
            best_by_species.values(),
            key=lambda row: (
                float(row["bitscore"]),
                int(row["length"]),
                float(row["pident"]),
            ),
            reverse=True,
        )
        if len(ranked) > 1 and float(ranked[0]["bitscore"]) <= float(
            ranked[1]["bitscore"]
        ):
            continue
        calls.append(
            {
                "biopsy": biopsy,
                "read_id": read_id,
                "marker": marker,
                "best_species": ranked[0]["species"],
            }
        )
    return calls


def biopsy_panel():
    samples = {row["biopsy"]: row["pcr_status"] for row in read("cohort_samples.tsv")}
    grouped = defaultdict(list)
    for row in biopsy_calls():
        grouped[(row["biopsy"], row["marker"])].append(row)
    panel = []
    for biopsy in BIOPSIES:
        for marker, label in MARKERS:
            rows = grouped.get((biopsy, marker), [])
            counts = Counter(row["best_species"] for row in rows)
            best_species, count = counts.most_common(1)[0] if counts else ("", 0)
            panel.append(
                {
                    "biopsy": biopsy,
                    "pcr_status": samples[biopsy],
                    "marker_label": label,
                    "best_species": best_species.replace("_", " "),
                    "independent_reads": count,
                }
            )
    return panel


def main():
    panel_a = public_panel()
    write("panel_a.tsv", list(panel_a[0]), panel_a)
    panel_b = biopsy_panel()
    write("panel_b.tsv", list(panel_b[0]), panel_b)


if __name__ == "__main__":
    main()
