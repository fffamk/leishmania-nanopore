#!/usr/bin/env python3
"""Recompute marker divergence from precomputed MAFFT alignments."""

import csv
import hashlib
import itertools
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "analysis_inputs"
OUTPUT = ROOT / "reproduced_source_data"
DNA = frozenset("ACGT")
MARKERS = [
    ("its1", "ITS1"),
    ("its2", "ITS2"),
    ("miniexon", "Miniexon"),
    ("hsp20", "HSP20"),
    ("hsp70", "HSP70 coding"),
    ("hsp70_3utr", "HSP70 3′ UTR"),
    ("cytb", "Cytochrome b"),
]
PRIMARY = {
    "its1": "einsi",
    "its2": "einsi",
    "miniexon": "einsi",
    "hsp20": "ginsi",
    "hsp70": "ginsi",
    "hsp70_3utr": "einsi",
    "cytb": "ginsi",
}
BENCHMARK_PAIRS = [
    ("Leishmania_infantum", "Leishmania_donovani"),
    ("Leishmania_panamensis", "Leishmania_guyanensis"),
    ("Leishmania_panamensis", "Leishmania_braziliensis"),
    ("Leishmania_guyanensis", "Leishmania_braziliensis"),
    ("Leishmania_mexicana", "Leishmania_infantum"),
    ("Leishmania_mexicana", "Leishmania_donovani"),
    ("Leishmania_mexicana", "Leishmania_panamensis"),
    ("Leishmania_mexicana", "Leishmania_guyanensis"),
    ("Leishmania_mexicana", "Leishmania_braziliensis"),
    ("Leishmania_infantum", "Leishmania_panamensis"),
    ("Leishmania_infantum", "Leishmania_guyanensis"),
    ("Leishmania_infantum", "Leishmania_braziliensis"),
    ("Leishmania_donovani", "Leishmania_panamensis"),
    ("Leishmania_donovani", "Leishmania_guyanensis"),
    ("Leishmania_donovani", "Leishmania_braziliensis"),
]
STUDY_SPECIES_PAIRS = [
    ("Leishmania_tropica", "Leishmania_aethiopica"),
    ("Leishmania_major", "Leishmania_tropica"),
    ("Leishmania_major", "Leishmania_aethiopica"),
]


def read_fasta(path):
    records = []
    name = None
    sequence = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(sequence).upper()))
                name = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)
    if name is not None:
        records.append((name, "".join(sequence).upper()))
    return records


def load(scope, marker, strategy):
    path = INPUT / f"{scope}_alignments" / f"{marker}.{strategy}.aligned.fa"
    records = []
    for header, sequence in read_fasta(path):
        _, species, accession, *_ = header.split("|")
        records.append(
            {"species": species, "accession": accession, "sequence": sequence}
        )
    return records


def pair_metrics(sequence_a, sequence_b):
    joint = [
        index
        for index, pair in enumerate(zip(sequence_a, sequence_b))
        if pair[0] in DNA and pair[1] in DNA
    ]
    if not joint:
        raise ValueError("Aligned sequences have no shared nucleotide columns")
    start, end = joint[0], joint[-1]
    overlap = list(zip(sequence_a[start : end + 1], sequence_b[start : end + 1]))
    union = sum(a in DNA or b in DNA for a, b in overlap)
    substitutions = sum(a in DNA and b in DNA and a != b for a, b in overlap)
    gaps = sum((a == "-" and b in DNA) or (b == "-" and a in DNA) for a, b in overlap)
    total = substitutions + gaps
    return total, 100 * total / union


def compare(sequence_a, sequence_b, marker, other_a, other_b):
    primary = pair_metrics(sequence_a, sequence_b)
    secondary = pair_metrics(other_a, other_b)
    stable = abs(primary[0] - secondary[0]) <= 2 and abs(primary[1] - secondary[1]) <= 1
    if PRIMARY[marker] == "ginsi":
        primary, secondary = secondary, primary
    return primary, stable


def short(species):
    return "L. " + species.split("_")[-1]


def format_count(value):
    return str(int(value)) if value == int(value) else f"{value:.1f}"


def write(name, fields, rows):
    OUTPUT.mkdir(exist_ok=True)
    with (OUTPUT / name).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def benchmark_panel():
    results = {}
    for marker, label in MARKERS:
        einsi = load("benchmark_reference", marker, "einsi")
        ginsi = {
            row["accession"]: row
            for row in load("benchmark_reference", marker, "ginsi")
        }
        for left, right in itertools.combinations(einsi, 2):
            primary, stable = compare(
                left["sequence"],
                right["sequence"],
                marker,
                ginsi[left["accession"]]["sequence"],
                ginsi[right["accession"]]["sequence"],
            )
            results[(frozenset((left["species"], right["species"])), marker)] = (
                primary,
                stable,
            )
    rows = []
    for species_a, species_b in BENCHMARK_PAIRS:
        for marker, label in MARKERS:
            (total, rate), stable = results[(frozenset((species_a, species_b)), marker)]
            rows.append(
                {
                    "species_pair": f"{short(species_a)} – {short(species_b)}",
                    "marker": label,
                    "differences": total,
                    "differences_per_100_columns": f"{rate:.6f}",
                    "alignment_stable": "yes" if stable else "no",
                }
            )
    return rows


def study_species_panel():
    rows = []
    for species_a, species_b in STUDY_SPECIES_PAIRS:
        for marker, label in MARKERS:
            einsi = load("study_species", marker, "einsi")
            ginsi = {
                row["accession"]: row for row in load("study_species", marker, "ginsi")
            }
            left = [row for row in einsi if row["species"] == species_a]
            right = [row for row in einsi if row["species"] == species_b]
            unique = {}
            for a in left:
                for b in right:
                    primary, stable = compare(
                        a["sequence"],
                        b["sequence"],
                        marker,
                        ginsi[a["accession"]]["sequence"],
                        ginsi[b["accession"]]["sequence"],
                    )
                    hash_a = hashlib.sha256(
                        a["sequence"].replace("-", "").encode()
                    ).hexdigest()
                    hash_b = hashlib.sha256(
                        b["sequence"].replace("-", "").encode()
                    ).hexdigest()
                    rounded = (primary[0], float(f"{primary[1]:.6f}"))
                    unique.setdefault((hash_a, hash_b), (rounded, stable))
            totals = [value[0][0] for value in unique.values()]
            rates = [value[0][1] for value in unique.values()]
            median_total = statistics.median(totals)
            median_rate = statistics.median(rates)
            sensitive = sum(not value[1] for value in unique.values())
            rows.append(
                {
                    "species_pair": f"{short(species_a)} – {short(species_b)}",
                    "marker": label,
                    "median_differences": format_count(median_total),
                    "minimum_differences": min(totals),
                    "maximum_differences": max(totals),
                    "unique_sequence_pairs": len(unique),
                    "alignment_sensitive_pairs": sensitive,
                    "median_differences_per_100_columns": f"{median_rate:.6f}",
                }
            )
    return rows


def main():
    benchmark = benchmark_panel()
    study_species = study_species_panel()
    write("panel_a.tsv", list(benchmark[0]), benchmark)
    write("panel_b.tsv", list(study_species[0]), study_species)


if __name__ == "__main__":
    main()
