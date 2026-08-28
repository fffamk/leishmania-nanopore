#!/usr/bin/env python3
"""Recompute cohort and held-out species-call composition."""

import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "analysis_inputs"
OUTPUT = ROOT / "reproduced_source_data"
LOCUS_WINDOW_BP = 10_000
REPORTED_SPECIES = [
    "Leishmania major",
    "Leishmania tropica",
    "Leishmania aethiopica",
    "Leishmania donovani",
    "Leishmania infantum",
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


def cohort_positive(rows):
    nuclear = [row for row in rows if row["is_nuclear"] == "true"]
    loci = {
        (row["best_target"], int(row["best_target_start"]) // LOCUS_WINDOW_BP)
        for row in nuclear
    }
    chromosomes = {row["best_chromosome"] for row in nuclear}
    return len(nuclear) >= 3 and len(loci) >= 3 and len(chromosomes) >= 2


def fraction(value, total):
    return "" if total == 0 else f"{value / total:.15g}"


def cohort_panel():
    assignments = read("cohort_target_reads.tsv.gz")
    for row in assignments:
        validate_assignment(row)
    by_biopsy = defaultdict(list)
    for row in assignments:
        by_biopsy[row["biopsy"]].append(row)

    kraken_rows = defaultdict(list)
    for row in read("cohort_kraken_reports.tsv.gz"):
        kraken_rows[row["biopsy"]].append(row)

    panel = []
    samples = read("sample_metadata.tsv")
    samples.sort(key=lambda row: row["pcr_status"] == "negative")
    for sample in samples:
        biopsy = sample["biopsy"]
        nuclear = [row for row in by_biopsy[biopsy] if row["is_nuclear"] == "true"]
        alignment_counts = Counter(row["best_species"] for row in nuclear)
        report = kraken_rows[biopsy]
        report_by_taxid = {row["taxid"]: row for row in report}
        kraken_counts = {
            "Leishmania_major": int(
                report_by_taxid.get("5664", {}).get("taxon_reads", 0)
            ),
            "Leishmania_tropica": int(
                report_by_taxid.get("5666", {}).get("taxon_reads", 0)
            ),
            "Leishmania_aethiopica": int(
                report_by_taxid.get("5667", {}).get("taxon_reads", 0)
            ),
        }
        leishmania = report_by_taxid.get("5658", {})
        method_rows = [
            ("alignment", alignment_counts, len(nuclear)),
            ("kraken2", kraken_counts, int(leishmania.get("clade_reads", 0))),
        ]
        positive = cohort_positive(by_biopsy[biopsy])
        for method, counts, method_total in method_rows:
            major = counts["Leishmania_major"]
            tropica = counts["Leishmania_tropica"]
            aethiopica = counts["Leishmania_aethiopica"]
            three_species_total = major + tropica + aethiopica
            panel.append(
                {
                    "biopsy": biopsy,
                    "pcr_status": sample["pcr_status"],
                    "sequencing_positive": "true" if positive else "false",
                    "method": method,
                    "composition_status": (
                        "reported" if positive else "suppressed_no_sequencing_call"
                    ),
                    "major_reads": major,
                    "tropica_reads": tropica,
                    "aethiopica_reads": aethiopica,
                    "three_species_reads": three_species_total,
                    "leishmania_reads": method_total,
                    "other_or_genus_only_reads": method_total - three_species_total,
                    "major_fraction_of_three_species_calls": fraction(
                        major, three_species_total
                    ),
                    "tropica_fraction_of_three_species_calls": fraction(
                        tropica, three_species_total
                    ),
                    "aethiopica_fraction_of_three_species_calls": fraction(
                        aethiopica, three_species_total
                    ),
                }
            )
    return panel


def parse_tags(fields):
    tags = {}
    for field in fields:
        pieces = field.split(":", 2)
        if len(pieces) == 3:
            key, kind, raw = pieces
            tags[key] = float(raw) if kind == "f" else int(raw) if kind == "i" else raw
    return tags


def parse_alignment(line, targets):
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 12:
        raise ValueError("PAF record has fewer than 12 fields")
    if fields[5] not in targets:
        raise ValueError(f"PAF target is absent from the manifest: {fields[5]}")
    target = targets[fields[5]]
    tags = parse_tags(fields[12:])
    return {
        "read_id": fields[0],
        "read_length": int(fields[1]),
        "query_start": int(fields[2]),
        "query_end": int(fields[3]),
        "target_name": fields[5],
        "target_start": min(int(fields[7]), int(fields[8])),
        "matches": int(fields[9]),
        "aligned_length": int(fields[10]),
        "mapq": int(fields[11]),
        "score": float(tags.get("AS", fields[9])),
        "target": target,
    }


def alignment_key(hit):
    return (
        -hit["score"],
        -hit["mapq"],
        -hit["matches"],
        -hit["aligned_length"],
        hit["target_name"],
        hit["target_start"],
    )


def reduce_alignment_group(hits):
    ordered = sorted(hits, key=alignment_key)
    best = ordered[0]
    best_by_genus = {}
    for hit in ordered:
        label = hit["target"]["genus"].casefold() or hit["target_name"].casefold()
        best_by_genus.setdefault(label, hit)
    genera = sorted(best_by_genus.values(), key=alignment_key)
    second_genus_score = genera[1]["score"] if len(genera) > 1 else 0.0
    genus_margin = genera[0]["score"] - second_genus_score
    identity = best["matches"] / best["aligned_length"] if best["aligned_length"] else 0
    query_coverage = (
        (best["query_end"] - best["query_start"]) / best["read_length"]
        if best["read_length"]
        else 0
    )
    qualifying = (
        best["read_length"] >= 500
        and best["mapq"] >= 20
        and best["aligned_length"] >= 500
        and identity >= 0.85
        and query_coverage >= 0.80
        and genus_margin > 0
        and best["target"]["genus"].casefold() == "leishmania"
        and best["target"]["compartment"].casefold() == "nuclear"
    )
    return best["read_id"], best["target"]["species"].replace("_", " "), qualifying


def competitive_holdout_calls(targets, truth):
    calls = defaultdict(Counter)
    seen = set()
    current_id = None
    group = []
    with gzip.open(INPUT / "holdout.paf.gz", "rt", encoding="utf-8") as handle:
        for line in handle:
            hit = parse_alignment(line, targets)
            if current_id is None:
                current_id = hit["read_id"]
            if hit["read_id"] != current_id:
                if current_id in seen:
                    raise RuntimeError(f"non-contiguous PAF query: {current_id}")
                seen.add(current_id)
                read_id, species, qualifying = reduce_alignment_group(group)
                unit = read_id.split("::", 1)[0]
                if unit not in truth:
                    raise RuntimeError(f"unknown holdout unit: {unit}")
                if qualifying:
                    calls[unit][species] += 1
                group = []
                current_id = hit["read_id"]
            group.append(hit)
    if group:
        read_id, species, qualifying = reduce_alignment_group(group)
        unit = read_id.split("::", 1)[0]
        if unit not in truth:
            raise RuntimeError(f"unknown holdout unit: {unit}")
        if qualifying:
            calls[unit][species] += 1
    return calls


def kraken_holdout_calls(taxonomy, truth):
    calls = defaultdict(Counter)
    with gzip.open(
        INPUT / "holdout.kraken.output.gz", "rt", encoding="utf-8"
    ) as handle:
        for line in handle:
            status, read_id, taxid, _length, _mapping = line.rstrip("\n").split("\t")
            unit = read_id.split("::", 1)[0]
            if unit not in truth:
                raise RuntimeError(f"unknown holdout unit: {unit}")
            if status == "C" and taxid in taxonomy:
                calls[unit][taxonomy[taxid]] += 1
    return calls


def holdout_panel():
    truth_rows = read("holdout_truth.tsv")
    truth = {row["unit_id"]: row for row in truth_rows}
    taxonomy = {row["taxid"]: row["species"] for row in read("holdout_taxonomy.tsv")}
    targets = {
        row["target_name"]: row for row in read("holdout_target_manifest.tsv.gz")
    }
    methods = [
        ("competitive_alignment", competitive_holdout_calls(targets, truth)),
        ("kraken2", kraken_holdout_calls(taxonomy, truth)),
    ]
    panel = []
    for method, calls in methods:
        for truth_row in truth_rows:
            unit = truth_row["unit_id"]
            total = sum(calls[unit].values())
            for species in REPORTED_SPECIES:
                count = calls[unit][species]
                panel.append(
                    {
                        "method": method,
                        "heldout_strain": truth_row["heldout_strain"],
                        "true_species": truth_row["true_species"],
                        "called_species": species,
                        "supporting_fragments": count,
                        "percent_among_species_calls": f"{100 * count / total:.15g}"
                        if total
                        else "",
                        "species_level_calls": total,
                    }
                )
    return panel


def main():
    panel_a = cohort_panel()
    fields = [
        "biopsy",
        "pcr_status",
        "sequencing_positive",
        "method",
        "composition_status",
        "major_reads",
        "tropica_reads",
        "aethiopica_reads",
        "three_species_reads",
        "leishmania_reads",
        "other_or_genus_only_reads",
        "major_fraction_of_three_species_calls",
        "tropica_fraction_of_three_species_calls",
        "aethiopica_fraction_of_three_species_calls",
    ]
    write("panel_a.tsv", fields, panel_a)
    panel_b = holdout_panel()
    write("panel_b.tsv", list(panel_b[0]), panel_b)


if __name__ == "__main__":
    main()
