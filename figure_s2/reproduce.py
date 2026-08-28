#!/usr/bin/env python3
"""Recompute the batch-size and platform cost results."""

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "analysis_inputs"
OUTPUT = ROOT / "reproduced_source_data"


def read_inputs():
    with (INPUT / "cost_inputs.tsv").open(newline="", encoding="utf-8") as handle:
        return {
            row["key"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")
        }


def number(values, key):
    return float(values[key])


def barcode_allocation(values, n):
    if n <= 24:
        libraries = number(values, "nbd96_libraries_at_24plex")
    elif n <= 48:
        libraries = number(values, "nbd96_libraries_at_48plex")
    else:
        libraries = number(values, "nbd96_libraries_at_96plex")
    return number(values, "nbd96_kit") / libraries


def consumables(values, n, flowcell_price):
    qubit = number(values, "qubit_hs_reagent") / number(
        values, "qubit_assays"
    ) + number(values, "qubit_tubes") / number(values, "qubit_tube_count")
    preparation = n * (
        number(values, "ffpe_repair_mix") / number(values, "ffpe_repair_specimens")
        + number(values, "endprep_module") / number(values, "endprep_specimens")
        + number(values, "blunt_ta_ligase") / number(values, "blunt_ta_specimens")
    )
    return (
        flowcell_price
        + barcode_allocation(values, n)
        + preparation
        + number(values, "quick_ligation_module")
        / number(values, "quick_ligation_libraries")
        + (n + 2) * qubit
        + n * number(values, "nanopore_minor_variable")
        + number(values, "nanopore_minor_fixed")
    )


def write(name, fields, rows):
    OUTPUT.mkdir(exist_ok=True)
    with (OUTPUT / name).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    values = read_inputs()
    minion_price = number(values, "flowcell_flo_min114")
    panel_a = []
    for n in range(1, 25):
        panel_a.append(
            {
                "method": "direct_nanopore",
                "strategy": "optimized_nbd96",
                "cost_scope": "consumables_only",
                "batch_size": n,
                "usd_per_specimen": f"{consumables(values, n, minion_price) / n:.10g}",
                "scenario_basis": "modeled 1–24-plex range",
                "observed_plex_reference": str(n == 16).lower(),
            }
        )
    for n in range(24, 97):
        panel_a.append(
            {
                "method": "direct_nanopore",
                "strategy": "optimized_nbd96",
                "cost_scope": "consumables_only",
                "batch_size": n,
                "usd_per_specimen": f"{consumables(values, n, minion_price) / n:.10g}",
                "scenario_basis": "modeled 24–96-plex range",
                "observed_plex_reference": "false",
            }
        )
    pcr_n = [int(value) for value in values["pcr_curve_n"].split(",")]
    pcr_cost = [float(value) for value in values["pcr_curve_usd"].split(",")]
    for n, cost in zip(pcr_n, pcr_cost):
        panel_a.append(
            {
                "method": "endpoint_ITS_PCR",
                "strategy": "base_throughput_model",
                "cost_scope": "fully_loaded",
                "batch_size": n,
                "usd_per_specimen": f"{cost:g}",
                "scenario_basis": "endpoint-PCR reference model",
                "observed_plex_reference": "false",
            }
        )
    write("panel_a.tsv", list(panel_a[0]), panel_a)

    channel_ratio = number(values, "promethion_channels") / number(
        values, "minion_channels"
    )
    minion_capacity = int(number(values, "observed_minion_capacity"))
    promethion_capacity = round(minion_capacity * channel_ratio)
    promethion_price = number(values, "flowcell_flo_pro114m_pack") / number(
        values, "flowcell_flo_pro114m_pack_count"
    )
    promethion_cost = (
        consumables(values, promethion_capacity, promethion_price) / promethion_capacity
    )
    promethion_relative_output = channel_ratio * minion_capacity / promethion_capacity
    minion_cost = consumables(values, minion_capacity, minion_price) / minion_capacity
    panel_b = [
        {
            "order": 1,
            "method": "endpoint_ITS_PCR",
            "scenario": "24-specimen fully loaded reference",
            "specimens_per_batch_or_flowcell": 24,
            "usd_per_specimen": f"{pcr_cost[-1]:g}",
            "cost_scope": "fully_loaded",
            "relative_output_per_specimen": "",
            "model_basis": "fully loaded endpoint-PCR reference model",
        },
        {
            "order": 2,
            "method": "PromethION",
            "scenario": "channel_scaled_depth_match",
            "specimens_per_batch_or_flowcell": promethion_capacity,
            "usd_per_specimen": f"{promethion_cost:.10g}",
            "cost_scope": "consumables_only",
            "relative_output_per_specimen": f"{promethion_relative_output:.10g}",
            "model_basis": "channel-normalized capacity model",
        },
        {
            "order": 3,
            "method": "MinION",
            "scenario": "observed_depth_reference",
            "specimens_per_batch_or_flowcell": minion_capacity,
            "usd_per_specimen": f"{minion_cost:.10g}",
            "cost_scope": "consumables_only",
            "relative_output_per_specimen": "1",
            "model_basis": "study cohort capacity and optimized SQK-NBD114.96 allocation",
        },
    ]
    write("panel_b.tsv", list(panel_b[0]), panel_b)


if __name__ == "__main__":
    main()
