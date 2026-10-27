#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import re
import argparse
import matplotlib.pyplot as plt
from collections import Counter

def parse_args():
    parser = argparse.ArgumentParser(
        description="Parse CSQ_stat.txt files and plot CSQ statistics."
    )

    parser.add_argument(
        "-i", "--input",
        default=".",
        help="Input directory containing *CSQ_stat.txt files [default: current directory]"
    )

    parser.add_argument(
        "-o", "--output",
        default="CSQ_stat",
        help="Output prefix [default: CSQ_stat]"
    )

    return parser.parse_args()


def parse_csq_stat(filename):
    """
    Parse one CSQ_stat.txt file.

    Example:
    label   count   frequency(%)
    Good    1161    97.07
    Moderate 33     2.76
    Poor    2       0.17
    Missing 9
    The database used in CSQ mode:embryophyta
    The fasta used in CSQ mode:Tcacao_523_v2.0.fa
    """

    result = {
        "Good": (0, 0.0),
        "Moderate": (0, 0.0),
        "Poor": (0, 0.0),
        "Missing": 0,
        "database": None,
        "fasta": None,
        "filename": os.path.basename(filename)
    }

    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            # Good / Moderate / Poor
            if line.startswith(("Good", "Moderate", "Poor")):
                fields = line.split()

                if len(fields) >= 3:
                    label = fields[0]
                    count = int(fields[1])
                    frequency = float(fields[2])

                    result[label] = (count, frequency)

            # Missing
            elif line.startswith("Missing"):
                fields = line.split()

                if len(fields) >= 2:
                    result["Missing"] = int(fields[1])

            # Database
            elif line.startswith("The database used in CSQ mode:"):
                result["database"] = line.split(":", 1)[1].strip()

            # FASTA
            elif line.startswith("The fasta used in CSQ mode:"):
                result["fasta"] = line.split(":", 1)[1].strip()

    return result

def read_all_csq_files(input_dir):

    pattern = os.path.join(input_dir, "*CSQ_stat.txt")
    files = sorted(glob.glob(pattern))

    print(f"[INFO] Found {len(files)} CSQ stat files.")

    all_records = []

    # -------------------------------------------------------
    # Step 1: Parse all files
    # -------------------------------------------------------

    for filename in files:

        record = parse_csq_stat(filename)
        all_records.append(record)

        print(
            f"[INFO] {record['filename']}: "
            f"database={record['database']}, "
            f"fasta={record['fasta']}, "
            f"Good={record['Good']}, "
            f"Moderate={record['Moderate']}, "
            f"Poor={record['Poor']}, "
            f"Missing={record['Missing']}"
        )

    if not all_records:
        return {}, None

    # -------------------------------------------------------
    # Step 2: Determine common database / fasta
    # -------------------------------------------------------

    databases = [
        x["database"]
        for x in all_records
        if x["database"]
    ]

    fastas = [
        x["fasta"]
        for x in all_records
        if x["fasta"]
    ]

    same_database = (
        len(databases) == len(all_records)
        and len(set(databases)) == 1
    )

    same_fasta = (
        len(fastas) == len(all_records)
        and len(set(fastas)) == 1
    )

    # -------------------------------------------------------
    # Step 3: Determine title and y-axis key
    # -------------------------------------------------------

    if same_database:
        # Common database -> title
        plot_title = databases[0]
        key_type = "fasta"

    elif same_fasta:
        # Common fasta -> title
        plot_title = fastas[0]
        key_type = "database"

    else:
        # Neither is common
        plot_title = None
        key_type = "database"

    print(f"[INFO] Plot title: {plot_title}")
    print(f"[INFO] Y-axis key: {key_type}")

    # -------------------------------------------------------
    # Step 4: Build dictionary
    # -------------------------------------------------------

    csqstat_dict = {}

    for record in all_records:

        if key_type == "fasta":
            key = record["fasta"]
        else:
            key = record["database"]

        # fallback
        if not key:
            key = record["filename"]

        if key in csqstat_dict:
            print(f"[WARN] Duplicate key '{key}' (from {record['filename']}); "
                  f"using the filename to keep both records.")
            key = record["filename"]

        csqstat_dict[key] = {
            "Good": record["Good"],
            "Moderate": record["Moderate"],
            "Poor": record["Poor"],
            "Missing": record["Missing"],
            "database": record["database"],
            "fasta": record["fasta"],
            "filename": record["filename"]
        }

    return csqstat_dict, plot_title


def flatten_records(csqstat_dict):

    records = []

    for key, item in csqstat_dict.items():

        x = item.copy()

        # y-axis label
        x["plot_label"] = key

        records.append(x)

    return records

def plot_csq_stats(records, output_prefix,plot_title=None):

    if not records:
        print("[WARNING] No CSQ records found.")
        return

    labels = [x["plot_label"] for x in records]

    good = [x["Good"][1] for x in records]
    moderate = [x["Moderate"][1] for x in records]
    poor = [x["Poor"][1] for x in records]
    missing = [x["Missing"] for x in records]

    n = len(records)

    # -------------------------------------------------------
    # Global font size
    # -------------------------------------------------------
    FONT_SIZE = 12

    plt.rcParams.update({
        "font.size": FONT_SIZE,
        "axes.labelsize": FONT_SIZE,
        "xtick.labelsize": FONT_SIZE,
        "ytick.labelsize": FONT_SIZE,
        "legend.fontsize": FONT_SIZE,
    })

    # -------------------------------------------------------
    # Colors
    # -------------------------------------------------------
    COLOR_GOOD = "#ef888a"
    COLOR_MODERATE = "#b8b8b8"
    COLOR_POOR = "#45465e"
    COLOR_MISSING = "#57606f"

    fig_height = max(5, n * 0.55 + 1.8)

    fig, (ax_missing, ax_csq) = plt.subplots(
        1,
        2,
        figsize=(11, fig_height),
        sharey=True,
        gridspec_kw={
            "width_ratios": [1.2, 4],
            "wspace": 0
        }
    )

    y = list(range(n))

    # =======================================================
    # Left: Missing count
    # =======================================================

    ax_missing.barh(
        y,
        missing,
        height=0.70,
        color=COLOR_MISSING,
        edgecolor="white",
        linewidth=0.6,
        hatch="///",
        label="Missing"
    )

    max_missing = max(missing) if max(missing) > 0 else 1

    ax_missing.set_xlim(
        max_missing * 1.20,
        0
    )

    ax_missing.set_xlabel(
        "Missing",
        fontsize=FONT_SIZE
    )

    ax_missing.set_yticks(y)
    ax_missing.set_yticklabels(
        labels,
        fontsize=FONT_SIZE
    )

    # Missing values
    for yi, value in zip(y, missing):

        if value > 0:
            ax_missing.text(
                value,
                yi,
                f"{value} ",
                va="center",
                ha="right",
                fontsize=FONT_SIZE
            )
        else:
            ax_missing.text(
                0,
                yi,
                "0 ",
                va="center",
                ha="right",
                fontsize=FONT_SIZE
            )

    # =======================================================
    # Right: Good / Moderate / Poor
    # =======================================================

    left_moderate = good

    left_poor = [
        g + m
        for g, m in zip(good, moderate)
    ]

    ax_csq.barh(
        y,
        good,
        height=0.70,
        color=COLOR_GOOD,
        edgecolor="white",
        linewidth=0.6,
        label="Good"
    )

    ax_csq.barh(
        y,
        moderate,
        left=left_moderate,
        height=0.70,
        color=COLOR_MODERATE,
        edgecolor="white",
        linewidth=0.6,
        label="Moderate"
    )

    ax_csq.barh(
        y,
        poor,
        left=left_poor,
        height=0.70,
        color=COLOR_POOR,
        edgecolor="white",
        linewidth=0.6,
        label="Poor"
    )

    ax_csq.set_xlim(0, 100)

    ax_csq.set_xlabel(
        "CSQ frequency (%)",
        fontsize=FONT_SIZE
    )

    # =======================================================
    # Add Good frequency (%) values
    # =======================================================

    for yi, value in zip(y, good):

        # Put label approximately at center of Good segment
        ax_csq.text(
            value / 2,
            yi,
            f"{value:.2f}%",
            ha="center",
            va="center",
            fontsize=FONT_SIZE
        )

    # =======================================================
    # Central dividing line
    # =======================================================

    ax_missing.spines["right"].set_visible(True)
    ax_missing.spines["right"].set_linewidth(1.0)
    ax_missing.spines["right"].set_color("black")

    ax_csq.spines["left"].set_visible(False)

    # =======================================================
    # Axis formatting
    # =======================================================

    # First record on top
    ax_missing.invert_yaxis()

    # y labels only on left
    ax_csq.tick_params(
        axis="y",
        left=False,
        labelleft=False
    )

    ax_missing.tick_params(
        axis="both",
        labelsize=FONT_SIZE
    )

    ax_csq.tick_params(
        axis="x",
        labelsize=FONT_SIZE
    )

    ax_missing.margins(y=0.02)
    ax_csq.margins(y=0.02)

    # =======================================================
    # Legend
    # =======================================================

    ax_csq.legend(
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.01),
        ncol=3,
        fontsize=FONT_SIZE
    )

    # =======================================================
    # Save
    # =======================================================
    if plot_title:
        fig.suptitle(
            plot_title,
            fontsize=12,
            x=0.0,
            y=0.98
        )
    pdf = output_prefix + ".pdf"
    png = output_prefix + ".png"
    svg = output_prefix + ".svg"

    plt.savefig(
        svg,
        bbox_inches="tight"
    )
    
    plt.savefig(
        pdf,
        bbox_inches="tight"
    )

    plt.savefig(
        png,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print("[INFO] Plot saved:")
    print(f"       {pdf}")
    print(f"       {png}")

def parse_args(args_list=None):

    parser = argparse.ArgumentParser(
        description="Merge and plot CSQ statistics."
    )

    parser.add_argument(
        "-i", "--input",
        default=".",
        help="Input directory containing *CSQ_stat.txt files"
    )

    parser.add_argument(
        "-o", "--output",
        default="CSQ_stat",
        help="Output prefix"
    )

    return parser.parse_args(args_list)

def main(args_list=None):
    args = parse_args(args_list)

    csqstat_dict, plot_title = read_all_csq_files(args.input)

    records = flatten_records(csqstat_dict)

    plot_csq_stats(
        records,
        args.output,
        plot_title
    )

if __name__ == "__main__":
    main()