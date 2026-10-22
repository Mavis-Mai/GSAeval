#!/usr/bin/env python3
"""Plot the top 25% CSQ observations for sampled Good/Check genes."""

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


CLASS_COLORS = {"Good": "#F28E2B", "Check": "#2E8B57"}


def find_column(table, candidates):
    """Find a column while ignoring case, spaces, hyphens, and underscores."""
    normalize = lambda value: "".join(
        character for character in str(value).lower() if character.isalnum()
    )
    available = {normalize(column): column for column in table.columns}
    for candidate in candidates:
        if normalize(candidate) in available:
            return available[normalize(candidate)]
    return None


def read_inputs(input_dir):
    input_dir = Path(input_dir).resolve()
    summary_file = input_dir / "CSQ_summary.tsv"
    assessment_file = input_dir / "gene_assessment.tsv"
    groups_file = input_dir / "species_groups.tsv"
    thresholds_file = input_dir / "assessment_thresholds.tsv"
    for filename in (summary_file, assessment_file, groups_file, thresholds_file):
        if not filename.is_file():
            raise FileNotFoundError(f"Required file not found: {filename}")

    summary = pd.read_csv(summary_file, sep="\t", index_col=0)
    summary = summary.apply(pd.to_numeric, errors="coerce")
    assessments = pd.read_csv(assessment_file, sep="\t")
    label_column = find_column(assessments, ("Gene_label", "gene label", "Assessment"))
    phylogenetic_column = find_column(
        assessments,
        ("Gene_phylogenetic_label", "Gene phylogenetic label"),
    )
    transcript_column = find_column(assessments, ("Transcript", "Gene"))
    gene_id_column = find_column(assessments, ("GeneID", "gene_id"))
    if transcript_column is None or label_column is None:
        raise ValueError("gene_assessment.tsv must contain Transcript and gene label columns")
    if gene_id_column is None:
        raise ValueError("gene_assessment.tsv must contain GeneID")
    assessments = assessments.rename(
        columns={transcript_column: "Transcript", gene_id_column: "GeneID"}
    )
    assessments = assessments.rename(columns={label_column: "Assessment"})
    if phylogenetic_column is None:
        raise ValueError(
            "gene_assessment.tsv must contain Gene_phylogenetic_label"
        )
    assessments = assessments.rename(
        columns={phylogenetic_column: "Phylogenetic_assessment"}
    )
    groups = pd.read_csv(groups_file, sep="\t", dtype=str)
    if not {"Species", "Relatedness_group"}.issubset(groups.columns):
        raise ValueError(
            "species_groups.tsv must contain Species and Relatedness_group columns"
        )
    group_names = groups["Relatedness_group"].str.strip().str.lower()
    invalid_groups = sorted(set(group_names.dropna()) - {"near", "far"})
    if invalid_groups:
        raise ValueError(
            "species_groups.tsv contains invalid groups: " + ", ".join(invalid_groups)
        )
    # Binary positions let the existing opacity function render Near as solid
    # and Far as transparent without requiring a detailed distance table.
    distances = pd.Series(
        np.where(group_names == "near", 0.0, 1.0),
        index=groups["Species"].astype(str),
        name="Relatedness_position",
    )

    # Current batch output is genes x species. Accept the older transposed form too.
    known_genes = set(assessments["Transcript"].astype(str))
    row_matches = len(known_genes.intersection(summary.index.astype(str)))
    column_matches = len(known_genes.intersection(summary.columns.astype(str)))
    if column_matches > row_matches:
        summary = summary.T
    summary.index = summary.index.astype(str)
    summary.columns = summary.columns.astype(str)
    guides = read_threshold_guides(thresholds_file)
    return summary, assessments, distances, guides


def read_threshold_guides(filename):
    """Read the overall Good threshold from wide or row-oriented TSV."""
    table = pd.read_csv(filename, sep="\t")
    aliases = {
        "overall_threshold": (
            "overall_threshold", "ajusted_threshold", "adjusted_threshold",
            "global_high_threshold",
        ),
    }
    guides = {}
    normalized_columns = {str(column).strip().lower(): column for column in table.columns}
    for output_name, names in aliases.items():
        for name in names:
            if name in normalized_columns:
                values = pd.to_numeric(table[normalized_columns[name]], errors="coerce").dropna()
                if len(values):
                    guides[output_name] = float(values.iloc[0])
                    break

    if {"threshold", "final_value"}.issubset(normalized_columns):
        name_column = normalized_columns["threshold"]
        value_column = normalized_columns["final_value"]
        row_names = table[name_column].astype(str).str.strip().str.lower()
        for output_name, names in aliases.items():
            if output_name in guides:
                continue
            matched = table.loc[row_names.isin(names), value_column]
            values = pd.to_numeric(matched, errors="coerce").dropna()
            if len(values):
                guides[output_name] = float(values.iloc[0])

    missing = [name for name in aliases if name not in guides]
    if missing:
        raise ValueError(
            "assessment_thresholds.tsv is missing threshold values: " + ", ".join(missing)
        )
    return guides


def sample_genes(summary, assessments, number_per_class, seed):
    rng = np.random.default_rng(seed)
    selected = []
    assessment_map = assessments.drop_duplicates("Transcript").set_index("Transcript")["Assessment"]
    assessment_map.index = assessment_map.index.astype(str)
    for label in ("Good", "Check"):
        candidates = [
            gene for gene in assessment_map.index[assessment_map == label]
            if gene in summary.index and summary.loc[gene].notna().any()
        ]
        if not candidates:
            print(f"[warning] No plottable {label} genes were found")
            continue
        count = min(number_per_class, len(candidates))
        chosen = rng.choice(candidates, size=count, replace=False).tolist()
        selected.extend((gene, label) for gene in chosen)
    if not selected:
        raise ValueError("No Good or Check genes are available for plotting")
    return selected


def opacity_by_relatedness(species, distances, near_alpha, far_alpha):
    available = distances.reindex(species)
    finite = available.dropna()
    if finite.empty:
        return pd.Series((near_alpha + far_alpha) / 2, index=species)
    if math.isclose(float(finite.min()), float(finite.max())):
        only_value = float(finite.iloc[0])
        fixed_alpha = near_alpha if only_value <= 0.0 else far_alpha
        return pd.Series(fixed_alpha, index=species)
    normalized = (available - finite.min()) / (finite.max() - finite.min())
    # Distance zero maps to near_alpha (solid by default); the most distant
    # species maps to far_alpha (high transparency by default).
    alpha = near_alpha + normalized.pow(1.5) * (far_alpha - near_alpha)
    return alpha.fillna((near_alpha + far_alpha) / 2).clip(0.0, 1.0)


def plot(summary, selected, distances, guides, output_prefix, seed, near_alpha, far_alpha):
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError(
            "Plotting requires matplotlib; install the project requirements first"
        ) from exc
    records = []
    for order, (gene, label) in enumerate(selected):
        for species, csq in summary.loc[gene].items():
            if pd.notna(csq):
                records.append({
                    "Gene": gene, "Assessment": label, "Species": species,
                    "CSQ": float(csq), "Gene_order": order,
                })
    data = pd.DataFrame(records)
    # Select exactly ceil(25%) of the available observations for each gene.
    # Stable ranking makes ties reproducible instead of coloring more than 25%.
    data["Top_25_percent"] = False
    for gene, indices in data.groupby("Gene", sort=False).groups.items():
        ranked = data.loc[list(indices)].sort_values(
            ["CSQ", "Species"], ascending=[False, True], kind="mergesort"
        )
        top_count = max(1, int(math.ceil(len(ranked) * 0.25)))
        data.loc[ranked.index[:top_count], "Top_25_percent"] = True
    gene_order = [gene for gene, _ in selected]
    species_order = list(dict.fromkeys(data["Species"]))
    alpha = opacity_by_relatedness(species_order, distances, near_alpha, far_alpha)
    # Sample the vivid turbo map in golden-ratio order so adjacent species in
    # the input do not receive neighboring, easily confused colors.
    color_positions = np.mod(np.arange(len(species_order)) * 0.61803398875, 1.0)
    palette = {
        species: plt.colormaps["turbo"](0.04 + 0.92 * position)
        for species, position in zip(species_order, color_positions)
    }

    width = max(10.0, 0.70 * len(gene_order) + 3.0)
    fig, ax = plt.subplots(figsize=(width, 7.0))

    rng = np.random.default_rng(seed + 1)
    x_lookup = {gene: index for index, gene in enumerate(gene_order)}
    # Draw background observations first so highlighted top-quartile points
    # remain visible. Relatedness controls opacity in both layers.
    for species in species_order:
        subset = data[data["Species"] == species]
        x = np.array([x_lookup[gene] for gene in subset["Gene"]], dtype=float)
        x += rng.uniform(-0.16, 0.16, size=len(x))
        background = ~subset["Top_25_percent"].to_numpy()
        highlighted = subset["Top_25_percent"].to_numpy()
        if background.any():
            ax.scatter(
                x[background], subset.loc[background, "CSQ"], s=38,
                color="#A9A9A9", marker="o", alpha=float(alpha.loc[species]),
                edgecolors="none", zorder=2,
            )
        if not highlighted.any():
            continue
        ax.scatter(
            x[highlighted], subset.loc[highlighted, "CSQ"], s=58,
            color=palette[species], marker="o",
            alpha=float(alpha.loc[species]), edgecolors="#222222", linewidths=0.55,
            label=species, zorder=3,
        )

    # Keep a small strip to identify which sampled genes are Good or Check.
    for index, (_, label) in enumerate(selected):
        ax.plot([index - 0.28, index + 0.28], [-0.035, -0.035],
                transform=ax.get_xaxis_transform(), color=CLASS_COLORS[label],
                linewidth=5, solid_capstyle="butt", clip_on=False)

    ax.set_ylim(bottom=max(0.0, min(0.8, float(data["CSQ"].min()) - 0.05)), top=1.015)
    ax.set_xlim(-0.5, len(gene_order) - 0.5)
    ax.set_xticks(range(len(gene_order)))
    ax.set_xticklabels(gene_order, rotation=55, ha="right")
    ax.set_xlabel("Randomly sampled target genes")
    ax.set_ylabel("CSQ")
    ax.set_title("Top 25% CSQ observations of sampled Good and Check genes")
    ax.grid(axis="y", linestyle=":", alpha=0.35)

    high_guide = guides["overall_threshold"]
    ax.axhline(high_guide, color="#C44E52", linestyle="--", linewidth=1.8,
               alpha=0.95, zorder=2)

    class_handles = [
        plt.Line2D([0], [0], color=color, linewidth=7, label=label)
        for label, color in CLASS_COLORS.items()
    ]
    class_handles.extend([
        plt.Line2D([0], [0], color="#C44E52", linestyle="--", linewidth=1.8,
                   label=f"overall threshold = {high_guide:.3f}"),
        plt.Line2D([0], [0], linestyle="", marker="o", markersize=7,
                   markerfacecolor="#A9A9A9", markeredgecolor="none",
                   label="Other 75%"),
    ])
    first_legend = ax.legend(handles=class_handles, title="Assessment",
                             loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ax.add_artist(first_legend)
    species_handles = [
        plt.Line2D([0], [0], linestyle="", markersize=7,
                   marker="o", markerfacecolor=palette[species],
                   markeredgecolor="#222222",
                   alpha=float(alpha.loc[species]), label=species)
        for species in species_order
    ]
    ax.legend(handles=species_handles, title="Top 25% (solid near → transparent far)",
              loc="upper left", bbox_to_anchor=(1.01, 0.78), ncol=max(1, math.ceil(len(species_order) / 20)))
    fig.tight_layout()

    output_prefix = Path(output_prefix).resolve()
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output_prefix) + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(str(output_prefix) + ".svg", bbox_inches="tight")
    plt.close(fig)

    sampled = pd.DataFrame(selected, columns=["Gene", "Assessment"])
    sampled.to_csv(str(output_prefix) + ".sampled_genes.tsv", sep="\t", index=False)
    opacity_table = pd.DataFrame({
        "Species": species_order,
        "Relatedness_group": distances.reindex(species_order).map(
            {0.0: "Near", 1.0: "Far"}
        ).fillna("Unknown").values,
        "Point_alpha": alpha.reindex(species_order).values,
    })
    opacity_table.to_csv(str(output_prefix) + ".species_opacity.tsv", sep="\t", index=False, na_rep="NA")
    return output_prefix


def plot_label_summary(assessments, output_prefix):
    """Draw one-count-per-GeneID pie charts for the two assessment labels."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError(
            "Plotting requires matplotlib; install the project requirements first"
        ) from exc

    loci = assessments.copy()
    loci["GeneID"] = loci["GeneID"].fillna("NA").astype(str).str.strip().replace("", "NA")
    # Keep unmapped transcripts distinct; otherwise all literal NA values
    # would incorrectly collapse into one artificial locus.
    missing = loci["GeneID"].str.upper().eq("NA")
    loci.loc[missing, "GeneID"] = "unmapped:" + loci.loc[missing, "Transcript"].astype(str)

    def conservative_label(values, priority):
        observed = {str(value).strip() if pd.notna(value) else "NA" for value in values}
        observed.discard("")
        return next((label for label in priority if label in observed), "NA")

    loci = loci.groupby("GeneID", sort=False, as_index=False).agg(
        Assessment=("Assessment", lambda values: conservative_label(
            values, ("Check", "Good", "NA")
        )),
        Phylogenetic_assessment=("Phylogenetic_assessment", lambda values: conservative_label(
            values, ("Phylogenetic discordance", "discordance",
                     "Phylogenetically consistent", "consistent", "NA")
        )),
    )
    gene_counts = (
        loci["Assessment"].fillna("NA").astype(str).str.strip()
        .replace("", "NA").value_counts()
    )
    phylogenetic_counts = (
        loci["Phylogenetic_assessment"].fillna("NA").astype(str).str.strip()
        .replace("", "NA").value_counts()
    )
    gene_colors = {
        "Good": CLASS_COLORS["Good"],
        "Check": CLASS_COLORS["Check"],
        "NA": "#BDBDBD",
    }
    phylogenetic_colors = {
        "Phylogenetically consistent": "#4E79A7",
        "Phylogenetic discordance": "#E15759",
        "consistent": "#4E79A7",
        "discordance": "#E15759",
        "NA": "#BDBDBD",
    }

    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    panels = (
        (axes[0], gene_counts, gene_colors, "Gene label"),
        (axes[1], phylogenetic_counts, phylogenetic_colors,
         "Gene phylogenetic label"),
    )
    for axis, counts, colors, title in panels:
        labels = [f"{label}\n(n={count})" for label, count in counts.items()]
        axis.pie(
            counts.values,
            labels=labels,
            colors=[colors.get(label, "#9C9C9C") for label in counts.index],
            autopct="%1.1f%%",
            startangle=90,
            counterclock=False,
            wedgeprops={"edgecolor": "white", "linewidth": 1.2},
            textprops={"fontsize": 10},
        )
        axis.set_title(title, fontsize=14, fontweight="bold")
        axis.axis("equal")
    fig.suptitle(
        f"Gene-locus assessment overview (n={len(loci)})",
        fontsize=16,
        fontweight="bold",
    )
    fig.tight_layout()

    output_prefix = Path(output_prefix).resolve()
    summary_prefix = Path(str(output_prefix) + ".label_summary")
    summary_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(summary_prefix) + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(str(summary_prefix) + ".svg", bbox_inches="tight")
    plt.close(fig)
    return summary_prefix


def parser():
    p = argparse.ArgumentParser(
        description="Draw CSQ scatter plots with only each gene's top 25% highlighted.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-i", "--input-dir", required=True, help="CSQ2spbatch output directory")
    p.add_argument("-o", "--output-prefix", default="CSQ_Good_Check_scatter")
    p.add_argument("-n", "--number-per-class", type=int, default=10)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--near-alpha", type=float, default=1.0,
                   help="Point opacity for the closest species (solid by default)")
    p.add_argument("--far-alpha", type=float, default=0.20,
                   help="Point opacity for the most distant species")
    return p


def main(args_list=None):
    args = parser().parse_args(args_list)
    if args.number_per_class < 1:
        raise ValueError("number-per-class must be at least 1")
    if not 0 <= args.near_alpha <= 1 or not 0 <= args.far_alpha <= 1:
        raise ValueError("alpha values must be between 0 and 1")
    summary, assessments, distances, guides = read_inputs(args.input_dir)
    selected = sample_genes(summary, assessments, args.number_per_class, args.seed)
    output = plot(summary, selected, distances, guides, args.output_prefix, args.seed,
                  args.near_alpha, args.far_alpha)
    summary_output = plot_label_summary(assessments, args.output_prefix)
    print(f"Scatter plots written to {output}.png and {output}.svg")
    print(f"Label summary written to {summary_output}.png and {summary_output}.svg")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
