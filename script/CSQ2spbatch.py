#!/usr/bin/env python3
"""Batch CSQ comparison, conserved-gene NJ trees, and gene assessment."""

import argparse
import math
import os
import re
import sys
from contextlib import contextmanager
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))


REQUIRED_MANIFEST_COLUMNS = {"species", "fasta", "gtf"}


@contextmanager
def working_directory(path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def safe_name(value):
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    if not name or name in {".", ".."}:
        raise ValueError(f"Invalid species name: {value!r}")
    return name


def read_manifest(filename,dir_path='.'):
    dir_path = Path(dir_path).resolve()
    manifest_path = dir_path / filename if (dir_path / filename).exists() else Path(filename).resolve()

    table = pd.read_csv(manifest_path, sep="\t", dtype=str).fillna("")
    missing = REQUIRED_MANIFEST_COLUMNS - set(table.columns)
    if missing:
        raise ValueError("Manifest is missing columns: " + ", ".join(sorted(missing)))
    if table.empty:
        raise ValueError("Manifest contains no comparison species")
    if table["species"].duplicated().any():
        duplicates = table.loc[table["species"].duplicated(), "species"].tolist()
        raise ValueError(f"Duplicate species names: {duplicates}")
    for column in ("fasta", "gtf", "genepairs"):
        if column not in table:
            table[column] = ""
        # 修改3：优先在dir_path下查找文件，否则保持原行为(使用文件内的源路径)
        table[column] = table[column].map(
            lambda x: str((dir_path / x).resolve()) if (x and (dir_path / x).exists()) 
            else str(Path(x).resolve()) if x else ""
        )
    for row in table.itertuples(index=False):
        for column in ("fasta", "gtf"):
            value = getattr(row, column)
            if not Path(value).is_file():
                raise FileNotFoundError(f"{row.species}: {column} file not found: {value}")
        if row.genepairs and not Path(row.genepairs).is_file():
            raise FileNotFoundError(f"{row.species}: genepairs file not found: {row.genepairs}")
    return table


def validate_target(fasta, gtf):
    fasta, gtf = Path(fasta).resolve(), Path(gtf).resolve()
    if not fasta.is_file():
        raise FileNotFoundError(f"Target FASTA not found: {fasta}")
    if not gtf.is_file():
        raise FileNotFoundError(f"Target GTF/GFF not found: {gtf}")
    return str(fasta), str(gtf)


def parse_gff_attributes(text):
    """Parse either GFF3 key=value or GTF key "value" attributes."""
    attributes = {}
    for field in text.strip().strip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            key, value = field.split("=", 1)
        elif " " in field:
            key, value = field.split(None, 1)
            value = value.strip().strip('"')
        else:
            continue
        attributes[key.strip()] = value.strip()
    return attributes


def load_transcript_gene_map(filename):
    """Return transcript-to-gene relationships from a target GFF3/GTF file."""
    mapping = {}
    transcript_types = {"mrna", "transcript"}
    with open(filename, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) != 9:
                print(
                    f"[warning] Skip malformed GFF/GTF line {line_number}: {filename}",
                    file=sys.stderr,
                )
                continue
            feature_type = fields[2].lower()
            attributes = parse_gff_attributes(fields[8])
            transcript_id = attributes.get("transcript_id")
            gene_id = attributes.get("gene_id")
            if feature_type in transcript_types:
                transcript_id = transcript_id or attributes.get("ID")
                gene_id = gene_id or attributes.get("Parent")
            if not transcript_id or not gene_id:
                continue
            # Multiple parents are legal GFF3. Preserve all of them rather
            # than silently assigning the transcript to an arbitrary locus.
            mapping[transcript_id] = ",".join(
                dict.fromkeys(part.strip() for part in gene_id.split(",") if part.strip())
            )
    return mapping


def load_target_gene_ids(filename):
    """Return every gene locus declared by a target GFF3/GTF file.

    Gene features are preferred, while transcript parent/gene_id attributes
    provide a useful fallback for annotation files that omit explicit gene
    rows.  Input order is retained for reproducible output.
    """
    gene_ids = []
    seen = set()
    transcript_types = {"mrna", "transcript"}
    with open(filename, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) != 9:
                print(
                    f"[warning] Skip malformed GFF/GTF line {line_number}: {filename}",
                    file=sys.stderr,
                )
                continue
            feature_type = fields[2].lower()
            attributes = parse_gff_attributes(fields[8])
            candidates = []
            if feature_type == "gene":
                candidates = [attributes.get("gene_id") or attributes.get("ID")]
            elif feature_type in transcript_types:
                parent = attributes.get("gene_id") or attributes.get("Parent")
                candidates = parent.split(",") if parent else []
            for gene_id in candidates:
                gene_id = str(gene_id or "").strip()
                if gene_id and gene_id not in seen:
                    seen.add(gene_id)
                    gene_ids.append(gene_id)
    return gene_ids


def run_pairwise(row, args, pair_dir):
    from csq2sp_oneStep import main as csq2sp_main

    output_prefix = "output"
    command = [
        "-f", args.target_fasta, "-g", args.target_gtf,
        "-F", row.fasta, "-G", row.gtf,
        "-w", args.weight_type, "-d", args.type,
        "-t", str(args.threads), "-o", output_prefix,
        "--no-reciprocal",
    ]
    if row.genepairs:
        command.extend(["-l", row.genepairs])
    pair_dir.mkdir(parents=True, exist_ok=True)
    print(f"[CSQ] {args.target_name} vs {row.species}", flush=True)
    with working_directory(pair_dir):
        csq2sp_main(command)
    output = pair_dir / f"{output_prefix}.csq.tsv"
    if not output.is_file():
        raise RuntimeError(f"CSQ2sp did not create {output}")
    return output


def load_pair_result(filename, species):
    data = pd.read_csv(filename, sep="\t")
    required = {"GeneA", "GeneB", "CSQ"}
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"{filename} is missing columns: {sorted(missing)}")
    data["CSQ"] = pd.to_numeric(data["CSQ"], errors="coerce")
    data = data.dropna(subset=["GeneA", "GeneB", "CSQ"])
    # One target transcript can have multiple hits; retain the highest-CSQ hit.
    best = data.sort_values("CSQ", ascending=False).drop_duplicates("GeneA")
    return best.set_index("GeneA")["CSQ"].rename(species), best.set_index("GeneA")["GeneB"]


def load_lib_sequences(filename):
    sequences = {}
    with open(filename, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 5:
                raise ValueError(f"Malformed library line {line_number}: {filename}")
            sequences[fields[0]] = "".join(fields[4].split("|"))
    return sequences


def protein(sequence):
    try:
        from Bio.Seq import Seq
    except ImportError as exc:
        raise RuntimeError("Biopython is required for sequence translation and NJ trees") from exc
    sequence = sequence.upper().replace("U", "T")
    sequence = sequence[: len(sequence) - len(sequence) % 3]
    if not sequence:
        return ""
    return str(Seq(sequence).translate(to_stop=False)).rstrip("*").replace("*", "X")


def normalized_alignment_distance(seq_a, seq_b, aligner):
    if not seq_a or not seq_b:
        return 1.0
    score = aligner.score(seq_a, seq_b)
    self_score = max(aligner.score(seq_a, seq_a), aligner.score(seq_b, seq_b))
    return max(0.0, min(1.0, 1.0 - score / self_score)) if self_score else 1.0


def make_nj_tree(names, sequences):
    try:
        from Bio.Align import PairwiseAligner
        from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
    except ImportError as exc:
        raise RuntimeError("Biopython is required for NJ tree construction") from exc
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1.0
    aligner.mismatch_score = 0.0
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1
    matrix = []
    for i, seq_a in enumerate(sequences):
        matrix.append([
            normalized_alignment_distance(seq_a, sequences[j], aligner)
            for j in range(i + 1)
        ])
    return DistanceTreeConstructor().nj(DistanceMatrix(names, matrix))


def select_conserved(summary, count, threshold, max_std):
    complete = summary.dropna(axis=1, how="any")
    if complete.empty:
        return []
    stats = pd.DataFrame({"minimum": complete.min(), "mean": complete.mean(), "std": complete.std(ddof=0)})
    selected = stats[(stats["minimum"] > threshold) & (stats["std"] <= max_std)]
    selected = selected.sort_values(["mean", "std"], ascending=[False, True])
    return selected.head(count).index.tolist()


def tree_distance(tree, first, second):
    try:
        return float(tree.distance(first, second))
    except (ValueError, AttributeError):
        return math.nan


def build_trees(genes, manifest, mappings, pair_dirs, target_name, tree_dir):
    if not genes:
        return []
    try:
        from Bio import Phylo
    except ImportError as exc:
        raise RuntimeError("Biopython is required for writing NJ trees") from exc
    tree_dir.mkdir(parents=True, exist_ok=True)
    target_lib = load_lib_sequences(pair_dirs[manifest.iloc[0]["species"]] / "query_lib.tsv")
    species_order = manifest["species"].tolist()
    trees = []
    for gene in genes:
        names = [target_name]
        sequences = [protein(target_lib.get(gene, ""))]
        for species in species_order:
            ref_lib = load_lib_sequences(pair_dirs[species] / "ref_lib.tsv")
            homolog = mappings[species].get(gene)
            names.append(species)
            sequences.append(protein(ref_lib.get(homolog, "")))
        if any(not sequence for sequence in sequences):
            print(f"[warning] Skip tree for {gene}: sequence missing", file=sys.stderr)
            continue
        tree = make_nj_tree(names, sequences)
        Phylo.write(tree, tree_dir / f"{safe_name(gene)}.nwk", "newick")
        trees.append((gene, tree))
    return trees


def consensus_target_distances(trees, target_name, species_order):
    records = []
    for gene, tree in trees:
        records.append({"gene": gene, **{
            species: tree_distance(tree, target_name, species) for species in species_order
        }})
    distances = pd.DataFrame(records).set_index("gene") if records else pd.DataFrame(columns=species_order)
    consensus = distances.mean(axis=0, skipna=True).sort_values()
    return distances, consensus


def write_consensus_tree(trees, names, output_dir):
    """Average patristic distances across gene trees and write an NJ consensus."""
    try:
        from Bio import Phylo
        from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
    except ImportError as exc:
        raise RuntimeError("Biopython is required for consensus NJ tree construction") from exc
    full = pd.DataFrame(0.0, index=names, columns=names)
    for i, first in enumerate(names):
        for second in names[:i]:
            values = [tree_distance(tree, first, second) for _, tree in trees]
            values = [value for value in values if not math.isnan(value)]
            distance = sum(values) / len(values) if values else math.nan
            full.loc[first, second] = full.loc[second, first] = distance
    if full.isna().any().any():
        print("[warning] Consensus tree skipped because some species distances are unavailable", file=sys.stderr)
        return None
    lower = [[float(full.iloc[i, j]) for j in range(i + 1)] for i in range(len(names))]
    tree = DistanceTreeConstructor().nj(DistanceMatrix(names, lower))
    Phylo.write(tree, output_dir / "species_consensus.nwk", "newick")
    return tree


def consensus_tree_target_distances(tree, target_name, species_order):
    """Calculate target-to-species patristic distances on the consensus NJ tree."""
    if tree is None:
        return pd.Series(index=species_order, dtype=float, name="Consensus_NJ_distance")
    values = {
        species: tree_distance(tree, target_name, species)
        for species in species_order
    }
    return pd.Series(values, dtype=float, name="Consensus_NJ_distance").sort_values()

def relatedness_groups(summary, nwk_file, target_name):
    """Select topology-consistent Near and Far groups from an NJ tree.

    Near species are collected from successive small sister branches along
    the path leading away from the target. The largest child branch at each
    node is treated as the continuation of the backbone.

    Far species are selected as a complete clade containing one of the
    topologically deepest species relative to the target.

    Species in the middle are excluded from both groups.
    """
    from Bio import Phylo
    import copy
    import math

    tree = Phylo.read(str(nwk_file), "newick")

    tree_species = {
        terminal.name
        for terminal in tree.get_terminals()
        if terminal.name
    }

    if target_name not in tree_species:
        raise ValueError(
            f"Target species {target_name!r} was not found "
            f"in NJ tree: {nwk_file}"
        )

    valid_species = {
        str(species)
        for species in summary.index
        if str(species) in tree_species
        and str(species) != target_name
    }

    species_count = len(valid_species)

    if species_count == 0:
        return [], []

    if species_count == 1:
        return sorted(valid_species), []

    quarter_count = max(
        1,
        int(math.ceil(species_count * 0.25)),
    )

    # Work on a copy so the original tree is not changed.
    rooted_tree = copy.deepcopy(tree)

    target_terminal = next(
        terminal
        for terminal in rooted_tree.get_terminals()
        if terminal.name == target_name
    )

    rooted_tree.root_with_outgroup(target_terminal)

    # Find the target again because rerooting changes the tree structure.
    target_terminal = next(
        terminal
        for terminal in rooted_tree.get_terminals()
        if terminal.name == target_name
    )

    def valid_terminals(clade):
        """Return valid comparison species below a clade."""
        return [
            terminal.name
            for terminal in clade.get_terminals()
            if terminal.name in valid_species
        ]

    def mean_target_distance(species_names):
        """Mean patristic distance from the target."""
        if not species_names:
            return math.inf

        return sum(
            float(rooted_tree.distance(target_terminal, species))
            for species in species_names
        ) / len(species_names)

    # Build child-to-parent relationships after rerooting.
    parent_map = {}

    for parent in rooted_tree.find_clades(order="level"):
        for child in parent.clades:
            parent_map[child] = parent

    # ---------------------------------------------------------------
    # Near group
    # ---------------------------------------------------------------

    # Rerooting on the target creates:
    #
    # root
    # ├── target
    # └── all remaining species
    #
    # Therefore, begin inside the non-target remainder rather than treating
    # the complete remainder as the first Near branch.
    root_children = [
        child
        for child in rooted_tree.root.clades
        if target_name not in {
            terminal.name
            for terminal in child.get_terminals()
        }
    ]

    if not root_children:
        return [], []

    # There should normally be one remainder child after rerooting.
    backbone = max(
        root_children,
        key=lambda clade: len(valid_terminals(clade)),
    )

    near_species = []
    near_set = set()

    while backbone is not None:
        child_branches = [
            child
            for child in backbone.clades
            if valid_terminals(child)
        ]

        if not child_branches:
            break

        if len(child_branches) == 1:
            backbone = child_branches[0]
            continue

        # The largest child is treated as the continuing backbone. Smaller
        # children represent successive target-side sister branches.
        continuation = max(
            child_branches,
            key=lambda child: (
                len(valid_terminals(child)),
                mean_target_distance(valid_terminals(child)),
            ),
        )

        sister_branches = [
            child
            for child in child_branches
            if child is not continuation
        ]

        # For multifurcating nodes, process closer sister branches first.
        sister_branches.sort(
            key=lambda child: (
                mean_target_distance(valid_terminals(child)),
                len(valid_terminals(child)),
            )
        )

        added_at_this_level = False

        for sister in sister_branches:
            sister_species = [
                species
                for species in valid_terminals(sister)
                if species not in near_set
            ]

            if not sister_species:
                continue

            proposed_size = len(near_species) + len(sister_species)

            current_difference = (
                abs(len(near_species) - quarter_count)
                if near_species
                else math.inf
            )
            proposed_difference = abs(
                proposed_size - quarter_count
            )

            # Always accept the first genuine sister branch. Afterwards,
            # accept a complete branch only if it makes the group size
            # closer to the desired quarter.
            if (
                not near_species
                or proposed_difference <= current_difference
            ):
                near_species.extend(sister_species)
                near_set.update(sister_species)
                added_at_this_level = True
            else:
                break

            if len(near_species) >= quarter_count:
                break

        if len(near_species) >= quarter_count:
            break

        # Continue moving away from the target through the largest branch.
        backbone = continuation

        # Stop if adding the sister branches at this level would move the
        # group farther from the target size.
        if sister_branches and not added_at_this_level:
            break

    # ---------------------------------------------------------------
    # Far group
    # ---------------------------------------------------------------

    available_terminals = [
        terminal
        for terminal in rooted_tree.get_terminals()
        if terminal.name in valid_species
        and terminal.name not in near_set
    ]

    if not available_terminals:
        return near_species, []

    def topology_depth(terminal):
        """Number of tree edges from the target-rooted root."""
        return len(rooted_tree.get_path(terminal))

    # Find a terminal at the deepest topology level. Use patristic distance
    # only to resolve ties.
    farthest_terminal = max(
        available_terminals,
        key=lambda terminal: (
            topology_depth(terminal),
            float(
                rooted_tree.distance(
                    target_terminal,
                    terminal,
                )
            ),
        ),
    )

    far_candidates = []
    current_clade = farthest_terminal

    while True:
        candidate_species = valid_terminals(current_clade)
        candidate_set = set(candidate_species)

        if (
            candidate_species
            and candidate_set.isdisjoint(near_set)
        ):
            far_candidates.append({
                "species": candidate_species,
                "size_difference": abs(
                    len(candidate_species) - quarter_count
                ),
                "mean_distance": mean_target_distance(
                    candidate_species
                ),
            })

        if current_clade not in parent_map:
            break

        parent = parent_map[current_clade]
        parent_species = set(valid_terminals(parent))

        # Do not expand a Far clade into the Near group.
        if not parent_species.isdisjoint(near_set):
            break

        current_clade = parent

    if not far_candidates:
        return near_species, []

    far_candidate = min(
        far_candidates,
        key=lambda candidate: (
            candidate["size_difference"],
            -candidate["mean_distance"],
        ),
    )

    near_species = sorted(
        near_species,
        key=lambda species: float(
            rooted_tree.distance(target_terminal, species)
        ),
    )

    far_species = sorted(
        far_candidate["species"],
        key=lambda species: float(
            rooted_tree.distance(target_terminal, species)
        ),
        reverse=True,
    )

    return near_species, far_species
    
 

def calculate_assessment_thresholds(summary, consensus, weight_type,nwk_file,target_name):
    """Return the fixed grading policy and phylogenetic species groups."""
    overall_threshold = 0.90 if weight_type == "plant" else 0.95
    thresholds = {
        "overall_threshold": overall_threshold,
        "missing_fraction": 0.60,
        "high_species_score": 0.98,
        "high_species_count": 2,
    }
    near_species, far_species = relatedness_groups(
            summary,
            nwk_file,
            target_name,
        )

    records = [
        {"threshold": "overall_threshold", "method": f"fixed_{weight_type}_policy",
         "final_value": overall_threshold},
        {"threshold": "missing_fraction", "method": "classification_policy",
         "final_value": 0.60},
        {"threshold": "high_species_score", "method": "classification_policy",
         "final_value": 0.98},
        {"threshold": "high_species_count", "method": "classification_policy",
         "final_value": 2},
    ]
    return thresholds, pd.DataFrame(records), near_species, far_species


def classify_genes(
    summary,
    consensus,
    args,
    thresholds=None,
    near_species=None,
    far_species=None,
    transcript_gene_map=None,
):
    """Classify genes and evaluate whether CSQ follows relatedness."""
    if thresholds is None:
        raise ValueError(
            "thresholds must be supplied by "
            "calculate_assessment_thresholds()"
        )

    if near_species is None or far_species is None:
        raise ValueError(
            "near_species and far_species must be supplied"
        )
    overall_threshold = thresholds["overall_threshold"]
    missing_threshold = thresholds["missing_fraction"]
    high_species_score = thresholds["high_species_score"]
    high_species_count = int(thresholds["high_species_count"])
    rows = []

    for gene in summary.columns:
        scores = pd.to_numeric(summary[gene], errors="coerce")
        valid = scores.dropna()
        missing_fraction = float(scores.isna().mean())
        gene_q1 = float(valid.quantile(0.25)) if len(valid) else math.nan
        gene_q2 = float(valid.quantile(0.50)) if len(valid) else math.nan
        gene_q3 = float(valid.quantile(0.75)) if len(valid) else math.nan
        high_species = [
            species for species, value in scores.items()
            if pd.notna(value) and float(value) >= high_species_score
        ]
        near_valid = scores.reindex(near_species).dropna()
        far_valid = scores.reindex(far_species).dropna()

        # Two independently near-perfect species are decisive even when much
        # of the batch is missing.  Otherwise, missingness is evaluated before
        # Q3 so a few remaining high observations cannot rescue a sparse gene.
        if len(high_species) >= high_species_count:
            label = "Good"
        elif missing_fraction >= missing_threshold:
            label = "NA"
        elif pd.notna(gene_q3) and gene_q3 > overall_threshold:
            label = "Good"
        else:
            label = "Check"

        phylogenetic_label = (
            "consistent"
            if len(near_valid) and len(far_valid)
            and float(near_valid.mean()) >= float(far_valid.mean())
            else "discordance"
        )

        rows.append({
            "Transcript": gene,
            "GeneID": (transcript_gene_map or {}).get(gene, "NA"),
            "CSQ-Q1": gene_q1,
            "CSQ-Q2": gene_q2,
            "CSQ-Q3": gene_q3,
            "Gene_label": label,
            "Gene_phylogenetic_label": phylogenetic_label,
        })
    return pd.DataFrame(rows)

def summarize_gene_loci(assessments, all_gene_ids=None):
    """Summarize transcript labels at the GeneID level.

    Rules:
      1. Any transcript is Good -> Good
      2. All transcripts are Check -> Check
      3. All transcripts are NA -> NA
      4. Mixed Check and NA -> NA
    """
    required_columns = {"GeneID", "Gene_label", "CSQ-Q3"}
    missing = required_columns - set(assessments.columns)
    if missing:
        raise ValueError(
            "Cannot summarize gene loci; missing columns: "
            + ", ".join(sorted(missing))
        )

    def summarize_labels(values):
        labels = {
            "NA" if pd.isna(value) or not str(value).strip()
            else str(value).strip()
            for value in values
        }

        if "Good" in labels:
            return "Good"

        if labels == {"Check"}:
            return "Check"

        if labels == {"NA"}:
            return "NA"

        # Mixed Check and NA.
        return "NA"

    def highest_q3(values):
        numeric = pd.to_numeric(values, errors="coerce").dropna()
        return float(numeric.max()) if len(numeric) else math.nan

    cleaned = assessments.assign(
        GeneID=assessments["GeneID"]
        .fillna("NA")
        .astype(str)
        .str.strip()
        .replace("", "NA"),
        Gene_label=assessments["Gene_label"]
        .fillna("NA")
        .astype(str)
        .str.strip()
        .replace("", "NA"),
    )

    locus_assessment = (
        cleaned
        .groupby("GeneID", as_index=False, sort=False)
        .agg(label=("Gene_label", summarize_labels))
    )

    # Select Q3 from transcripts that support the final locus label, then use
    # the highest available value.  This preserves the highest Good
    # transcript and also retains a measurable Q3 for Check/NA loci.
    grouped_assessments = {
        gene_id: group for gene_id, group in cleaned.groupby("GeneID", sort=False)
    }
    q3_values = []
    for row in locus_assessment.itertuples(index=False):
        group = grouped_assessments[row.GeneID]
        matching = group[group["Gene_label"] == row.label]
        q3_values.append(highest_q3(matching["CSQ-Q3"]))
    locus_assessment["CSQ-Q3"] = q3_values

    if all_gene_ids is not None:
        complete = pd.DataFrame({"GeneID": list(dict.fromkeys(all_gene_ids))})
        locus_assessment = complete.merge(locus_assessment, on="GeneID", how="left")
        locus_assessment["label"] = locus_assessment["label"].fillna(
            "Missing homologous genes"
        )

    return locus_assessment


def plot_locus_outputs(locus_assessment, assessments, threshold, output_dir):
    """Write Gene_label proportions, a pie chart, and a locus Q3 rank plot."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError(
            "Plotting requires matplotlib; install the project requirements first"
        ) from exc

    output_dir = Path(output_dir)
    # Match plotCSQ2batch.py: count each GeneID once and give Good priority
    # when a locus has multiple transcript assessments.  Keep transcripts
    # without a GeneID distinct instead of collapsing every NA ID together.
    gene_rows = assessments.copy()
    gene_rows["GeneID"] = (
        gene_rows["GeneID"].fillna("NA").astype(str).str.strip().replace("", "NA")
    )
    missing_gene_id = gene_rows["GeneID"].str.upper().eq("NA")
    gene_rows.loc[missing_gene_id, "GeneID"] = (
        "unmapped:" + gene_rows.loc[missing_gene_id, "Transcript"].astype(str)
    )

    def gene_label(values):
        observed = {
            "NA" if pd.isna(value) or not str(value).strip() else str(value).strip()
            for value in values
        }
        return next((label for label in ("Good", "Check", "NA") if label in observed), "NA")

    gene_labels = gene_rows.groupby("GeneID", sort=False)["Gene_label"].agg(gene_label)
    counts = gene_labels.value_counts().reindex(["Good", "Check", "NA"], fill_value=0)
    total = int(counts.sum())
    proportions = pd.DataFrame({
        "label": counts.index,
        "count": counts.values,
        "percentage": (counts.values / total * 100.0) if total else 0.0,
    })
    proportions.to_csv(
        output_dir / "gene_locus_proportions.tsv", sep="\t", index=False,
        float_format="%.4f",
    )

    if total:
        fig, ax = plt.subplots(figsize=(7, 6))
        nonzero = counts[counts > 0]
        colors = {"Good": "#a4d8f4", "Check": "#f7c78c", "NA": "#BDBDBD"}
        ax.pie(
            nonzero.values,
            labels=[f"{label}\n(n={count})" for label, count in nonzero.items()],
            colors=[colors[label] for label in nonzero.index],
            autopct="%1.1f%%", startangle=90, counterclock=False,
            wedgeprops={"edgecolor": "white", "linewidth": 1.2},
        )
        ax.set_title(f"Gene_label classification (n={total})")
        ax.axis("equal")
        fig.tight_layout()
        fig.savefig(output_dir / "gene_locus_proportions.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    ranked = locus_assessment.copy()
    ranked["CSQ-Q3"] = pd.to_numeric(ranked["CSQ-Q3"], errors="coerce")
    ranked = ranked.dropna(subset=["CSQ-Q3"]).sort_values(
        "CSQ-Q3", ascending=False, kind="mergesort"
    ).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(max(8.0, min(18.0, len(ranked) / 150.0)), 6))
    if len(ranked):
        ranks = ranked.index.to_numpy() + 1
        above = ranked["CSQ-Q3"].to_numpy() > threshold
        ax.scatter(ranks[~above], ranked.loc[~above, "CSQ-Q3"], s=9,
                   color="#9E9E9E", edgecolors="none", label="Below threshold")
        ax.scatter(ranks[above], ranked.loc[above, "CSQ-Q3"], s=11,
                   color="#D62728", edgecolors="none", label="Above threshold")
    ax.axhline(threshold, color="#D62728", linestyle="--", linewidth=1.4,
               label=f"Threshold = {threshold:.3f}")
    ax.set_xlabel("Gene rank (descending CSQ-Q3)")
    ax.set_ylabel("CSQ-Q3")
    ax.set_title("Gene-locus CSQ-Q3 rank")
    ax.grid(axis="y", linestyle=":", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "generank.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print("Gene_label proportions:")
    for row in proportions.itertuples(index=False):
        print(f"  {row.label}: {row.count} ({row.percentage:.2f}%)")

def filter_suspected_same_species(summary, fraction_threshold=0.95, perfect_csq=0.999):
    """Remove species whose near-perfect CSQ covers almost all target genes."""
    if not 0 < fraction_threshold <= 1:
        raise ValueError("duplicate-species-perfect-fraction must be in (0, 1]")
    if not 0 <= perfect_csq <= 1:
        raise ValueError("duplicate-species-perfect-csq must be between 0 and 1")
    total_genes = len(summary.columns)
    if total_genes == 0:
        raise ValueError("CSQ summary contains no target genes")

    records, excluded = [], []
    for species in summary.index:
        scores = pd.to_numeric(summary.loc[species], errors="coerce")
        valid_count = int(scores.notna().sum())
        perfect_count = int((scores >= perfect_csq).sum())
        # Use all target genes as denominator. Missing values are not perfect.
        perfect_fraction = perfect_count / total_genes
        is_excluded = perfect_fraction >= fraction_threshold
        if is_excluded:
            excluded.append(species)
        records.append({
            "Species": species,
            "Total_target_genes": total_genes,
            "Valid_CSQ_count": valid_count,
            "Valid_CSQ_fraction": valid_count / total_genes,
            "Near_perfect_CSQ_count": perfect_count,
            "Near_perfect_CSQ_fraction": perfect_fraction,
            "Perfect_CSQ_cutoff": perfect_csq,
            "Exclusion_fraction_threshold": fraction_threshold,
            "Status": "Excluded_suspected_same_species" if is_excluded else "Retained",
        })
    return summary.drop(index=excluded), excluded, pd.DataFrame(records)


def exclude_manifest_species(manifest, requested_values):
    """Remove explicitly requested species before any pairwise calculation."""
    requested = []
    for value in requested_values or []:
        requested.extend(
            item.strip() for item in str(value).split(",") if item.strip()
        )
    # Preserve user order while removing duplicates.
    requested = list(dict.fromkeys(requested))
    if not requested:
        return manifest, []

    available = manifest["species"].astype(str).tolist()
    unknown = [species for species in requested if species not in available]
    if unknown:
        raise ValueError(
            "Requested exclusion species not found in manifest: "
            + ", ".join(unknown)
            + ". Available species: "
            + ", ".join(available)
        )
    filtered = manifest.loc[
        ~manifest["species"].isin(requested)
    ].reset_index(drop=True)
    if filtered.empty:
        raise ValueError("All manifest species were removed by --excusion-species")
    return filtered, requested


def parser():
    p = argparse.ArgumentParser(
        description="Run CSQ2sp against multiple species, build conserved-gene NJ trees, and assess genes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-f", "--target-fasta", required=True)
    p.add_argument("-g", "--target-gtf", required=True)
    p.add_argument("-w", "--weight-type", required=True, choices=["plant", "animal"])
    p.add_argument("-n", "--target-name", default="target")
    p.add_argument("-m", "--manifest", help="TSV columns: species, fasta, gtf, optional genepairs")
    p.add_argument("-o", "--output-dir", default="CSQ2spbatch_output")
    p.add_argument("-t", "--threads", type=int, default=8)
    p.add_argument("-d", "--type", choices=["CDS", "exon"], default="CDS")
    p.add_argument("--reuse", action="store_true", help="Reuse existing pairwise output.csq.tsv files")
    p.add_argument(
        "--excusion-species", "--exclude-species", action="append", default=[],
        metavar="SPECIES",
        help=(
            "Remove species from the manifest before CSQ2sp calculation. "
            "Repeat the option or provide comma-separated names."
        ),
    )
    p.add_argument(
        "--duplicate-species-perfect-fraction", type=float, default=0.95,
        help="Exclude a species when this fraction of all target genes has near-perfect CSQ.",
    )
    p.add_argument(
        "--duplicate-species-perfect-csq", type=float, default=0.999,
        help="CSQ treated as effectively equal to 1 for duplicate-species detection.",
    )
    return p

def main(args_list=None):
    args = parser().parse_args(args_list)
    args.target_fasta, args.target_gtf = validate_target(args.target_fasta, args.target_gtf)
    args.target_name = safe_name(args.target_name)
    if args.manifest:
        manifest = read_manifest(args.manifest)
    else:
        dirpath = Path(__file__).parent.parent
        datasetpath=os.path.join(dirpath,"Genome_dataset",args.weight_type)
        manifestfile=os.path.join(datasetpath,"species_manifest.tsv")
        manifest = read_manifest(manifestfile,datasetpath)

    manifest, manually_excluded_species = exclude_manifest_species(
        manifest, args.excusion_species
    )
    if args.target_name in set(manifest["species"]):
        raise ValueError("Target name must differ from every comparison species name")
    output_dir = Path(args.output_dir).resolve()
    pair_root = output_dir / "pairwise"
    pair_root.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({
        "Species": manually_excluded_species,
        "Reason": "Excluded_by_command_line",
    }).to_csv(
        output_dir / "manually_excluded_species.tsv", sep="\t", index=False
    )
    if manually_excluded_species:
        print(
            "[info] Species excluded before CSQ2sp calculation: "
            + ", ".join(manually_excluded_species)
        )
    # These detailed distance tables were produced by older versions.  The
    # consensus tree plus the concise near/far group table supersede them.
    for obsolete_name in ("species_relatedness.tsv", "species_distance_matrix.tsv"):
        obsolete_file = output_dir / obsolete_name
        if obsolete_file.is_file():
            obsolete_file.unlink()

    scores, mappings, pair_dirs = [], {}, {}
    for row in manifest.itertuples(index=False):
        pair_dir = pair_root / safe_name(row.species)
        pair_dirs[row.species] = pair_dir
        result = pair_dir / "output.csq.tsv"
        if not (args.reuse and result.is_file()):
            result = run_pairwise(row, args, pair_dir)
        score, mapping = load_pair_result(result, row.species)
        scores.append(score)
        mappings[row.species] = mapping

    summary = pd.concat(scores, axis=1).T
    summary.index.name = "Species"
    summary.columns.name = "Target_gene"
    summary, excluded_species, similarity_check = filter_suspected_same_species(
        summary,
        args.duplicate_species_perfect_fraction,
        args.duplicate_species_perfect_csq,
    )
    similarity_check.to_csv(
        output_dir / "species_similarity_check.tsv", sep="\t", index=False
    )
    if excluded_species:
        print(
            "[warning] Excluded suspected same-species comparisons: "
            + ", ".join(excluded_species),
            file=sys.stderr,
        )
        manifest = manifest.loc[
            ~manifest["species"].isin(excluded_species)
        ].reset_index(drop=True)
    if summary.empty:
        raise ValueError(
            "All comparison species were excluded as suspected same-species matches; "
            "inspect species_similarity_check.tsv or adjust the duplicate-species thresholds"
        )
    summary_tmp = summary.T
    summary_tmp.to_csv(output_dir / "CSQ_summary.tsv", sep="\t", na_rep="NA")
    
    conserved_count, conserved_threshold, conserved_max_std = 20, 0.95, 0.05
    genes = select_conserved(
        summary, conserved_count, conserved_threshold, conserved_max_std
    )
    pd.DataFrame({"Gene": genes}).to_csv(output_dir / "conserved_genes.tsv", sep="\t", index=False)
    if len(genes) < conserved_count:
        print(f"[warning] Only {len(genes)} genes satisfy the conserved-gene criteria", file=sys.stderr)

    trees = build_trees(genes, manifest, mappings, pair_dirs, args.target_name, output_dir / "NJ_trees")
    consensus_tree = None
    if trees:
        consensus_tree = write_consensus_tree(
            trees, [args.target_name] + manifest["species"].tolist(), output_dir
        )
    tree_distances, mean_gene_tree_distances = consensus_target_distances(
        trees, args.target_name, manifest["species"].tolist()
    )
    tree_distances.to_csv(output_dir / "NJ_target_distances.tsv", sep="\t", na_rep="NA")
    consensus = consensus_tree_target_distances(
        consensus_tree, args.target_name, manifest["species"].tolist()
    )
    if consensus.dropna().empty:
        print(
            "[warning] Consensus-tree distances are unavailable; using mean "
            "gene-tree distances for relatedness grouping",
            file=sys.stderr,
        )
        consensus = mean_gene_tree_distances.rename("Mean_NJ_distance")
    nwk_file = output_dir / "species_consensus.nwk"
    if not nwk_file.is_file():
        raise FileNotFoundError(
            "Consensus NJ tree is required for topology-based Near/Far "
            f"grouping, but it was not generated: {nwk_file}"
        )
    thresholds, threshold_table, near_species, far_species = calculate_assessment_thresholds(
        summary,
        consensus,
        args.weight_type,
        nwk_file,
        args.target_name,
    )
    species_groups = pd.DataFrame(
        [(species, "Near") for species in near_species]
        + [(species, "Far") for species in far_species],
        columns=["Species", "Relatedness_group"],
    )
    species_groups.to_csv(
        output_dir / "species_groups.tsv", sep="\t", index=False
    )
    threshold_table["near_species"] = ",".join(near_species)
    threshold_table["far_species"] = ",".join(far_species)
    threshold_table.to_csv(
        output_dir / "assessment_thresholds.tsv", sep="\t", index=False, na_rep="NA"
    )
    transcript_gene_map = load_transcript_gene_map(args.target_gtf)
    assessments = classify_genes(
        summary,
        consensus,
        args,
        thresholds=thresholds,
        near_species=near_species,
        far_species=far_species,
        transcript_gene_map=transcript_gene_map,
    )
    missing_gene_ids = int(assessments["GeneID"].eq("NA").sum())
    if missing_gene_ids:
        print(
            f"[warning] No GeneID mapping found for {missing_gene_ids} target transcripts",
            file=sys.stderr,
        )
    assessments.to_csv(output_dir / "gene_assessment.tsv", sep="\t", index=False, na_rep="NA")
    all_gene_ids = load_target_gene_ids(args.target_gtf)
    if not all_gene_ids:
        raise ValueError(
            "No gene loci could be read from the target GFF3/GTF file: "
            f"{args.target_gtf}"
        )
    locus_assessment = summarize_gene_loci(assessments, all_gene_ids)

    locus_assessment.to_csv(
        output_dir / "gene_locus_classify.tsv",
        sep="\t",
        index=False,
        na_rep="",
    )
    plot_locus_outputs(
        locus_assessment,
        assessments,
        thresholds["overall_threshold"],
        output_dir,
    )
    print(f"Done. Results: {output_dir}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
