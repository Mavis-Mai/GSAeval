#!/usr/bin/env python3
"""CSQ2sp batch v2: shared junction libraries and two-column manifests.

This version reuses the reporting/classification pipeline from CSQ2spbatch,
but replaces library construction and pairwise execution so the target library
is built once per batch and reference libraries are built or loaded once.
"""

import argparse
import os
import sys
from pathlib import Path
from urllib.parse import unquote

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

import CSQ2batch as legacy


_DATA_CACHE = {}
_PROTEIN_CACHE = {}
_LIBRARY_PATH_CACHE = {}


def _annotation_attributes(text):
    """Parse attributes and provide case-insensitive key lookup."""
    parsed = legacy.parse_gff_attributes(text)
    return {
        str(key).strip().lower(): unquote(str(value).strip().strip('"'))
        for key, value in parsed.items()
    }


def _first_attribute(attributes, names):
    for name in names:
        value = attributes.get(name.lower(), "").strip()
        if value:
            return value
    return ""


def _split_annotation_ids(value):
    return list(dict.fromkeys(
        item.strip() for item in str(value or "").split(",") if item.strip()
    ))


def _transcript_and_gene_ids(feature_type, attributes):
    """Extract IDs from GFF3/GTF transcript records and GTF child records."""
    transcript_keys = ("transcript_id", "transcriptid", "transcript")
    if feature_type in {"transcript", "mrna"}:
        transcript_keys += ("id",)
    transcript_id = _first_attribute(attributes, transcript_keys)
    gene_id = _first_attribute(
        attributes,
        ("gene_id", "geneid", "gene", "locus_tag", "parent"),
    )
    # On exon/CDS rows, a GFF3 Parent normally points to a transcript rather
    # than a gene. Only use Parent as gene ID for transcript/mRNA features.
    if feature_type not in {"transcript", "mrna"}:
        gene_id = _first_attribute(
            attributes, ("gene_id", "geneid", "gene", "locus_tag")
        )
    return transcript_id, _split_annotation_ids(gene_id)


def load_transcript_gene_map(filename):
    """Map transcripts to genes even when explicit gene features are absent."""
    mapping = {}
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
            feature_type = fields[2].strip().lower()
            attributes = _annotation_attributes(fields[8])
            transcript_id, gene_ids = _transcript_and_gene_ids(
                feature_type, attributes
            )
            if not transcript_id or not gene_ids:
                continue
            previous = _split_annotation_ids(mapping.get(transcript_id, ""))
            mapping[transcript_id] = ",".join(
                dict.fromkeys(previous + gene_ids)
            )
    if not mapping:
        print(
            f"[warning] No transcript-to-gene relationships found in {filename}",
            file=sys.stderr,
        )
    return mapping


def load_target_gene_ids(filename):
    """Collect explicit and transcript-derived gene IDs in annotation order."""
    explicit_gene_ids = []
    transcript_gene_ids = []
    explicit_seen, transcript_seen = set(), set()
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
            feature_type = fields[2].strip().lower()
            attributes = _annotation_attributes(fields[8])
            if feature_type == "gene":
                value = _first_attribute(
                    attributes,
                    ("gene_id", "geneid", "gene", "locus_tag", "id"),
                )
                for gene_id in _split_annotation_ids(value):
                    if gene_id not in explicit_seen:
                        explicit_seen.add(gene_id)
                        explicit_gene_ids.append(gene_id)
                continue
            _, gene_ids = _transcript_and_gene_ids(feature_type, attributes)
            for gene_id in gene_ids:
                if gene_id not in transcript_seen:
                    transcript_seen.add(gene_id)
                    transcript_gene_ids.append(gene_id)

    # Keep declared gene features first, then append loci recovered from
    # transcript/mRNA (or GTF child) rows that were not explicitly declared.
    result = list(explicit_gene_ids)
    result.extend(gene_id for gene_id in transcript_gene_ids if gene_id not in explicit_seen)
    return result


def _resolve_input(value, base_dir):
    if not value:
        return ""
    path = Path(value)
    if not path.is_absolute():
        path = base_dir / path
    return str(path.resolve())


def read_manifest(filename, dir_path="."):
    """Read either species+lib or species+fasta+gtf manifest layouts."""
    base_hint = Path(dir_path).resolve()
    candidate = base_hint / filename
    manifest_path = candidate if candidate.is_file() else Path(filename).resolve()
    table = pd.read_csv(manifest_path, sep="\t", dtype=str).fillna("")
    if "species" not in table.columns:
        raise ValueError("Manifest must contain a species column")
    if table.empty:
        raise ValueError("Manifest contains no comparison species")
    if table["species"].duplicated().any():
        duplicates = table.loc[table["species"].duplicated(), "species"].tolist()
        raise ValueError(f"Duplicate species names: {duplicates}")

    for column in ("lib", "fasta", "gtf", "genepairs"):
        if column not in table.columns:
            table[column] = ""
    base_dir = manifest_path.parent
    for column in ("lib", "fasta", "gtf", "genepairs"):
        table[column] = table[column].map(lambda value: _resolve_input(value, base_dir))

    for row in table.itertuples(index=False):
        if row.lib:
            if not Path(row.lib).is_file():
                raise FileNotFoundError(f"{row.species}: lib file not found: {row.lib}")
        elif row.fasta and row.gtf:
            for label, value in (("fasta", row.fasta), ("gtf", row.gtf)):
                if not Path(value).is_file():
                    raise FileNotFoundError(f"{row.species}: {label} file not found: {value}")
        else:
            raise ValueError(
                f"{row.species}: provide either lib, or both fasta and gtf"
            )
        if row.genepairs and not Path(row.genepairs).is_file():
            raise FileNotFoundError(
                f"{row.species}: genepairs file not found: {row.genepairs}"
            )
    return table


def _valid_file(path):
    path = Path(path)
    return path.is_file() and path.stat().st_size > 0


def _build_library(fasta, gtf, output_stem, args):
    from csq_oneStep import build_junciton_lib

    # build_junciton_lib -> parse_gtf -> print_debug expects this attribute.
    # Older batch parsers did not define it, so keep the direct-call path safe.
    if not hasattr(args, "debug"):
        args.debug = False
    output_stem = Path(output_stem).resolve()
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    output_file = Path(str(output_stem) + "_lib.tsv")
    if _valid_file(output_file) and not getattr(args, "rebuild_libs", False):
        print(f"[library] Reuse {output_file}")
        return output_file
    print(f"[library] Build {output_file}")
    result = build_junciton_lib(
        fasta, gtf, args.type, str(output_stem), args
    )
    result = Path(result).resolve()
    if not _valid_file(result):
        raise RuntimeError(f"Library build did not produce a valid file: {result}")
    return result


def _batch_root_from_pair_dir(pair_dir):
    return Path(pair_dir).resolve().parent.parent


def ensure_target_library(args, batch_root):
    key = ("target", str(Path(batch_root).resolve()))
    if key in _LIBRARY_PATH_CACHE:
        return _LIBRARY_PATH_CACHE[key]
    result = _build_library(
        args.target_fasta,
        args.target_gtf,
        batch_root / "libraries" / "target",
        args,
    )
    _LIBRARY_PATH_CACHE[key] = result
    return result


def ensure_reference_library(row, args, batch_root):
    if row.lib:
        return Path(row.lib).resolve()
    key = ("reference", str(Path(batch_root).resolve()), str(row.species))
    if key in _LIBRARY_PATH_CACHE:
        return _LIBRARY_PATH_CACHE[key]
    result = _build_library(
        row.fasta,
        row.gtf,
        batch_root / "libraries" / legacy.safe_name(row.species),
        args,
    )
    _LIBRARY_PATH_CACHE[key] = result
    return result


def _load_data(path):
    from loaddata import load_junction_data

    key = str(Path(path).resolve())
    if key not in _DATA_CACHE:
        _DATA_CACHE[key] = load_junction_data(key)
    return _DATA_CACHE[key]


def _protein_file(lib_path, output_stem):
    from csq_oneStep import getProt

    key = str(Path(lib_path).resolve())
    output = Path(str(output_stem) + ".fa").resolve()
    if key not in _PROTEIN_CACHE or not _valid_file(output):
        output.parent.mkdir(parents=True, exist_ok=True)
        generated = getProt(_load_data(lib_path), str(output.with_suffix("")), False)
        output = Path(generated).resolve()
        _PROTEIN_CACHE[key] = output
    return _PROTEIN_CACHE.get(key, output)


def run_pairwise(row, args, pair_dir):
    """Calculate one comparison without rebuilding either junction library."""
    from csq2sp_oneStep import csq_calculation2db, read_table, single_diamond_blast

    pair_dir = Path(pair_dir).resolve()
    pair_dir.mkdir(parents=True, exist_ok=True)
    batch_root = _batch_root_from_pair_dir(pair_dir)
    target_lib = ensure_target_library(args, batch_root)
    reference_lib = ensure_reference_library(row, args, batch_root)
    target_data = _load_data(target_lib)
    reference_data = _load_data(reference_lib)

    print(f"[CSQ v2] {args.target_name} vs {row.species}", flush=True)
    with legacy.working_directory(pair_dir):
        if row.genepairs:
            gene_pairs = read_table(row.genepairs)
        else:
            target_protein = _protein_file(
                target_lib, batch_root / "libraries" / "target_protein"
            )
            reference_protein = _protein_file(
                reference_lib,
                batch_root / "libraries" / f"{legacy.safe_name(row.species)}_protein",
            )
            gene_pairs = single_diamond_blast(
                str(target_protein), str(reference_protein), args.threads
            )
        output = csq_calculation2db(
            target_data, reference_data, gene_pairs, "output", args.weight_type
        )
    output = pair_dir / output
    if not _valid_file(output):
        raise RuntimeError(f"CSQ calculation did not produce a valid file: {output}")
    return output


def build_trees(genes, manifest, mappings, pair_dirs, target_name, tree_dir):
    """Build trees directly from shared libraries instead of per-pair copies."""
    if not genes:
        return []
    try:
        from Bio import Phylo
    except ImportError as exc:
        raise RuntimeError("Biopython is required for writing NJ trees") from exc
    # args and batch root are registered immediately before legacy.main runs.
    args = _ACTIVE_ARGS["args"]
    batch_root = Path(args.output_dir).resolve()
    target_sequences = legacy.load_lib_sequences(
        ensure_target_library(args, batch_root)
    )
    reference_sequences = {
        row.species: legacy.load_lib_sequences(
            ensure_reference_library(row, args, batch_root)
        )
        for row in manifest.itertuples(index=False)
    }
    tree_dir = Path(tree_dir)
    tree_dir.mkdir(parents=True, exist_ok=True)
    trees = []
    for gene in genes:
        names = [target_name]
        sequences = [legacy.protein(target_sequences.get(gene, ""))]
        for species in manifest["species"]:
            homolog = mappings[species].get(gene)
            names.append(species)
            sequences.append(
                legacy.protein(reference_sequences[species].get(homolog, ""))
            )
        if any(not sequence for sequence in sequences):
            print(f"[warning] Skip tree for {gene}: sequence missing", file=sys.stderr)
            continue
        tree = legacy.make_nj_tree(names, sequences)
        Phylo.write(tree, tree_dir / f"{legacy.safe_name(gene)}.nwk", "newick")
        trees.append((gene, tree))
    return trees


_ACTIVE_ARGS = {}
_legacy_parser = legacy.parser


def parser():
    result = _legacy_parser()
    for action in result._actions:
        if action.dest == "manifest":
            action.help = (
                "TSV with species+lib, or species+fasta+gtf; optional genepairs"
            )
    if not any(action.dest == "debug" for action in result._actions):
        result.add_argument(
            "--debug", action="store_true",
            help="Enable detailed library-building diagnostics.",
        )
    result.add_argument(
        "--rebuild-libs", action="store_true",
        help="Rebuild generated target/reference libraries instead of reusing cache.",
    )
    original_parse_args = result.parse_args

    def parse_args(args=None, namespace=None):
        parsed = original_parse_args(args, namespace)
        _ACTIVE_ARGS["args"] = parsed
        return parsed

    result.parse_args = parse_args
    return result


def main(args_list=None):
    legacy.read_manifest = read_manifest
    legacy.run_pairwise = run_pairwise
    legacy.build_trees = build_trees
    legacy.load_transcript_gene_map = load_transcript_gene_map
    legacy.load_target_gene_ids = load_target_gene_ids
    legacy.parser = parser
    return legacy.main(args_list)


if __name__ == "__main__":
    raise SystemExit(main())
