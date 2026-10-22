"""Introduce exon errors in GTF and label their measured coding consequences.

Modified transcripts receive a newly predicted longest ATG-initiated ORF.
Unmodified transcripts retain their original CDS annotation.
"""
import argparse
import copy
import csv
import random
import re
import sys
import time
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

try:
    import numpy as np
except ImportError:
    np = None


# Standard nuclear genetic code, codons ordered T,C,A,G at each position.
_CODONS = [a+b+c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG']
_CODE = dict(zip(_CODONS, 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))


@lru_cache(maxsize=4096)
def translate(dna):
    peptide = ''.join(_CODE.get(dna[i:i+3], 'X') for i in range(0, len(dna)-2, 3))
    return peptide.split('*', 1)[0]


def _protein_alignment_python(a, b):
    """Global alignment: match +1, mismatch/gap -1; identity includes gap columns.

    Ties prefer diagonal, then deletion, then insertion. Linear memory, O(n*m) time.
    """
    if a == b:
        return 100.0, 100.0 if a else 0.0, 100.0 if b else 0.0
    if not a or not b:
        return 0.0, 0.0, 0.0
    previous = [(-j, 0, j, 0) for j in range(len(b)+1)]
    for i, aa in enumerate(a, 1):
        current = [(-i, 0, i, 0)]
        for j, bb in enumerate(b, 1):
            s, matches, length, paired = previous[j-1]
            diagonal = (s + (1 if aa == bb else -1), matches + (aa == bb), length+1, paired+1)
            s, matches, length, paired = previous[j]
            deletion = (s-1, matches, length+1, paired)
            s, matches, length, paired = current[-1]
            insertion = (s-1, matches, length+1, paired)
            current.append(max((diagonal, deletion, insertion), key=lambda x: x[0]))
        previous = current
    return (100.0 * previous[-1][1] / previous[-1][2],
            100.0 * previous[-1][3] / len(a), 100.0 * previous[-1][3] / len(b))


def _protein_alignment_numpy(a, b):
    """Exact same recurrence/ties, vectorized one row at a time, O(m) memory."""
    width = len(b)+1
    indices = np.arange(width, dtype=np.int64)
    letters = np.frombuffer(b.encode('ascii'), dtype=np.uint8)
    scores, matches = -indices, np.zeros(width, dtype=np.int64)
    lengths, paired = indices.copy(), np.zeros(width, dtype=np.int64)
    for i, aa in enumerate(a, 1):
        equal = letters == ord(aa)
        diagonal = scores[:-1] + np.where(equal, 1, -1)
        deletion = scores[1:]-1
        choose_diagonal = diagonal >= deletion
        # An insertion from column k to j costs j-k. The running maximum
        # solves that dependency; the most recent equal maximum preserves
        # diagonal > deletion > insertion tie-breaking at each cell.
        candidates = np.empty(width, dtype=np.int64)
        candidates[0] = -i
        candidates[1:] = np.maximum(diagonal, deletion)+indices[1:]
        envelope = np.maximum.accumulate(candidates)
        anchors = np.maximum.accumulate(np.where(candidates == envelope, indices, 0))
        base_matches = np.empty(width, dtype=np.int64)
        base_lengths = np.empty(width, dtype=np.int64)
        base_paired = np.empty(width, dtype=np.int64)
        base_matches[0], base_lengths[0], base_paired[0] = 0, i, 0
        base_matches[1:] = np.where(choose_diagonal, matches[:-1]+equal, matches[1:])
        base_lengths[1:] = np.where(choose_diagonal, lengths[:-1], lengths[1:])+1
        base_paired[1:] = np.where(choose_diagonal, paired[:-1]+1, paired[1:])
        scores = envelope-indices
        matches, paired = base_matches[anchors], base_paired[anchors]
        lengths = base_lengths[anchors]+indices-anchors
    return (100.0*int(matches[-1])/int(lengths[-1]),
            100.0*int(paired[-1])/len(a), 100.0*int(paired[-1])/len(b))


@lru_cache(maxsize=4096)
def protein_alignment(a, b, backend='auto'):
    if backend not in ('auto', 'numpy', 'python'):
        raise ValueError('Unknown alignment backend: ' + backend)
    if backend == 'numpy' and np is None:
        raise ValueError('numpy backend requires NumPy: python -m pip install numpy')
    if a == b:
        return 100.0, 100.0 if a else 0.0, 100.0 if b else 0.0
    if not a or not b:
        return 0.0, 0.0, 0.0
    # A subsequence alignment reaches the score upper bound with all residues
    # of the shorter protein matched and the minimum possible number of gaps.
    shorter, longer = (a, b) if len(a) <= len(b) else (b, a)
    iterator = iter(longer)
    if all(any(char == other for other in iterator) for char in shorter):
        count = len(shorter)
        return 100.0*count/len(longer), 100.0*count/len(a), 100.0*count/len(b)
    if np is not None and backend != 'python' and len(a)*len(b) >= 4096:
        return _protein_alignment_numpy(a, b)
    return _protein_alignment_python(a, b)


def protein_identity(a, b):
    return protein_alignment(a, b)[0]


def read_fasta(path):
    sequences = {}
    name = None
    with open(path, encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('>'):
                header = line[1:].split()
                if not header:
                    raise ValueError('Empty FASTA header')
                name = header[0]
                if name in sequences:
                    raise ValueError('Duplicate FASTA sequence: ' + name)
                sequences[name] = []
            elif line.strip():
                if name is None:
                    raise ValueError('FASTA sequence before header')
                sequence = line.strip().upper()
                if set(sequence) - set('ACGTRYMKSWBDHVN'):
                    raise ValueError('Invalid DNA characters in FASTA: ' + name)
                sequences[name].append(sequence)
    result = {k: ''.join(v) for k, v in sequences.items()}
    if not result or any(not sequence for sequence in result.values()):
        raise ValueError('FASTA contains no sequences or an empty sequence')
    return result


class Feature:
    def __init__(self, line):
        self.fields = line.rstrip('\n').split('\t')
        if len(self.fields) != 9:
            raise ValueError('Expected nine tab-separated GTF columns')
        self.attrs = dict(re.findall(r'(\S+)\s+"([^"]*)"\s*;?', self.fields[8]))
        if not self.attrs and self.fields[8] != '.':
            raise ValueError('This script accepts GTF attributes; convert GFF3 to GTF first')
        self.start, self.end = int(self.fields[3]), int(self.fields[4])
        if self.start < 1 or self.end < self.start:
            raise ValueError('Invalid feature coordinates')

    @property
    def kind(self):
        return self.fields[2]

    def line(self):
        fields = self.fields[:]
        fields[3:5] = [str(self.start), str(self.end)]
        fields[8] = ' '.join(f'{k} "{v}";' for k, v in self.attrs.items())
        return '\t'.join(fields) + '\n'


def ordered(features):
    return sorted(features, key=lambda f: f.start, reverse=bool(features and features[0].fields[6] == '-'))


def coding_dna(cds, genome):
    if not cds:
        return ''
    parts = []
    for feature in ordered(cds):
        sequence = genome[feature.fields[0]]
        if feature.end > len(sequence):
            raise ValueError('CDS exceeds FASTA sequence bounds')
        part = sequence[feature.start-1:feature.end]
        if feature.fields[6] == '-':
            part = part.translate(str.maketrans('ACGTRYMKBDHVN', 'TGCAYRKMVHDBN'))[::-1]
        parts.append(part)
    first_phase = ordered(cds)[0].fields[7]
    if first_phase not in ('0', '1', '2'):
        raise ValueError('Missing or invalid CDS phase')
    return ''.join(parts)[int(first_phase):]


class FastGTFModifier:
    EXON_ERRORS = ('exon_deletion', 'exon_truncation', 'exon_shift', 'add_exon', 'merge_exon', 'split_exon')
    ERROR_TYPES = EXON_ERRORS + ('gene_fusion', 'gene_split')

    def __init__(self, error_rate=0.1, genome=None, seed=None, error_types=None,
                 max_fusion_gap=10000, min_orf_aa=1, alignment_backend='auto', progress_every=100,
                 require_protein_change=True, max_attempts=20):
        if not 0 <= error_rate <= 1:
            raise ValueError('error_rate must be between 0 and 1')
        self.error_rate = error_rate
        if alignment_backend not in ('auto', 'numpy', 'python') or progress_every < 0:
            raise ValueError('Invalid alignment backend or progress interval')
        if alignment_backend == 'numpy' and np is None:
            raise ValueError('Install NumPy to use --alignment-backend numpy')
        self.alignment_backend, self.progress_every = alignment_backend, progress_every
        if max_attempts < 1:
            raise ValueError('max_attempts must be >= 1')
        self.require_protein_change, self.max_attempts = require_protein_change, max_attempts
        self.genome_path = Path(genome).resolve() if genome else None
        self.genome = read_fasta(genome) if genome else None
        self.rng = random.Random(seed)
        self.error_types = list(error_types or self.ERROR_TYPES)
        if max_fusion_gap < 2 or min_orf_aa < 1:
            raise ValueError('max_fusion_gap must be >= 2; min_orf_aa must be >= 1')
        self.max_fusion_gap, self.min_orf_aa = max_fusion_gap, min_orf_aa
        if set(self.error_types) - set(self.ERROR_TYPES):
            raise ValueError('Unsupported error type')

    def _mutate(self, exons, error):
        if error not in self.EXON_ERRORS:
            return None
        result = sorted(copy.deepcopy(exons), key=lambda f: f.start)
        if not result:
            return None
        target = self.rng.choice(result)
        if error == 'exon_deletion':
            if len(result) < 2:
                return None
            result.remove(target)
        elif error == 'exon_truncation':
            if target.end <= target.start:
                return None
            amount = self.rng.randint(1, max(1, (target.end-target.start)//2))
            if self.rng.choice([True, False]):
                target.start += amount
            else:
                target.end -= amount
        elif error == 'exon_shift':
            amount = self.rng.choice([-1, 1]) * self.rng.randint(1, 100)
            target.start += amount
            target.end += amount
        elif error == 'add_exon':
            # GTF coordinates are inclusive: reserve two intronic bases on
            # each side, plus at least three bases for the inserted exon.
            min_gap, min_length = 2, 3
            gaps = [(a.end+min_gap+1, b.start-min_gap-1)
                    for a, b in zip(result, result[1:])
                    if b.start-a.end-1 >= 2*min_gap+min_length]
            if not gaps:
                return None
            lo, hi = self.rng.choice(gaps)
            target = copy.deepcopy(target)
            length = self.rng.randint(min_length, min(90, hi-lo+1))
            target.start = self.rng.randint(lo, hi-length+1)
            target.end = target.start+length-1
            result.append(target)
        elif error == 'merge_exon':
            if len(result) < 2:
                return None
            i = self.rng.randrange(len(result)-1)
            result[i].end = result[i+1].end
            del result[i+1]
        elif error == 'split_exon':
            candidates = [f for f in result if f.end-f.start >= 5]
            if not candidates:
                return None
            target = self.rng.choice(candidates)
            right = copy.deepcopy(target)
            cut = self.rng.randint(target.start+1, target.end-2)
            right.start = cut+2  # Remove an actual base: artificial intron.
            target.end = cut
            result.append(right)
        result.sort(key=lambda f: f.start)
        if any(f.start < 1 or f.end < f.start for f in result):
            return None
        if any(a.end >= b.start for a, b in zip(result, result[1:])):
            return None
        if self.genome and any(f.end > len(self.genome[f.fields[0]]) for f in result):
            return None
        return result

    def _new_id(self, base):
        candidate, number = base, 1
        while candidate in self._used_ids:
            number += 1
            candidate = f'{base}_{number}'
        self._used_ids.add(candidate)
        if hasattr(self, '_pending_ids'):
            self._pending_ids.append(candidate)
        return candidate

    def _discard_pending_ids(self):
        self._used_ids.difference_update(self._pending_ids)
        self._pending_ids = []

    def _protein_changed(self, output, sources, groups):
        """Cheap strict acceptance check: no alignment needed for retries."""
        new = translate(coding_dna([f for f in output if f.kind == 'CDS'], self.genome))
        if 'X' in new:
            return False
        for source in sources:
            cds = [f for f in groups[source] if f.kind == 'CDS']
            if not cds:
                return False
            old = translate(coding_dna(cds, self.genome))
            if not old or 'X' in old or old == new:
                return False
        return True

    def _has_evaluable_protein(self, features):
        cds = [f for f in features if f.kind == 'CDS']
        if not cds:
            return False
        peptide = translate(coding_dna(cds, self.genome))
        return bool(peptide) and 'X' not in peptide

    def _fusion_fragment(self, features, keep_prefix):
        """Keep a transcript prefix/suffix cut INSIDE its annotated CDS.

        This creates a chimeric exon chain, rather than appending a second
        complete ORF behind the first gene's original stop codon.
        """
        cds = ordered([f for f in features if f.kind == 'CDS'])
        total = sum(f.end-f.start+1 for f in cds)
        if total < 3:
            return None
        position = self.rng.randint(1, total-2)
        for child in cds:
            length = child.end-child.start+1
            if position < length:
                cut = child.start+position if child.fields[6] == '+' else child.end-position
                break
            position -= length
        keep_low = keep_prefix == (cds[0].fields[6] == '+')
        result = []
        for exon in [f for f in features if f.kind == 'exon']:
            child = copy.deepcopy(exon)
            if keep_low:
                child.end = min(child.end, cut)
            else:
                child.start = max(child.start, cut)
            if child.start <= child.end:
                result.append(child)
        return result or None

    def _progress(self, message):
        if self.progress_every:
            print(message, file=sys.stderr, flush=True)

    def _predict(self, exons, gid, tid, source_ids, event_id):
        """Longest ATG ORF; first in transcript order wins length ties.

        CDS excludes the stop codon, per GTF convention. An open-ended ORF
        ends at the final complete codon; no three_prime_utr is inferred without
        a stop codon.
        """
        exons = ordered(copy.deepcopy(exons))
        pieces, mapping, offset = [], [], 0
        for exon in exons:
            sequence = self.genome[exon.fields[0]][exon.start-1:exon.end]
            if exon.fields[6] == '-':
                sequence = sequence.translate(str.maketrans('ACGTRYMKBDHVN', 'TGCAYRKMVHDBN'))[::-1]
            pieces.append(sequence)
            mapping.append((offset, offset+len(sequence), exon))
            offset += len(sequence)
        sequence = ''.join(pieces)
        best = None
        # Each frame is scanned once: starts within an existing ORF cannot
        # produce a longer ORF ending at the same stop.
        for frame in range(3):
            start = None
            for pos in range(frame, len(sequence)-2, 3):
                codon = sequence[pos:pos+3]
                if start is None and codon == 'ATG':
                    start = pos
                if codon in ('TAA', 'TAG', 'TGA') and start is not None:
                    candidate = (start, pos, True)
                    if best is None or (pos-start, -start) > (best[1]-best[0], -best[0]):
                        best = candidate
                    start = None
            if start is not None:
                end = len(sequence)-(len(sequence)-frame) % 3
                candidate = (start, end, False)
                if best is None or (end-start, -start) > (best[1]-best[0], -best[0]):
                    best = candidate
        if best is not None and (best[1]-best[0])//3 < self.min_orf_aa:
            best = None
        attrs = dict(gene_id=gid, transcript_id=tid,
                     source_gene_ids=','.join(dict.fromkeys(g for g, _ in source_ids)),
                     source_transcript_ids=','.join(t for _, t in source_ids), event_id=event_id,
                     cds_prediction='longest_ATG_ORF',
                     orf_status=('complete' if best[2] else 'partial_3prime') if best else 'no_orf')

        def make(kind, start, end, template):
            child = copy.deepcopy(template)
            child.fields[1:3] = ['annotation_error_simulator', kind]
            child.fields[5], child.fields[7] = '.', '.'
            child.start, child.end = start, end
            child.attrs = attrs.copy()  # Drop stale protein/exon/Parent/Name identifiers.
            return child

        tx = make('transcript', min(e.start for e in exons), max(e.end for e in exons), exons[0])
        result = [tx]
        for number, exon in enumerate(exons, 1):
            child = make('exon', exon.start, exon.end, exon)
            child.attrs.update(exon_number=str(number), exon_id=self._new_id(f'{tid}.exon{number}'))
            result.append(child)

        def emit(kind, lo, hi):
            features = []
            for start, end, exon in mapping:
                left, right = max(lo, start), min(hi, end)
                if left >= right:
                    continue
                if exon.fields[6] == '-':
                    a, b = exon.end-(right-start)+1, exon.end-(left-start)
                else:
                    a, b = exon.start+left-start, exon.start+right-start-1
                features.append(make(kind, a, b, exon))
            if kind in ('start_codon', 'stop_codon'):
                consumed = 0
                for child in features:
                    child.fields[7] = str((-consumed) % 3)
                    consumed += child.end-child.start+1
            return features

        if best:
            start, end, complete = best
            cds = emit('CDS', start, end)
            consumed = 0
            protein_id = self._new_id(f'{tid}.protein')
            for child in cds:
                child.fields[7] = str((-consumed) % 3)
                child.attrs['protein_id'] = protein_id
                consumed += child.end-child.start+1
            result.extend(cds)
            result.extend(emit('start_codon', start, start+3))
            if complete:
                result.extend(emit('stop_codon', end, end+3))
            result.extend(emit('five_prime_utr', 0, start))
            if complete:
                result.extend(emit('three_prime_utr', end+3, len(sequence)))
        return result

    def _plan_events(self, groups):
        """One selection per source gene; a consumed gene cannot be edited twice."""
        genes = defaultdict(list)
        for tid, features in groups.items():
            keys = {(f.fields[0], f.fields[6], f.attrs.get('gene_id')) for f in features}
            if len(keys) != 1 or not features[0].attrs.get('gene_id'):
                raise ValueError('Inconsistent or missing gene_id in transcript: ' + tid)
            exons = sorted([f for f in features if f.kind == 'exon'], key=lambda f: f.start)
            if any(b.start <= a.end for a, b in zip(exons, exons[1:])):
                raise ValueError('Overlapping exons in transcript: ' + tid)
            genes[features[0].attrs['gene_id']].append(tid)
        bounds = {}
        for gid, tids in genes.items():
            children = [f for tid in tids for f in groups[tid]]
            if len({(f.fields[0], f.fields[6]) for f in children}) != 1:
                raise ValueError('Gene spans multiple sequences/strands: ' + gid)
            bounds[gid] = (children[0].fields[0], min(f.start for f in children),
                           max(f.end for f in children), children[0].fields[6])
        gene_order = sorted(genes, key=lambda g: bounds[g][:3])
        handled, plans = set(), []
        event_number = 0
        self._pending_ids = []

        def unchanged(gid, status):
            for tid in genes[gid]:
                plans.append((copy.deepcopy(groups[tid]), [tid], 'none', status, 'NA'))

        for index, gid in enumerate(gene_order):
            if self.progress_every and index % self.progress_every == 0:
                self._progress(f'[structure] genes {index}/{len(gene_order)}; products {len(plans)}')
            if gid in handled:
                continue
            if self.rng.random() >= self.error_rate:
                unchanged(gid, 'not_selected')
                handled.add(gid)
                continue
            if self.require_protein_change and not any(self._has_evaluable_protein(groups[t]) for t in genes[gid]):
                unchanged(gid, 'no_evaluable_source_protein')
                handled.add(gid)
                continue
            choices = self.error_types[:]
            self.rng.shuffle(choices)
            choices *= self.max_attempts if self.require_protein_change else 1
            selected = False
            attempted_products = False
            for error in choices:
                self._discard_pending_ids()
                # Each product is (exons, source transcript IDs, new gene ID).
                products, partner = [], None
                if error == 'gene_fusion':
                    if index+1 >= len(gene_order):
                        continue
                    other = gene_order[index+1]
                    a, b = bounds[gid], bounds[other]
                    gap = b[1]-a[2]-1
                    if (other in handled or a[0] != b[0] or a[3] != b[3] or
                            a[3] not in ('+', '-') or not 2 <= gap <= self.max_fusion_gap):
                        continue
                    if any(not any(f.kind == 'exon' for f in groups[t]) for t in genes[gid]+genes[other]):
                        continue
                    if self.require_protein_change and not all(
                            self._has_evaluable_protein(groups[t]) for t in genes[gid]+genes[other]):
                        continue
                    partner = other
                    new_gid = self._new_id(f'{gid}__{other}.fusion')
                    for left in genes[gid]:
                        for right in genes[other]:
                            sources = sorted([left, right], key=lambda t: min(f.start for f in groups[t]),
                                             reverse=a[3] == '-')
                            if self.require_protein_change:
                                prefix = self._fusion_fragment(groups[sources[0]], True)
                                suffix = self._fusion_fragment(groups[sources[1]], False)
                                if prefix is None or suffix is None:
                                    products = []
                                    break
                                exons = prefix+suffix
                            else:
                                exons = [f for t in sources for f in groups[t] if f.kind == 'exon']
                            products.append((exons, sources, new_gid))
                        if not products:
                            break
                elif error == 'gene_split':
                    # Split only in a common intronic gap across ALL isoforms.
                    if any(not any(f.kind == 'exon' for f in groups[t]) for t in genes[gid]):
                        continue
                    if self.require_protein_change and not all(self._has_evaluable_protein(groups[t]) for t in genes[gid]):
                        continue
                    all_exons = sorted([f for t in genes[gid] for f in groups[t] if f.kind == 'exon'],
                                       key=lambda f: f.start)
                    merged = []
                    for exon in all_exons:
                        if merged and exon.start <= merged[-1][1]+1:
                            merged[-1][1] = max(merged[-1][1], exon.end)
                        else:
                            merged.append([exon.start, exon.end])
                    gaps = [(a[1], b[0]) for a, b in zip(merged, merged[1:]) if b[0]-a[1]-1 >= 2]
                    if not gaps or bounds[gid][3] not in ('+', '-'):
                        continue
                    left_end, right_start = self.rng.choice(gaps)
                    cut = (left_end+right_start)//2
                    gids = [self._new_id(f'{gid}.split1'), self._new_id(f'{gid}.split2')]
                    for t in genes[gid]:
                        if self.require_protein_change and not self._has_evaluable_protein(groups[t]):
                            continue
                        exons = [f for f in groups[t] if f.kind == 'exon']
                        partitions = [[f for f in exons if f.end <= cut],
                                      [f for f in exons if f.start > cut]]
                        if bounds[gid][3] == '-':
                            partitions.reverse()
                        for part, new_gid in zip(partitions, gids):
                            if part:
                                products.append((part, [t], new_gid))
                else:
                    for t in genes[gid]:
                        if self.require_protein_change and not self._has_evaluable_protein(groups[t]):
                            continue
                        exons = [f for f in groups[t] if f.kind == 'exon']
                        mutated = self._mutate(exons, error)
                        if mutated is not None and exons[0].fields[6] in ('+', '-'):
                            products.append((mutated, [t], gid))
                if not products:
                    continue
                attempted_products = True
                event_id = f'event_{event_number+1:06d}'
                product_sources = set()
                candidate_plans = []
                for number, (exons, sources, new_gid) in enumerate(products, 1):
                    new_tid = self._new_id(f'{new_gid}.{error}.t{number}')
                    source_ids = [(groups[t][0].attrs['gene_id'], t) for t in sources]
                    output = self._predict(exons, new_gid, new_tid, source_ids, event_id)
                    if self.require_protein_change and not self._protein_changed(output, sources, groups):
                        if error in ('gene_fusion', 'gene_split'):
                            # Structural events must pass atomically for every
                            # product/source pair; do not silently drop isoforms.
                            candidate_plans = []
                            break
                        continue
                    candidate_plans.append((output, sources, error, 'applied', event_id))
                    product_sources.update(sources)
                if not candidate_plans:
                    continue
                event_number += 1
                plans.extend(candidate_plans)
                self._pending_ids = []  # Commit accepted IDs.
                if error in self.EXON_ERRORS:
                    for t in genes[gid]:
                        if t not in product_sources:
                            status = 'no_protein_changing_error' if self.require_protein_change else 'no_applicable_error'
                            plans.append((copy.deepcopy(groups[t]), [t], 'none', status, 'NA'))
                handled.add(gid)
                if partner:
                    handled.add(partner)
                selected = True
                break
            if not selected:
                self._discard_pending_ids()
                status = 'no_protein_changing_error' if attempted_products and self.require_protein_change else 'no_applicable_error'
                unchanged(gid, status)
                handled.add(gid)
        return plans

    def _label(self, before, after):
        old = [f for f in before if f.kind == 'CDS']
        new = [f for f in after if f.kind == 'CDS']
        result = dict(coding_impact='unknown', protein_identity_pct='NA', protein_identity_loss_pct='NA',
                      original_protein_length='NA', modified_protein_length='NA', cds_length_delta='NA',
                      frameshift='NA', comparison_status='missing_genome',
                      original_aligned_pct='NA', modified_aligned_pct='NA')
        if not old:
            result.update(coding_impact='protein_gained' if new else 'noncoding',
                          comparison_status='no_original_CDS')
            if self.genome is not None:
                result['modified_protein_length'] = len(translate(coding_dna(new, self.genome)))
            return result
        result['cds_length_delta'] = sum(f.end-f.start+1 for f in new)-sum(f.end-f.start+1 for f in old)
        result['frameshift'] = str(result['cds_length_delta'] % 3 != 0).lower()
        if any(f.attrs.get('cds_prediction') == 'longest_ATG_ORF' for f in after):
            # A newly selected ORF has no defined inherited reading-frame anchor.
            result['frameshift'] = 'NA'
        if self.genome is None:
            return result
        dna_old, dna_new = coding_dna(old, self.genome), coding_dna(new, self.genome)
        a, b = translate(dna_old), translate(dna_new)
        result.update(original_protein_length=len(a), modified_protein_length=len(b))
        if 'X' in a+b or not a:
            result['comparison_status'] = 'ambiguous_sequence_or_empty_original_protein'
            return result
        identity, original_aligned, modified_aligned = protein_alignment(a, b, self.alignment_backend)
        result.update(coding_impact='protein_unchanged' if a == b else 'protein_changed',
                      protein_identity_pct=f'{identity:.6f}', protein_identity_loss_pct=f'{100-identity:.6f}',
                      original_protein_length=len(a), modified_protein_length=len(b), comparison_status='ok',
                      original_aligned_pct=f'{original_aligned:.6f}', modified_aligned_pct=f'{modified_aligned:.6f}')
        return result

    def modify_gtf(self, input_gtf, output_gtf):
        started = time.perf_counter()
        if self.genome is None:
            raise ValueError('Reference genome FASTA is required; supply genome= or -g')
        outputs = [Path(output_gtf), Path(str(output_gtf)+'.error_report.tsv'),
                   Path(str(output_gtf)+'.original_proteins.fa'),
                   Path(str(output_gtf)+'.modified_proteins.fa')]
        if Path(input_gtf).resolve() in [p.resolve() for p in outputs]:
            raise ValueError('Input and output paths must differ')
        if self.genome_path in [p.resolve() for p in outputs]:
            raise ValueError('Output paths must not overwrite the reference FASTA')
        groups = defaultdict(list)
        backend = 'numpy' if np is not None and self.alignment_backend != 'python' else 'python'
        self._progress(f'[input] reading GTF; alignment backend={backend}')
        if backend == 'python':
            self._progress('[hint] NumPy is unavailable or disabled; long protein comparisons may be slow')
        standalone, comments = [], []
        with open(input_gtf, encoding='utf-8') as handle:
            for line in handle:
                if not line.strip() or line.startswith('#'):
                    comments.append(line)
                    continue
                feature = Feature(line)
                seqid = feature.fields[0]
                if seqid not in self.genome:
                    raise ValueError(f'GTF sequence {seqid!r} is missing from FASTA')
                if feature.end > len(self.genome[seqid]):
                    raise ValueError(f'GTF coordinates exceed FASTA bounds: {seqid}:{feature.start}-{feature.end}')
                if feature.kind == 'CDS':
                    if feature.fields[7] not in ('0', '1', '2'):
                        raise ValueError(f'Invalid CDS phase at {seqid}:{feature.start}')
                    if feature.fields[6] not in ('+', '-'):
                        raise ValueError(f'CDS requires + or - strand at {seqid}:{feature.start}')
                tid = feature.attrs.get('transcript_id')
                if tid:
                    groups[tid].append(feature)
                else:
                    standalone.append(feature)
        self._used_ids = {value for features in groups.values() for f in features
                          for key, value in f.attrs.items()
                          if key in ('gene_id', 'transcript_id', 'exon_id', 'protein_id')}
        self._used_ids.update(f.attrs['gene_id'] for f in standalone if 'gene_id' in f.attrs)
        results, report = [], []
        modified_proteins = {}
        self._progress(f'[structure] input transcripts={len(groups)}')
        self._progress(f'[structure] require_protein_change={self.require_protein_change}; max_rounds={self.max_attempts}')
        plans = self._plan_events(groups)
        self._progress(f'[structure] applied_products={sum(p[3] == "applied" for p in plans)}; '
                       f'preserved_transcripts={sum(p[3] != "applied" for p in plans)}')
        self._progress(f'[protein] comparing {len(plans)} output transcripts')
        last_progress = time.perf_counter()
        for index, (output, sources, error, status, event_id) in enumerate(plans):
            now = time.perf_counter()
            if self.progress_every and (index % self.progress_every == 0 or now-last_progress >= 30):
                self._progress(f'[protein] {index}/{len(plans)}; elapsed={now-started:.1f}s')
                last_progress = now
            gid, tid = output[0].attrs['gene_id'], output[0].attrs['transcript_id']
            comparisons = []
            for source in sources:
                labels = self._label(groups[source], output)
                comparisons.append(labels)
                report.append(dict(gene_id=gid, transcript_id=tid,
                                   source_gene_id=groups[source][0].attrs['gene_id'],
                                   source_transcript_id=source, event_id=event_id,
                                   error_type=error, error_status=status,
                                   cds_prediction=output[0].attrs.get('cds_prediction', 'original_annotation'),
                                   orf_status=output[0].attrs.get('orf_status', 'original_annotation'), **labels))
            # Fusion has two baselines: keep each numerical comparison in TSV,
            # rather than assigning one source's identity to the fused transcript.
            summary = comparisons[0].copy()
            if len(comparisons) > 1:
                for key in summary:
                    summary[key] = 'NA'
                impacts = {c['coding_impact'] for c in comparisons}
                summary['coding_impact'] = ('protein_changed' if 'protein_changed' in impacts else
                                            'protein_gained' if 'protein_gained' in impacts else
                                            'unknown' if 'unknown' in impacts else
                                            'protein_unchanged' if 'protein_unchanged' in impacts else 'noncoding')
                summary['comparison_status'] = 'multiple_sources_see_report'
            for f in output:
                f.attrs.update(error_type=error, error_status=status,
                               **{key: str(value) for key, value in summary.items()})
            cds = [f for f in output if f.kind == 'CDS']
            if cds or status == 'applied' or any(f.kind == 'CDS' for t in sources for f in groups[t]):
                modified_proteins[tid] = translate(coding_dna(cds, self.genome))
            results.extend(output)
        # Rebuild exactly one gene row for each surviving/new gene. Remove all
        # old gene-level rows of consumed genes, including stale gene metadata.
        spans = defaultdict(list)
        for f in results:
            spans[f.attrs['gene_id']].append(f)
        original_genes = {f.attrs.get('gene_id'): f for f in standalone if f.kind == 'gene'}
        source_genes = {features[0].attrs['gene_id'] for features in groups.values()}
        preserved = [f for f in standalone if f.kind != 'gene' and
                     (not f.attrs.get('gene_id') or f.attrs.get('gene_id') in spans)]
        preserved.extend(f for gid, f in original_genes.items() if gid not in source_genes)
        for gid, children in spans.items():
            if gid in original_genes:
                gene = copy.deepcopy(original_genes[gid])
            else:
                gene = copy.deepcopy(children[0])
                gene.fields[2], gene.fields[7] = 'gene', '.'
                gene.attrs = {'gene_id': gid}
            gene.start, gene.end = min(f.start for f in children), max(f.end for f in children)
            for key in ('transcript_id', 'protein_id', 'exon_id', 'exon_number'):
                gene.attrs.pop(key, None)
            for key in ('source_gene_ids', 'event_id', 'error_type', 'error_status'):
                values = list(dict.fromkeys(f.attrs[key] for f in children if key in f.attrs))
                if values:
                    gene.attrs[key] = ','.join(values)
            preserved.append(gene)
        self._progress('[output] writing GTF, report and protein FASTA')
        with open(output_gtf, 'w', encoding='utf-8') as handle:
            handle.writelines(comments)
            rank = {'gene': 0, 'transcript': 1, 'exon': 2}
            for f in sorted(preserved+results, key=lambda f: (f.fields[0], f.start,
                            rank.get(f.kind, 3), f.end, f.attrs.get('transcript_id', ''))):
                handle.write(f.line())
        with open(str(output_gtf)+'.error_report.tsv', 'w', encoding='utf-8', newline='') as handle:
            columns = ['gene_id', 'transcript_id', 'source_gene_id', 'source_transcript_id',
                       'event_id', 'error_type', 'error_status', 'cds_prediction', 'orf_status', 'coding_impact',
                       'protein_identity_pct', 'protein_identity_loss_pct', 'original_aligned_pct',
                       'modified_aligned_pct', 'original_protein_length', 'modified_protein_length',
                       'cds_length_delta', 'frameshift', 'comparison_status']
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t')
            writer.writeheader()
            writer.writerows(report)
        original_proteins = {tid: translate(coding_dna([f for f in features if f.kind == 'CDS'], self.genome))
                             for tid, features in groups.items() if any(f.kind == 'CDS' for f in features)}
        for proteins, suffix in ((original_proteins, '.original_proteins.fa'),
                                 (modified_proteins, '.modified_proteins.fa')):
            with open(str(output_gtf)+suffix, 'w', encoding='utf-8') as handle:
                for tid, peptide in proteins.items():
                    handle.write('>' + tid + '\n')
                    for start in range(0, len(peptide), 60):
                        handle.write(peptide[start:start+60] + '\n')
        self._progress(f'[done] elapsed={time.perf_counter()-started:.1f}s; output={output_gtf}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input', required=True, help='Input GTF (not GFF3)')
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-e', '--error_rate', type=float, default=0.1)
    parser.add_argument('-g', '--genome', '--fasta', required=True, help='Reference genome DNA FASTA')
    parser.add_argument('--seed', type=int)
    parser.add_argument('--max-fusion-gap', type=int, default=10000, help='Maximum intergenic gap in bp')
    parser.add_argument('--min-orf-aa', type=int, default=1, help='Minimum predicted ORF length in amino acids')
    parser.add_argument('--alignment-backend', choices=('auto', 'numpy', 'python'), default='auto',
                        help='Exact alignment backend; auto uses NumPy when available')
    parser.add_argument('--progress-every', type=int, default=100, help='Log every N genes/products; 0 disables logs')
    parser.add_argument('--max-attempts', type=int, default=20, help='Maximum candidate rounds per selected gene')
    parser.add_argument('--allow-protein-unchanged', action='store_true',
                        help='Disable the default requirement that every applied error changes its source proteins')
    parser.add_argument('--error-types', nargs='+', choices=FastGTFModifier.ERROR_TYPES)
    args = parser.parse_args()
    FastGTFModifier(args.error_rate, args.genome, args.seed, args.error_types,
                    args.max_fusion_gap, args.min_orf_aa, args.alignment_backend,
                    args.progress_every, not args.allow_protein_unchanged,
                    args.max_attempts).modify_gtf(args.input, args.output)


if __name__ == '__main__':
    main()
