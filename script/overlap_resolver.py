from dataclasses import dataclass, field
from typing import List, Dict, Optional

@dataclass
class Region:
    Chr: str
    Start: int
    End: int
    Strand: str
    cov: float
    Length: int
    flag: int = field(default=0)  # 1 if affected per rules (see description)

    def copy(self) -> "Region":
        return Region(
            self.Chr, self.Start, self.End, self.Strand,
            self.cov, self.Length, self.flag
        )

    def update_length(self):
        self.Length = max(0, self.End - self.Start + 1)

def overlap_len(a: Region, b: Region) -> int:
    """
    Calculate the overlap length between two regions.
    """
    left = max(a.Start, b.Start)
    right = min(a.End, b.End)
    return max(0, right - left + 1)

def resolve_pair(a: Region, b: Region, gap: int) -> (Optional[Region], Optional[Region], bool):
    """
    Resolve a pair of overlapping regions on the same chromosome.

    Rules:
    1. If overlap > 50% of the smaller region -> drop smaller, set flag=1 on longer (kept).
    2. Else trim the smaller region to remove overlap and ensure 'gap' bp separation; 
       set flag=1 on the trimmed one. If trimming collapses the region, drop it.
    
    Also handles cases where regions are too close (<gap) but not overlapping.
    
    Returns:
        Tuple of (modified_a, modified_b, changed_flag)
    """
    ovl = overlap_len(a, b)

    # If no overlap, enforce minimum gap if too close
    if ovl == 0:
        changed = False
        if a.End + gap >= b.Start:
            smaller, larger = (a, b) if a.Length <= b.Length else (b, a)
            if smaller is b:
                new_start = larger.End + gap + 1
                if new_start <= b.End:
                    b.Start = new_start
                    b.update_length()
                    b.flag = 1
                    changed = True
                else:
                    b = None
                    changed = True
            else:
                new_end = larger.Start - gap - 1
                if new_end >= a.Start:
                    a.End = new_end
                    a.update_length()
                    a.flag = 1
                    changed = True
                else:
                    a = None
                    changed = True
        return a, b, changed

    # There is overlap: determine smaller vs larger by current length
    smaller, larger = (a, b) if a.Length <= b.Length else (b, a)
    if ovl > 0.5 * smaller.Length:
        # Rule 1: drop smaller; mark larger
        larger.flag = 1
        if smaller is a:
            return None, larger, True
        else:
            return larger, None, True

    # Rule 2: trim the smaller to remove overlap and ensure gap; mark smaller
    changed = False
    if smaller.Start <= larger.Start:
        # overlap on smaller's right side -> trim End
        new_end = larger.Start - gap - 1
        if new_end < smaller.Start:
            # collapsed -> drop smaller
            if smaller is a:
                return None, larger, True
            else:
                return larger, None, True
        smaller.End = new_end
        smaller.update_length()
        smaller.flag = 1
        changed = True
    else:
        # overlap on smaller's left side -> trim Start
        new_start = larger.End + gap + 1
        if new_start > smaller.End:
            if smaller is a:
                return None, larger, True
            else:
                return larger, None, True
        smaller.Start = new_start
        smaller.update_length()
        smaller.flag = 1
        changed = True

    # Return in original identity order
    if smaller is a:
        return a, larger, changed
    else:
        return larger, b, changed

def resolve_overlaps_one_chr(regs: List[Region], gap: int, max_passes: int = 1000) -> List[Region]:
    """
    Process regions on one chromosome to resolve overlaps.
    
    Args:
        regs: List of regions (must be on same chromosome)
        gap: Minimum gap to enforce between regions
        max_passes: Maximum iterations to attempt
        
    Returns:
        List of processed regions with overlaps resolved
    """
    regs = sorted(regs, key=lambda r: (r.Start, r.End))
    passes = 0
    changed_any = True
    while changed_any and passes < max_passes:
        passes += 1
        changed_any = False
        i = 0
        while i < len(regs) - 1:
            a, b = regs[i], regs[i+1]
            a2, b2, changed = resolve_pair(a.copy(), b.copy(), gap)
            if not changed:
                i += 1
                continue
            changed_any = True
            new_items = []
            if a2 is not None: new_items.append(a2)
            if b2 is not None: new_items.append(b2)
            regs = regs[:i] + new_items + regs[i+2:]
            regs.sort(key=lambda r: (r.Start, r.End))
            i = max(0, i - 1)
        regs.sort(key=lambda r: (r.Start, r.End))
    return regs

def process_regions(all_regions: List[Region], gap: int = 50) -> List[Region]:
    """
    Main function to process regions across all chromosomes.
    
    Args:
        all_regions: List of all regions to process
        gap: Minimum gap to enforce between regions
        
    Returns:
        List of processed regions with overlaps resolved
    """
    # Group by chromosome (strand-agnostic). For strand-aware processing, key by (Chr, Strand).
    by_chr: Dict[str, List[Region]] = {}
    for r in all_regions:
        by_chr.setdefault(r.Chr, []).append(r)

    result: List[Region] = []
    for chrom, lst in by_chr.items():
        cleaned = resolve_overlaps_one_chr(lst, gap=gap)
        result.extend(cleaned)

    return sorted(result, key=lambda r: (r.Chr, r.Start, r.End))
