from dataclasses import replace, dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from utils.config import (
    MIN_REGION_SEPARATION_TO_KEEP_SEPARATE,
    DISORDERED_DSSP_CODES,
    GROUP_DISORDERED,
)
from utils.classes import ComparedSegment, ModelComparison, Residue, OpenSegment, ExtractedModel
from utils.dssp import (
    classify_foldswitch,
    get_sequence,
    run_dssp,
    normalize_dssp_code,
)
from utils.utils import compact_reason
from utils.outputs import log_input_too_disordered
from pathlib import Path

TOO_DISORDERED_THRESHOLD_PERCENT = 55.0

def input_is_too_disordered(
    jobname: str,
    summarydir: Path,
    dominant_dssp_code: Sequence,
) -> bool:
    from utils.dssp import percent_unassigned
    pct = percent_unassigned(dominant_dssp_code) # checks only for - or space
    if pct > TOO_DISORDERED_THRESHOLD_PERCENT:
        log_input_too_disordered(jobname, summarydir, pct)
        return True
    return False

def is_good_physical_prediction(
    candidate_model: ExtractedModel,
    dominant_residues: Sequence,
    verbose: bool = False,
) -> bool:
    """
    Returns True if the candidate model passes all physical prediction filters.
    Filters are applied in order and are cumulative — all must pass.

    Step 1a: Phenix filter for good residues. Predictive must be >= 15% OR Unpacked high pLDDT >= 30%
    Step 1b: Unphysical must be <= 5%.
    Step 2: Unassigned DSSP filter — alternative unassigned percent must not exceed
            dominant by more than UNASSIGNED_MAX_ALTERNATIVE_PERCENT.
    """
    from utils.phenix_filter import get_mode_fractions
    from utils.dssp import run_dssp, percent_unassigned
    from utils.config import UNASSIGNED_MAX_ALTERNATIVE_PERCENT

    # Step 1: Phenix filter
    fractions = get_mode_fractions(
        candidate_model.pdb_path, 
        candidate_model.dssp_path.with_suffix(".phenix"),
    )

    if fractions is None:
        return True #Phenix failed but keep going

    if fractions is not None:
        predictive_pct = fractions.get("Predictive") * 100.0
        unpacked_highPLDDt_pct = fractions.get("Unpacked high pLDDT") * 100.0
        unphysical_pct = fractions.get("Unphysical") * 100.0

        if predictive_pct < 15.0 and unpacked_highPLDDt_pct < 30.0:
            if verbose:
                print(
                    f"  {candidate_model.model_name}: rejected, not enough good residues"
                    f"(Predictive={predictive_pct:.1f}% < 15% and Unpacked high pLDDT={unpacked_highPLDDt_pct:.1f}% < 30%).\n"
                )
            return False

        if unphysical_pct > 5.0:
            if verbose:
                print(
                    f"  {candidate_model.model_name}: rejected, too many bad residues"
                    f"(Unphysical={unphysical_pct:.1f}% > 5%)."
                )
            return False

    # Step 2: Unassigned DSSP filter (only reached if phenix passed)
    try:
        alt_residues = run_dssp(candidate_model)
        dom_pct = percent_unassigned(dominant_residues)
        alt_pct = percent_unassigned(alt_residues)
        if alt_pct > dom_pct + UNASSIGNED_MAX_ALTERNATIVE_PERCENT:
            if verbose:
                print(
                    f"  {candidate_model.model_name}: rejected "
                    f"(alt_unassigned={alt_pct:.1f}% > dom={dom_pct:.1f}% "
                    f"+ {UNASSIGNED_MAX_ALTERNATIVE_PERCENT:.1f}%)."
                )
            return False
    except Exception:
        if verbose:
            print(f"  {candidate_model.model_name}: DSSP failed, skipping unassigned filter.")

    return True

    
# Helper of helper, processes raw DSSP into list with structures (Beta, Helix, Kappa, Disordered) and bridges short disordered gaps
def _build_strand_groups(
    residues: Sequence[Residue],
    min_separation_to_keep_separate: int,
) -> List[str]:
    """
    Returns a flat per-position group assignment list (HELIX, BETA, KAPPA, DISORDERED)
    for a single strand, with short all-disordered gaps bridged where the same exact
    DSSP code appears on both sides.
    """
    # Pass 1: assign each position its group directly
    groups = [r.group for r in residues]

    # Pass 2: bridge short disordered gaps between same-code ordered regions
    n = len(groups)
    i = 0
    while i < n:
        if groups[i] == GROUP_DISORDERED:
            i += 1
            continue

        # Find the end of this ordered run
        run_end = i
        while run_end + 1 < n and groups[run_end + 1] == groups[i]:
            run_end += 1

        # Find the start of the next ordered run of the same group
        gap_start = run_end + 1
        gap_end = gap_start
        while gap_end < n and groups[gap_end] == GROUP_DISORDERED:
            gap_end += 1

        if gap_end >= n:
            break

        gap_length = gap_end - gap_start
        next_run_group = groups[gap_end]
        same_group = next_run_group == groups[i]
        same_code = (
            normalize_dssp_code(residues[i].dssp_code)
            == normalize_dssp_code(residues[gap_end].dssp_code)
        )
        short_enough = gap_length < min_separation_to_keep_separate

        if same_group and same_code and short_enough:
            for j in range(gap_start, gap_end):
                groups[j] = groups[i]

        i = run_end + 1

    return groups

# Helper, main logic for assigning FS classification
def compare_dssp_codes(
    dominant: Sequence[Residue],
    alternative: Sequence[Residue],
    min_run_length: int,
    min_separation_to_keep_separate: int,
) -> Tuple[List[ComparedSegment], List[ComparedSegment], float]:
    """
    Independently builds group segments for each DSSP sequence, then compares
    position by position to find fold-switch and disorder-change regions.
    Also computes same_ss_preserved_percent from the merged group assignments.

    Only called after RMSD gate and sequence identity are already confirmed
    in compare_candidate_models.

    Returns (switch_segments, disordered_change_segments, same_ss_preserved_percent).
    """
    dom_groups = _build_strand_groups(dominant, min_separation_to_keep_separate)
    alt_groups = _build_strand_groups(alternative, min_separation_to_keep_separate)

    # Walk both per-position group arrays to find contiguous change segments
    raw_segments: List[ComparedSegment] = []
    open_segment: Optional[OpenSegment] = None

    for idx, (dom_group, alt_group) in enumerate(zip(dom_groups, alt_groups), start=1):
        category = classify_foldswitch(dom_group, alt_group)

        is_same_run = (
            open_segment is not None
            and open_segment.dominant_group == dom_group
            and open_segment.alternative_group == alt_group
            and open_segment.category == category
        )

        if open_segment is not None and not is_same_run:
            raw_segments.append(ComparedSegment(
                start=open_segment.start,
                end=idx - 1,
                dominant_group=open_segment.dominant_group,
                alternative_group=open_segment.alternative_group,
                category=open_segment.category,
            ))
            open_segment = None

        if category is not None and open_segment is None:
            open_segment = OpenSegment(
                start=idx,
                dominant_group=dom_group,
                alternative_group=alt_group,
                category=category,
            )

    #close last run
    if open_segment is not None:
        raw_segments.append(ComparedSegment(
            start=open_segment.start,
            end=len(dominant),
            dominant_group=open_segment.dominant_group,
            alternative_group=open_segment.alternative_group,
            category=open_segment.category,
        ))

    qualifying = [s for s in raw_segments if s.end - s.start + 1 >= min_run_length]
    switch_segments = [s for s in qualifying if s.category == "switch"]
    disordered_change_segments = [s for s in qualifying if s.category == "disordered"]

    same_ss = sum(
        1 for d, a in zip(dom_groups, alt_groups)
        if d == a and d != GROUP_DISORDERED
    )
    total_ordered = sum(1 for d in dom_groups if d != GROUP_DISORDERED)
    same_ss_preserved_percent = 100.0 * same_ss / total_ordered if total_ordered else 0.0
    return switch_segments, disordered_change_segments, same_ss_preserved_percent


# Main function called at pdb/pse file level called by pdb_helpers.py helper file
def compare_candidate_models(
    dominant_model,
    alternative_models: Sequence,
    rmsd_results: Dict[str, Tuple[Optional[float], str]],
    min_run_length: int,
    rmsd_threshold: float,
) -> Tuple[List[ModelComparison], Dict[str, List[Residue]]]:
    dominant_residues = run_dssp(dominant_model)
    dominant_sequence = get_sequence(dominant_residues)
    assignments: Dict[str, List[Residue]] = {}
    comparisons: List[ModelComparison] = []
    assignments[dominant_model.model_name] = dominant_residues

    from utils.phenix_filter import get_mode_fractions
    from utils.phenix_filter import PHENIX_FILTERED_MODES
    for encounter_index, model in enumerate(alternative_models):
        fractions = get_mode_fractions(model.pdb_path, model.dssp_path.with_suffix(".phenix"),)
        rmsd, rmsd_error = rmsd_results.get(
            model.model_name, (None, "missing_rmsd_result")
        )
        comparison = ModelComparison(
            alternative_model=model.model_name,
            encounter_index=encounter_index,
            rmsd_to_dominant=rmsd,
            switch_segments=[],
            disordered_change_segments=[],
            # Phenix filter is cached (already calculated in pipelines.py)
            # needs to be stored in Comparison for eventual display per alternative model
            barbed_wire_fraction=sum(fractions.get(mode, 0.0) for mode in PHENIX_FILTERED_MODES),
        )
        if rmsd_error or rmsd is None:
            continue
        if rmsd < rmsd_threshold:
            continue

        alternative_residues = run_dssp(model)
        assignments[model.model_name] = alternative_residues
        alternative_sequence = get_sequence(alternative_residues)

        sequence_match = dominant_sequence == alternative_sequence
        comparison = replace(comparison, sequences_identical=sequence_match)
        if not sequence_match:
            raise ValueError(
                f"DSSP-derived sequences must be exact and pre-aligned; sequence mismatch between "
                f"{dominant_model.model_name} and {model.model_name}: "
                f"{len(dominant_sequence)} != {len(alternative_sequence)}"
            )

        try:
            switch_segments, disordered_change_segments, same_ss_preserved = compare_dssp_codes(
                dominant_residues,
                alternative_residues,
                min_run_length,
                MIN_REGION_SEPARATION_TO_KEEP_SEPARATE,
            )
        except Exception as exc:
            comparisons.append(
                add_invalid_reason(comparison, f"comparison_failed:{compact_reason(exc)}")
            )
            continue

        comparisons.append(replace(
            comparison,
            switch_segments=switch_segments,
            disordered_change_segments=disordered_change_segments,
            same_ss_preserved_percent=same_ss_preserved,
        ))

    return comparisons, assignments


# ---------------------------------------------------------------------------
# Hit selection
# ---------------------------------------------------------------------------

def add_invalid_reason(comparison: ModelComparison, reason: str) -> ModelComparison:
    reason = compact_reason(reason)
    existing = [part for part in comparison.invalid_reason.split(";") if part]
    if reason and reason not in existing:
        existing.append(reason)
    return replace(comparison, invalid_reason=";".join(existing))


def no_hit_comparison(n_searched: int) -> ModelComparison:
    # Only called if all possible structures are traversed
    return ModelComparison(
        alternative_model="none",
        encounter_index=0,
        switch_segments=[],
        disordered_change_segments=[],
        invalid_reason=compact_reason(
            f"no hits found after searching {n_searched} structures"
        ),
    )

def is_genuine_hit(
    #Also can be implemented later to filter out any wrongly classified hits and reduce false positives
    comparison: ModelComparison,
    assignments: Dict[str, List[Residue]],
    dominant_model_name: str,
) -> bool:
    if comparison.classification != "likely":
        return False
    return True

# Called by pipelines.py to identify the first hit in a PSE file
def first_hit(
    comparisons: Sequence[ModelComparison],
    assignments: Dict[str, List[Residue]],
    dominant_model_name: str,
) -> Optional[ModelComparison]:
    for comparison in comparisons:
        if is_genuine_hit(comparison, assignments, dominant_model_name):
            return comparison
        # recovery_regions = recovery_by_model.get(comparison.alternative_model, [])
        # if any(region.average_distance >= RECOVERY_REGION_AVERAGE_DISTANCE_THRESHOLD for region in recovery_regions):
        #     return comparison
    return None


def record_assignments(
    destination: Dict[str, List[Residue]],
    order: List[str],
    new_assignments: Dict[str, List[Residue]],
) -> None:
    for model_name, residues in new_assignments.items():
        if model_name not in destination:
            order.append(model_name)
        destination[model_name] = residues


def keep_going_summary_fields(
    comparisons: Sequence[ModelComparison],
) -> Dict[str, object]:
    likely_hits = [
        comparison.alternative_model
        for comparison in comparisons
        if comparison.classification == "likely"
    ]
    return {
        "keep_going_n_likely_hit_structures": len(likely_hits),
        "keep_going_likely_hit_pdb_files": ";".join(likely_hits),
    }