import re
from typing import Dict, List, Optional, Sequence, Tuple

from .recovery_classes import (
    BridgeReshuffledRegion, ChangedRegion, RigidBodyRegion, ResidueSS, OpenDisorderSegment
)

from utils.dssp import normalize_dssp_code
from utils.config import CFG, MIN_REGION_SEPARATION_TO_KEEP_SEPARATE

HELIX_DSSP_CODES      = set(CFG["dssp_codes"]["helix"])
BETA_DSSP_CODES       = set(CFG["dssp_codes"]["beta"])
DISORDERED_DSSP_CODES = set(CFG["dssp_codes"]["disordered"])

MIN_BETA_BRIDGE_RESHUFFLED_RESIDUES = 2
BETA_BRIDGE_PARTNER_SHIFT_TOLERANCE = 5
RIGID_BODY_WINDOW_SIZE              = 10
RIGID_BODY_CA_DISTANCE_THRESHOLD    = 3000 #in Angstroms, Basically turning it off for debugging

DETECTOR_BRIDGE_RESHUFFLE = "beta_bridge_reshuffle"
DETECTOR_DISORDER_BETA    = "disorder_to_beta"
DETECTOR_DISORDER_HELIX   = "disorder_to_helix"
DETECTOR_RIGID_BODY       = "rigid_body_helix"


def _bp_map(residues: Sequence[ResidueSS]) -> Dict[int, str]:
    return {r.position_index: r.residue_number for r in residues}


def _resi_int(resi: Optional[str]) -> Optional[int]:
    if resi is None:
        return None
    try:
        return int(re.sub(r"[^0-9\-]", "", resi))
    except (ValueError, TypeError):
        return None


def _partner_reshuffled(
    dom: ResidueSS, alt: ResidueSS,
    dom_map: Dict[int, str], alt_map: Dict[int, str],
    tolerance: int,
) -> bool:
    for di, ai in [(dom.bp1, alt.bp1), (dom.bp2, alt.bp2)]:
        dp = dom_map.get(di) if di != 0 else None
        ap = alt_map.get(ai) if ai != 0 else None
        if (dp is None) != (ap is None):
            return True
        if dp and ap:
            dint, aint = _resi_int(dp), _resi_int(ap)
            if dint is not None and aint is not None and abs(dint - aint) > tolerance:
                return True
    return False


def detect_bridge_reshuffle(
    dominant: Sequence[ResidueSS],
    alternative: Sequence[ResidueSS],
) -> List[BridgeReshuffledRegion]:
    if len(dominant) != len(alternative):
        return []
    dom_map = _bp_map(dominant)
    alt_map = _bp_map(alternative)
    raw: List[BridgeReshuffledRegion] = []
    current: List[int] = []

    for idx, (dr, ar) in enumerate(zip(dominant, alternative), start=1):
        both_beta = (
            normalize_dssp_code(dr.dssp_code) in BETA_DSSP_CODES and
            normalize_dssp_code(ar.dssp_code) in BETA_DSSP_CODES
        )
        if both_beta and _partner_reshuffled(dr, ar, dom_map, alt_map, BETA_BRIDGE_PARTNER_SHIFT_TOLERANCE):
            current.append(idx)
        else:
            if len(current) >= MIN_BETA_BRIDGE_RESHUFFLED_RESIDUES:
                raw.append(BridgeReshuffledRegion(current[0], current[-1], list(current)))
            current = []

    # Close any run still open at end of sequence
    if len(current) >= MIN_BETA_BRIDGE_RESHUFFLED_RESIDUES:
        raw.append(BridgeReshuffledRegion(current[0], current[-1], list(current)))

    if not raw:
        return []
    merged: List[BridgeReshuffledRegion] = [raw[0]]
    for region in raw[1:]:
        previous = merged[-1]
        if region.start - previous.end - 1 < MIN_REGION_SEPARATION_TO_KEEP_SEPARATE:
            merged[-1] = BridgeReshuffledRegion(
                start=previous.start,
                end=region.end,
                positions=sorted(set(previous.positions + region.positions)),
            )
        else:
            merged.append(region)
    return merged
def _detect_disorder_transition(
    dominant: Sequence[ResidueSS],
    alternative: Sequence[ResidueSS],
    ordered_codes: set,
    ordered_group_name: str,
    min_run: int,
) -> List[ChangedRegion]:
    if len(dominant) != len(alternative):
        return []

    raw: List[ChangedRegion] = []
    open_segment: Optional[OpenDisorderSegment] = None

    for idx, (dr, ar) in enumerate(zip(dominant, alternative), start=1):
        dc = normalize_dssp_code(dr.dssp_code)
        ac = normalize_dssp_code(ar.dssp_code)

        # Check if this position is a disorder-to-ordered transition in either direction
        case_a = dc in DISORDERED_DSSP_CODES and ac in ordered_codes
        case_b = dc in ordered_codes and ac in DISORDERED_DSSP_CODES
        is_transition = case_a or case_b

        new_dom = "COIL" if dc in DISORDERED_DSSP_CODES else ordered_group_name
        new_alt = "COIL" if ac in DISORDERED_DSSP_CODES else ordered_group_name
        is_same_run = (
            open_segment is not None
            and open_segment.dominant_code_group == new_dom
            and open_segment.alternative_code_group == new_alt
        )

        # Close the open segment if this position doesn't continue it
        if open_segment is not None and (not is_transition or not is_same_run):
            if idx - open_segment.start >= min_run:
                raw.append(ChangedRegion(
                    start=open_segment.start,
                    end=idx - 1,
                    dominant_code_group=open_segment.dominant_code_group,
                    alternative_code_group=open_segment.alternative_code_group,
                ))
            open_segment = None

        # Open a new segment if this position is a transition
        if is_transition and open_segment is None:
            open_segment = OpenDisorderSegment(
                start=idx,
                dominant_code_group=new_dom,
                alternative_code_group=new_alt,
            )

    # Close any segment still open at end of sequence
    if open_segment is not None and len(dominant) - open_segment.start + 1 >= min_run:
        raw.append(ChangedRegion(
            start=open_segment.start,
            end=len(dominant),
            dominant_code_group=open_segment.dominant_code_group,
            alternative_code_group=open_segment.alternative_code_group,
        ))

    # Merge adjacent regions separated to account for turns between beta strands/helices
    if not raw:
        return []
    merged: List[ChangedRegion] = [raw[0]]
    for region in raw[1:]:
        previous = merged[-1]
        if region.start - previous.end - 1 < MIN_REGION_SEPARATION_TO_KEEP_SEPARATE:
            merged[-1] = ChangedRegion(
                start=previous.start,
                end=region.end,
                dominant_code_group=previous.dominant_code_group,
                alternative_code_group=(
                    previous.alternative_code_group
                    if previous.alternative_code_group == region.alternative_code_group
                    else f"{previous.alternative_code_group}/{region.alternative_code_group}"
                ),
            )
        else:
            merged.append(region)
    return merged


def detect_disorder_to_beta(dominant, alternative, min_run):
    return _detect_disorder_transition(
        dominant, alternative, BETA_DSSP_CODES, "BETA", min_run
    )


def detect_disorder_to_helix(dominant, alternative, min_run):
    return _detect_disorder_transition(
        dominant, alternative, HELIX_DSSP_CODES, "HELIX", min_run
    )


def _sliding_window_max(
    distances: Dict[int, float],
    positions: List[int],
    window_size: int,
    ca_threshold: float,
) -> Tuple[float, bool]:
    max_mean = 0.0
    exceeded = False
    for i in range(max(1, len(positions) - window_size + 1)):
        window = positions[i: i + window_size]
        dists = [distances[p] for p in window if p in distances]
        if not dists:
            continue
        mean = sum(dists) / len(dists)
        max_mean = max(max_mean, mean)
        if mean >= ca_threshold:
            exceeded = True
    return max_mean, exceeded


def detect_rigid_body_helix(
    dominant: Sequence[ResidueSS],
    alternative: Sequence[ResidueSS],
    distances: Dict[int, float],
    min_run: int,
) -> List[RigidBodyRegion]:
    helix_runs: List[List[int]] = []
    current: List[int] = []
    for idx, (dr, ar) in enumerate(zip(dominant, alternative), start=1):
        if (normalize_dssp_code(dr.dssp_code) in HELIX_DSSP_CODES and
                normalize_dssp_code(ar.dssp_code) in HELIX_DSSP_CODES):
            current.append(idx)
        else:
            if len(current) >= min_run:
                helix_runs.append(current)
            current = []
    if len(current) >= min_run:
        helix_runs.append(current)

    regions: List[RigidBodyRegion] = []
    for run in helix_runs:
        max_win, exceeded = _sliding_window_max(
            distances, run, RIGID_BODY_WINDOW_SIZE, RIGID_BODY_CA_DISTANCE_THRESHOLD
        )
        if exceeded:
            codes = [normalize_dssp_code(dominant[i - 1].dssp_code) for i in run]
            regions.append(RigidBodyRegion(
                run[0], run[-1],
                max(set(codes), key=codes.count),
                max_win,
            ))
    return regions