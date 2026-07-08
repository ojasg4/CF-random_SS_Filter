import math
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from utils.classes import ExtractedModel
from utils.dssp import normalize_dssp_code
from utils.structure import (
    extract_models_from_pse,
    find_dominant_model,
    find_keep_going_pdb_files,
    extracted_models_from_keep_going_pdb_files,
    pymol_session,
)
from .recovery_classes import RecoveryResult, ResidueSS
from .recovery_logic import (
    detect_bridge_reshuffle,
    detect_disorder_to_beta,
    detect_disorder_to_helix,
    detect_rigid_body_helix,
    DETECTOR_BRIDGE_RESHUFFLE,
    DETECTOR_DISORDER_BETA,
    DETECTOR_DISORDER_HELIX,
    DETECTOR_RIGID_BODY,
)
from .recovery_io import write_recovery_log, result_to_row, load_residues, dssp_string, sequence_from_residues


def _ca_distances(
    cmd,
    dominant_pdb: Path,
    alternative_pdb: Path,
    dominant: Sequence[ResidueSS],
    alternative: Sequence[ResidueSS],
    obj_idx: int,
) -> Dict[int, float]:
    dom_obj, alt_obj = f"rb_dom_{obj_idx}", f"rb_alt_{obj_idx}"
    distances: Dict[int, float] = {}
    try:
        cmd.load(str(dominant_pdb), dom_obj)
        cmd.load(str(alternative_pdb), alt_obj)
        cmd.cealign(f"{dom_obj} and name CA", f"{alt_obj} and name CA", quiet=1)

        def coords(obj):
            return {
                (a.chain.strip(), str(a.resi).strip()):
                (float(a.coord[0]), float(a.coord[1]), float(a.coord[2]))
                for a in cmd.get_model(f"{obj} and name CA").atom
            }

        dc, ac = coords(dom_obj), coords(alt_obj)
        for idx, (dr, ar) in enumerate(zip(dominant, alternative), start=1):
            dk = (dr.chain_id.strip(), dr.residue_number.strip())
            ak = (ar.chain_id.strip(), ar.residue_number.strip())
            if dk in dc and ak in ac:
                distances[idx] = math.sqrt(
                    sum((a - b) ** 2 for a, b in zip(dc[dk], ac[ak]))
                )
    except Exception:
        pass
    finally:
        for obj in [dom_obj, alt_obj]:
            try: cmd.delete(obj)
            except Exception: pass
    return distances


def _analyze_candidate(
    dominant_model: ExtractedModel,
    alternative_model: ExtractedModel,
    source: str,
    encounter_index: int,
    assignments: Dict[str, List[ResidueSS]],
    args,
    cmd=None,
    obj_idx: int = 0,
) -> RecoveryResult:
    if dominant_model.model_name not in assignments:
        assignments[dominant_model.model_name] = load_residues(dominant_model)
    dom = assignments[dominant_model.model_name]

    try:
        alt = load_residues(alternative_model)
    except Exception as exc:
        return RecoveryResult(
            alternative_model=alternative_model.model_name, source=source,
            encounter_index=encounter_index, sequences_identical=None,
            detector_hit="", invalid_reason=str(exc),
        )
    assignments[alternative_model.model_name] = alt

    if sequence_from_residues(dom) != sequence_from_residues(alt):
        return RecoveryResult(
            alternative_model=alternative_model.model_name, source=source,
            encounter_index=encounter_index, sequences_identical=False,
            detector_hit="", invalid_reason="sequence_mismatch",
        )

    dom_dssp, alt_dssp = dssp_string(dom), dssp_string(alt)

    # DETECTOR 1: Beta bridge-partner reshuffle (toggleable)
    if args.bridge_reshuffling:
        regions = detect_bridge_reshuffle(dom, alt)
        if regions:
            return RecoveryResult(
                alternative_model=alternative_model.model_name, source=source,
                encounter_index=encounter_index, sequences_identical=True,
                detector_hit=DETECTOR_BRIDGE_RESHUFFLE,
                bridge_reshuffled_regions=regions,
                dominant_dssp=dom_dssp, alternative_dssp=alt_dssp,
            )

    # DETECTOR 2: Disorder to or from Beta
    regions = detect_disorder_to_beta(dom, alt, args.min_run_length)
    if regions:
        return RecoveryResult(
            alternative_model=alternative_model.model_name, source=source,
            encounter_index=encounter_index, sequences_identical=True,
            detector_hit=DETECTOR_DISORDER_BETA,
            disorder_change_regions=regions,
            dominant_dssp=dom_dssp, alternative_dssp=alt_dssp,
        )

    # DETECTOR 3: Disorder to or from Helix
    regions = detect_disorder_to_helix(dom, alt, args.min_run_length)
    if regions:
        return RecoveryResult(
            alternative_model=alternative_model.model_name, source=source,
            encounter_index=encounter_index, sequences_identical=True,
            detector_hit=DETECTOR_DISORDER_HELIX,
            disorder_change_regions=regions,
            dominant_dssp=dom_dssp, alternative_dssp=alt_dssp,
        )

    # DETECTOR 4: Rigid-body helix displacement
    if cmd is not None:
        distances = _ca_distances(
            cmd, dominant_model.pdb_path, alternative_model.pdb_path,
            dom, alt, obj_idx,
        )
        regions = detect_rigid_body_helix(dom, alt, distances, args.min_run_length)
        if regions:
            return RecoveryResult(
                alternative_model=alternative_model.model_name, source=source,
                encounter_index=encounter_index, sequences_identical=True,
                detector_hit=DETECTOR_RIGID_BODY,
                rigid_body_regions=regions,
                dominant_dssp=dom_dssp, alternative_dssp=alt_dssp,
            )

    return RecoveryResult(
        alternative_model=alternative_model.model_name, source=source,
        encounter_index=encounter_index, sequences_identical=True,
        detector_hit="", dominant_dssp=dom_dssp, alternative_dssp=alt_dssp,
    )


def run_recovery_pipeline(
    args,
    job_name: str,
    input_pse: Path,
    work_dir: Path,
    raw_dssp_dir: Path,
    summary_base_dir: Path,
) -> Tuple[Dict, bool]:
    log_file = work_dir / "recovery.txt"
    assignments: Dict[str, List[ResidueSS]] = {}
    results_by_model: Dict[str, RecoveryResult] = {}
    selected: Optional[RecoveryResult] = None
    n_analyzed = 0

    try:
        # 1. Extract structures from PSE
        models = extract_models_from_pse(input_pse, work_dir, raw_dssp_dir)
        dominant_model = find_dominant_model(models, args.dominant_label_patterns)
        main_alternatives = [
            m for m in models if m.model_name != dominant_model.model_name
        ]

        with pymol_session() as cmd:
            # 2. Search PSE clusters
            for enc_idx, alt in enumerate(main_alternatives):
                result = _analyze_candidate(
                    dominant_model, alt, "main_pse", enc_idx,
                    assignments, args, cmd, enc_idx,
                )
                results_by_model[result.alternative_model] = result
                n_analyzed += 1
                if result.recovered and selected is None:
                    selected = result

            # 3. If no hit, walk keep-going PDB directories
            if selected is None:
                keep_going_models = extracted_models_from_keep_going_pdb_files(
                    find_keep_going_pdb_files(input_pse),
                    raw_dssp_dir,
                    input_pse.resolve().parent,
                )
                enc_base = len(main_alternatives)
                for enc_idx, alt in enumerate(keep_going_models):
                    result = _analyze_candidate(
                        dominant_model, alt, "keep_going",
                        enc_base + enc_idx, assignments, args,
                        cmd, enc_base + enc_idx,
                    )
                    results_by_model[result.alternative_model] = result
                    n_analyzed += 1
                    if result.recovered:
                        selected = result
                        break

        # 4. Write log and return row
        write_recovery_log(
            log_file, job_name, input_pse,
            dominant_model.model_name, assignments,
            results_by_model, selected, args, summary_base_dir,
        )

        if selected and selected.recovered:
            return result_to_row(
                selected, input_pse, "recovered", n_analyzed,
                log_file, raw_dssp_dir, summary_base_dir,
            ), True

        best_failed = next(
            (r for r in sorted(
                results_by_model.values(), key=lambda r: r.encounter_index
            ) if not r.invalid_reason),
            None,
        )
        return result_to_row(
            best_failed, input_pse, "failed", n_analyzed,
            log_file, raw_dssp_dir, summary_base_dir,
        ), False

    except Exception as exc:
        write_recovery_log(
            log_file, job_name, input_pse, "Dominant",
            assignments, results_by_model, None, args, summary_base_dir,
        )
        return result_to_row(
            None, input_pse, "failed", n_analyzed,
            log_file, raw_dssp_dir, summary_base_dir,
            invalid_reason=str(exc),
        ), False