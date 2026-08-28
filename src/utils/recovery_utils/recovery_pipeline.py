import math
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from utils.classes import ExtractedModel
from utils.dssp import normalize_dssp_code
from utils.pdb_helpers import (
    find_dominant_model,
    pymol_session,
)

from utils.chainsaw_filter import parse_chainsaw_tsv

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
from utils.chainsaw_filter import (
    score_work_directory,
    passes_chainsaw_filter,
    chainsaw_summary_fields,
)
from utils.config import CFG
HELIX_DSSP_CODES      = set(CFG["dssp_codes"]["helix"])
BETA_DSSP_CODES       = set(CFG["dssp_codes"]["beta"])
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
    if args.strict: # if a case has the same exact ordered residue counts within 1, ignore
        print("Make sure to enable the final filter in recovery logic.py as well!")
        alt_helix_n  = sum(1 for r in alt_dssp if r in HELIX_DSSP_CODES)
        dom_helix_n  = sum(1 for r in dom_dssp if r in HELIX_DSSP_CODES)
        alt_strand_n = sum(1 for r in alt_dssp if r in BETA_DSSP_CODES)
        dom_strand_n = sum(1 for r in dom_dssp if r in BETA_DSSP_CODES)
        composition_shift = abs(alt_helix_n - dom_helix_n) + abs(alt_strand_n - dom_strand_n)
        if composition_shift < 1:
            return RecoveryResult(
                alternative_model=alternative_model.model_name, source=source,
                encounter_index=encounter_index, sequences_identical=True,
                detector_hit="", dominant_dssp=dom_dssp, alternative_dssp=alt_dssp,
            )

    # DETECTOR 1: Disorder to or from Beta
    regions = detect_disorder_to_beta(dom, alt, args.min_run_length)
    if regions:
        return RecoveryResult(
            alternative_model=alternative_model.model_name, source=source,
            encounter_index=encounter_index, sequences_identical=True,
            detector_hit=DETECTOR_DISORDER_BETA,
            disorder_change_regions=regions,
            dominant_dssp=dom_dssp, alternative_dssp=alt_dssp,
        )

    # DETECTOR 2: Disorder to or from Helix
    regions = detect_disorder_to_helix(dom, alt, args.min_run_length)
    if regions:
        return RecoveryResult(
            alternative_model=alternative_model.model_name, source=source,
            encounter_index=encounter_index, sequences_identical=True,
            detector_hit=DETECTOR_DISORDER_HELIX,
            disorder_change_regions=regions,
            dominant_dssp=dom_dssp, alternative_dssp=alt_dssp,
        )

    # DETECTOR 3: Rigid-body helix displacement
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

    # DETECTOR 4: Beta bridge-partner reshuffle (toggleable)
    if args.beta_bridges:
        regions = detect_bridge_reshuffle(dom, alt)
        if regions:
            return RecoveryResult(
                alternative_model=alternative_model.model_name, source=source,
                encounter_index=encounter_index, sequences_identical=True,
                detector_hit=DETECTOR_BRIDGE_RESHUFFLE,
                bridge_reshuffled_regions=regions,
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
        # 1. Load all cached structures from the work directory.
        #    The main pipeline (run with --keep-work) already extracted every
        #    good-prediction PDB here, alongside its .dssp and .phenix cache.
        #    Nothing is re-extracted or recomputed.
        pdbs = sorted(work_dir.glob("*.pdb"))
        if not pdbs:
            raise ValueError(
                f"No cached PDBs in {work_dir}. "
                f"Run the main pipeline with --keep-work first."
            )
        models = [
            ExtractedModel(
                model_name=p.stem,
                pdb_path=p,
                dssp_path=p.with_suffix(".dssp"),
            )
            for p in pdbs
        ]

        dominant_model = find_dominant_model(models, args.dominant_label_patterns)

        main_alternatives = [
            m for m in models if m.model_name != dominant_model.model_name
        ]

        with pymol_session() as cmd:
            # Scan every cached structure (PSE-cluster and keep-going alike —
            # both already live in the work dir).
            for enc_idx, alt in enumerate(main_alternatives):
                result = _analyze_candidate(
                    dominant_model, alt, "work_dir", enc_idx,
                    assignments, args, cmd, enc_idx,
                )
                results_by_model[result.alternative_model] = result
                n_analyzed += 1
                if result.recovered and selected is None:
                    selected = result
            
            # Chainsaw
            chainsaw_rows = {}
            # Chainsaw — read the cached TSV written by the main pipeline
            if getattr(args, "chainsaw", False):
                tsv_path = work_dir / "chainsaw_results.tsv"
                if not tsv_path.exists():
                    raise ValueError(
                        f"Missing cached chainsaw_results.tsv at {tsv_path}. "
                        f"Run the main pipeline with --chainsaw --keep-work first."
                    )
                chainsaw_rows = parse_chainsaw_tsv(tsv_path)

                if selected is not None:
                    recovered_hits = sorted(
                        (r for r in results_by_model.values() if r.recovered),
                        key=lambda r: r.encounter_index,
                    )
                    selected = None
                    for candidate in recovered_hits:
                        row = chainsaw_rows.get(candidate.alternative_model)
                        if passes_chainsaw_filter(row):
                            selected = candidate
                            break

        # 4. Write log and return row
        write_recovery_log(
            log_file, job_name, input_pse,
            dominant_model.model_name, assignments,
            results_by_model, selected, args, summary_base_dir,
        )

        if selected and selected.recovered:
            chainsaw_fields = {}
            if getattr(args, "chainsaw", False):
                chainsaw_fields.update(
                    chainsaw_summary_fields("dom", chainsaw_rows.get(dominant_model.model_name))
                )
                chainsaw_fields.update(
                    chainsaw_summary_fields("alt", chainsaw_rows.get(selected.alternative_model))
                )
            return result_to_row(
                selected, input_pse, "recovered", n_analyzed,
                log_file, raw_dssp_dir, summary_base_dir,
                chainsaw_fields=chainsaw_fields,
            ), True

        best_failed = next(
            (r for r in sorted(
                results_by_model.values(), key=lambda r: r.encounter_index
            ) if not r.invalid_reason),
            None,
        )
        chainsaw_fields = {}
        if getattr(args, "chainsaw", False):
            chainsaw_fields.update(
                chainsaw_summary_fields("dom", chainsaw_rows.get(dominant_model.model_name))
            )
            if best_failed is not None:
                chainsaw_fields.update(
                    chainsaw_summary_fields("alt", chainsaw_rows.get(best_failed.alternative_model))
                )
        return result_to_row(
            best_failed, input_pse, "failed", n_analyzed,
            log_file, raw_dssp_dir, summary_base_dir,
            chainsaw_fields=chainsaw_fields,
        ), False

    except Exception as exc:
        raise
        write_recovery_log(
            log_file, job_name, input_pse, "Dominant",
            assignments, results_by_model, None, args, summary_base_dir,
        )
        return result_to_row(
            None, input_pse, "failed", n_analyzed,
            log_file, raw_dssp_dir, summary_base_dir,
            invalid_reason=str(exc),
        ), False