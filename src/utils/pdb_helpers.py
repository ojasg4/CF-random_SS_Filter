import math
import re
from contextlib import contextmanager
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from utils.config import MIN_REGION_SEPARATION_TO_KEEP_SEPARATE
from utils.classes import ExtractedModel, ModelComparison, Residue
from utils.dssp import run_dssp, get_sequence
from utils.utils import compact_reason, relative_path_for_summary, safe_path_component


# ---------------------------------------------------------------------------
# PyMOL session utils
# ---------------------------------------------------------------------------

@contextmanager
def pymol_session():
    try:
        import pymol2

        session = pymol2.PyMOL()
        session.start()
        cmd = session.cmd
    except Exception:
        import pymol

        pymol.finish_launching(["pymol", "-cq"])
        session, cmd = None, pymol.cmd
    try:
        yield cmd
    finally:
        if session is not None:
            session.stop()


def model_cache_path(work_dir: Path, model_name: str, suffix: str) -> Path:
    """
    Return the work-dir cache path for a model.

    The model_name is the canonical spreadsheet/log key.

    Examples:
      model_name="Dominant", suffix=".dssp"
        -> <work_dir>/Dominant.dssp

      model_name="WP_000590537_predicted_models_rand_54_max_2_ext_4/WP_000590537_unrelaxed_rank_001_alphafold2_ptm_model_2_seed_003.pdb",
      suffix=".dssp"
        -> <work_dir>/WP_000590537_predicted_models_rand_54_max_2_ext_4/
           WP_000590537_unrelaxed_rank_001_alphafold2_ptm_model_2_seed_003.dssp
    """
    rel = Path(model_name)

    if rel.is_absolute() or ".." in rel.parts:
        raise ValueError(
            f"Model name must be a relative cache key, not an absolute/escaping path: {model_name}"
        )

    if rel.suffix.lower() == ".pdb":
        rel = rel.with_suffix(suffix)
    else:
        rel = Path(str(rel) + suffix)

    output_path = work_dir / rel
    output_path.parent.mkdir(parents=True, exist_ok=True)
    return output_path


def pymol_residue_selection(
    object_name: str, residues: Sequence[Residue], positions: Iterable[int]
) -> str:
    terms = []
    residues_by_position = dict(enumerate(residues, start=1))
    for position in sorted(set(positions)):
        residue = residues_by_position.get(position)
        if residue is None or not residue.residue_number:
            continue
        resi = f"{residue.residue_number}"
        term = f"({object_name} and resi {resi}"
        if residue.chain_id:
            term += f" and chain {residue.chain_id}"
        term += ")"
        terms.append(term)
    return " or ".join(terms) if terms else "none"


def pymol_ca_rows(
    cmd, selection: str
) -> List[Tuple[str, str, Tuple[float, float, float]]]:
    model = cmd.get_model(selection)
    return [
        (
            atom.chain.strip(),
            str(atom.resi).strip(),
            (float(atom.coord[0]), float(atom.coord[1]), float(atom.coord[2])),
        )
        for atom in model.atom
    ]


# ---------------------------------------------------------------------------
# Extracting structures from a PSE file (or from keep-going PDB siblings)
# ---------------------------------------------------------------------------

def extract_models_from_pse(
    pse_path: Path, extracted_pdb_dir: Path, raw_dssp_dir: Path, cmd
) -> List[ExtractedModel]:
    """
    Extract PSE objects into work_dir.

    Here raw_dssp_dir is treated as the DSSP/Phenix cache root. In the patched
    pipeline, raw_dssp_dir will be work_dir.
    """
    extracted_pdb_dir.mkdir(parents=True, exist_ok=True)
    raw_dssp_dir.mkdir(parents=True, exist_ok=True)

    models: List[ExtractedModel] = []

    cmd.reinitialize()
    cmd.load(str(pse_path))
    objects = cmd.get_object_list("all")
    print(f"PSE objects found: {objects}")

    for obj in objects:
        try:
            atom_count = cmd.count_atoms(obj)
        except Exception:
            atom_count = 0

        if atom_count == 0:
            continue

        # Use one canonical model key everywhere.
        model_name = safe_path_component(obj)

        pdb_path = extracted_pdb_dir / f"{model_name}.pdb"
        dssp_path = model_cache_path(raw_dssp_dir, model_name, ".dssp")

        cmd.save(str(pdb_path), selection=obj, state=1)

        if pdb_path.exists() and pdb_path.stat().st_size > 0:
            models.append(
                ExtractedModel(
                    model_name=model_name,
                    pdb_path=pdb_path,
                    dssp_path=dssp_path,
                )
            )

    if not models:
        raise ValueError(f"No molecular objects extracted from {pse_path}")

    return models


def find_dominant_model(
    models: Sequence[ExtractedModel], dominant_labels: Sequence[str]
) -> ExtractedModel:
    matches = [model for model in models if model.model_name in dominant_labels]
    if len(matches) > 1:
        names = "\n".join(f"  - {model.model_name}" for model in models)
        raise ValueError(
            f"Expected exactly one dominant/reference model with an exact case-sensitive label, but found {len(matches)}.\n"
            f"Dominant labels: {', '.join(dominant_labels)}\n"
            f"Available models:\n{names}"
        )
    if not matches:
        #picks the first model if there are none named "Dominant"/[whatever the override CLI input was] in the given pse file
        return models[0]
    return matches[0]

def keep_going_directory_max_value(path: Path) -> int:
    # ignores the full MSA sampled models. Looks for "_max_" in pattern directory to get all of the subsampled MSA directories
    match = re.search(r"_max_(\d+)(?:_|$)", path.name)
    return int(match.group(1)) if match else -1


def find_keep_going_pdb_files(pse_file: Path) -> List[Path]:
    parent_dir = pse_file.resolve().parent
    model_dirs = sorted(
        (
            path
            for path in parent_dir.iterdir()
            if path.is_dir()
            and "_predicted_models_rand_" in path.name
        ),
        key=lambda path: (keep_going_directory_max_value(path), path.name),
        reverse=True,
    )
    if not model_dirs:
        raise ValueError(
            f"No keep-going directories matching *_predicted_models_rand_* were found in {parent_dir}"
        )
    pdb_files: List[Path] = []
    for model_dir in model_dirs:
        pdbs = sorted(path for path in model_dir.rglob("*.pdb") if path.is_file())
        if len(pdbs) != 25:
            raise ValueError(
                f"Expected exactly 25 PDB files in keep-going directory {model_dir}, but found {len(pdbs)}"
            )
        pdb_files.extend(pdbs)
    return pdb_files


def extracted_models_from_keep_going_pdb_files(
    pdb_files: Sequence[Path],
    raw_dssp_dir: Path,
    pse_parent: Path,
) -> List[ExtractedModel]:
    """
    Build ExtractedModel objects for keep-going PDBs.

    The source PDB remains wherever it already lives, usually under the
    original CF-random directory.

    The canonical model_name is the relative PDB path from the PSE parent.
    This is the value that should appear in:
      - models_by_name
      - ModelComparison.alternative_model
      - final_comparison.txt
      - summary CSV alt_model
      - spreadsheet

    The DSSP cache is written under raw_dssp_dir, which the patched pipeline
    passes as work_dir.
    """
    models: List[ExtractedModel] = []
    raw_dssp_dir.mkdir(parents=True, exist_ok=True)

    for pdb_path in pdb_files:
        model_name = relative_path_for_summary(pdb_path, pse_parent)
        dssp_path = model_cache_path(raw_dssp_dir, model_name, ".dssp")

        models.append(
            ExtractedModel(
                model_name=model_name,
                pdb_path=pdb_path,
                dssp_path=dssp_path,
            )
        )

    return models
# ---------------------------------------------------------------------------
# RMSD to dominant
# ---------------------------------------------------------------------------

def calculate_rmsds_to_dominant(
    cmd,
    dominant_model: ExtractedModel,
    alternative_models: Sequence[ExtractedModel],
) -> Dict[str, Tuple[Optional[float], str]]:
    output: Dict[str, Tuple[Optional[float], str]] = {}
    cmd.reinitialize()
    dominant_object = "rmsd_dominant"
    cmd.load(str(dominant_model.pdb_path), dominant_object)
    target_selection = f"{dominant_object} and name CA"
    if cmd.count_atoms(target_selection) == 0:
        target_selection = dominant_object
    for idx, model in enumerate(alternative_models, start=1):
        object_name = f"rmsd_alt_{idx}"
        try:
            cmd.load(str(model.pdb_path), object_name)
            mobile_selection = f"{object_name} and name CA"
            if cmd.count_atoms(mobile_selection) == 0:
                mobile_selection = object_name
            if (
                cmd.count_atoms(mobile_selection) == 0
                or cmd.count_atoms(target_selection) == 0
            ):
                output[model.model_name] = (None, "no_atoms_for_rmsd")
            else:
                result = cmd.cealign(target_selection, mobile_selection, quiet=1)
                rmsd = float(result["RMSD"])
                if not math.isfinite(rmsd) or rmsd >= 100000.0:
                    output[model.model_name] = (
                        None,
                        "rmsd_calculation_failed:cealign_unusable_result",
                    )
                else:
                    output[model.model_name] = (rmsd, "")
        except Exception as exc:
            output[model.model_name] = (
                None,
                f"rmsd_calculation_failed:{compact_reason(exc)}",
            )
        finally:
            try:
                cmd.delete(object_name)
            except Exception:
                pass
    return output


def analyze_alternatives_with_dssp(
    cmd,
    dominant_model: ExtractedModel,
    alternative_models: Sequence[ExtractedModel],
    min_run_length: int,
    rmsd_threshold: float,
) -> Tuple[List[ModelComparison], Dict[str, List[Residue]]]:
    from utils.comparison import compare_candidate_models
    rmsd_results = calculate_rmsds_to_dominant(cmd, dominant_model, alternative_models)
    # skipped for now to save time (3-4 seconds per pdb)
    # rmsd_results = {
    #     m.model_name: (999.0, "")
    #     for m in alternative_models
    # }

    return compare_candidate_models(
        dominant_model,
        alternative_models,
        rmsd_results,
        min_run_length,
        rmsd_threshold,
    )


# ---------------------------------------------------------------------------
# Recovery regions: same-DSSP-code residues that still moved in space
# ---------------------------------------------------------------------------

def residue_key(residue: Residue) -> Tuple[str, str]:
    return residue.chain_id.strip(), f"{residue.residue_number}".strip()


def distance_between_points(
    left: Tuple[float, float, float], right: Tuple[float, float, float]
) -> float:
    return math.sqrt(
        (left[0] - right[0]) ** 2
        + (left[1] - right[1]) ** 2
        + (left[2] - right[2]) ** 2
    )


def ca_distances_by_position_after_alignment(
    cmd,
    dominant_model: ExtractedModel,
    alternative_model: ExtractedModel,
    dominant_residues: Sequence[Residue],
    alternative_residues: Sequence[Residue],
    object_index: int,
) -> Dict[int, float]:
    dominant_object = f"recover_dominant_{object_index}"
    alternative_object = f"recover_alternative_{object_index}"
    try:
        cmd.load(str(dominant_model.pdb_path), dominant_object)
        cmd.load(str(alternative_model.pdb_path), alternative_object)
        dominant_selection = f"{dominant_object} and name CA"
        alternative_selection = f"{alternative_object} and name CA"
        if cmd.count_atoms(dominant_selection) == 0:
            dominant_selection = dominant_object
        if cmd.count_atoms(alternative_selection) == 0:
            alternative_selection = alternative_object
        cmd.cealign(dominant_selection, alternative_selection, quiet=1)
        dominant_rows = pymol_ca_rows(cmd, f"{dominant_object} and name CA")
        alternative_rows = pymol_ca_rows(cmd, f"{alternative_object} and name CA")
        dominant_by_key = {(chain, resi): coord for chain, resi, coord in dominant_rows}
        alternative_by_key = {
            (chain, resi): coord for chain, resi, coord in alternative_rows
        }
        dominant_by_order = [coord for _, _, coord in dominant_rows]
        alternative_by_order = [coord for _, _, coord in alternative_rows]
        distances: Dict[int, float] = {}
        for idx, (dominant_residue, alternative_residue) in enumerate(
            zip(dominant_residues, alternative_residues), start=1
        ):
            dominant_coord = dominant_by_key.get(residue_key(dominant_residue))
            alternative_coord = alternative_by_key.get(residue_key(alternative_residue))
            if dominant_coord is None and idx <= len(dominant_by_order):
                dominant_coord = dominant_by_order[idx - 1]
            if alternative_coord is None and idx <= len(alternative_by_order):
                alternative_coord = alternative_by_order[idx - 1]
            if dominant_coord is not None and alternative_coord is not None:
                distances[idx] = distance_between_points(
                    dominant_coord, alternative_coord
                )
        return distances
    except Exception:
        return {}
    finally:
        try:
            cmd.delete(dominant_object)
        except Exception:
            pass
        try:
            cmd.delete(alternative_object)
        except Exception:
            pass


def close_recovery_region(
    output: List,
    positions: List[int],
    dssp_code: Optional[str],
    distances_by_position: Dict[int, float],
    min_run_length: int,
) -> None:
    from utils.classes import RecoveryRegion
    if len(positions) < min_run_length or dssp_code is None:
        return
    distances = [
        distances_by_position[position]
        for position in positions
        if position in distances_by_position
    ]
    if not distances:
        return
    output.append(
        RecoveryRegion(
            start=positions[0],
            end=positions[-1],
            dssp_code=dssp_code,
            positions=list(positions),
            average_distance=sum(distances) / len(distances),
        )
    )


def same_dssp_recovery_regions(
    dominant_residues: Sequence[Residue],
    alternative_residues: Sequence[Residue],
    distances_by_position: Dict[int, float],
    min_run_length: int,
    average_distance_threshold: float,
) -> List:
    from utils.dssp import normalize_dssp_code
    raw_regions = []
    current_positions: List[int] = []
    current_code: Optional[str] = None

    for position, (dominant_residue, alternative_residue) in enumerate(
        zip(dominant_residues, alternative_residues), start=1
    ):
        dominant_code = normalize_dssp_code(dominant_residue.dssp_code)
        alternative_code = normalize_dssp_code(alternative_residue.dssp_code)
        if dominant_code != alternative_code:
            close_recovery_region(
                raw_regions,
                current_positions,
                current_code,
                distances_by_position,
                min_run_length,
            )
            current_positions = []
            current_code = None
            continue
        if current_positions and current_code != dominant_code:
            close_recovery_region(
                raw_regions,
                current_positions,
                current_code,
                distances_by_position,
                min_run_length,
            )
            current_positions = []
        current_code = dominant_code
        current_positions.append(position)

    close_recovery_region(
        raw_regions,
        current_positions,
        current_code,
        distances_by_position,
        min_run_length,
    )
    return [
        region
        for region in raw_regions
        if region.average_distance >= average_distance_threshold
    ]


def calculate_recovery_regions_by_model(
    cmd,
    dominant_model: ExtractedModel,
    models_by_name: Dict[str, ExtractedModel],
    assignments: Dict[str, List[Residue]],
    comparisons_by_model: Dict[str, ModelComparison],
    rmsd_threshold: float,
    min_run_length: int,
) -> Dict[str, List]:
    dominant_residues = assignments.get(dominant_model.model_name)
    if dominant_residues is None:
        return {}
    output: Dict[str, List] = {}
    cmd.reinitialize()
    for object_index, (model_name, comparison) in enumerate(
        comparisons_by_model.items(), start=1
    ):
        if comparison.switch_regions:
            continue
        alternative_model = models_by_name.get(model_name)
        alternative_residues = assignments.get(model_name)
        if alternative_model is None or alternative_residues is None:
            continue
        distances = ca_distances_by_position_after_alignment(
            cmd,
            dominant_model,
            alternative_model,
            dominant_residues,
            alternative_residues,
            object_index,
        )
        output[model_name] = same_dssp_recovery_regions(
            dominant_residues,
            alternative_residues,
            distances,
            min_run_length,
            rmsd_threshold,
        )
    return output


def recovery_region_text(regions: Sequence) -> str:
    if not regions:
        return "-"
    return ",".join(
        f"{region.start}-{region.end}:{region.dssp_code}:{region.average_distance:.3f}"
        for region in regions
    )

def recovery_summary_fields_for_selected(
    cmd,
    selected_comparison: ModelComparison,
    dominant_model: ExtractedModel,
    alternative_model: ExtractedModel,
    dominant_residues: Sequence[Residue],
    alternative_residues: Sequence[Residue],
    rmsd_threshold: float,
    min_run_length: int,
) -> dict:
    fields: dict = {}
    cmd.reinitialize()
    distances = ca_distances_by_position_after_alignment(
        cmd,
        dominant_model,
        alternative_model,
        dominant_residues,
        alternative_residues,
        1,
    )
    regions = same_dssp_recovery_regions(
        dominant_residues,
        alternative_residues,
        distances,
        min_run_length,
        rmsd_threshold,
    )
    if regions:
        fields.update(
            {
                "same_ss_recovery_region_count": len(regions),
                "same_ss_recovery_regions": recovery_region_text(regions),
                "same_ss_recovery_max_avg_distance": f"{max(region.average_distance for region in regions):.3f}",
            }
        )
    return fields
