import csv
import fcntl
import shutil
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from utils.config import (
    SUMMARY_FILES,
    SUMMARY_BASE_PREFIX_COLUMNS,
    SUMMARY_BASE_SUFFIX_COLUMNS,
    KEEP_GOING_COLUMNS,
    RECOVERY_METRIC_COLUMNS,
    PHENIX_MODE_COLUMNS,
    UNASSIGNED_METRIC_COLUMNS,
    CHAINSAW_COLUMNS,
)
from utils.classes import ComparedSegment, ModelComparison, Residue, segment_annotation
from utils.dssp import (
    get_sequence,
    get_dssp,
    unassigned_metrics_for_residues
)
from utils.pdb_helpers import (
    pymol_session,
    pymol_residue_selection,
    recovery_region_text,
    recovery_summary_fields_for_selected,
    calculate_recovery_regions_by_model,
)
from utils.utils import (
    compact_reason,
    relative_path_for_summary,
    short_model_name,
    wrap_fasta,
)

from utils.phenix_filter import phenix_fields_for_model, ALL_MODES


# ---------------------------------------------------------------------------
# Colored output PSE for the selected pair
# ---------------------------------------------------------------------------

def log_input_too_disordered(job_name: str, summary_dir: Path, percent: float) -> None:
    summary_dir.mkdir(parents=True, exist_ok=True)
    log_path = summary_dir / "too_disordered.txt"
    log_path.touch(exist_ok=True)
    print(f"{job_name} was too disordered, not analyzed\t{percent:.2f}\n")
    with log_path.open("a") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        try:
            fh.write(f"{job_name}\t{percent:.2f}\n")
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)


def create_selected_colored_pse(
    cmd,
    dominant_model,
    alternative_model,
    dominant_residues: Sequence[Residue],
    alternative_residues: Sequence[Residue],
    selected_comparison: ModelComparison,
    output_pse: Path,
) -> None:
    output_pse.parent.mkdir(parents=True, exist_ok=True)
    cmd.reinitialize()
    cmd.load(str(dominant_model.pdb_path), "dominant")
    cmd.load(str(alternative_model.pdb_path), "alternative")
    try:
        cmd.cealign("dominant and name CA", "alternative and name CA", quiet=1)
    except Exception:
        cmd.cealign("dominant", "alternative", quiet=1)
    cmd.bg_color("white")
    cmd.hide("everything")
    cmd.show("cartoon", "dominant or alternative")
    cmd.color("gray70", "dominant or alternative")

    switch_positions = []
    for segment in selected_comparison.switch_regions:
        switch_positions.extend(range(segment.start, segment.end + 1))

    disordered_positions = []
    for segment in selected_comparison.disordered_regions:
        disordered_positions.extend(range(segment.start, segment.end + 1))

    disordered_selection_dominant = pymol_residue_selection(
        "dominant", dominant_residues, disordered_positions
    )
    disordered_selection_alternative = pymol_residue_selection(
        "alternative", alternative_residues, disordered_positions
    )
    switch_selection_dominant = pymol_residue_selection(
        "dominant", dominant_residues, switch_positions
    )
    switch_selection_alternative = pymol_residue_selection(
        "alternative", alternative_residues, switch_positions
    )

    if disordered_positions:
        cmd.color(
            "orange",
            f"({disordered_selection_dominant}) or ({disordered_selection_alternative})",
        )
    if switch_positions:
        cmd.color(
            "red",
            f"({switch_selection_dominant}) or ({switch_selection_alternative})",
        )

    cmd.set("ray_opaque_background", 0)
    cmd.orient("dominant or alternative")
    cmd.save(str(output_pse))


# ---------------------------------------------------------------------------
# Saving the selected pair as a fold package
# ---------------------------------------------------------------------------

def save_fold_package(
    fold_root: Path,
    job_name: str,
    dominant_sequence: str,
    alternative_sequence: str,
    dominant_model,
    alternative_model,
    selected_comparison: ModelComparison,
    rmsd_threshold: float,
) -> Path:
    fold_dir = fold_root / job_name
    fold_dir.mkdir(parents=True, exist_ok=True)
    with (fold_dir / "conf1.fasta").open("w") as handle:
        handle.write(
            f">{job_name}_conf1_{dominant_model.model_name}\n{wrap_fasta(dominant_sequence)}\n"
        )
    with (fold_dir / "conf2.fasta").open("w") as handle:
        handle.write(
            f">{job_name}_conf2_{alternative_model.model_name}\n{wrap_fasta(alternative_sequence)}\n"
        )
    if dominant_sequence == alternative_sequence:
        with (fold_dir / "sequence.fasta").open("w") as handle:
            handle.write(f">{job_name}\n{wrap_fasta(dominant_sequence)}\n")
    shutil.copyfile(dominant_model.pdb_path, fold_dir / "conf1.pdb")
    shutil.copyfile(alternative_model.pdb_path, fold_dir / "conf2.pdb")
    switch_range = ";".join(
        f"{segment.start}-{segment.end}"
        for segment in selected_comparison.switch_regions
    )
    with (fold_dir / "metadata.tsv").open("w") as handle:
        handle.write("category\t" + selected_comparison.classification + "\n")
        handle.write("switch_range\t" + switch_range + "\n")
        handle.write("alternative_model\t" + alternative_model.model_name + "\n")
        handle.write(
            "rmsd_to_dominant\t"
            + (
                f"{selected_comparison.rmsd_to_dominant:.3f}"
                if selected_comparison.rmsd_to_dominant is not None
                else "N/A"
            )
            + "\n"
        )
        handle.write("rmsd_threshold\t" + f"{rmsd_threshold:.3f}" + "\n")
        handle.write(
            "sequences_identical\t"
            + str(dominant_sequence == alternative_sequence).lower()
            + "\n"
        )
        handle.write("invalid_reason\t" + selected_comparison.invalid_reason + "\n")
    return fold_dir


# ---------------------------------------------------------------------------
# Per-run debug dump: final_comparison.txt
# ---------------------------------------------------------------------------

def compact_comparison_annotation(
    comparison: Optional[ModelComparison],
    model_name: str,
    base_dir: Path,
) -> str:
    if comparison is None:
        return f"model={model_name} role=dominant"

    annotations = ",".join(comparison.annotations) if comparison.annotations else "-"
    rmsd = (
        f"{comparison.rmsd_to_dominant:.3f}"
        if comparison.rmsd_to_dominant is not None
        else "NA"
    )
    reason = comparison.invalid_reason or "-"

    return (
        f"model={model_name} "
        f"class={comparison.classification} "
        f"rmsd={rmsd} "
        f"sw={sum(segment.end - segment.start + 1 for segment in comparison.switch_regions)}/{len(comparison.switch_regions)} "
        f"dischg={comparison.n_disordered_change_residues}/{len(comparison.disordered_regions)} "
        f"ann={annotations} "
        f"reason={reason}"
    )


def write_final_comparison_file(
    cmd,
    output_path: Path,
    job_name: str,
    run_mode: str,
    input_pse: Path,
    dominant_model_name: str,
    assignments: Dict[str, List[Residue]],
    assignment_order: Sequence[str],
    comparisons_by_model: Dict[str, ModelComparison],
    models_by_name: dict,
    selected_comparison: ModelComparison,
    rmsd_threshold: float,
    min_run_length: int,
    base_dir: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ordered_names: List[str] = []
    if dominant_model_name in assignments:
        ordered_names.append(dominant_model_name)
    for model_name in assignment_order:
        if model_name in assignments and model_name not in ordered_names:
            ordered_names.append(model_name)

    with output_path.open("w") as handle:
        handle.write(f"# job={job_name}\n")
        handle.write(f"# mode={run_mode}\n")
        handle.write(f"# input={relative_path_for_summary(input_pse, base_dir)}\n")
        handle.write(f"# dominant={dominant_model_name}\n")
        handle.write(f"# selected={selected_comparison.alternative_model}\n")
        handle.write(f"# selected_class={selected_comparison.classification}\n")
        handle.write(
            f"# rmsd_threshold={rmsd_threshold:.3f} min_run_length={min_run_length}\n"
        )
        handle.write("# disorder=percent/n_regions/avg_region_length\n")
        handle.write("# same_ss_preserved=percent_of_dominant_dssp_codes_preserved\n")
        handle.write("# recovery=same_dssp_region:avg_ca_distance_after_alignment\n")
        if not ordered_names:
            handle.write("# no_successful_dssp_assignments\n")
            return

        recovery_by_model = calculate_recovery_regions_by_model(
            cmd,
            models_by_name[dominant_model_name],
            models_by_name,
            assignments,
            comparisons_by_model,
            rmsd_threshold,
            min_run_length,
        )
        shared_sequence = get_sequence(assignments[ordered_names[0]])
        handle.write(shared_sequence + "\n")
        for model_name in ordered_names:
            residues = assignments[model_name]
            structure = get_dssp(residues)
            comparison = (
                None
                if model_name == dominant_model_name
                else comparisons_by_model.get(model_name)
            )
            annotation = compact_comparison_annotation(comparison, model_name, base_dir)
            recovery = (
                "-"
                if model_name == dominant_model_name
                else recovery_region_text(recovery_by_model.get(model_name, []))
            )
            dominant_residues = assignments[dominant_model_name]
            same_ss_preserved = (
                100.0
                if model_name == dominant_model_name
                else comparisons_by_model[model_name].same_ss_preserved_percent or 0.0
                if model_name in comparisons_by_model
                else 0.0
            )
            handle.write(
                f"{structure}\t{annotation} same_ss_preserved={same_ss_preserved:.1f}% rec={recovery}\n"
            )


# ---------------------------------------------------------------------------
# Summary CSV utils
# ---------------------------------------------------------------------------

def read_summary_rows(path: Path) -> List[Dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return []
        rows = []
        for row in reader:
            if not row:
                continue
            row = dict(row)
            row.pop("n_invalid_candidate_structures", None)
            rows.append(row)
        return rows


def write_summary_rows(
    path: Path,
    rows: Sequence[Dict[str, object]],
    fieldnames: Sequence[str],
) -> None:
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    with tmp_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in fieldnames})
    tmp_path.replace(path)


def summary_fieldnames_for_rows(rows: Sequence[Dict[str, object]]) -> List[str]:
    max_change_idx = 2
    extras = set()
    reserved = set(SUMMARY_BASE_PREFIX_COLUMNS)
    reserved.update(SUMMARY_BASE_SUFFIX_COLUMNS)
    reserved.update(KEEP_GOING_COLUMNS)
    reserved.update(RECOVERY_METRIC_COLUMNS)
    reserved.update(PHENIX_MODE_COLUMNS)
    reserved.update(UNASSIGNED_METRIC_COLUMNS)
    reserved.update(CHAINSAW_COLUMNS)
    for row in rows:
        for key in row.keys():
            if key.startswith("change_"):
                try:
                    max_change_idx = max(max_change_idx, int(key.split("_", 1)[1]))
                except ValueError:
                    extras.add(key)
            elif key not in reserved:
                extras.add(key)
    fieldnames = list(SUMMARY_BASE_PREFIX_COLUMNS)
    fieldnames.extend(CHAINSAW_COLUMNS)
    fieldnames.extend(f"change_{idx}" for idx in range(1, max_change_idx + 1))
    if any(any(column in row for column in KEEP_GOING_COLUMNS) for row in rows):
        fieldnames.extend(KEEP_GOING_COLUMNS)
    fieldnames.extend(RECOVERY_METRIC_COLUMNS)
    fieldnames.extend(PHENIX_MODE_COLUMNS)
    fieldnames.extend(UNASSIGNED_METRIC_COLUMNS)
    fieldnames.extend(SUMMARY_BASE_SUFFIX_COLUMNS)
    fieldnames.extend(sorted(extras.difference(fieldnames)))
    return fieldnames
    
def append_or_update_summary(
    summary_dir: Path,
    selected_comparison: ModelComparison,
    output_comparison: ModelComparison, # Output comparison is used to display the metrics from the best possible failure. Automatically
    # set to either the success or the best possible failure when called by assemble_outputs called by pipelines.py
    # Original output of selected_comparison still needs to be passed for the failure reason if applicable
    # assigned in assemble_outputs
    input_pse: Path,
    run_mode: str,
    job_name: str,
    summary_base_dir: Path,
    selected_colored_pse: Optional[Path],
    raw_dssp_dir: Path,
    saved_folds_dir: Optional[Path],
    dominant_model_name: str,
    rmsd_threshold: float,
    extra_summary_fields: Optional[Dict[str, object]] = None,
) -> None:
    summary_dir.mkdir(parents=True, exist_ok=True)
    lock_path = summary_dir / ".summary.lock"
    classification = selected_comparison.classification
    tier_file = SUMMARY_FILES[classification]
    change_columns = {
        f"change_{idx}": value
        for idx, value in enumerate(output_comparison.annotations, start=1)
    }
    row = {
        "input_pse": relative_path_for_summary(input_pse, summary_base_dir),
        "alt_model": output_comparison.alternative_model,
        "classification": classification,
        "rmsd_to_dominant": (
            f"{output_comparison.rmsd_to_dominant:.3f}"
            if output_comparison.rmsd_to_dominant is not None
            else "N/A"
        ),
        "n_switch_residues": sum(
            segment.end - segment.start + 1
            for segment in output_comparison.switch_regions
        ),
        "n_switch_regions": len(output_comparison.switch_regions),
        "n_disordered_change_residues": output_comparison.n_disordered_change_residues,
        "n_disordered_change_regions": len(output_comparison.disordered_regions),
        "invalid_reason": selected_comparison.invalid_reason, # This is the only one that needs to come from original
        # because the reason should be failed if that is accurate but everything else should load from the best possible model
        "run_mode": run_mode,
        "sequences_identical": (
            ""
            if output_comparison.sequences_identical is None
            else str(output_comparison.sequences_identical).lower()
        ),
        "selected_colored_pse": relative_path_for_summary(
            selected_colored_pse, summary_base_dir
        ),
        "raw_dssp_dir": relative_path_for_summary(raw_dssp_dir, summary_base_dir),
        "saved_folds_dir": (
            relative_path_for_summary(saved_folds_dir, summary_base_dir)
            if saved_folds_dir
            else ""
        ),
        **change_columns,
    }
    if extra_summary_fields:
        row.update(extra_summary_fields)
    with lock_path.open("w") as lock_handle:
        fcntl.flock(lock_handle, fcntl.LOCK_EX)
        try:
            rows_by_file: Dict[str, List[Dict[str, object]]] = {}
            for filename in SUMMARY_FILES.values():
                path = summary_dir / filename
                rows = read_summary_rows(path)
                rows = [
                    existing
                    for existing in rows
                    if existing.get("job_name") != job_name
                ]
                if filename == tier_file:
                    rows.append(row)
                rows_by_file[filename] = rows
            for filename, rows in rows_by_file.items():
                fieldnames = summary_fieldnames_for_rows(rows)
                write_summary_rows(summary_dir / filename, rows, fieldnames)
        finally:
            fcntl.flock(lock_handle, fcntl.LOCK_UN)


# ---------------------------------------------------------------------------
# Shared output assembly (used by both pipelines)
# ---------------------------------------------------------------------------

# main function called by pipelines.py
def assemble_outputs(
    args,
    job_name: str,
    run_mode: str,
    work_dir: Path,
    protein_log_dir: Path,
    raw_dssp_dir: Path,
    dominant_model,
    models_by_name: dict,
    all_assignments: Dict[str, List[Residue]],
    assignment_order: List[str],
    comparisons_by_model: Dict[str, ModelComparison],
    selected_comparison: ModelComparison,
    extra_summary_fields: Dict[str, object],
    cmd,
    best_failed_comparison: Optional[ModelComparison] = None,
) -> None:
    # For output purposes, use the best failed comparison when no hit was found
    output_comparison = (
        best_failed_comparison
        if selected_comparison.alternative_model == "none" and best_failed_comparison is not None
        else selected_comparison
    )

    alternative_model = models_by_name.get(output_comparison.alternative_model)
    dominant_residues = all_assignments.get(dominant_model.model_name)
    alternative_residues = (
        all_assignments.get(alternative_model.model_name)
        if alternative_model is not None
        else None
    )
    selected_pse_path = args.summary_dir / f"{job_name}_selected_colored.pse"
    final_comparison_path = work_dir / "final_comparison.txt"

    selected_colored_pse: Optional[Path] = None
    saved_folds_dir = None

    have_both_residue_sets = (
        dominant_residues is not None and alternative_residues is not None
    )
    dominant_residues = all_assignments.get(dominant_model.model_name)
    if dominant_residues is not None:
        extra_summary_fields.update(unassigned_metrics_for_residues("dom_poor_structure", dominant_residues))
    if alternative_residues is not None:
        extra_summary_fields.update(unassigned_metrics_for_residues("alt_poor_structure", alternative_residues))
    extra_summary_fields.update(
        phenix_fields_for_model(
            "dom",
            dominant_model.pdb_path,
            dominant_model.dssp_path.with_suffix(".phenix"),
        )
    )

    if alternative_model is not None:
        extra_summary_fields.update(
            phenix_fields_for_model(
                "alt",
                alternative_model.pdb_path,
                alternative_model.dssp_path.with_suffix(".phenix"),
            )
        )
    if output_comparison.same_ss_preserved_percent is not None:
        extra_summary_fields["same_ss_preserved"] = f"{output_comparison.same_ss_preserved_percent:.2f}"


    if alternative_model is not None and have_both_residue_sets:
        extra_summary_fields.update(
            recovery_summary_fields_for_selected(
                cmd,
                output_comparison,
                dominant_model,
                alternative_model,
                dominant_residues,
                alternative_residues,
                args.rmsd_threshold,
                args.min_run_length,
            )
        )

        dominant_sequence = get_sequence(dominant_residues)
        alternative_sequence = get_sequence(alternative_residues)
        selected_colored_pse = selected_pse_path
        create_selected_colored_pse(
            cmd,
            dominant_model=dominant_model,
            alternative_model=alternative_model,
            dominant_residues=dominant_residues,
            alternative_residues=alternative_residues,
            selected_comparison=selected_comparison,
            output_pse=selected_colored_pse,
        )

        if args.save_folds:
            saved_folds_dir = save_fold_package(
                fold_root=args.folds_root,
                job_name=job_name,
                dominant_sequence=dominant_sequence,
                alternative_sequence=alternative_sequence,
                dominant_model=dominant_model,
                alternative_model=alternative_model,
                selected_comparison=selected_comparison,
                rmsd_threshold=args.rmsd_threshold,
            )

    if args.keep_work:
        write_final_comparison_file(
            cmd,
            output_path=final_comparison_path,
            job_name=job_name,
            run_mode=run_mode,
            input_pse=args.pse_file,
            dominant_model_name=dominant_model.model_name,
            assignments=all_assignments,
            assignment_order=assignment_order,
            comparisons_by_model=comparisons_by_model,
            models_by_name=models_by_name,
            selected_comparison=output_comparison,
            rmsd_threshold=args.rmsd_threshold,
            min_run_length=args.min_run_length,
            base_dir=args.pse_file.resolve().parent,
        )

    append_or_update_summary(
        summary_dir=args.summary_dir,
        selected_comparison=selected_comparison,
        output_comparison=output_comparison,
        input_pse=args.pse_file,
        run_mode=run_mode,
        job_name=job_name,
        summary_base_dir=args.summary_dir.resolve().parent,
        selected_colored_pse=selected_colored_pse,
        raw_dssp_dir=raw_dssp_dir,
        saved_folds_dir=saved_folds_dir,
        dominant_model_name=dominant_model.model_name,
        rmsd_threshold=args.rmsd_threshold,
        extra_summary_fields=extra_summary_fields,
    )

    if args.verbose:
        _print_verbose_summary(
            job_name, run_mode, dominant_model, selected_comparison, output_comparison,
            args, extra_summary_fields, raw_dssp_dir,
            selected_colored_pse, saved_folds_dir,
            work_dir if args.keep_work else None,
            final_comparison_path if args.keep_work else None,
        )


def _print_verbose_summary(
    job_name: str,
    run_mode: str,
    dominant_model,
    selected_comparison: ModelComparison,
    output_comparison: ModelComparison,
    args,
    extra_summary_fields: Dict[str, object],
    raw_dssp_dir: Path,
    selected_colored_pse: Optional[Path],
    saved_folds_dir: Optional[Path],
    work_dir: Optional[Path],
    final_comparison_path: Optional[Path],
) -> None:
    print(f"Job name: {job_name}")
    print(f"Run mode: {run_mode}")
    print(f"Dominant model: {dominant_model.model_name}")
    print(f"Alternative model: {output_comparison.alternative_model}")
    print(
        f"RMSD to dominant: {output_comparison.rmsd_to_dominant if output_comparison.rmsd_to_dominant is not None else 'NA'}"
    )
    print(f"RMSD threshold: {args.rmsd_threshold}")
    print(f"Classification: {selected_comparison.classification}")
    print(f"Switch residues: {sum(segment.end - segment.start + 1 for segment in output_comparison.switch_regions)}")
    print(f"Switch regions: {len(output_comparison.switch_regions)}")
    print(f"Disordered-change residues: {output_comparison.n_disordered_change_residues}")
    print(f"Disordered-change regions: {len(output_comparison.disordered_regions)}")

    print("Dominant phenix modes:")
    for mode in ALL_MODES:
        key = f"dom_phenix_{mode}"
        print(f"  {mode}: {extra_summary_fields.get(key, 'unavailable')}")

    print("Alternative phenix modes:")
    for mode in ALL_MODES:
        key = f"alt_phenix_{mode}"
        print(f"  {mode}: {extra_summary_fields.get(key, 'unavailable')}")
    print(f"Bad structures skipped: {extra_summary_fields.get('n_skipped', 0)}")

    if selected_comparison.invalid_reason:
        print(f"Invalid reason: {selected_comparison.invalid_reason}")
    if extra_summary_fields:
        for key, value in extra_summary_fields.items():
            if "_phenix_" in key:
                continue # Already printed before in a neatly formatted manner, but other fields in extra_summary exist
            print(f"{key}: {value}")
    print(f"Raw DSSP directory: {raw_dssp_dir}")
    print(
        f"Selected colored PSE: {selected_colored_pse if selected_colored_pse else 'not written'}"
    )
    print(f"Summary directory: {args.summary_dir}")
    if saved_folds_dir:
        print(f"Saved folds directory: {saved_folds_dir}")
    if work_dir:
        print(f"Kept work directory: {work_dir}")
    if final_comparison_path:
        print(f"Final comparison file: {final_comparison_path}")
    print("--------------------------------------------------")
