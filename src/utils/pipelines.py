from dataclasses import replace
from pathlib import Path
from typing import Dict, List

from utils.classes import ModelComparison, Residue
from utils.structure import (
    extract_models_from_pse,
    find_dominant_model,
    find_keep_going_pdb_files,
    extracted_models_from_keep_going_pdb_files,
    analyze_alternatives_with_dssp,
)
from utils.comparison import (
    first_hit,
    is_genuine_hit,
    no_hit_comparison,
    record_assignments,
    keep_going_summary_fields,
)
from utils.outputs import assemble_outputs

from utils.phenix_filter import get_mode_fractions
from utils.phenix_filter import PHENIX_FILTERED_MODES, PHENIX_THRESHOLDS

def run_pipeline_default(
    args,
    job_name: str,
    protein_log_dir: Path,
    work_dir: Path,
) -> None:
    raw_dssp_dir = protein_log_dir / "raw_dssp_runs"
    raw_dssp_dir.mkdir(parents=True, exist_ok=True)

    #  1. Extract structures from PSE
    extracted_dir = work_dir
    models = extract_models_from_pse(args.pse_file, extracted_dir, raw_dssp_dir)
    dominant_model = find_dominant_model(models, args.dominant_label_patterns)
    candidate_models = [m for m in models if m.model_name != dominant_model.model_name]
    print(f"Models: {', '.join(m.model_name for m in models)}, Alternatives: {', '.join(m.model_name for m in candidate_models)}")

    if not candidate_models:
        raise ValueError("No alternative models were available for comparison")

    models_by_name = {model.model_name: model for model in models}

    #  1b. Disorder filter: check barbed_wire_fractions
    n_phenix_skipped = 0
    dominant_fractions = get_mode_fractions(dominant_model.pdb_path)
    dominant_summary = (
                ", ".join(f"{mode}={dominant_fractions.get(mode, 0.0):.3f}" for mode in PHENIX_FILTERED_MODES)
            )
    alternative_models = []
    for m in candidate_models:
        fractions = get_mode_fractions(m.pdb_path)
        if fractions is not None and any(
            fractions.get(mode, 0.0) > PHENIX_THRESHOLDS.get(mode, 1.0)
            for mode in PHENIX_FILTERED_MODES
        ):
            alt_summary = ", ".join(f"{mode}={fractions.get(mode, 0.0):.3f}" for mode in PHENIX_FILTERED_MODES)
            print(f"Alternative model {m.model_name} was skipped: ({alt_summary}). Dominant was: {dominant_summary}.")
            n_phenix_skipped += 1
        else:
            alternative_models.append(m)

    #  2. RMSD gate + DSSP comparison
    comparisons, assignments = analyze_alternatives_with_dssp(
        dominant_model,
        alternative_models,
        args.min_run_length,
        args.rmsd_threshold,
    )

    all_assignments: Dict[str, List[Residue]] = {}
    assignment_order: List[str] = []
    record_assignments(all_assignments, assignment_order, assignments)

    comparisons_by_model = {c.alternative_model: c for c in comparisons}

    #  3. Select best hit 
    selected_comparison = first_hit(comparisons, assignments, dominant_model.model_name)
    best_failed = None
    if selected_comparison is None:
        all_comparisons = list(comparisons_by_model.values())
        best_failed = max(
            (c for c in all_comparisons if c.rmsd_to_dominant is not None),
            key=lambda c: c.rmsd_to_dominant,
            default=None,
        )
        selected_comparison = no_hit_comparison(len(alternative_models))


    #  4. Write outputs 
    extra_fields = {"n_phenix_skipped": n_phenix_skipped}
    assemble_outputs(
        args, job_name, "default", work_dir, protein_log_dir, raw_dssp_dir,
        dominant_model, models_by_name, all_assignments, assignment_order,
        comparisons_by_model, selected_comparison, extra_fields,
        best_failed_comparison=best_failed,
    )



def run_pipeline_keep_going(
    args,
    job_name: str,
    protein_log_dir: Path,
    work_dir: Path,
) -> None:
    raw_dssp_dir = protein_log_dir / "raw_dssp_runs"
    raw_dssp_dir.mkdir(parents=True, exist_ok=True)

    #  1. Extract structures from PSE 
    extracted_dir = work_dir
    models = extract_models_from_pse(args.pse_file, extracted_dir, raw_dssp_dir)
    dominant_model = find_dominant_model(models, args.dominant_label_patterns)
    candidate_models = [m for m in models if m.model_name != dominant_model.model_name]

    models_by_name = {model.model_name: model for model in models}

    #  1b. Disorder filter: check barbed_wire_fractions
    n_phenix_skipped = 0
    dominant_fractions = get_mode_fractions(dominant_model.pdb_path)
    dominant_summary = (
                ", ".join(f"{mode}={dominant_fractions.get(mode, 0.0):.3f}" for mode in PHENIX_FILTERED_MODES)
            )
    main_alternative_models = []
    for m in candidate_models:
        fractions = get_mode_fractions(m.pdb_path)
        if fractions is not None and any(
            fractions.get(mode, 0.0) > PHENIX_THRESHOLDS.get(mode, 1.0)
            for mode in PHENIX_FILTERED_MODES
        ):
            alt_summary = ", ".join(f"{mode}={fractions.get(mode, 0.0):.3f}" for mode in PHENIX_FILTERED_MODES)
            print(f"Alternative model {m.model_name} was skipped: ({alt_summary}). Dominant was: {dominant_summary}.")
            n_phenix_skipped += 1
        else:
            main_alternative_models.append(m)

    all_assignments: Dict[str, List[Residue]] = {}
    assignment_order: List[str] = []
    comparisons_by_model: Dict[str, ModelComparison] = {}
    n_searched = 0

    #  2. Search PSE clusters first 
    if main_alternative_models:
        main_comparisons, main_assignments = analyze_alternatives_with_dssp(
            dominant_model,
            main_alternative_models,
            args.min_run_length,
            args.rmsd_threshold,
        )
        record_assignments(all_assignments, assignment_order, main_assignments)
        comparisons_by_model.update({c.alternative_model: c for c in main_comparisons})
        n_searched += len(main_alternative_models)
        selected_comparison = first_hit(main_comparisons, main_assignments, dominant_model.model_name)
        # Can also return none if secondary filters are implemented in comparison.py
    else:
        selected_comparison = None

    #  3. If no hit yet, walk the raw PDB directories one-by-one 
    if selected_comparison is None:
        keep_going_models = extracted_models_from_keep_going_pdb_files(
            find_keep_going_pdb_files(args.pse_file),
            raw_dssp_dir,
            args.pse_file.resolve().parent,
        )
        models.extend(keep_going_models)
        models_by_name.update({m.model_name: m for m in keep_going_models})

        encounter_index = len(main_alternative_models) if main_alternative_models else 0
        for alternative_model in keep_going_models:
            fractions = get_mode_fractions(alternative_model.pdb_path)
            if fractions is not None and any(
                fractions.get(mode, 0.0) > PHENIX_THRESHOLDS.get(mode, 1.0)
                for mode in PHENIX_FILTERED_MODES
            ):
                alt_summary = ", ".join(f"{mode}={fractions.get(mode, 0.0):.3f}" for mode in PHENIX_FILTERED_MODES)
                print(f"Alternative model {alternative_model.model_name} was skipped: ({alt_summary}). Dominant was: {dominant_summary}.")
                n_phenix_skipped += 1
                continue

            one_comparisons, one_assignments = analyze_alternatives_with_dssp(
                dominant_model,
                [alternative_model],
                args.min_run_length,
                args.rmsd_threshold,
            )
            n_searched += 1
            if one_assignments:
                record_assignments(all_assignments, assignment_order, one_assignments)
            if not one_comparisons:
                continue
            current_comparison = replace(one_comparisons[-1], encounter_index=encounter_index)
            encounter_index += 1
            comparisons_by_model[current_comparison.alternative_model] = current_comparison
            if is_genuine_hit(current_comparison, one_assignments, dominant_model.model_name):
                selected_comparison = current_comparison
                break

    best_failed = None
    if selected_comparison is None:
        all_comparisons = list(comparisons_by_model.values())
        best_failed = max(
            (c for c in all_comparisons if c.rmsd_to_dominant is not None),
            key=lambda c: c.rmsd_to_dominant,
            default=None,
        )
        selected_comparison = no_hit_comparison(len(main_alternative_models))


    #  4. Write outputs 
    extra_fields = keep_going_summary_fields([selected_comparison])
    extra_fields["n_phenix_skipped"] = n_phenix_skipped
    assemble_outputs(
        args, job_name, "keep_going", work_dir, protein_log_dir, raw_dssp_dir,
        dominant_model, models_by_name, all_assignments, assignment_order,
        comparisons_by_model, selected_comparison, extra_fields,
        best_failed_comparison=best_failed,
    )

