import json
from pathlib import Path
from typing import Optional
from typing import Dict


def load_config(config_path: Optional[Path] = None) -> dict:
    if config_path is None:
        config_path = Path(__file__).parent.parent.parent / "config" / "config.json"
    with open(config_path) as f:
        return json.load(f)


CFG = load_config()

MIN_CONSECUTIVE_CHANGED_RESIDUES       = CFG["thresholds"]["min_consecutive_changed_residues"]
MIN_REGION_SEPARATION_TO_KEEP_SEPARATE = CFG["thresholds"]["min_region_separation_to_keep_separate"]
RMSD_THRESHOLD                         = CFG["thresholds"]["rmsd_threshold"]
PHENIX_MODE_COLUMNS = CFG["csv_columns"]["phenix_mode_columns"]

HELIX_DSSP_CODES      = set(CFG["dssp_codes"]["helix"])
BETA_DSSP_CODES       = set(CFG["dssp_codes"]["beta"])
KAPPA_DSSP_CODES      = set(CFG["dssp_codes"]["kappa"])
DISORDERED_DSSP_CODES = set(CFG["dssp_codes"]["disordered"])

GROUP_HELIX      = CFG["groups"]["helix"]
GROUP_BETA       = CFG["groups"]["beta"]
GROUP_KAPPA      = CFG["groups"]["kappa"]
GROUP_DISORDERED = CFG["groups"]["disordered"]

GROUP_SINGLE_LETTER  = CFG["group_single_letter"]
NONDISORDERED_GROUPS = {GROUP_HELIX, GROUP_BETA, GROUP_KAPPA}

SUMMARY_FILES               = CFG["summary_files"]
SUMMARY_BASE_PREFIX_COLUMNS = CFG["csv_columns"]["base_prefix"]
SUMMARY_BASE_SUFFIX_COLUMNS = CFG["csv_columns"]["base_suffix"]
KEEP_GOING_COLUMNS          = CFG["csv_columns"]["keep_going"]
RECOVERY_METRIC_COLUMNS     = CFG["csv_columns"]["recovery_metrics"]
UNASSIGNED_MAX_ALTERNATIVE_PERCENT = CFG["unassigned_filter"]["max_alternative_percent"]

PHENIX_MODE_COLUMNS   = CFG["csv_columns"]["phenix_mode_columns"]
UNASSIGNED_METRIC_COLUMNS = CFG["csv_columns"]["unassigned_metrics"]

CHAINSAW_COLUMNS = CFG["csv_columns"].get("chainsaw_columns", [])