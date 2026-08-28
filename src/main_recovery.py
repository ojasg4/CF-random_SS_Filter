#!/usr/bin/env python3
"""
Reads unlikely_predictions.csv, runs four ordered
detectors per protein, writes recovered.csv and recovery_failed.csv.

Detector order (short-circuits on first hit):
  1. beta_bridge_reshuffle  (optional, --no-bridge-reshuffle to skip)
  2. disorder_to_beta
  3. disorder_to_helix
  4. rigid_body_helix
"""

import argparse
import csv
import sys
import tempfile
from pathlib import Path
from typing import Dict, List

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from utils.config import MIN_CONSECUTIVE_CHANGED_RESIDUES
from utils.recovery_utils.recovery_pipeline import run_recovery_pipeline
from utils.recovery_utils.recovery_io import write_output_csvs

DOMINANT_LABEL_PATTERNS = ["Dominant"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fold-switch recovery screen."
    )
    parser.add_argument(
        "--summary-dir", type=Path, default=Path("fold_switch_dssp_summary")
    )
    parser.add_argument("--unlikely-csv", type=Path, default=None)
    parser.add_argument("--logs-root",    type=Path, default=Path("logs"))
    parser.add_argument("--work-root",    type=Path, default=Path("fold_switch_dssp_work"))
    parser.add_argument(
        "--dominant-label-patterns", nargs="+", default=DOMINANT_LABEL_PATTERNS
    )
    parser.add_argument(
        "--beta-bridges", action="store_true",
        help="Enable beta bridge reshuffling detection.",
    )
    parser.add_argument(
        "--strict", action="store_true",
        help="Enables stricter parameters for recovery centered around rejecting similar secondary structures. Might exclude rigid body translations."
    )
    parser.add_argument(
        "--min-run-length", type=int, default=MIN_CONSECUTIVE_CHANGED_RESIDUES
    )
    parser.add_argument("--verbose",   action="store_true")
    parser.add_argument("--chainsaw", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        summary_dir      = args.summary_dir
        summary_base_dir = summary_dir.resolve().parent
        unlikely_csv     = args.unlikely_csv or (summary_dir / "unlikely_predictions.csv")

        if not unlikely_csv.exists():
            raise ValueError(f"No unlikely predictions CSV found at {unlikely_csv}")
        rows = [dict(r) for r in csv.DictReader(unlikely_csv.open())]
        if not rows:
            raise ValueError(f"No rows found in {unlikely_csv}")

        recovered_rows: List[Dict] = []
        failed_rows:    List[Dict] = []

        for idx, row in enumerate(rows, start=1):
            if args.verbose:
                print(f"[{idx}/{len(rows)}] {row.get('input_pse', '')}")

            input_pse = summary_base_dir / row.get("input_pse", "")
            raw = row.get("raw_dssp_dir", "").strip()
            work_dir = summary_base_dir / raw if raw else args.work_root / input_pse.stem
            job_name = work_dir.name
            raw_dssp_dir = args.logs_root / job_name / "raw_dssp_runs"

            work_dir = args.work_root / job_name
            work_dir.mkdir(parents=True, exist_ok=True)
            out_row, recovered = run_recovery_pipeline(
                args, job_name, input_pse,
                work_dir, raw_dssp_dir, summary_base_dir,
            )


            if recovered:
                recovered_rows.append(out_row)
                if args.verbose:
                    print(f"  → recovered via {out_row.get('recovery_detector', '?')}")
            else:
                failed_rows.append(out_row)
                if args.verbose:
                    print("  → failed")

        write_output_csvs(summary_dir, recovered_rows, failed_rows)

        if args.verbose:
            print(f"\nRecovered: {len(recovered_rows)}  Failed: {len(failed_rows)}")
        return 0

    except Exception as exc:
        if args.verbose:
            print(f"ERROR: {exc}", file=sys.stderr)
            raise
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
