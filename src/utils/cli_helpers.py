import re
from pathlib import Path
from typing import List


def derive_job_name(pse_file: Path) -> str:
    """
    Derive a job name from a PSE path.
    Strips the -structures_of_interest suffix, uses the stem if informative,
    otherwise walks the path from the top for the first WP_ component so that
    sibling variants (e.g. WP_000617148 vs WP_000617148_confirm) stay distinct.
    """
    stem = pse_file.stem
    for suffix in ("-structures_of_interest", "_structures_of_interest"):
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    if stem and stem != "structures_of_interest":
        return stem
    for part in pse_file.resolve().parts:
        if re.match(r"^WP_\d+", part):
            return part
    return pse_file.parent.name


def collect_pse_files(inputs: List[Path], pse_glob: str) -> List[Path]:
    """
    Expand each input into concrete PSE files.
      - A .pse file is used directly.
      - A directory is searched recursively for pse_glob.
    """
    pse_files: List[Path] = []
    for p in inputs:
        if p.is_dir():
            pse_files.extend(sorted(p.rglob(pse_glob)))
        elif p.suffix == ".pse":
            pse_files.append(p)
    return pse_files