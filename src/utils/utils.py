import os
import re
from pathlib import Path
from typing import Optional


def compact_reason(text: str) -> str:
    text = str(text).replace(";", ",")
    text = re.sub(r"\s+", "_", text.strip())
    return text[:400]


def safe_path_component(text: str) -> str:
    text = str(text).strip()
    text = text.replace(os.sep, "_")
    if os.altsep:
        text = text.replace(os.altsep, "_")
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text)
    return text.strip("_") or "job"


def wrap_fasta(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def relative_path_for_summary(path: Optional[Path], base_dir: Path) -> str:
    if path is None:
        return ""
    path = Path(path)
    try:
        return str(path.resolve().relative_to(base_dir.resolve()))
    except Exception:
        try:
            return os.path.relpath(str(path.resolve()), str(base_dir.resolve()))
        except Exception:
            return str(path)


def short_model_name(model_name: str, base_dir: Path) -> str:
    text = str(model_name)
    path = Path(text)
    if not path.is_absolute():
        return text
    return relative_path_for_summary(path, base_dir)
