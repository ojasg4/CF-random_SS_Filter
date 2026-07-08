from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

from utils.config import GROUP_SINGLE_LETTER, GROUP_DISORDERED


@dataclass
class ExtractedModel:
    model_name: str
    pdb_path: Path
    dssp_path: Path


@dataclass
class Residue:
    amino_acid: str
    dssp_code: str
    group: str
    chain_id: str
    residue_number: str


@dataclass
class ComparedSegment:
    start: int
    end: int
    dominant_group: str
    alternative_group: str
    category: str

@dataclass
class OpenSegment:
    start: int
    dominant_group: str
    alternative_group: str
    category: str

@dataclass
class RecoveryRegion:
    start: int
    end: int
    dssp_code: str
    positions: List[int]
    average_distance: float


@dataclass
class ModelComparison:
    alternative_model: str
    encounter_index: int
    rmsd_to_dominant: Optional[float] = None
    switch_segments: Optional[List[ComparedSegment]] = None
    disordered_change_segments: Optional[List[ComparedSegment]] = None
    invalid_reason: str = ""
    sequences_identical: Optional[bool] = None
    same_ss_preserved_percent: Optional[float] = None
    barbed_wire_fraction: Optional[float] = None

    @property
    def switch_regions(self) -> List[ComparedSegment]:
        return self.switch_segments or []

    @property
    def disordered_regions(self) -> List[ComparedSegment]:
        return self.disordered_change_segments or []

    @property
    def n_disordered_change_residues(self) -> int:
        return sum(
            segment.end - segment.start + 1 for segment in self.disordered_regions
        )

    @property
    def classification(self) -> str:
        return "likely" if self.switch_regions else "unlikely"

    @property
    def annotations(self) -> List[str]:
        return [
            segment_annotation(segment)
            for segment in sorted(
                self.switch_regions + self.disordered_regions,
                key=lambda segment: (segment.start, segment.end, segment.category),
            )
        ]


def segment_annotation(segment: ComparedSegment) -> str:
    left = GROUP_SINGLE_LETTER[segment.dominant_group]
    right = GROUP_SINGLE_LETTER[segment.alternative_group]
    return f"{segment.start}-{segment.end} {left}->{right}"
