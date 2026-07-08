from dataclasses import dataclass, field
from typing import List, Optional


@dataclass
class ResidueSS:
    position_index: int
    amino_acid: str
    dssp_code: str
    chain_id: str
    residue_number: str
    bp1: int = 0
    bp2: int = 0

@dataclass
class OpenDisorderSegment:
    start: int
    dominant_code_group: str
    alternative_code_group: str

@dataclass
class ChangedRegion:
    start: int
    end: int
    dominant_code_group: str
    alternative_code_group: str

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    @property
    def annotation(self) -> str:
        return f"{self.start}-{self.end}:{self.dominant_code_group}->{self.alternative_code_group}"


@dataclass
class BridgeReshuffledRegion:
    start: int
    end: int
    positions: List[int]

    @property
    def length(self) -> int:
        return len(self.positions)

    @property
    def annotation(self) -> str:
        return f"{self.start}-{self.end}"


@dataclass
class RigidBodyRegion:
    start: int
    end: int
    dssp_code_group: str
    max_window_distance: float

    @property
    def annotation(self) -> str:
        return f"{self.start}-{self.end}:{self.dssp_code_group}:{self.max_window_distance:.2f}A"


@dataclass
class RecoveryResult:
    alternative_model: str
    source: str
    encounter_index: int
    sequences_identical: Optional[bool]
    detector_hit: str
    bridge_reshuffled_regions: List[BridgeReshuffledRegion] = field(default_factory=list)
    disorder_change_regions:   List[ChangedRegion]          = field(default_factory=list)
    rigid_body_regions:        List[RigidBodyRegion]        = field(default_factory=list)
    dominant_dssp:    str = ""
    alternative_dssp: str = ""
    invalid_reason:   str = ""

    @property
    def recovered(self) -> bool:
        return bool(self.detector_hit)

    @property
    def bridge_reshuffled_residue_count(self) -> int:
        return sum(r.length for r in self.bridge_reshuffled_regions)

    @property
    def disorder_change_residue_count(self) -> int:
        return sum(r.length for r in self.disorder_change_regions)

    @property
    def rigid_body_max_window_distance(self) -> float:
        if not self.rigid_body_regions:
            return 0.0
        return max(r.max_window_distance for r in self.rigid_body_regions)