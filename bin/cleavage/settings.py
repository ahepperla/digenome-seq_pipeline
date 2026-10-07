"""Caller settings, read from one JSON object written by main.nf.

Every name matches a pipeline parameter. main.nf has already applied the
--keep_multimappers overrides (both minimum MAPQs and the support mean-MAPQ
filter become 0), and nf-schema has already checked types and ranges, so this
module only checks that the object has exactly the expected keys.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, fields
from pathlib import Path

ANALYSES = ("digenome", "ndigenome")


@dataclass(frozen=True)
class CallerSettings:
    analysis: str
    keep_multimappers: bool

    digenome_overhang: int
    digenome_pair_window: int
    digenome_min_mapq: int
    digenome_forward_cutoff: int
    digenome_reverse_cutoff: int
    digenome_depth_cutoff: int
    digenome_fraction_cutoff: float
    digenome_pair_score_cutoff: float

    ndigenome_min_count: int
    ndigenome_min_fraction: float
    ndigenome_min_mapq: int
    ndigenome_opposite_window: int
    ndigenome_ambiguous_min_count: int
    ndigenome_ambiguous_min_fraction: float

    cleavage_artifact_window: int
    cleavage_max_softclip_fraction: float
    cleavage_max_indel_fraction: float
    cleavage_min_support_mean_mapq: float
    cleavage_control_min_depth: int
    cleavage_control_max_fraction: float
    cleavage_control_min_fold: float
    cleavage_control_max_q: float

    @property
    def min_mapq(self) -> int:
        """The minimum MAPQ of the selected mode."""
        if self.analysis == "ndigenome":
            return self.ndigenome_min_mapq
        return self.digenome_min_mapq

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def settings_from_dict(values: dict[str, object]) -> CallerSettings:
    expected = {field.name for field in fields(CallerSettings)}
    missing = sorted(expected - set(values))
    unknown = sorted(set(values) - expected)
    if missing or unknown:
        raise ValueError(
            f"Caller settings mismatch: missing {missing}, unknown {unknown}"
        )
    if values["analysis"] not in ANALYSES:
        raise ValueError(f"Unknown analysis: {values['analysis']}")
    return CallerSettings(**values)


def load_settings(path: Path) -> CallerSettings:
    return settings_from_dict(json.loads(Path(path).read_text()))
