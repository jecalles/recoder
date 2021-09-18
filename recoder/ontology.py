from typing import List
from dataclasses import dataclass, field
import pathlib
from datetime import datetime

import pandas as pd

from . import cli


__all__ = [
    "Recoding", "Log", "AminoSet",
]


@dataclass
class Recoding:
    original_seq: str
    recoded_seq: str
    score: float


@dataclass
class AminoSet:
    on_frame: List[str]
    off_frame: List[str]

    def flatten(self):
        return self.on_frame + self.off_frame


@dataclass
class Log:
    original_seq: List[str] = field(default_factory=list)
    recoded_seq: List[str] = field(default_factory=list)
    score: List[float] = field(default_factory=list)
    gene: List[str] = field(default_factory=list)
    position: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    def write_log(self, filename: str) -> None:
        outpath = pathlib.Path(f"{filename} @ {datetime.now()}.csv")
        df = pd.DataFrame.from_dict(vars(self))
        df.to_csv(outpath)
