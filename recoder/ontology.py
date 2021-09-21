import pathlib
from typing import List
from dataclasses import dataclass, field
from datetime import datetime

import pandas as pd

from synbio.polymers import DNA

__all__ = [
    # Dataclasses
    "Site", "Recoding", "Log"
]

@dataclass
class Site:
    seq_str: str
    on_frame_dna: str
    on_frame_prot: str
    off_frame_dna: str
    off_frame_prot: str

    @classmethod
    def from_seq(cls, seq:str):
        split_seq = seq.split('|')
        on_frame_DNA = DNA(''.join(split_seq))
        off_frame_DNA = DNA(split_seq[1])

        new_site = cls(
            seq_str=seq,
            on_frame_dna=str(on_frame_DNA),
            off_frame_dna = str(off_frame_DNA),
            on_frame_prot=str(on_frame_DNA.translate()),
            off_frame_prot = str(off_frame_DNA.translate())
        )
        return new_site

    def get_all_aminos(self) -> str:
        return self.on_frame_prot + self.off_frame_prot

@dataclass
class Recoding:
    original_site: Site
    recoded_site: Site
    score: float


@dataclass
class Log:
    name_of_recoding: List[str] = field(default_factory=list)
    original_seq: List[str] = field(default_factory=list)
    recoded_seq: List[str] = field(default_factory=list)
    score: List[float] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    def write_log(self, filename: str) -> None:
        outpath = pathlib.Path(f"{filename} @ {datetime.now()}.csv")
        df = pd.DataFrame.from_dict(vars(self))
        df.to_csv(outpath)

