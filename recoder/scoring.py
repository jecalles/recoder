from pathlib import Path

import pandas as pd

from ontology import Site

__all__ = [
    # matrices
    'blosum62',
    # scoring functions
    'sum_score', 'square_sum_score'
]


PARENT_FILEPATH = Path(__file__).absolute().parent
blosum62 = pd.read_csv(
    Path(PARENT_FILEPATH, "data/formatted_blosum62.csv"),
    index_col=0
)


def sum_score(
        matrix: pd.DataFrame,
        original_aminos: Site,
        recoded_aminos: Site
) -> float:

    return sum(
        matrix[old_aa][old_aa] - matrix[old_aa][new_aa]
        for old_aa, new_aa in zip(
            original_aminos.get_all_aminos(), recoded_aminos.get_all_aminos()
        )
    )


def distribute_changes_score(
        matrix: pd.DataFrame,
        original_aminos: Site,
        recoded_aminos: Site
) -> float:
    """ A scoring function that prefers many small amino acid changes
    distributed across both reading frames rather than few big changes in one
    particular reading frame"""
    on_frame_score = sum(
        (matrix[old_aa][old_aa] - matrix[old_aa][new_aa])**1.1
        for old_aa, new_aa in zip(
            original_aminos.on_frame_prot, recoded_aminos.on_frame_prot
        )
    )
    off_frame_score = sum(
        (matrix[old_aa][old_aa] - matrix[old_aa][new_aa])**1.1
        for old_aa, new_aa in zip(
            original_aminos.off_frame_prot, recoded_aminos.off_frame_prot
        )
    )
    return on_frame_score**1.5 + off_frame_score**1.5
