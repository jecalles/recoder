from typing import List

import pandas as pd

from ontology import *

__all__ = [
    # matrices
    'blosum62',
    # scoring functions
    'sum_score'
]

blosum62 = pd.read_csv("formatted_blosum62.csv", index_col=0)


def sum_score(
        matrix: pd.DataFrame,
        original_aminos: AminoSet,
        recoded_aminos: AminoSet
) -> float:

    return sum(
        matrix[old_aa][new_aa]
        for old_aa, new_aa in zip(
            original_aminos.flatten(), recoded_aminos.flatten()
        )
    )
