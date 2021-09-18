import itertools
from functools import partial
from typing import Callable, List, Optional, Set

import pandas as pd
from synbio.polymers import DNA, Protein
from synbio.utils import dNTPs

import cli
import scoring
from ontology import *

__all__ = [
    # Classes
    "Recoder"
]


def get_aminos(split_seq: List[str]) -> AminoSet:
    on_frame = DNA(''.join(split_seq))
    off_frame = DNA(split_seq[1])

    def to_str(prot: Protein) -> List[str]:
        return [str(aa) for aa in prot]

    return AminoSet(
        to_str(on_frame.translate()),
        to_str(off_frame.translate())
    )


class Recoder:
    SEQ_INPUT_PROMPT = "input sequence; highlight serine codon with brackets \
        (" "e.g., 'g|tcg|ct'): "
    RECODING_CHOICE_PROMPT = "which recoding would you like to choose? " \
                             "(please input int): "

    def __init__(
            self,
            codons_to_remove: Optional[Set[str]] = None,
            dist_mat: Optional[pd.DataFrame] = None,
            score: Optional[
                Callable[
                    [pd.DataFrame, List[str], List[str]],  # takes 3 args
                    float]  # outputs 1 arg
            ] = None,
            verbose: Optional[bool] = None
    ):
        self.codons_to_remove = codons_to_remove \
            if codons_to_remove is not None \
            else {'TCG', 'TCA', 'TGA'}

        self.matrix = dist_mat if dist_mat is not None else scoring.blosum62

        scorefxn = score if score is not None else scoring.sum_score
        self.score = partial(scorefxn, self.matrix)
        self.verbose = verbose if verbose is not None else False
        self.log = Log()

        self.__call__()

    def __call__(self):
        while True:
            split_seq = self.ask_for_seq()
            recodings = self.generate_recodings(split_seq)
            choice = self.choose_recodings(recodings)
            self.update_log(choice)

            if split_seq is None:
                break

        self.write_log_prompt()

    @cli.userinterface(SEQ_INPUT_PROMPT)
    def ask_for_seq(self, user_input: str) -> List[str]:
        split_seq = user_input.split("|")

        # tests
        assert len(split_seq) == 3
        assert all(DNA()._seq_check(subseq) for subseq in split_seq)

        recode_frame = DNA(''.join(split_seq))
        assert len(recode_frame) == 9
        _ = recode_frame.translate()

        off_frame = DNA(split_seq[1])
        assert len(off_frame) == 6
        _ = off_frame.translate()

        return split_seq

    @cli.pipeline
    def generate_recodings(self, split_seq: List[str]) -> List[Recoding]:
        possible_nt_changes = [
            ''.join(tup)
            for tup in itertools.product(dNTPs, repeat=len(split_seq[1]))
        ]

        all_possible_sequences = [
            (split_seq[0], recode, split_seq[2])
            for recode in possible_nt_changes
        ]
        recoded_sequences = [
            list(tup) for tup in all_possible_sequences
            if ''.join(tup)[3:6].upper() not in self.codons_to_remove
        ]

        original_aminos = get_aminos(split_seq)
        print(f"original_aminos: {original_aminos}")
        recoded_aminos_list = [get_aminos(seq) for seq in recoded_sequences]

        scores = [
            self.score(original_aminos, recoded_aminos)
            for recoded_aminos in recoded_aminos_list
        ]

        recodings = [
            Recoding(
                original_seq="|".join(split_seq),
                recoded_seq="|".join(tup[0]),
                score=tup[1]
            ) for tup in zip(recoded_sequences, scores)
        ]
        ranked_recodings = sorted(
            recodings, key=lambda recoding: -recoding.score
        )

        return ranked_recodings

    @staticmethod
    def report_recodings(ranked_recodings: List[Recoding]) -> None:
        def report_(rec: Recoding):
            print(f"recoded sequence: {rec.recoded_seq} | score: "
                  f"{rec.score}")
            aminos = get_aminos(rec.recoded_seq.split("|"))
            print(f"recode frame aminos: {aminos[0:3]} | frame 2 aminos:"
                  f" {aminos[3:]}")

        print(f"I found {len(ranked_recodings)} recodings")
        print(f"The first ten recodings are: ")
        for i, recoding in enumerate(ranked_recodings[:10]):
            print(f"{i}: {'/' * 40}")
            report_(recoding)

    @cli.pipeline
    @cli.userinterface(RECODING_CHOICE_PROMPT)
    def choose_recodings(user_input: str,
                         recodings: List[Recoding]) -> Recoding:
        return recodings[int(user_input)]

    @cli.pipeline
    def update_log(self, choice: Recoding) -> Log:

        @cli.userinterface("gene name: ")
        def ask_for_gene(user_input):
            return user_input

        @cli.userinterface("position of recoding (aa of TCG/TCA codon): ")
        def ask_for_position(user_input):
            return int(user_input)

        @cli.userinterface("recoding notes: ")
        def ask_for_notes(user_input):
            return user_input

        self.log.original_seq.append(choice.original_seq)
        self.log.recoded_seq.append(choice.recoded_seq)
        self.log.score.append(choice.score)
        self.log.gene.append(ask_for_gene())
        self.log.position.append(ask_for_position())
        self.log.notes.append(ask_for_notes())

        return self.log

    @cli.pipeline
    @cli.userinterface("write log to filepath? ")
    def write_log_prompt(self, user_input: str) -> None:
        self.log.write_log(user_input)
