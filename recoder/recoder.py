import itertools
from functools import partial
from typing import Callable, List, Optional, Set

import cli
import pandas as pd
from synbio.polymers import DNA
from synbio.utils import dNTPs

import scoring
from ontology import *


# FUNCTIONS
def input_seq_validator(seq):
    split_seq = seq.split("|")

    bools = []
    # tests
    try:
        bools.append(len(split_seq) == 3)
        bools.append(all(DNA()._seq_check(subseq) for subseq in split_seq))

        recode_frame = DNA(''.join(split_seq))
        bools.append(len(recode_frame) == 9)
        _ = recode_frame.translate()

        off_frame = DNA(split_seq[1])
        bools.append(len(off_frame) == 6)
        _ = off_frame.translate()

        bools.append(True)
    except:
        bools.append(False)

    return all(bools)


def generate_recodings(
        user_input: str,
        codons_to_remove: Optional[Set[str]] = None,
        dist_mat: Optional[pd.DataFrame] = None,
        scoring_fxn: Optional[
            Callable[
                [pd.DataFrame, Site, Site],  # takes 3 args
                float  # outputs 1 arg
            ]
        ] = None,
) -> List[Recoding]:
    # set optional parameters
    if codons_to_remove is None:
        codons_to_remove = {'TCG', 'TCA', 'TGA'}

    if dist_mat is None:
        dist_mat = scoring.blosum62

    if scoring_fxn is None:
        scoring_fxn = scoring.sum_score
    scoring_fxn = partial(scoring_fxn, dist_mat)

    # generate site for original sequence
    split_seq = user_input.split('|')
    original_site = Site.from_seq(user_input)
    print(f"Original aminos: on frame = {original_site.on_frame_prot}; "
          f"off_frame = {original_site.off_frame_prot}")

    # generate all possible recoded sites
    possible_nt_changes = [
        ''.join(tup)
        for tup in itertools.product(dNTPs, repeat=len(split_seq[1]))
    ]
    all_possible_recoded_sites = [
        (split_seq[0], recode, split_seq[2])
        for recode in possible_nt_changes
    ]
    recoded_sites = [
        Site.from_seq('|'.join(tup)) for tup in all_possible_recoded_sites
        if ''.join(tup)[3:6].upper() not in codons_to_remove
    ]

    # calculate scores for each recoding
    scores = [
        scoring_fxn(original_site, recoded_site)
        for recoded_site in recoded_sites
    ]

    # generate all recodings and rank them
    recodings = [
        Recoding(
            original_site=original_site,
            recoded_site=recoded_site,
            score=score
        ) for recoded_site, score in zip(recoded_sites, scores)
    ]
    ranked_recodings = sorted(
        recodings, key=lambda recoding: -recoding.score
    )

    return ranked_recodings


def report_recodings(ranked_recodings: List[Recoding]):
    def report_(rec: Recoding):
        print(f"recoded sequence: {rec.recoded_site.seq_str} | "
              f"score: {rec.score}")
        print(f"recode frame aminos: {rec.recoded_site.on_frame_prot} | "
              f"off frame aminos: {rec.recoded_site.off_frame_prot}")

    print(f"I found {len(ranked_recodings)} recodings")
    print(f"The first ten recodings are: ")

    for i, recoding in enumerate(ranked_recodings[:10]):
        print(f"{i}: {'/' * 40}")
        report_(recoding)


def update_log(
        log: Log, rec: Recoding,
        gene_name: str, position: str, notes: str
) -> Log:
    log.original_seq.append(rec.original_site.seq_str)
    log.recoded_seq.append(rec.recoded_site.seq_str)
    log.score.append(rec.score)
    log.gene.append(gene_name)
    log.position.append(position)
    log.notes.append(notes)

    return log


if __name__ == "__main__":
    # script parameters
    SCORING_FXN = None

    # long prompts
    seq_input_prompt = "input sequence; highlight serine codon with " \
                       "brackets (e.g., 'ct|ttcgga|t'): "
    choice_prompt = f"which recoding would you like to choose? (please input " \
                    f"int): "

    # initialize log and loop variable
    log = Log()
    looping = True

    # main event loop
    while looping:
        # loop: input sequence to recode
        input_seq = cli.prompt_user(seq_input_prompt, input_seq_validator)
        if looping and input_seq is not None:
            # generate recodings for input sequence
            ranked_recodings = generate_recodings(
                input_seq, scoring_fxn=SCORING_FXN
            )
            report_recodings(ranked_recodings)

            # choose particular recoding
            choice_validator = cli.choice_validator_factory(ranked_recodings)
            chosen_recoding_ix = cli.prompt_user(choice_prompt,
                                                 choice_validator)
        else:
            looping = False
            chosen_recoding_ix = None
            ranked_recodings = None

        # loop: choose and log recoding
        if looping and chosen_recoding_ix is not None:
            chosen_recoding = ranked_recodings[int(chosen_recoding_ix)]

            # update log
            gene_name = cli.prompt_user("gene name: ")
            position = cli.prompt_user("position of recoding (aa of TCG/TCA "
                                       "codon): ")
            notes = cli.prompt_user("notes for this recoding?: ")
            log = update_log(log, chosen_recoding, gene_name, position, notes)
        else:
            looping = False

    # write log
    filename = cli.prompt_user("write log to file? (input filename): ")
    if filename is not None:
        log.write_log(filename)
