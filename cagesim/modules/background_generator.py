'''Class to generate background sequences
author: Kata Ferenc
email: katalitf@uio.no
last modified: 27th Feb 2023
'''

import random
import numpy as np

class BackgroundSeq:
    '''Class to generate background sequences
    '''

    def __init__(
            self,
            vocabulary: list[str],
            l_sequence: int,
            sequence_prior: list[float],
            seed: int = 43) -> None:
        '''Initialize generator
        Parameters
        ----------
        vocabulary:
            allowed characters in the sequence (eg [A, C, T, G] or 'ACTG')
        l_sequence:
            the length of the sequences
        sequence_prior:
            prior probabilities of each letter in the sequence
            (assume independence, dinucleotide combinations may be considered later)
        seed:
            value of the seed (only used if set_seed is True)
            (default = 43)
        '''
        self.seed = seed
        self.vocabulary = vocabulary
        self.l_sequence = l_sequence
        self.sequence_prior = sequence_prior


    # TODO: add prior to sequence composition, eg GC-content
    def generate_background_seq(
            self,
            n_sequences: int,
            set_seed: bool = True) -> str:
        '''Generates random sequences
        Parameters
        ----------
        n_sequences:
            the number of sequences to create
        set_seed:
            whether to use seed when generating the motifs
            (default = True)
        Return
        ------
        random_sequences:
            generated sequences to return
        '''
        if set_seed:
            random.seed(self.seed)
            np.random.seed(self.seed)

        for _ in range(n_sequences):
            random_sequence = ''.join(
                random.choice(self.vocabulary) for _ in range(self.l_sequence)
                )
            yield random_sequence
