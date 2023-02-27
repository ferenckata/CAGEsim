'''Class to add motifs to sequences
author: Kata Ferenc
email: katalitf@uio.no
last modified: 27th Feb 2023
'''

import random
import numpy as np

class InsertMotif:
    '''Class to add motifs to sequences
    '''

    def __init__(
            self,
            seed: int = 43) -> None:
        '''Initialize inserter
        Parameters
        ----------
        seed:
            value of the seed (only used if set_seed is True)
            (default = 43)
        '''
        self.seed = seed
        self.num_positions = None
        self.motif_instances = None


    def sequence_partition(self, sequence: str, parts: int) -> list[str]:
        '''Cut the sequence into as many parts as many potential insertions can happen
        Parameters
        ----------
        sequence:
            character string
        parts:
            number of parts to cut into
        Return
        ------
        seq_parts:
            list of subsequences
        '''
        len_part = int(len(sequence)/parts)
        seq_parts = [sequence[i : (i + len_part)] for i in range(0, len(sequence), len_part)]
        return seq_parts
    

    def select_positions(
            self,
            k_positions: int,
            l_motif: int,
            l_sequence: int,
            set_seed: bool = True) -> list[int]:
        '''Generate positions within sequences to insert motifs
        Parameters
        ----------
        k_positions:
            number of positions to generate
        l_motif:
            length of motif to insert
        l_sequence:
            length of sequence where the motif is to be inserted
        set_seed:
            whether to use seed when generating the motifs
            (default = True)
        Return
        ------
        motif_indeces:
            list of indeces in each subsequence for each motif to be inserted
        '''
        if set_seed:
            random.seed(self.seed)
            np.random.seed(self.seed)

        motif_indeces = [random.randint(0, l_sequence - l_motif) for _ in range(k_positions)]
        return motif_indeces


    def insert_motif(
            self,
            random_sequence: str,
            motif: np.ndarray,
            position: int,
            vocabulary: list(str),
            set_seed: bool = True) -> str:
        '''Adds a given motif in the sequences
        Parameters
        ----------
        random_sequence:
            input string of a random sequence
        motif:
            PWM, ie the probability of each character to be inserted (len_motif x len_vocab),
            where row sum = 1
        position:
            where should the sequence be inserted within the random sequence
        vocabulary:
            allowed characters in the sequence (eg [A, C, T, G] or 'ACTG')
            order is considered corresponding to the order of probabilities within the motif
        set_seed:
            whether to use seed when generating the motifs
            (default = True)
        seed:
            value of the seed (only used if set_seed is True)
            (default = 43)
        Return
        ------
        mod_random_sequence:
            sequence with motif inserted
        '''
        if set_seed:
            random.seed(self.seed)
            np.random.seed(self.seed)

        if (position + motif.shape[0]) > len(random_sequence):
            print(f'generated insertion position is {position}')
            print(f'motif shape is {motif.shape[0]}')
            print(f'length of generated sequence is {random_sequence}')
            raise IndexError("Motif does not fit into sequence. \
                Provide shorter motif, smaller position index or longer sequence.")

        if len(vocabulary) != motif.shape[1]:
            print(f'length of vocabulary is {len(vocabulary)}')
            print(f'motif length is {motif.shape[1]}')
            raise Exception("Vocabulary length does not match motif shape. \
                Please provide probability for all letters.")

        str_seq = list(random_sequence)
        for i in range(motif.shape[0]):
            str_seq[position + i] = random.choices(vocabulary, weights=motif[i,], k=1)[0]
        mod_random_sequence = ''.join(str_seq)

        return mod_random_sequence