'''Class to assign motifs and cell types to topics
author: Kata Ferenc
email: katalitf@uio.no
last modified: 27th Feb 2023
'''

import random
import numpy as np

class TopicSelector:
    '''Class to assign motifs and cell types to topics'''

    def __init__(self) -> None:
        self.n_topics = None
        self.motif_topic_posterior = None
        self.topic_cell_posterior = None
        self.num_positions = None
        self.selected_motif_indeces = None
        self.sequence_signal_per_cell = None

    def select_num_positions(self):
        '''Select number of motifs to be inserted per sequence'''
        pass

    def select_num_topics(self):
        '''Select number of topics to be selected from per sequence'''
        pass

    def select_topic(self):
        '''Select topic(s) per sequence'''
        pass

    def select_motif(self):
        '''Select motif (index) from topic'''
        pass

    def calculate_cell_signal(self):
        '''Calculate overall expected cell signal given number of motifs per topic'''
        pass

    # TODO: select from topic
    def select_motifs(
            self,
            motif_list: list(np.ndarray),
            num_motifs: int,
            cell_type_matrix: np.ndarray,
            motif_prob: float = 0.1,
            set_seed: bool = True
        ) -> list(list(str), list(int), list(np.ndarray)):
        '''Randomly select motifs and their respective locations
        Parameters
        ----------
        motif_list:
            a list of PWM-style motifs with probabilities for each letter
        num_motifs:
            number of motifs to be inserted
        cell_type_matrix:
            a 2D matrix containing the probability of each cell type for each motif
            (n_motif x n_cell_type)
        motif_prob:
            probability (threshold) of a motif to be inserted [0,1]
        set_seed:
            whether to use seed when generating the motifs
            (default = True)
        Return
        ------
        motifs:
            list of selected motifs (PWM-style)
        passing_indeces:
            index of substrings and corresponding motifs to be used
        cell_types:
            list of cell types corresponding to selected motifs
        '''
        if set_seed:
            random.seed(self.seed)
            np.random.seed(self.seed)
        # assign probability to each motif
        k_probabilities = [random.random() for _ in range(num_motifs)]
        # select motif from list (with replacement)
        motif_idxs = random.choices(range(len(motif_list)), k = num_motifs)
        motifs = list(motif_list[i] for i in motif_idxs)
        # fetch the cell annotation of the sequence
        cell_types = list(cell_type_matrix[i,] for i in motif_idxs)
        # check if any probability is below insertion threshold
        passing_indeces = [idx for idx, val in enumerate(k_probabilities) if val < motif_prob]
        return passing_indeces, motifs, cell_types
