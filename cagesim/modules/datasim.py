'''Class for simulating sequence data with inserted patterns'''
# code modified from Tim / Ieva

import random
import numpy as np

from cagesim.utils.inputoutput import InputOutput

class MotivedSequenceSimulator:
    '''Class to support data simulation tasks
    '''


    @classmethod
    def generate_patterned_seq(
            cls,
            total_num_motifs: int,
            n_sequences: int,
            l_sequence: int,
            num_inserted_motifs: int,
            vocabulary: list[str],
            l_motif: int = 5,
            motif_dirichlet_prior: tuple[float] = (1,1,1,1),
            cell_dirichlet_prior: tuple[float] = (1,1,1,1,1,1,1,1),
            seq_proportion: float = 0.9,
            motif_probability_threshold: float = 0.1,
            set_seed: bool = True,
            seed: int = 43
        ) -> tuple[str, np.ndarray]:
        '''
        Parameters
        ----------
        total_num_motifs:
            the total number of motifs to be generated
        n_sequences:
            the number of sequences to create
        l_sequence:
            the length of the sequences
        num_inserted_motifs:
            number of motifs to be inserted in each sequence
        vocabulary:
            allowed characters in the sequence (eg [A, C, T, G] or 'ACTG')
        l_motif:
            length of motif to insert
            (default = 5)
        motif_dirichlet_prior:
            prior weights of each letter as described in the np.random.dirichlet manual
            if all values are uniformly high, each letter gets about the same
            probability in the PWM
            if all values are uniformly low, one (randomly selected) letter
            gets the highest probability
            tuple length should be same as vocabulary size
            (default = (1,1,1,1))
        cell_dirichlet_prior:
            prior weights of each cell type as described in the np.random.dirichlet manual
            if all values are uniformly high, each cell type gets about the same
            probability for a motif
            if all values are uniformly low, one (randomly selected) cell type
            gets the highest probability
            tuple length should be the number of cell types to be simulated
            (default = (1,1,1,1,1,1,1,1))
        seq_proportion:
            proportion of sequences that should have the motif inserted (uniform distribution)
            (default = 0.9)
        motif_probability_threshold:
            probability (threshold) of a motif to be inserted [0,1]
            (default = 0.1)
        set_seed:
            whether to use seed when generating the motifs
            (default = True)
        seed:
            value of the seed (only used if set_seed is True)
            (default = 43)
        Return
        ------
        final_sequence:
            random sequence with inserted motifs
        p_cell_seq:
            probability of cell type given sequence
        '''
        if set_seed:
            random.seed(seed)
            np.random.seed(seed)

        if num_inserted_motifs is None:
            # select number of motifs
            num_inserted_motifs = random.randint(0, l_sequence/l_motif - 1)

        # get motifs to insert
        motif_list, cell_type_matrix = cls.generate_motif(
            n_motif = total_num_motifs,
            l_motif = l_motif,
            motif_dirichlet_prior = motif_dirichlet_prior,
            cell_dirichlet_prior = cell_dirichlet_prior,
            set_seed = set_seed,
            seed = seed)

        random_state = random.getstate()
        numpy_random_state = np.random.get_state()
        randomness_lines = []
        randomness_lines.append(f"random state: {random_state}\n")
        randomness_lines.append(f"numpy random state: {numpy_random_state}\n")
        # save random state to file
        InputOutput.list_to_file(randomness_lines, 'random_state')
        # save motifs to file
        InputOutput.list_to_file(motif_list, 'motifs')
        # save cell type matrix to file
        InputOutput.list_to_file(cell_type_matrix, 'motif_cell_type_matrix')

        sequence_generator = cls.generate_random_seq(
                n_sequences = n_sequences,
                l_sequence = l_sequence,
                vocabulary = vocabulary,
                set_seed = set_seed,
                seed = seed)

        for rand_seq in sequence_generator:
            # decide whether motif should be inserted into the sequence
            if random.random() < seq_proportion:
                # select motifs to insert
                passing_indeces, motifs, cell_types = cls.select_motifs(
                    motif_list = motif_list,
                    num_motifs = num_inserted_motifs,
                    cell_type_matrix = cell_type_matrix,
                    motif_prob = motif_probability_threshold,
                    set_seed = False)

                # if any motif passed the insertion threshold, continue
                if len(passing_indeces) > 0:

                    # cut sequence to as many as the maximal possible number of motifs would be
                    seq_parts = cls.sequence_partition(
                        sequence = rand_seq,
                        parts = num_inserted_motifs)

                    # randomly select positions for each motif within subsequence
                    motif_indeces = cls.select_positions(
                        k_positions = num_inserted_motifs,
                        l_motif = l_motif,
                        l_sequence = int(l_sequence/num_inserted_motifs),
                        set_seed = False)

                    # insert motifs into parts
                    motif_sequences = []

                    # initialize the unnormalized cell type probability array
                    # for each sequence with a small value
                    p_seq_cell = np.ones(cell_types[0].shape)

                    for k in range(num_inserted_motifs):
                        if k in passing_indeces:
                            sequence_with_motif = cls.insert_motif(
                                random_sequence = seq_parts[k],
                                motif = motifs[k],
                                position = motif_indeces[k],
                                vocabulary = vocabulary,
                                set_seed = False)
                            motif_sequences.append(sequence_with_motif)
                            # if the motif is in the sequence, multiply with cell type probability
                            # this here is p(seq | cell_type)
                            p_seq_cell *= np.array(cell_types[k])
                        else:
                            motif_sequences.append(seq_parts[k])
                    final_sequence = ''.join(motif_sequences)
                    # p (cell_type | seq) = p(seq | cell_type) * p(cell_type) / p(seq),
                    # where p(seq) = sum_i(p(seq | t_i) * p(t_i))
                    p_cell = 1/len(cell_dirichlet_prior)
                    p_seqcell = p_seq_cell * p_cell
                    p_seq = np.sum(p_seqcell)
                    p_cell_seq = p_seqcell / p_seq
                else:
                    final_sequence = rand_seq
                    p_cell_seq = np.zeros(len(cell_dirichlet_prior))
            else:
                final_sequence = rand_seq
                p_cell_seq = np.zeros(len(cell_dirichlet_prior))

            yield final_sequence, p_cell_seq
