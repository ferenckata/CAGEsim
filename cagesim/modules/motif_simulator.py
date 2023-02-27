'''Class to generate and select sequence motifs
author: Kata Ferenc
email: katalitf@uio.no
last modified: 27th Feb 2023
'''

import random
import numpy as np

class MotifSimulator:
    '''Class to generate and select sequence motifs
    '''

    def __init__(
            self,
            n_motif: int,
            l_motif: int = 5,
            seed: int = 43) -> None:
        '''Initialize simulator
        Parameters
        ----------
        n_motif:
            the number of motifs to be generated
        l_motif:
            the length of each motif (assume fixed lenght, may change it later)
            (default = 5)
        seed:
            value of the seed (only used if set_seed is True)
            (default = 43)
        '''
        self.n_motif = n_motif
        self.l_motif = l_motif
        self.seed = seed
        self.archetypes = None
        self.motif_instances = None


    # TODO: define topics first with motifs and assign topics to cell types (another function maybe)
    def generate_motif_archetypes(
            self,
            motif_dirichlet_prior: tuple[float] = (1,1,1,1),
            cell_dirichlet_prior: tuple[float] = (1,1,1,1,1,1,1,1),
            set_seed: bool = True
        ) -> tuple[list[np.ndarray], np.ndarray]:
        '''Generate motif PFMs
        Parameters
        ----------
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
        set_seed:
            whether to use seed when generating the motifs
            (default = True)
        Return
        ------
        motif_list:
            a list of PWM-style motifs with probabilities for each letter
        cell_type_matrix:
            a 2D matrix containing the probability of each cell type for each motif
            (n_motif x n_cell_type)
        '''
        if set_seed:
            random.seed(self.seed)
            np.random.seed(self.seed)
        # generate priors for each motif
        motif_priors = np.random.dirichlet(motif_dirichlet_prior, self.n_motif)
        motif_list = []
        cell_type_matrix = np.zeros((self.n_motif, len(cell_dirichlet_prior)))
        for m_index, motif in enumerate(motif_priors):
            # select letters from the generated Dirichlet prior for each position of each motif
            motif_list.append(np.random.multinomial(1, motif, self.l_motif))
            # probabilistically assign each motif to cell types
            cell_type_matrix[m_index, ] = np.random.dirichlet(cell_dirichlet_prior, 1)

        return motif_list, cell_type_matrix


    def generate_motif_instances(self):
        '''Simulate motifs from archetype PFMs
        (assume position independence)'''
        pass