'''Define simulation methods for nodes'''
import yaml
from yaml.loader import SafeLoader

from cagesim.modules.background_generator import BackgroundSeq


class SimMethods:

    def __init__(self, config_file_name) -> None:
        with open(config_file_name) as config_file:
            self.config = yaml.load(config_file, Loader=SafeLoader)

    def topic_selection(self, num_topic, cell_type):
        raise NotImplementedError


    def simulate_pfm(self, motif_length):
        raise NotImplementedError


    def simulate_motif(self, motif_prior, topic, orientation):
        raise NotImplementedError


    def simulate_activity(self, topic, num_motifs):
        raise NotImplementedError


    def sim_background(self) -> str:
        vocabulary = self.config["background_vocab"]
        l_sequence = self.config["background_len"]
        sequence_prior = self.config["background_prior"]
        n_sequences = self.config["num_sequences"]
        backseq = BackgroundSeq(vocabulary, l_sequence, sequence_prior)
        return backseq.generate_background_seq(n_sequences)


    def simulate_cage(self, num_motifs, motif, background_seq):
        raise NotImplementedError

