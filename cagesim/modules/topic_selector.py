'''Select motifs in topics and topics in cells
author: Kata Ferenc
email: katalitf@uio.no
last modified: 27th Feb 2023
'''

class TopicSelector:
    '''Class to select motif-topic and topic-cell correspondance'''

    def __init__(self) -> None:
        self.n_motifs = None
        self.n_topics = None
        self.n_cell_types = None
        self.motif_topic_prior = None
        self.topic_cell_prior = None
        self.motif_topic_posterior = None
        self.topic_cell_posterior = None

    
    def assign_motif_archetype_to_topic(self):
        '''Assign each motif archetype to a topic with a probability'''
        pass

    def assign_topic_to_cell(self):
        '''Assign topics to each cell type with certain probability'''
        pass
