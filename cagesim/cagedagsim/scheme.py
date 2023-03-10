import dagsim.base as ds
import numpy as np
from cagesim.cagedagsim.sim_methods import SimMethods

num_topics = ds.Node(
    name="num_topics",
    function=np.random.poisson,
    args=[2.0]
)

# 10 cell types with equal probability for now
# later add to config TODO
cell_type = ds.Node(
    name="cell_type",
    function=np.random.randint,
    args = [0,10]
)

topic = ds.Node(
    name="topic",
    function=SimMethods.topic_selection,
    args=[num_topics, cell_type]
)

motif_length = ds.Node(
    name="motif_length",
    function=np.random.normal,
    args=[7.0,1.0]
)

motif_prior = ds.Node(
    name="PFM",
    function=SimMethods.simulate_pfm,
    args=[motif_length]
)

num_motifs = ds.Node(
    name="num_motif",
    function=np.random.poisson,
    args=[3.0]
)

orientation = ds.Node(
    name="orientation",
    function=np.random.binomial,
    args=[1,.5]
)

motif = ds.Node(
    name="motif",
    function=SimMethods.simulate_motif,
    args=[motif_prior, topic, num_motifs, orientation]
)

activity_score = ds.Node(
    name="activity_score",
    function=SimMethods.simulate_activity,
    args=[topic]
)

background_seq = ds.Node(
    name="background",
    function=SimMethods.sim_background
)

motived_seq = ds.Node(
    name="cageseq",
    function=SimMethods.simulate_cage,
    args=[motif, background_seq]
)