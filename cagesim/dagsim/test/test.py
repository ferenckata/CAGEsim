# Testing dagsim

#%%
import dagsim.base as ds
import numpy as np

#%%
from random import choices

def simulate_sequence(seq_len, p_head):
    return "".join(choices(
        ["H", "T"],
        [p_head, 1-p_head],
        k=seq_len
    ))
#%%
sequence_length = ds.Node(
    name="seq_len",
    function=np.random.randint,
    args=[10,20]
)
p_head = ds.Node(
    name="p_head",
    function=np.random.uniform
)
sequence = ds.Node(
    name="sequence",
    function=simulate_sequence,
    args=[sequence_length, p_head]
)
#%%
listNodes = [sequence_length, p_head, sequence]
my_graph = ds.Graph(listNodes, "Graph1")
#%%
my_graph.draw()
#%%
data = my_graph.simulate(
    num_samples=100,
    csv_name="hello_world"
)

#%%
import numpy as np
#%%
np.random.binomial(1,0.5,10)


