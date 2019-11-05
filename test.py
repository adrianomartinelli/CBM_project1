#%%
from suffix_trees import STree
st = STree.STree("abcdefghab")

# %%
start = st.root
leaves =start._get_leaves()
junctions = []