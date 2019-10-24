from os import path as osp

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import math
import random
from collections import defaultdict

from nxviz.plots import ArcPlot, CircosPlot, BasePlot
import nxviz

curr_path = osp.dirname(osp.abspath(__file__))

tf_fams_df = pd.read_csv('source_data/annotation/motifs_prefinal.tsv', sep='\t', header=None)
tf_fams = {row[2]: row[5] for idx,row in tf_fams_df.iterrows()}

def load_motifs_network():
    df = pd.read_csv("source_data/best_motifs.tsv",sep="\t")
    G = nx.Graph()
    
    recognize_counts = defaultdict(int)


    # tfs = set([rec[0] for idx, rec in df.iterrows()]).union([rec[3] for idx, rec in df.iterrows()])
    # for tf in tfs:
    #     fam = tf_fams.get(tf,'unknown')
    #     G.add_node(tf, tf_fam=fam)

        # target_tf, motif, auc, recognizing_tf = rec[0:4]
    for idx, rec in df.iterrows():
        target_tf, motif, auc, recognizing_tf = rec[0:4]
        target_fam = tf_fams.get(target_tf,'unknown')
        recognizing_fam = tf_fams.get(recognizing_tf,'unknown')

        if not isinstance(target_fam, str):
            target_fam = 'unknown'
        if not isinstance(recognizing_fam, str):
            recognizing_fam = 'unknown'

        if target_fam == 'unknown':# or recognizing_fam == 'unknown':
            continue
        recognize_counts[recognizing_tf] += 1

        G.add_node(target_tf, bipartite="target_tf", tf_fam=target_fam)
        G.add_node(recognizing_tf, bipartite="recognizing_tf", tf_fam=recognizing_fam)
        if target_tf != recognizing_tf:
            G.add_edge(target_tf, recognizing_tf, target_fam=target_fam, recognizing_fam=recognizing_fam)

    return G, recognize_counts

G,recognize_counts = load_motifs_network()
m = CircosPlot(
    G, node_grouping="tf_fam", node_color="tf_fam",
    group_label_position='middle', group_label_color=True,
    #, node_order="tf_fam"
    edge_color='target_fam',
    node_labels= lambda node: recognize_counts[node] >= 4,
    node_label_layout='rotation',
    fontsize=8,

)
m.compute_node_colors()
m.draw()
plt.show()


# # # Make an ArcPlot.
# a = ArcPlot(
#     G, node_color="tf_fam", node_grouping="tf_fam"#, node_order="connectivity"
# )
# a.draw()
# plt.show()
