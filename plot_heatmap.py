import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['figure.figsize'] = 14,12
matplotlib.rcParams['font.weight'] = "medium"
matplotlib.rcParams['axes.labelweight'] = 'medium'
matplotlib.rcParams['figure.titleweight'] = 'medium'
matplotlib.rcParams['axes.titleweight'] = 'medium'

sns.set(font_scale=1.1, style="ticks", font="Lato")
df = pd.read_csv('tf_classes_heatmap.tsv', sep='\t',index_col=0)
g = sns.heatmap(df,vmin=0.5, vmax = 1, linewidths=.25, linecolor='darkgray',cmap='YlOrRd')
g = g.get_figure()
g.tight_layout()
g.savefig(f"heatmap_classes.svg", dpi=300)