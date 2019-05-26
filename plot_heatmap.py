import sys
import re
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['figure.figsize'] = 18,16
matplotlib.rcParams['font.weight'] = "medium"
matplotlib.rcParams['axes.labelweight'] = 'medium'
matplotlib.rcParams['figure.titleweight'] = 'medium'
matplotlib.rcParams['axes.titleweight'] = 'medium'

heatmap_matrix = sys.argv[1]
heatmap_svg = sys.argv[2]

sns.set(font_scale=0.8, style="ticks", font="Lato")
df = pd.read_csv(heatmap_matrix, sep='\t',index_col=0)
labels = [re.sub(r'.*--skip-me--.*','', k) for k in df.keys()]
g = sns.heatmap(df,
    vmin=0.5, vmax = 1, 
    cmap='YlOrRd',
    linewidths=.25, linecolor='darkgray',
    xticklabels=labels,
    yticklabels=labels,
)
g = g.get_figure()
g.tight_layout()
g.savefig(heatmap_svg, dpi=300)
