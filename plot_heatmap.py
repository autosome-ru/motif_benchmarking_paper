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

no_labels = False
if '--no-labels' in sys.argv:
    sys.argv.remove('--no-labels')
    no_labels = True

if '--correlation' in sys.argv:
    # correlation mode
    sys.argv.remove('--correlation')
    heatmap_args = dict(vmin= -1, vmax = 1, cmap='coolwarm')
else:
    # ROC mode
    heatmap_args = dict(vmin=0.5, vmax = 1, cmap='YlOrRd')

heatmap_matrix = sys.argv[1]
heatmap_svg = sys.argv[2]

sns.set(font_scale=2.5, style="ticks", font="Lato")
df = pd.read_csv(heatmap_matrix, sep='\t',index_col=0)
labels = [re.sub(r'.*--skip-me--.*','', k) for k in df.keys()]

if no_labels:
    labels_cfg = {'xticklabels': False, 'yticklabels': False}
else:
    labels_cfg = {'xticklabels': labels, 'yticklabels': labels}

g = sns.heatmap(df,
    **heatmap_args,
    square=True,
    linewidths=.25, linecolor='darkgray',
    **labels_cfg,
)

g = g.get_figure()
g.tight_layout()
g.savefig(heatmap_svg, dpi=300)
