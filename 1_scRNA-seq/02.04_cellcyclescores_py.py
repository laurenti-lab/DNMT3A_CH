# %% [markdown]
# # 02.04 cell cycle analysis
#
# Author:GM <br>
# useful links: https://satijalab.org/seurat/archive/v2.4/cell_cycle_vignette.html

# %%
# %matplotlib widget
# %load_ext watermark
 
import warnings
warnings.filterwarnings('ignore')
 
import os, sys, json, operator, getpass
from pathlib import Path
from datetime import datetime
 
import numpy as np
import pandas as pd
import scanpy as sc
 
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
 
import ipywidgets as widgets
 
 
sc.settings.verbosity = 3             
sc.settings.file_format_figs = 'svg'  
sc.settings.savefigs = False
 
# plt.rcParams['font.family'] = 'sans-serif'
# plt.rcParams['font.sans-serif'] = 'Arial'
# plt.rc('font', size=14)
 
home = str(Path.home())
user = getpass.getuser()
 
# %watermark --iversions

# %%
basedir = os.path.join(laurenti, user, project)
sc.settings.writedir = os.path.join(basedir, 'output/objects/')
basename = '_CHIP27_'
os.makedirs(os.path.join(basedir, 'output/figures', '02.05'), exist_ok=True)
figdir = os.path.join(basedir, 'output/figures', '02.05')

# %% [markdown]
# ### Load data

# %%
# %%time
data= sc.read('20241010_CHIP27_after_lognorm_BBKNN_Kwok_Zeng')

# %% [markdown]
# ### Cell cycle scores

# %%
cell_cycle_genes = [x.strip() for x in open('/home/gm686/giovanna/CHIP27/references/regev_lab_cell_cycle_genes.txt')]

# %%
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]

# %%
# %%time
sc.tl.score_genes_cell_cycle(data, s_genes=s_genes, g2m_genes=g2m_genes)

# %%
data.uns['phase_colors'] = ['#1f77b4', '#ff7f0e', '#2ca02c']


# %%
sc.settings._vector_friendly = True


# %%
sc.pl.umap(data, color='phase', size=10)

# %%
figname = os.path.join(figdir, prefix + '_cellcyclescore.pdf')
plt.gcf().savefig(figname, dpi=200, bbox_inches='tight')

# %% [markdown]
# #### export anndata

# %%
# %%time

filename = prefix + basename + 'after_lognorm_BBKNN_Kwok_Zeng_cellcycle'
print(filename)
print('\n')

sc.write(filename, data)
