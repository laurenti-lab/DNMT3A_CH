# %% [markdown]
# # 01.02_Cell and genes Filtering
#
# Author: GM <br>

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

from scipy.stats import median_abs_deviation
 
 
sc.settings.verbosity = 3             
sc.settings.file_format_figs = 'svg'  # set this to 'svg' (notebook) or 'pdf' (files) if you want vector graphics
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

# %%
basename = '_CHIP27_'

# %%
figdir = os.path.join(basedir, 'output/figures', '01.02')

if not os.path.exists(figdir):
    os.makedirs(figdir)
    print(figdir, "created successfully")
else:
    print(figdir, "already exists")


# %% [markdown]
# ---

# %% [markdown]
# ### Import data

# %%
data = sc.read('20240918_CHIP27_filtered_gene_bc_expression_minus_797_putative_doublets')

# %% [markdown]
# #### annotate data:

# %%
genotype_dict={'SITTA2': 'MUT', 'SITTD2': 'MUT', 'SITTC2': 'MUT', 'SITTB2': 'WT', 'SITTE2': 'WT', 'SITTF2': 'WT'}
treatment_dict={'SITTA2':'CTR', 'SITTB2':'CTR', 'SITTD2':'CTR', 'SITTE2':'CTR', 'SITTC2':'IFNy', 'SITTF2':'IFNy'}
sort_dict={'SITTA2':'not_enriched', 'SITTB2':'not_enriched', 'SITTD2':'enriched', 'SITTE2':'enriched', 'SITTC2':'enriched', 'SITTF2':'enriched'}

# %%
data.obs['DNMT3A_genotype'] = data.obs['library'].map(genotype_dict)
data.obs['treatment'] = data.obs['library'].map(treatment_dict)
data.obs['FACS'] = data.obs['library'].map(sort_dict)

# %% [markdown]
# ---

# %% [markdown]
# ### Calculate QC metrix

# %%
data.var["mt"] = data.var_names.str.startswith("MT-")
data.var["ribo"] = data.var_names.str.startswith(("RPS", "RPL"))
#data.var["ercc"] = adata.var_names.str.contains("ERCC")
data.var["hb"] = data.var_names.str.contains(("^HB[^(P)]"))

# %%
sc.pp.calculate_qc_metrics(
    data, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=None, log1p=True)


# %% [markdown]
# #### General overview QC metrics and top 20 genes

# %%
sc.pl.violin(data, ['n_genes_by_counts','total_counts','log1p_total_counts', 'pct_counts_mt', 'pct_counts_hb', 'pct_counts_ribo'],
             jitter=0.4, multi_panel=True)

# %%
fig, ax = plt.subplots(1,1, figsize=(8,8),dpi=150 )
sc.pl.highest_expr_genes(data, n_top=20, ax=ax, show=False)

# %% [markdown]
# ---

# %% [markdown]
# ### filtering step 1: genes at least in 3 cells and cells with at least 500 genes

# %% [markdown]
# #### genes at least in 3 cells

# %%
sc.pp.filter_genes(data, min_cells=3)

# %% [markdown]
# #### cells with at least 500 genes

# %%
sc.pp.filter_cells(data, min_genes=500)

# %% [markdown]
# ### filtering step 2: total counts outliers
#

# %%
p2 = sns.displot(data.obs["total_counts"], bins=100, kde=False)
ax=p2.ax
ax.set_ylim(0, 6000)

# %%
p3 = sns.displot(data.obs["log1p_total_counts"], bins=100, kde=False)
ax=p3.ax

# %% [markdown]
# **MADs** (median absolute deviations) method <br>
# see :  https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html

# %%
# upper limit = UL
log_total_counts_median = np.median(data.obs['log1p_total_counts'])
log_total_counts_mad = np.median(np.abs(data.obs['log1p_total_counts'] - log_total_counts_median))

log_total_counts_4mad = log_total_counts_median - (log_total_counts_mad * 4)
log_total_counts_4mad_right = log_total_counts_median + (log_total_counts_mad * 4)
print('4 MAD (log1p):', log_total_counts_4mad, 'UL:', log_total_counts_4mad_right)
print('4 MAD (transformed):', np.expm1(log_total_counts_4mad),'UL:', np.expm1(log_total_counts_4mad_right))

log_total_counts_35mad = log_total_counts_median - (log_total_counts_mad * 3.5)
log_total_counts_35mad_right = log_total_counts_median + (log_total_counts_mad * 3.5)
print('3.5 MAD (log1p):', log_total_counts_35mad, 'UL:', log_total_counts_35mad_right)
print('3.5 MAD (transformed):', np.expm1(log_total_counts_35mad), 'UL:', np.expm1(log_total_counts_35mad_right))

log_total_counts_3mad = log_total_counts_median - (log_total_counts_mad * 3)
log_total_counts_3mad_right = log_total_counts_median + (log_total_counts_mad * 3)
print('3 MAD (log1p):', log_total_counts_3mad, 'UL:', log_total_counts_3mad_right)
print('3 MAD (transformed):', np.expm1(log_total_counts_3mad), 'UL:', np.expm1(log_total_counts_3mad_right))


# %%
ax.axvline(log_total_counts_4mad, color='orange', linewidth=1)
ax.axvline(log_total_counts_4mad_right, color='orange', linewidth=1)

#ax.axvline(log_total_counts_3mad, color='purple', linewidth=1)
#ax.axvline(log_total_counts_3mad_right, color='purple', linewidth=1)

ax.axvline(log_total_counts_35mad, color='green', linewidth=1)
ax.axvline(log_total_counts_35mad_right, color='green', linewidth=1)

plt.draw()

# %% [markdown]
# **THRESHOLDS SET MANUALLY** <br>
# MADs only visualised to have an idea. Then threshold set manually:we have clearly two distributions of counts, and doing the MAD method skews it towards the most abundant population 
#
#
# upper theshold picked: 85000 counts. (I calculated 4MADs and went a bit more permissive to include more cells)
# lower threshold picked: 1300 counts. (It is 3.5MADS rounded up to the closest thousand)

# %%

ax.axvline(np.log1p(85000), color='red', linewidth=2) #this is the best one, it removes nicely the end of the tail
ax.axvline(np.log1p(1300), color='red', linewidth=2) 
plt.draw()

# %% [markdown]
# **check that the thresholds are equally applicable to all of our libraries**

# %%
samples = data.obs.groupby('library')

#initialise list
plots = []

for sample, sample_data in samples:
    
    p = sns.displot(sample_data['log1p_total_counts'], bins=100, kde=False, log_scale=False)
    
    # Access the Matplotlib axis object
    ax = p.ax
    
    ax.set(title=f'library {sample} - log1p Total Counts Distribution')
    
    # ax.axvline(log_total_counts_4mad, color='orange', linewidth=1)
    # ax.axvline(log_total_counts_4mad_right, color='orange', linewidth=1)
    
    # ax.axvline(log_total_counts_3mad, color='green', linewidth=1)
    # ax.axvline(log_total_counts_3mad_right, color='green', linewidth=1)
    
    # ax.axvline(log_total_counts_35mad, color='red', linewidth=1)
    # ax.axvline(log_total_counts_35mad_right, color='red', linewidth=1)

    ax.axvline(np.log1p(85000), color='red', linewidth=2)
    ax.axvline(np.log1p(1300), color='red', linewidth=2)
    
    p.fig.set_size_inches(8, 8)
    plt.show()

# %% [markdown]
# #### filtering:

# %%
LT=1300
UT=85000

# %%
percentage_saved= (
    data[(data.obs.total_counts > LT) & (data.obs.total_counts < UT)].obs.groupby(['library']).size() / data.obs.groupby(['library']).size()
) * 100

percentage_df = percentage_saved.reset_index(name='Percentage_COUNTS_passed')

print(percentage_df)

# %%
# %%time
data=data[(data.obs.total_counts > LT) & (data.obs.total_counts < UT)].copy()

# %% [markdown]
# ### filtering step 3: mitochondrial percentage
#

# %%
sc.pl.violin(data, ['pct_counts_mt', 'pct_counts_hb', 'pct_counts_ribo'],
             jitter=0.4, multi_panel=True)

# %%
#make colour palette
norm = plt.Normalize(data.obs['pct_counts_mt'].min(), data.obs['pct_counts_mt'].max())
sm = plt.cm.ScalarMappable(cmap=sns.cubehelix_palette(as_cmap=True), norm=norm)
sm.set_array([])

# %%
fig, ax = plt.subplots(1,2, figsize=(14,6) )
sns.scatterplot(x='total_counts', y='n_genes_by_counts', 
                data=data.obs[['total_counts', 'n_genes_by_counts', 'pct_counts_mt']] , 
                hue='pct_counts_mt', 
                ec=None, ax=ax[0] )

ax[1].scatter(data.obs['total_counts'], data.obs['pct_counts_mt'], c='lightgrey', s=1.5)

# %%
for a in ax[0].collections: 
    a.set_sizes([5])

# Remove the legend and add a colorbar
#ax[0].get_legend().remove()
fig.colorbar(sm, ax=ax[0], shrink=0.95)

ax[0].axes.set_ylabel('Number of genes')
ax[0].axhline(500, color='green', linewidth=0.75)
ax[0].axvline(1300, color='red', linewidth=0.75)#remember to change this if you change totcoutns threshold above

ax[0].axes.set_xlabel('')
ax[1].axes.set_ylabel('percentage of mitochondrial genes')

ax[1].set_title('Fraction of mitochondrial genes')
fig.text(0.5, 0.01, 'total counts', ha='center', va='center')

# %%
ax[1].axhline(20, color='yellow', linewidth=0.75)
ax[1].axhline(15, color='orange', linewidth=0.75)
ax[1].axhline(10, color='red', linewidth=0.75)


plt.text(60000, 10, sum(data.obs['pct_counts_mt'] >= 10), fontsize=8, va='center',  backgroundcolor='w')
plt.text(60000, 15, sum(data.obs['pct_counts_mt'] >= 15), fontsize=8, va='center',  backgroundcolor='w')
plt.text(60000, 20, sum(data.obs['pct_counts_mt'] >= 20), fontsize=8, va='center',  backgroundcolor='w')

# %%
plt.gcf()

# %% [markdown]
# **COMMENT** <br>
# 15 % works fine

# %%
print('percentage discarded with MT cut at 15')
data[ data.obs['pct_counts_mt'] >= 15 ].obs.groupby(['library']).size()/data.obs.groupby(['library']).size()*100

# %%
percentage_saved= (
    data[data.obs.pct_counts_mt < 15].obs.groupby(['library']).size() / data.obs.groupby(['library']).size()
) * 100

percentage_df = percentage_saved.reset_index(name='Percentage_MITO_passed')

print(percentage_df)

# %% [markdown]
# #### filtering:

# %%
# %%time
data=data[data.obs.pct_counts_mt < 15].copy()

# %% [markdown]
# ### export filtered anndata

# %%
# %%time

filename = prefix + basename + 'postQC_beforenorm'
print(filename)
print('\n')

sc.write(filename, data)
