import celltypemark.resource as resource
import celltypemark.plotting as pl
import scanpy as sc
from anndata import AnnData
from collections import defaultdict
import pandas as pd
import os
marker_genes = resource.load_resource()

def score(adata: AnnData, marker_genes: defaultdict) -> AnnData:
    for k, v in marker_genes.items():
        if k not in adata.obs:
            try:
                sc.tl.score_genes(adata, v, ctrl_size=len(v), score_name=k)
            except:
                next
        else:
            print(f"{k} already exists in adata.obs")
    return adata

def mark(adata: AnnData, marker_genes: defaultdict, by: str='leiden', save: str=None, plot: bool=True) -> AnnData:
    marker_exist = set(adata.obs.columns) & set(marker_genes.keys())
    marker_df = pd.concat([adata.obs.pop(x) for x in marker_exist], axis=1)
    auto_annote = marker_df.idxmax(axis=1)
    adata.obs['celltypemark'] = auto_annote
    ## verify "by" is in adata.obs
    if by not in adata.obs:
        raise ValueError(f"{by} is not in adata.obs")
    ## if by not None
    if by is not None:
        auto_annot_df = pd.crosstab(adata.obs[by], adata.obs.celltypemark, normalize='columns')*100
        auto_annote_by = auto_annot_df.idxmax(axis=1)
        adata.obs[f'celltypemark_{by}'] = adata.obs[by].map(auto_annote_by)
        adata.obs[f'celltypemark_{by}'] = pd.Categorical(adata.obs[f'celltypemark_{by}'].astype(str), categories=adata.obs[f'celltypemark_{by}'].astype(str).dropna().unique())
        if plot:
            if save is not None:
                ## create directory if not exists
                os.makedirs('celltypemark_out/', exist_ok=True)
                pl.heatmap(auto_annot_df, f'celltypemark_out/{save}_{by}_heatmap.pdf')
            else:
                pl.heatmap(auto_annot_df)
    if save:
        auto_annot_df.to_csv(f'celltypemark_out/{save}.csv')
    return adata