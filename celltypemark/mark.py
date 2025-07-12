import celltypemark.resource as resource
import celltypemark.plotting as pl
import scanpy as sc
from anndata import AnnData
from collections import defaultdict
import pandas as pd
import os
from scipy.stats import fisher_exact
import numpy as np
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

def fisher(deg: pd.DataFrame, marker_genes: defaultdict, top_n: int = 20) -> pd.DataFrame:
    """
    Perform Fisher's exact test to assess enrichment of top cluster markers 
    against predefined marker gene sets.
    
    Parameters:
    -----------
    deg : pd.DataFrame
        DataFrame with cluster marker genes. Expected columns: 'cluster', 'gene', 
        and optionally ranking columns like 'logfoldchanges', 'pvals_adj'
    marker_genes : defaultdict
        Dictionary of cell type marker genes from HuBMAP resource
    top_n : int
        Number of top genes per cluster to use for Fisher's test (default: 20)
    
    Returns:
    --------
    pd.DataFrame
        Results with columns: cluster, cell_type, odds_ratio, p_value, 
        overlap_genes, cluster_genes_in_celltype, total_cluster_genes, 
        celltype_genes_total
    """
    
    results = []
    
    # Get all unique genes from marker_genes for background
    all_marker_genes = set()
    for genes in marker_genes.values():
        all_marker_genes.update(genes)
    
    # Group by cluster and get top N genes
    for cluster in deg['cluster'].unique():
        cluster_data = deg[deg['cluster'] == cluster]
        
        # Sort by significance or fold change if available
        if 'pvals_adj' in cluster_data.columns:
            cluster_data = cluster_data.sort_values('pvals_adj')
        elif 'logfoldchanges' in cluster_data.columns:
            cluster_data = cluster_data.sort_values('logfoldchanges', ascending=False)
        
        # Get top N genes for this cluster
        top_genes = set(cluster_data.head(top_n)['gene'].tolist())
        
        # Test enrichment against each cell type
        for cell_type, cell_type_genes in marker_genes.items():
            cell_type_genes = set(cell_type_genes)
            
            # Create 2x2 contingency table
            # Overlap between top cluster genes and cell type genes
            overlap = top_genes & cell_type_genes & all_marker_genes
            a = len(overlap)  # top cluster genes that are also cell type markers
            
            # Top cluster genes that are not cell type markers (but in marker universe)
            b = len((top_genes & all_marker_genes) - cell_type_genes)
            
            # Cell type genes that are not in top cluster genes
            c = len(cell_type_genes - top_genes)
            
            # All other genes in marker universe
            d = len(all_marker_genes - top_genes - cell_type_genes)
            
            # Skip if no overlap or invalid contingency table
            if a == 0 or (a + b) == 0 or (a + c) == 0:
                continue
                
            # Perform Fisher's exact test
            try:
                odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')
                
                results.append({
                    'cluster': cluster,
                    'cell_type': cell_type,
                    'odds_ratio': odds_ratio,
                    'p_value': p_value,
                    'overlap_genes': len(overlap),
                    'overlap_gene_names': ';'.join(sorted(overlap)),
                    'cluster_genes_in_celltype': a,
                    'total_cluster_genes': len(top_genes & all_marker_genes),
                    'celltype_genes_total': len(cell_type_genes),
                    'contingency_a': a,
                    'contingency_b': b,
                    'contingency_c': c,
                    'contingency_d': d
                })
            except (ValueError, ZeroDivisionError):
                continue
    
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        # Adjust p-values for multiple testing
        from statsmodels.stats.multitest import multipletests
        _, results_df['p_value_adj'], _, _ = multipletests(results_df['p_value'], method='fdr_bh')
        
        # Sort by cluster and p-value
        results_df = results_df.sort_values(['cluster', 'p_value'])
    
    return results_df

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