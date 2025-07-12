# CellTypeMark

CellTypeMark is a Python package for cell type annotation in single-cell RNA sequencing data using marker genes. It provides tools for scoring and annotating cell types based on predefined marker gene sets.

## Features

- Score cell types using marker genes (using scanpy score_genes)
- Automatic cell type annotation based on marker gene expression
- Visualization of cell type annotations
- Integration with Scanpy and AnnData objects
- Support for custom marker gene sets

## Installation

```bash
uv pip install git+https://github.com/brainfo/celltypemark.git
```

## Usage

### Load your AnnData object
```python
import scanpy as sc
import celltypemark as ctm
adata = sc.read_h5ad("your_data.h5ad")
```

### Score cell types using default marker genes
```python
adata = ctm.score(adata, ctm.marker_genes)
```

### Score cell types using specified marker genes

- Example, any existing local file

```python
from pathlib import Path
resource_path = Path("your gene marker file")
marker_genes = ctm.load_resource(resource_path)
```

- Example, if your resource file is not yet on disk, download from resource_url

```python
from pathlib import Path
resource_path = Path("your gene marker file")
marker_genes = ctm.load_resource(resource_path, resource_url='https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Tabula_Muris')
```

### Annotate cell types

1. by scoring marker gene expression activities

```python
# For using marker sets from enrichr, uppercase var_names
adata.var_names = adata.var_names.str.upper()
adata = ctm.score(adata, ctm.marker_genes)
adata = ctm.mark(adata, ctm.marker_genes, by='leiden', save='results', plot=True)
```

2. by fisher's exact test for marker gene overlaps

Perform Fisher's exact test to assess enrichment of cluster marker genes against predefined cell type markers:

```python
import pandas as pd

# Load cluster marker genes (Excel file with 'cluster' and 'gene' columns)
degs = pd.read_excel('cluster_markers.xlsx')

# Perform Fisher's exact test using default HuBMAP markers
fisher_results = ctm.fisher(degs, marker_genes, top_n=20)

significant_results = fisher_results[fisher_results['p_value_adj'] < 0.05]

print(significant_results[['cluster', 'cell_type', 'odds_ratio', 'p_value_adj', 'overlap_genes']])

     cluster   cell_type  odds_ratio   p_value_adj  overlap_genes
0          0    B.Bcells  268.516854  8.313690e-17             14
6          0      Plasma  257.827027  1.290801e-16             14
1          0   Cycling B   24.675651  2.044974e-07              8
5          0  I.Lymphoid   17.833333  5.915358e-06              7
2          0  Follicular   14.984749  5.096367e-05              6
..       ...         ...         ...           ...            ...
148        9  CD69- Mast  142.581081  1.433895e-16             12
149        9        ILCs    8.585106  1.444172e-02              3
156       10      Plasma         inf  4.033077e-14             11
152       10    B.Bcells   49.449275  1.942378e-08              8
153       10   Cycling B    7.847561  2.107963e-02              3

[91 rows x 5 columns]
```

### Output

Add in the adata.obs:
- the scores of the keys in the gene sets for each observation
- the predicted key for each observation
- (if by) the predicted key for each by group, e.g., leiden
- (if save) save the scores of the keys for each by group in a txt file and the heatmap of the data under cellmarkoutput/

## Dependencies

- scanpy

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
