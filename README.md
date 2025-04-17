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
pip install celltypemark
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

```python
# For using marker sets from enrichr, uppercase var_names
adata.var_names = adata.var_names.str.upper()
adata = ctm.score(adata, ctm.marker_genes)
adata = ctm.mark(adata, ctm.marker_genes, by='leiden', save='results', plot=True)
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
