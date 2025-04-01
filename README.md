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

```python
import scanpy as sc
import celltypemark as ctm

# Load your AnnData object
adata = sc.read_h5ad("your_data.h5ad")

# Score cell types using marker genes
adata = ctm.score(adata, ctm.marker_genes)

# Annotate cell types
adata = ctm.mark(adata, ctm.marker_genes, by='leiden', save='results', plot=True)
```

## Dependencies

- scanpy

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
