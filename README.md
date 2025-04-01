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

## Citation

If you use CellTypeMark in your research, please cite:

```
@software{celltypemark2024,
  author = {brainfo},
  title = {CellTypeMark: A Python package for cell type annotation in single-cell RNA sequencing data},
  year = {2024},
  publisher = {GitHub},
  url = {https://github.com/yourusername/celltypemark}
}
```

---

Copyright (c) 2024

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. 