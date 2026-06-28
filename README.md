# MX-ONMTF

Community detection in multiplex networks based on Orthogonal Nonnegative Matrix Tri-Factorization.

**Authors:** Meiby Ortiz-Bouza and Selin Aviyente\
Department of Electrical and Computer Engineering, Michigan State University, MI\
https://ieeexplore.ieee.org/abstract/document/10385068

## Overview

MX-ONMTF is a multiplex community detection method that identifies communities that are **common** across layers as well as those that are **unique** (private) to each layer. The method represents the adjacency matrix of each layer as the sum of two low-rank matrix factorizations corresponding to the common and private communities, respectively. Unlike most existing methods which require the number of communities to be pre-determined, MX-ONMTF also introduces a two-stage method to determine the number of common and private communities.

![image](https://github.com/ortizbou/MX-ONMTF/assets/92049169/7c7631c5-eb38-4673-bc7d-c7d40d746bed)

## Repository Structure

```
MX-ONMTF/
├── matlab/                        # MATLAB implementation
│   ├── MX_ONMTF.m                  # Main algorithm (unified: synthetic + real)
│   ├── ONMTF.m                     # Single-layer baseline (unified: NMI or ModDen)
│   ├── findingk.m                  # Automatic community count detection
│   ├── assigncomm.m                # Community assignment post-processing
│   ├── patchmult.m                 # Common community detection per layer
│   ├── examples/
│   │   ├── demo_synthetic.m         # Example with synthetic data + ground truth
│   │   └── demo_real.m              # Example with real data (no ground truth)
│   └── helpers/
│       ├── EigGap.m                 # Eigenvalue gap method
│       ├── ModDen.m                 # Modularity Density metric
│       ├── getNMI.m                 # Normalized Mutual Information
│       ├── normadj.m                # Normalized adjacency matrix
│       ├── rnorm.m                  # Row normalization
│       └── vec.m                    # Vectorization utility
├── python/                        # Python implementation
│   ├── pyproject.toml               # Package metadata and dependencies
│   ├── src/mxonmtf/
│   │   ├── factorization.py         # Core MX-ONMTF algorithm
│   │   ├── single_layer.py          # Single-layer ONMTF baseline
│   │   ├── community.py             # Community assignment post-processing
│   │   ├── parameter_selection.py   # Automatic parameter selection (findingk)
│   │   ├── metrics.py               # NMI and Modularity Density
│   │   └── utils.py                 # Graph utilities
│   ├── tests/                       # pytest test suite
│   └── examples/
│       └── tutorial.py              # End-to-end tutorial script
└── data/                          # Shared example data
```

## Installation

### Python

```bash
cd python
pip install -e .          # basic install
pip install -e ".[dev]"   # with dev tools (pytest, ruff)
```

Requires Python >= 3.9, NumPy >= 1.22, SciPy >= 1.8.

### MATLAB

No installation needed. Add `matlab/` and `matlab/helpers/` to your MATLAB path.

## Quick Start

### Python

```python
import numpy as np
from mxonmtf import mx_onmtf, findingk

# Al is a list of adjacency matrices (one per layer)
# Step 1: Find community counts (or set manually)
kc, k, kpl = findingk(Al)

# Step 2: Run MX-ONMTF
# Real data (trace minimization)
result = mx_onmtf(Al, k=k, kc=kc, kpl=kpl)

# Synthetic data (NMI selection)
result = mx_onmtf(Al, k=k, kc=kc, kpl=kpl, ground_truth=labels)

# Access results
result.H_best       # common community factor matrices
result.Hl_best      # private community factor matrices
result.Clusters     # community labels
result.averagedNMI  # NMI (when ground_truth provided)
```

See `python/examples/tutorial.py` for a complete walkthrough.

### MATLAB

```matlab
% Step 1: Find community counts
[kc, k, kpl] = findingk(Al);

% Step 2: Run MX-ONMTF
results = MX_ONMTF(Al, k, kc, kpl);                                        % real data
results = MX_ONMTF(Al, k, kc, kpl, 'ground_truth', GTl, 'realizations', 100); % synthetic
results = MX_ONMTF(Al, k, kc, kpl, 'eta', 0.3, 'max_iter', 5000);          % custom params
```

See `matlab/examples/demo_synthetic.m` and `matlab/examples/demo_real.m` for complete walkthroughs.

## Testing

```bash
cd python
pytest -v    # 49 tests
```

## Citation

```bibtex
@article{ortiz2024community,
  title={Community Detection in Multiplex Networks Based on Orthogonal Nonnegative Matrix Tri-Factorization},
  author={Ortiz-Bouza, Meiby and Aviyente, Selin},
  journal={IEEE Access},
  volume={12},
  pages={6423--6436},
  year={2024},
  publisher={IEEE}
}
```
