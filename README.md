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
├── matlab/                 # MATLAB implementation
│   ├── MX_ONMTF.m          # Main algorithm (unified: synthetic + real)
│   ├── ONMTF.m             # Single-layer baseline (unified: NMI or ModDen)
│   ├── findingk.m          # Automatic community count detection
│   ├── assigncomm.m        # Community assignment post-processing
│   ├── patchmult.m         # Common community detection per layer
│   ├── Steps.mlx           # Example workflow (MATLAB Live Script)
│   └── helpers/
│       ├── EigGap.m         # Eigenvalue gap method
│       ├── ModDen.m         # Modularity Density metric
│       ├── getNMI.m         # Normalized Mutual Information
│       ├── normadj.m        # Normalized adjacency matrix
│       ├── rnorm.m          # Row normalization
│       └── vec.m            # Vectorization utility
├── python/                 # Python implementation (coming soon)
└── data/                   # Shared example data
```

## Usage (MATLAB)

**Step 1.** Use `findingk.m` to find the number of common communities and the number of private and total communities per layer. If the number of communities is known, this step can be skipped.

```matlab
[kc, k, kpl] = findingk(Al);
```

**Step 2.** Run MX-ONMTF:

```matlab
% For real data (no ground truth) — uses trace minimization
results = MX_ONMTF(Al, k, kc, kpl);

% For synthetic data with ground truth — uses NMI for solution selection
results = MX_ONMTF(Alr, k, kc, kpl, 'ground_truth', GTlr, 'realizations', 100);

% Override any parameter
results = MX_ONMTF(Al, k, kc, kpl, 'eta', 0.3, 'max_iter', 5000);
```

See `matlab/Steps.mlx` for a complete walkthrough.

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
