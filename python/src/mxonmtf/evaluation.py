"""Community quality evaluation without ground truth.

Provides per-layer and overall quality metrics for assessing community
detection results when ground truth labels are unavailable.
"""

from __future__ import annotations

import numpy as np

from mxonmtf.metrics import modularity_density


def evaluate(
    Al: list[np.ndarray],
    labels_per_layer: list[np.ndarray],
    k: np.ndarray | list[int],
) -> dict:
    """Evaluate community detection quality per layer.

    Computes Modularity Density for each layer and overall.

    Parameters
    ----------
    Al : list of np.ndarray
        Adjacency matrices per layer.
    labels_per_layer : list of np.ndarray
        Community labels per layer (1-indexed).
    k : array-like of int
        Total communities per layer.

    Returns
    -------
    dict
        Dictionary with keys:
        - 'per_layer': list of dicts with 'modularity_density' and
          'modularity_density_norm' per layer
        - 'mean_modularity_density': average across layers
        - 'cluster_sizes': list of dicts mapping community id to count
    """
    k = np.asarray(k, dtype=int)

    if len(Al) != len(labels_per_layer):
        raise ValueError("Al and labels_per_layer must have the same length.")

    per_layer = []
    cluster_sizes = []

    for l, (A, labels) in enumerate(zip(Al, labels_per_layer)):
        md, md_norm = modularity_density(int(k[l]), A, labels)
        per_layer.append({
            "modularity_density": md,
            "modularity_density_norm": md_norm,
        })

        unique, counts = np.unique(labels[labels > 0], return_counts=True)
        cluster_sizes.append({int(u): int(c) for u, c in zip(unique, counts)})

    mean_md = np.mean([p["modularity_density"] for p in per_layer])

    return {
        "per_layer": per_layer,
        "mean_modularity_density": float(mean_md),
        "cluster_sizes": cluster_sizes,
    }
