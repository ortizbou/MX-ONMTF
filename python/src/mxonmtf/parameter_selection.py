"""Automatic parameter selection.

Determines the number of common communities (kc), total communities per
layer (k), and private communities per layer (kpl) using eigenvalue gap
analysis and hierarchical clustering.

Reference:
    Algorithm 1 in M. Ortiz-Bouza and S. Aviyente, "Community Detection in
    Multiplex Networks Based on Orthogonal Nonnegative Matrix Tri-Factorization,"
    IEEE Access, vol. 12, pp. 6423-6436, 2024.
"""

from __future__ import annotations

import numpy as np
from scipy.cluster.hierarchy import linkage

from mxonmtf.utils import normalized_adjacency


def eiggap(
    A: np.ndarray,
    eiggap_thresh: float = 0.95,
    seed: int = 100,
) -> int:
    """Determine optimal number of communities using eigenvalue gap method.

    Generates a random null-model graph with the same density, computes the
    eigenvalue gap threshold from its normalized Laplacian, then finds
    eigenvalue gaps in the actual network that exceed this threshold.

    Parameters
    ----------
    A : np.ndarray
        Square adjacency matrix (n x n).
    eiggap_thresh : float
        Fraction of max null-model gap used as threshold. Default 0.95.
    seed : int
        Random seed for null model generation. Default 100.

    Returns
    -------
    int
        Estimated number of communities.
    """
    A = np.asarray(A, dtype=float)
    n = A.shape[0]
    density = 2.0 * A.sum() / (n * (n - 1))

    # Generate random null-model graph with same density
    rng = np.random.RandomState(seed)
    G = (rng.rand(n, n) < density).astype(float)
    G = np.triu(G, k=1)
    G = G + G.T

    # Normalized Laplacian of null model
    NM = np.eye(n) - normalized_adjacency(G)
    dn = np.abs(np.linalg.eigvalsh(NM))
    dn = np.sort(dn)[::-1]

    dif = np.zeros(n - 1)
    for i in range(1, n - 1):  # skip index 0 (set to 0 like MATLAB)
        dif[i] = dn[i] - dn[i + 1]
    delta = eiggap_thresh * np.max(dif)

    # Eigenvalue decomposition of actual network
    d = np.abs(np.linalg.eigvalsh(A.astype(float)))
    d = np.sort(d)[::-1]

    gap = np.zeros(n)
    for i in range(n - 1):
        if d[i] > 0 and d[i + 1] > 0:
            gap[i] = abs(d[i]) - abs(d[i + 1])

    # Find last eigenvalue gap exceeding threshold
    k = 1
    for j in range(n):
        if gap[j] > delta:
            k = j + 1  # 1-indexed community count

    return k


def findingk(
    Al: list[np.ndarray],
    eiggap_thresh: float = 0.95,
    onmtf_fn=None,
    **onmtf_kwargs,
) -> tuple[int, np.ndarray, np.ndarray]:
    """Find the number of common and private communities per layer.

    Implements Algorithm 1 from the paper: for each layer, estimate the number
    of communities via eigenvalue gap, run single-layer ONMTF to get factor
    matrices, then use hierarchical clustering to identify which communities
    are common across layers.

    Parameters
    ----------
    Al : list of np.ndarray
        List of adjacency matrices, one per layer. All must be square and
        same size (n x n).
    eiggap_thresh : float
        Threshold for eigenvalue gap method. Default 0.95.
    onmtf_fn : callable, optional
        Single-layer ONMTF function with signature ``(A, k, **kwargs) -> (U, metric)``.
        If None, imports from ``mxonmtf.single_layer``.
    **onmtf_kwargs
        Additional keyword arguments passed to ``onmtf_fn``.

    Returns
    -------
    tuple[int, np.ndarray, np.ndarray]
        (kc, k, kpl) where:
        - kc: number of common communities
        - k: array of total communities per layer
        - kpl: array of private communities per layer
    """
    if onmtf_fn is None:
        from mxonmtf.single_layer import onmtf

        onmtf_fn = onmtf

    L = len(Al)
    n = Al[0].shape[0]

    # Validate inputs
    for idx, A in enumerate(Al):
        if A.shape[0] != A.shape[1]:
            raise ValueError(f"Adjacency matrix for layer {idx} is not square.")
        if A.shape[0] != n:
            raise ValueError("All layers must have the same number of nodes.")

    # Step 1: Find k_l and U_l for each layer
    k = np.zeros(L, dtype=int)
    U = []
    for l in range(L):
        k[l] = eiggap(Al[l], eiggap_thresh=eiggap_thresh)
        U_l, _ = onmtf_fn(Al[l], k[l], **onmtf_kwargs)
        U.append(U_l)

    # Step 2: Hierarchical clustering on concatenated factor matrices
    X = np.vstack(U)
    Z = linkage(X)

    # Step 3: Count common communities
    kc = 0
    n_rows = X.shape[0]
    cut = n_rows // 2 - 1

    for i in range(1, n_rows - 1):  # i=1 corresponds to MATLAB i=2
        if max(Z[i - 1, 0], Z[i - 1, 1]) <= n_rows:
            kc += 1

        d_val = (Z[i, 2] - Z[i - 1, 2]) / Z[i - 1, 2] if Z[i - 1, 2] != 0 else 0
        if d_val >= 0.5:
            cut = i - 1
            break

    # Step 4: Compute private communities per layer
    kpl = np.zeros(L, dtype=int)
    k_cumsum = np.cumsum(k)
    for l in range(L):
        start = 1 + (k_cumsum[l - 1] if l > 0 else 0)
        end = k_cumsum[l]

        # Count communities from this layer that appear in the clustering cut
        Z_indices = Z[:cut + 1, :2]
        count = np.sum((Z_indices >= start) & (Z_indices <= end))
        kpl[l] = k[l] - count

        if (k[l] - kc) > kpl[l]:
            kpl[l] = k[l] - kc

    return kc, k, kpl
