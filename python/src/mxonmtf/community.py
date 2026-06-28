"""Community assignment post-processing.

Assigns nodes to common or private communities based on factor matrices,
including detection of which common communities exist per layer.

Reference:
    Algorithm 3 in M. Ortiz-Bouza and S. Aviyente, "Community Detection in
    Multiplex Networks Based on Orthogonal Nonnegative Matrix Tri-Factorization,"
    IEEE Access, vol. 12, pp. 6423-6436, 2024.
"""

from __future__ import annotations

import numpy as np


def patchmult(
    A: np.ndarray, H: np.ndarray, kc: int
) -> tuple[np.ndarray, np.ndarray]:
    """Determine which common communities are present in a layer.

    Computes the ratio of within-community to between-community density for
    each common community, and sorts by this ratio to identify which common
    communities are least present (lowest ratio).

    Parameters
    ----------
    A : np.ndarray
        Adjacency matrix (n x n) for one layer.
    H : np.ndarray
        Common community factor matrix (n x kc).
    kc : int
        Number of common communities.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        (pos, q): sorted positions and density ratios.
    """
    n = A.shape[1]
    I1 = np.argmax(H, axis=1)  # community assignment per node

    # Binary membership matrix
    h = np.zeros_like(H)
    for i in range(kc):
        h[I1 == i, i] = 1.0

    di = np.zeros(kc)
    d0 = np.zeros(kc)

    for k_idx in range(kc):
        members = np.where(h[:, k_idx] == 1)[0]
        edges = len(members)

        if edges <= 1:
            di[k_idx] = 0.0
            d0[k_idx] = 1.0  # avoid division by zero
            continue

        # Within-community density
        mask = np.outer(h[:, k_idx], h[:, k_idx])
        di[k_idx] = np.sum(A * mask) / (edges * (edges - 1))

        # Between-community edges
        total_edges_from_members = A[members, :].sum()
        within_edges = np.sum(A * mask)
        between_edges = total_edges_from_members - within_edges

        if edges * (n - edges) > 0:
            d0[k_idx] = between_edges / (edges * (n - edges))
        else:
            d0[k_idx] = 1.0

    # Ratio and sort
    q = np.where(d0 != 0, di / d0, 0.0)
    pos = np.argsort(q)  # ascending order (least present first)

    return pos, q


def assign_communities(
    Al: list[np.ndarray],
    H: np.ndarray,
    Hl: list[np.ndarray],
    kc: int,
    k: np.ndarray,
    kpl: np.ndarray,
    L: int,
    mode: str = "same",
    perc: float = 0.8,
) -> tuple[list[np.ndarray], np.ndarray]:
    """Assign nodes to common or private communities.

    Parameters
    ----------
    Al : list of np.ndarray
        Adjacency matrices per layer.
    H : np.ndarray
        Common community factor matrix (n x kc).
    Hl : list of np.ndarray
        Private community factor matrices per layer (each n x kpl[l]).
    kc : int
        Number of common communities.
    k : array-like of int
        Total communities per layer.
    kpl : array-like of int
        Private communities per layer.
    L : int
        Number of layers.
    mode : str
        Assignment mode: 'flex' or 'same'. Default 'same'.
    perc : float
        Threshold for common community assignment in 'same' mode. Default 0.8.

    Returns
    -------
    tuple[list[np.ndarray], np.ndarray]
        (Il, ClustersSupra): per-layer labels and supralayer labels.
    """
    k = np.asarray(k, dtype=int)
    kpl = np.asarray(kpl, dtype=int)
    n = H.shape[0]

    # Determine which common communities exist in each layer
    Hc = []
    Il = [np.zeros(n, dtype=int) for _ in range(L)]

    for l in range(L):
        Hc_l = H.copy()
        m = kc - (k[l] - kpl[l])
        if m > 0:
            pos, _ = patchmult(Al[l], H, kc)
            # Zero out least-present common communities
            Hc_l[:, pos[:m]] = 0.0
        Hc.append(Hc_l)

    if mode == "flex":
        for l in range(L):
            combined = np.hstack([Hc[l], Hl[l]])
            Il[l] = np.argmax(combined, axis=1) + 1  # 1-indexed

        for l in range(L):
            for m_idx in range(1, kpl[l] + 1):
                mask = Il[l] == (kc + m_idx)
                Il[l][mask] += int(np.sum(k[:l]))

        clusters_supra = np.concatenate(Il)

    elif mode == "same":
        for node in range(n):
            for g in range(kc):
                # Check which layers have this common community
                H0l = np.column_stack([Hc[l][:, g] for l in range(L)])
                layers_present = np.where(H0l.sum(axis=0) != 0)[0]

                if len(layers_present) == 0:
                    continue

                Hall = np.hstack([Hl[l] for l in layers_present])

                for l in layers_present:
                    count_greater = np.sum(Hc[l][node, g] > Hall[node, :])
                    if count_greater > perc * Hall.shape[1]:
                        Il[l][node] = g + 1  # 1-indexed

        # Assign unassigned nodes to private communities
        for l in range(L):
            unassigned = np.where(Il[l] == 0)[0]
            if len(unassigned) > 0 and Hl[l].shape[1] > 0:
                private_labels = np.argmax(Hl[l][unassigned, :], axis=1) + 1
                Il[l][unassigned] = private_labels + kc + int(np.sum(k[:l]))

        clusters_supra = np.concatenate(Il)

    else:
        raise ValueError("mode must be 'flex' or 'same'.")

    return Il, clusters_supra
