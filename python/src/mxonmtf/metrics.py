"""Community detection quality metrics.

Normalized Mutual Information (NMI) and Modularity Density for evaluating
community detection results.
"""

from __future__ import annotations

import numpy as np


def normalized_mutual_information(A: np.ndarray, B: np.ndarray) -> float:
    """Compute the Normalized Mutual Information between two partitions.

    NMI measures similarity between two community partitions. If A and B are
    identical, NMI = 1. If A and B are independent, NMI approaches 0.

    Parameters
    ----------
    A : np.ndarray
        Community labels for partition A (length n).
    B : np.ndarray
        Community labels for partition B (length n).

    Returns
    -------
    float
        Normalized Mutual Information in [0, 1].

    References
    ----------
    [1] Kuncheva & Hadjitodorov, 2004, IEEE Intern. Conf. on Syst. Man and Cybern.
    [2] Alexander-Bloch et al., 2012, NeuroImage.
    """
    A = np.asarray(A).ravel()
    B = np.asarray(B).ravel()

    if len(A) != len(B):
        raise ValueError("Partitions A and B must have the same length.")

    N = len(A)
    Ca = np.unique(A)
    Cb = np.unique(B)

    numerator = 0.0
    D1 = 0.0
    D2 = 0.0

    for i in Ca:
        N_idot = np.sum(A == i)
        D1 += N_idot * np.log(N_idot / N)

        for j in Cb:
            N_dotj = np.sum(B == j)
            N_ij = np.sum((A == i) & (B == j))

            if N_ij > 0:
                numerator += N_ij * np.log((N_ij * N) / (N_idot * N_dotj))

    for j in Cb:
        N_dotj = np.sum(B == j)
        D2 += N_dotj * np.log(N_dotj / N)

    denominator = D1 + D2

    if denominator == 0:
        return 1.0  # both partitions are trivial (single community)

    return -2.0 * numerator / denominator


def modularity_density(
    k: int, A: np.ndarray, labels: np.ndarray
) -> tuple[float, float]:
    """Compute Modularity Density and Normalized Modularity Density.

    Modularity Density quantifies community structure quality when ground truth
    is unavailable. It balances internal and external edge density.

    Parameters
    ----------
    k : int
        Number of communities.
    A : np.ndarray
        Square adjacency matrix (n x n).
    labels : np.ndarray
        Community label for each node (length n).

    Returns
    -------
    tuple[float, float]
        (ModularityDensity, ModularityDensityNorm).
    """
    A = np.asarray(A, dtype=float)
    labels = np.asarray(labels).ravel()

    mod_density = 0.0
    mod_density_norm = 0.0

    for c in range(1, k + 1):
        members = np.where(labels == c)[0]
        if len(members) == 0:
            continue

        abs_v = len(members)
        non_members = np.where(labels != c)[0]

        # Internal edges
        LV = A[np.ix_(members, members)].sum()

        # External edges
        if len(non_members) == 0:
            LVn = 0.0
            abs_vn = 1  # avoid division by zero
        else:
            LVn = A[np.ix_(members, non_members)].sum()
            abs_vn = len(non_members)

        mod_density += (LV - LVn) / abs_v

        # Normalized version
        internal_denom = abs_v * (abs_v - 1) if abs_v > 1 else 1
        mod_density_norm += LV / internal_denom - LVn / (abs_v * abs_vn)

    return mod_density, mod_density_norm
