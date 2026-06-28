"""Single-layer ONMTF baseline.

Applies Orthogonal Nonnegative Matrix Tri-Factorization to individual
adjacency matrices for single-layer community detection.

When ground_truth is provided, uses NMI for solution selection.
When omitted, uses Modularity Density.

Reference:
    M. Ortiz-Bouza and S. Aviyente, "Community Detection in Multiplex Networks
    Based on Orthogonal Nonnegative Matrix Tri-Factorization," IEEE Access,
    vol. 12, pp. 6423-6436, 2024.
"""

from __future__ import annotations

import numpy as np

from mxonmtf.metrics import modularity_density, normalized_mutual_information


def onmtf(
    A: np.ndarray,
    k: int,
    ground_truth: np.ndarray | None = None,
    eta: float = 0.5,
    runs: int = 50,
    max_iter: int = 1000,
    conv_thresh: float = 1e-3,
    rng: np.random.Generator | None = None,
) -> tuple[np.ndarray, float]:
    """Single-layer Orthogonal Nonnegative Matrix Tri-Factorization.

    Factorizes A ≈ U1 × C1 × U1' with orthogonality constraint on U1.

    Update rules:
        U1 ← U1 × ((A U1 C1') / (U1 U1' A U1 C1'))^eta
        C1 ← C1 × ((U1' A U1) / (U1' U1 C1 U1' U1))^eta

    Parameters
    ----------
    A : np.ndarray
        Square adjacency matrix (n x n).
    k : int
        Number of communities.
    ground_truth : np.ndarray, optional
        Ground truth community labels (length n). When provided, NMI is used
        for solution selection. When None, Modularity Density is used.
    eta : float
        Learning rate in (0, 1]. Default 0.5.
    runs : int
        Number of independent runs. Default 50.
    max_iter : int
        Maximum iterations per run. Default 1000.
    conv_thresh : float
        Convergence threshold. Default 1e-3.
    rng : np.random.Generator, optional
        Random number generator for reproducibility.

    Returns
    -------
    tuple[np.ndarray, float]
        (U1_best, best_metric): best factor matrix and corresponding metric.
    """
    A = np.asarray(A, dtype=float)
    n = A.shape[0]

    if A.shape[0] != A.shape[1]:
        raise ValueError("A must be a square matrix.")
    if k < 1:
        raise ValueError("k must be a positive integer.")
    if not (0 < eta <= 1):
        raise ValueError("eta must be in (0, 1].")

    if rng is None:
        rng = np.random.default_rng()

    use_nmi = ground_truth is not None

    U1_best = None
    best_metric = -np.inf

    for j in range(runs):
        # Initialize U1, C1
        U1 = rng.random((n, k))
        C1 = np.diag(rng.random(k))

        # Iterative update
        for _ in range(max_iter):
            # U1 update
            num_U1 = A @ U1 @ C1.T
            den_U1 = U1 @ U1.T @ A @ U1 @ C1.T
            den_U1 = np.maximum(den_U1, 1e-16)  # avoid division by zero
            U1_new = U1 * (num_U1 / den_U1) ** eta

            # C1 update
            num_C1 = U1.T @ A @ U1
            den_C1 = U1.T @ U1 @ C1 @ (U1.T @ U1)
            den_C1 = np.maximum(den_C1, 1e-16)
            C1_new = C1 * (num_C1 / den_C1) ** eta

            # NaN handling
            if np.all(np.isnan(U1_new)):
                U1_new = rng.random((n, k))
            if np.all(np.isnan(C1_new)):
                C1_new = rng.random((k, k))

            # Convergence check
            if (
                np.linalg.norm(U1 - U1_new) < conv_thresh
                and np.linalg.norm(C1 - C1_new) < conv_thresh
            ):
                U1 = U1_new
                C1 = C1_new
                break

            U1 = U1_new
            C1 = C1_new

        # Assign communities and compute metric
        labels = np.argmax(U1, axis=1) + 1  # 1-indexed to match MATLAB

        if use_nmi:
            metric_j = normalized_mutual_information(ground_truth, labels)
        else:
            metric_j, _ = modularity_density(k, A, labels)

        if metric_j > best_metric:
            best_metric = metric_j
            U1_best = U1.copy()

        # Early stopping for NMI mode
        if use_nmi and abs(best_metric - 1.0) < 1e-3:
            break

    return U1_best, best_metric
