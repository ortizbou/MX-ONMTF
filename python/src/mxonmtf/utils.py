"""Graph utility functions.

Normalized adjacency matrix computation and row normalization helpers.
"""

from __future__ import annotations

import numpy as np
from scipy import sparse


def normalized_adjacency(A: np.ndarray) -> np.ndarray:
    """Compute the normalized adjacency matrix D^{-1/2} L D^{-1/2}.

    Computes the normalized Laplacian-style matrix where L = D - A and D is
    the degree matrix. Self-loops are handled separately.

    Based on the SGWT toolbox (David K. Hammond, 2010).

    Parameters
    ----------
    A : np.ndarray
        Square adjacency matrix (n x n).

    Returns
    -------
    np.ndarray
        Normalized adjacency matrix (n x n), dense.
    """
    A = np.asarray(A, dtype=float)
    n = A.shape[0]
    degrees = A.sum(axis=1)
    diag_w = np.diag(A).copy()

    # Extract non-diagonal entries
    A_sparse = sparse.coo_matrix(A)
    mask = A_sparse.row != A_sparse.col
    ni = A_sparse.row[mask]
    nj = A_sparse.col[mask]
    w = A_sparse.data[mask]

    # Diagonal entries of normalized matrix
    dL = np.zeros(n)
    nonzero_deg = degrees != 0
    dL[nonzero_deg] = diag_w[nonzero_deg] / degrees[nonzero_deg]

    # Non-diagonal entries
    ndL = w / np.sqrt(degrees[ni] * degrees[nj])

    # Build sparse matrix and convert to dense
    rows = np.concatenate([ni, np.arange(n)])
    cols = np.concatenate([nj, np.arange(n)])
    data = np.concatenate([ndL, dL])
    L = sparse.coo_matrix((data, (rows, cols)), shape=(n, n))
    return L.toarray()


def row_normalize(M: np.ndarray) -> np.ndarray:
    """Normalize each row of a matrix by its Euclidean norm.

    Rows with zero norm are left unchanged.

    Parameters
    ----------
    M : np.ndarray
        Input matrix.

    Returns
    -------
    np.ndarray
        Row-normalized matrix.
    """
    M = np.asarray(M, dtype=float)
    norms = np.linalg.norm(M, axis=1, keepdims=True)
    norms[norms == 0] = 1.0  # avoid division by zero
    return M / norms
