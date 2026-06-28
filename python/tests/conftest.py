"""Shared test fixtures for MX-ONMTF tests."""

import numpy as np
import pytest


@pytest.fixture
def rng():
    """Seeded random number generator for reproducible tests."""
    return np.random.default_rng(42)


@pytest.fixture
def planted_partition_graph(rng):
    """Factory fixture for generating planted partition model graphs.

    Returns a function that creates a symmetric adjacency matrix given
    community labels and edge probabilities.
    """

    def _make(labels: np.ndarray, p_in: float = 0.4, p_out: float = 0.05) -> np.ndarray:
        n = len(labels)
        A = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                p = p_in if labels[i] == labels[j] else p_out
                if rng.random() < p:
                    A[i, j] = 1.0
                    A[j, i] = 1.0
        return A

    return _make
