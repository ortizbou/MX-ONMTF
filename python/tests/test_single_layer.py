"""Tests for mxonmtf.single_layer."""

import numpy as np
import pytest
from mxonmtf.single_layer import onmtf


class TestONMTF:
    def test_returns_factor_and_metric(self, planted_partition_graph):
        """Should return (U1_best, best_metric) tuple."""
        labels = np.array([1] * 20 + [2] * 20)
        A = planted_partition_graph(labels, p_in=0.5, p_out=0.02)
        U1, metric = onmtf(A, k=2, runs=3, max_iter=100, rng=np.random.default_rng(42))
        assert U1.shape == (40, 2)
        assert isinstance(metric, float)

    def test_output_shape(self, planted_partition_graph):
        """Factor matrix should be (n, k)."""
        labels = np.array([1] * 15 + [2] * 15 + [3] * 15)
        A = planted_partition_graph(labels)
        U1, _ = onmtf(A, k=3, runs=2, max_iter=50, rng=np.random.default_rng(0))
        assert U1.shape == (45, 3)

    def test_nmi_mode_with_ground_truth(self, planted_partition_graph):
        """When ground_truth is provided, should use NMI and return value in [0, 1]."""
        labels = np.array([1] * 20 + [2] * 20)
        A = planted_partition_graph(labels, p_in=0.5, p_out=0.02)
        U1, nmi = onmtf(
            A, k=2, ground_truth=labels,
            runs=5, max_iter=200, rng=np.random.default_rng(42)
        )
        assert 0.0 <= nmi <= 1.0

    def test_modularity_mode_without_ground_truth(self, planted_partition_graph):
        """When ground_truth is None, should use Modularity Density."""
        labels = np.array([1] * 20 + [2] * 20)
        A = planted_partition_graph(labels, p_in=0.5, p_out=0.02)
        U1, mod_den = onmtf(A, k=2, runs=3, max_iter=100, rng=np.random.default_rng(42))
        assert isinstance(mod_den, float)

    def test_reproducible_with_seed(self, planted_partition_graph):
        """Same rng seed should produce same result."""
        labels = np.array([1] * 20 + [2] * 20)
        A = planted_partition_graph(labels)
        U1_a, m_a = onmtf(A, k=2, runs=3, max_iter=50, rng=np.random.default_rng(99))
        U1_b, m_b = onmtf(A, k=2, runs=3, max_iter=50, rng=np.random.default_rng(99))
        np.testing.assert_array_equal(U1_a, U1_b)
        assert m_a == m_b

    def test_nonnegative_factors(self, planted_partition_graph):
        """Factor matrix should be nonnegative (initialized from rand)."""
        labels = np.array([1] * 20 + [2] * 20)
        A = planted_partition_graph(labels)
        U1, _ = onmtf(A, k=2, runs=2, max_iter=50, rng=np.random.default_rng(42))
        # May have some near-zero negatives from numerical issues, but should be ~nonneg
        assert np.all(U1 >= -1e-10)

    def test_validation_non_square(self):
        """Should raise ValueError for non-square matrix."""
        with pytest.raises(ValueError, match="square"):
            onmtf(np.ones((3, 4)), k=2)

    def test_validation_bad_k(self):
        """Should raise ValueError for k < 1."""
        with pytest.raises(ValueError, match="positive"):
            onmtf(np.eye(5), k=0)

    def test_validation_bad_eta(self):
        """Should raise ValueError for eta out of range."""
        with pytest.raises(ValueError, match="eta"):
            onmtf(np.eye(5), k=2, eta=0.0)
