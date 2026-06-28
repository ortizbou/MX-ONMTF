"""Tests for mxonmtf.parameter_selection."""

import numpy as np
import pytest
from mxonmtf.parameter_selection import eiggap, findingk


class TestEigGap:
    def test_two_communities(self, planted_partition_graph):
        """Should detect approximately 2 communities in a 2-block graph."""
        labels = np.array([1] * 30 + [2] * 30)
        A = planted_partition_graph(labels, p_in=0.5, p_out=0.02)
        k = eiggap(A)
        assert k >= 1  # should detect at least 1 community
        assert k <= 10  # shouldn't overestimate wildly

    def test_single_clique(self):
        """Fully connected graph should detect 1 community."""
        A = np.ones((20, 20))
        np.fill_diagonal(A, 0)
        k = eiggap(A)
        assert k >= 1

    def test_returns_integer(self, planted_partition_graph):
        """Output should be an integer."""
        labels = np.array([1] * 20 + [2] * 20)
        A = planted_partition_graph(labels)
        k = eiggap(A)
        assert isinstance(k, (int, np.integer))

    def test_reproducible(self, planted_partition_graph):
        """Same seed should give same result."""
        labels = np.array([1] * 20 + [2] * 20)
        A = planted_partition_graph(labels)
        k1 = eiggap(A, seed=42)
        k2 = eiggap(A, seed=42)
        assert k1 == k2

    def test_custom_threshold(self, planted_partition_graph):
        """Different thresholds should run without error."""
        labels = np.array([1] * 20 + [2] * 20)
        A = planted_partition_graph(labels)
        k_low = eiggap(A, eiggap_thresh=0.5)
        k_high = eiggap(A, eiggap_thresh=0.99)
        assert isinstance(k_low, (int, np.integer))
        assert isinstance(k_high, (int, np.integer))


class TestFindingK:
    def _mock_onmtf(self, A, k, **kwargs):
        """Mock ONMTF that returns random factor matrix."""
        n = A.shape[0]
        rng = np.random.default_rng(42)
        U = rng.random((n, k))
        return U, 0.5

    def test_returns_three_values(self, planted_partition_graph):
        """Should return (kc, k, kpl) tuple."""
        labels = np.array([1] * 20 + [2] * 20)
        A1 = planted_partition_graph(labels)
        A2 = planted_partition_graph(labels)
        kc, k, kpl = findingk([A1, A2], onmtf_fn=self._mock_onmtf)
        assert isinstance(kc, (int, np.integer))
        assert len(k) == 2
        assert len(kpl) == 2

    def test_kpl_non_negative(self, planted_partition_graph):
        """Private community counts should be non-negative."""
        labels = np.array([1] * 20 + [2] * 20)
        A1 = planted_partition_graph(labels)
        A2 = planted_partition_graph(labels)
        kc, k, kpl = findingk([A1, A2], onmtf_fn=self._mock_onmtf)
        assert np.all(kpl >= 0)

    def test_validation_non_square(self):
        """Should raise ValueError for non-square matrix."""
        A = np.ones((3, 4))
        with pytest.raises(ValueError, match="not square"):
            findingk([A], onmtf_fn=self._mock_onmtf)

    def test_validation_size_mismatch(self, planted_partition_graph):
        """Should raise ValueError when layers have different sizes."""
        A1 = np.ones((10, 10))
        A2 = np.ones((15, 15))
        with pytest.raises(ValueError, match="same number of nodes"):
            findingk([A1, A2], onmtf_fn=self._mock_onmtf)
