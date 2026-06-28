"""Tests for mxonmtf.factorization."""

import numpy as np
import pytest
from mxonmtf.factorization import MXONMTFResult, mx_onmtf


class TestMXONMTF:
    def _make_multiplex(self, planted_partition_graph):
        """Helper: create a 2-layer multiplex with known structure."""
        labels1 = np.array([1] * 20 + [2] * 20)
        labels2 = np.array([1] * 20 + [2] * 20)
        A1 = planted_partition_graph(labels1, p_in=0.5, p_out=0.02)
        A2 = planted_partition_graph(labels2, p_in=0.5, p_out=0.02)
        return [A1, A2], [labels1, labels2]

    def test_returns_result_object(self, planted_partition_graph):
        """Should return an MXONMTFResult."""
        Al, _ = self._make_multiplex(planted_partition_graph)
        result = mx_onmtf(
            Al, k=np.array([2, 2]), kc=2, kpl=np.array([0, 0]),
            runs=2, max_iter=50, rng=np.random.default_rng(42)
        )
        assert isinstance(result, MXONMTFResult)

    def test_trace_mode_has_clusters(self, planted_partition_graph):
        """Without ground truth, should produce clusters via trace minimization."""
        Al, _ = self._make_multiplex(planted_partition_graph)
        result = mx_onmtf(
            Al, k=np.array([2, 2]), kc=2, kpl=np.array([0, 0]),
            runs=2, max_iter=50, rng=np.random.default_rng(42)
        )
        assert len(result.Clusters) > 0
        assert len(result.H_best) > 0

    def test_nmi_mode_with_ground_truth(self, planted_partition_graph):
        """With ground truth, should use NMI and return NMI stats."""
        Al, GTl = self._make_multiplex(planted_partition_graph)
        result = mx_onmtf(
            Al, k=np.array([2, 2]), kc=2, kpl=np.array([0, 0]),
            ground_truth=GTl,
            runs=3, max_iter=100, rng=np.random.default_rng(42)
        )
        assert result.averagedNMI is not None
        assert 0.0 <= result.averagedNMI <= 1.0

    def test_with_private_communities(self, planted_partition_graph):
        """Should handle layers with private communities."""
        labels1 = np.array([1] * 15 + [2] * 15 + [3] * 15)
        labels2 = np.array([1] * 15 + [2] * 15 + [4] * 15)
        A1 = planted_partition_graph(labels1, p_in=0.4, p_out=0.03)
        A2 = planted_partition_graph(labels2, p_in=0.4, p_out=0.03)
        result = mx_onmtf(
            [A1, A2], k=np.array([3, 3]), kc=2, kpl=np.array([1, 1]),
            runs=2, max_iter=50, rng=np.random.default_rng(42)
        )
        assert len(result.Clusters) > 0

    def test_reproducible(self, planted_partition_graph):
        """Same rng should produce same result."""
        Al, _ = self._make_multiplex(planted_partition_graph)
        kwargs = dict(
            k=np.array([2, 2]), kc=2, kpl=np.array([0, 0]),
            runs=2, max_iter=50,
        )
        r1 = mx_onmtf(Al, **kwargs, rng=np.random.default_rng(99))
        r2 = mx_onmtf(Al, **kwargs, rng=np.random.default_rng(99))
        np.testing.assert_array_equal(r1.H_best[0], r2.H_best[0])

    def test_validation_bad_eta(self):
        """Should raise ValueError for eta out of range."""
        with pytest.raises(ValueError, match="eta"):
            mx_onmtf([np.eye(5)], k=np.array([2]), kc=1, kpl=np.array([1]), eta=0.0)

    def test_validation_bad_kc(self):
        """Should raise ValueError for kc < 1."""
        with pytest.raises(ValueError, match="kc"):
            mx_onmtf([np.eye(5)], k=np.array([2]), kc=0, kpl=np.array([1]))

    def test_validation_k_kpl_mismatch(self):
        """Should raise ValueError when k and kpl have different lengths."""
        with pytest.raises(ValueError, match="same length"):
            mx_onmtf([np.eye(5)], k=np.array([2, 3]), kc=1, kpl=np.array([1]))
