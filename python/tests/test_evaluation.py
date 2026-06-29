"""Tests for mxonmtf.evaluation."""

import numpy as np
import pytest
from mxonmtf.evaluation import evaluate


class TestEvaluate:
    def _make_data(self):
        """Two-layer network with clear communities."""
        # Two disjoint cliques per layer
        A = np.zeros((6, 6))
        A[:3, :3] = 1.0
        A[3:, 3:] = 1.0
        np.fill_diagonal(A, 0)
        labels = np.array([1, 1, 1, 2, 2, 2])
        return [A, A], [labels, labels], np.array([2, 2])

    def test_returns_dict(self):
        """Should return a dict with expected keys."""
        Al, labels, k = self._make_data()
        result = evaluate(Al, labels, k)
        assert "per_layer" in result
        assert "mean_modularity_density" in result
        assert "cluster_sizes" in result

    def test_per_layer_count(self):
        """Should have one entry per layer."""
        Al, labels, k = self._make_data()
        result = evaluate(Al, labels, k)
        assert len(result["per_layer"]) == 2

    def test_modularity_positive_for_clear_structure(self):
        """Modularity density should be positive for clean communities."""
        Al, labels, k = self._make_data()
        result = evaluate(Al, labels, k)
        assert result["mean_modularity_density"] > 0

    def test_cluster_sizes(self):
        """Should correctly report cluster sizes."""
        Al, labels, k = self._make_data()
        result = evaluate(Al, labels, k)
        assert result["cluster_sizes"][0] == {1: 3, 2: 3}

    def test_length_mismatch(self):
        """Should raise ValueError if Al and labels have different lengths."""
        A = np.eye(5)
        with pytest.raises(ValueError, match="same length"):
            evaluate([A, A], [np.ones(5)], np.array([1, 1]))
