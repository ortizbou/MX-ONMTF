"""Tests for mxonmtf.utils."""

import numpy as np
import pytest
from mxonmtf.utils import normalized_adjacency, row_normalize


class TestNormalizedAdjacency:
    def test_identity_like(self):
        """Isolated nodes (identity-like) should return zeros."""
        A = np.eye(3)
        result = normalized_adjacency(A)
        # Diagonal should be diag_w / degrees = 1/1 = 1
        np.testing.assert_array_almost_equal(np.diag(result), [1.0, 1.0, 1.0])

    def test_simple_graph(self):
        """Triangle graph should produce a valid normalized matrix."""
        A = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=float)
        result = normalized_adjacency(A)
        assert result.shape == (3, 3)
        # Symmetric input should produce symmetric output
        np.testing.assert_array_almost_equal(result, result.T)

    def test_disconnected_node(self):
        """Node with zero degree should not produce NaN."""
        A = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=float)
        result = normalized_adjacency(A)
        assert not np.any(np.isnan(result))
        # Node 3 has degree 0, its diagonal should be 0
        assert result[2, 2] == 0.0

    def test_output_shape(self):
        """Output shape should match input."""
        A = np.zeros((5, 5))
        result = normalized_adjacency(A)
        assert result.shape == (5, 5)


class TestRowNormalize:
    def test_basic(self):
        """Each row should have unit norm after normalization."""
        M = np.array([[3.0, 4.0], [1.0, 0.0], [0.0, 2.0]])
        result = row_normalize(M)
        norms = np.linalg.norm(result, axis=1)
        np.testing.assert_array_almost_equal(norms, [1.0, 1.0, 1.0])

    def test_zero_row(self):
        """Zero rows should remain zero (not produce NaN)."""
        M = np.array([[1.0, 0.0], [0.0, 0.0]])
        result = row_normalize(M)
        assert not np.any(np.isnan(result))
        np.testing.assert_array_equal(result[1], [0.0, 0.0])

    def test_already_normalized(self):
        """Already-normalized rows should be unchanged."""
        M = np.array([[1.0, 0.0], [0.0, 1.0]])
        result = row_normalize(M)
        np.testing.assert_array_almost_equal(result, M)
