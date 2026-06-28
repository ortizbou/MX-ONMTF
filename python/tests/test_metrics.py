"""Tests for mxonmtf.metrics."""

import numpy as np
import pytest
from mxonmtf.metrics import modularity_density, normalized_mutual_information


class TestNormalizedMutualInformation:
    def test_identical_partitions(self):
        """NMI of identical partitions should be 1."""
        A = np.array([1, 1, 2, 2, 3, 3])
        assert normalized_mutual_information(A, A) == pytest.approx(1.0)

    def test_different_partitions(self):
        """NMI of different partitions should be between 0 and 1."""
        A = np.array([1, 1, 2, 2, 3, 3])
        B = np.array([1, 2, 1, 2, 1, 2])
        nmi = normalized_mutual_information(A, B)
        assert 0.0 <= nmi <= 1.0

    def test_permuted_labels(self):
        """NMI should be 1 for relabeled but equivalent partitions."""
        A = np.array([1, 1, 2, 2, 3, 3])
        B = np.array([3, 3, 1, 1, 2, 2])  # same partition, different labels
        assert normalized_mutual_information(A, B) == pytest.approx(1.0)

    def test_single_community(self):
        """NMI with all nodes in one community should be 1 when both match."""
        A = np.array([1, 1, 1, 1])
        assert normalized_mutual_information(A, A) == pytest.approx(1.0)

    def test_length_mismatch(self):
        """Should raise ValueError for different-length inputs."""
        with pytest.raises(ValueError):
            normalized_mutual_information(np.array([1, 2]), np.array([1, 2, 3]))

    def test_symmetry(self):
        """NMI(A, B) should equal NMI(B, A)."""
        A = np.array([1, 1, 2, 2, 3, 3])
        B = np.array([1, 2, 2, 3, 3, 1])
        assert normalized_mutual_information(A, B) == pytest.approx(
            normalized_mutual_information(B, A)
        )


class TestModularityDensity:
    def test_perfect_communities(self):
        """Clear community structure should have positive modularity density."""
        # Two disjoint cliques
        A = np.zeros((6, 6))
        A[:3, :3] = 1.0
        A[3:, 3:] = 1.0
        np.fill_diagonal(A, 0)
        labels = np.array([1, 1, 1, 2, 2, 2])
        md, md_norm = modularity_density(2, A, labels)
        assert md > 0

    def test_single_community(self):
        """Single community: no external edges, LVn=0."""
        A = np.ones((4, 4))
        np.fill_diagonal(A, 0)
        labels = np.array([1, 1, 1, 1])
        md, md_norm = modularity_density(1, A, labels)
        # All edges are internal
        assert md > 0

    def test_empty_community(self):
        """Empty communities should be skipped without error."""
        A = np.ones((4, 4))
        np.fill_diagonal(A, 0)
        labels = np.array([1, 1, 1, 1])  # community 2 is empty
        md, md_norm = modularity_density(2, A, labels)
        assert md > 0

    def test_returns_two_values(self):
        """Should return both density and normalized density."""
        A = np.eye(4)
        labels = np.array([1, 1, 2, 2])
        result = modularity_density(2, A, labels)
        assert len(result) == 2
