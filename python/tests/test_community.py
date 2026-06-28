"""Tests for mxonmtf.community."""

import numpy as np
import pytest
from mxonmtf.community import assign_communities, patchmult


class TestPatchmult:
    def test_returns_pos_and_q(self):
        """Should return sorted positions and density ratios."""
        A = np.ones((6, 6))
        np.fill_diagonal(A, 0)
        H = np.zeros((6, 2))
        H[:3, 0] = 1.0
        H[3:, 1] = 1.0
        pos, q = patchmult(A, H, kc=2)
        assert len(pos) == 2
        assert len(q) == 2

    def test_pos_is_permutation(self):
        """Positions should be a permutation of range(kc)."""
        A = np.ones((6, 6))
        np.fill_diagonal(A, 0)
        H = np.zeros((6, 3))
        H[:2, 0] = 1.0
        H[2:4, 1] = 1.0
        H[4:, 2] = 1.0
        pos, _ = patchmult(A, H, kc=3)
        assert sorted(pos) == [0, 1, 2]


class TestAssignCommunities:
    def test_flex_mode(self):
        """Flex mode should produce labels for all nodes."""
        n = 20
        kc = 2
        L = 2
        k = np.array([3, 3])
        kpl = np.array([1, 1])

        rng = np.random.default_rng(42)
        A1 = np.ones((n, n)); np.fill_diagonal(A1, 0)
        A2 = np.ones((n, n)); np.fill_diagonal(A2, 0)
        H = rng.random((n, kc))
        Hl = [rng.random((n, kpl[l])) for l in range(L)]

        Il, supra = assign_communities(
            [A1, A2], H, Hl, kc, k, kpl, L, mode='flex'
        )
        assert len(Il) == L
        assert all(len(Il[l]) == n for l in range(L))
        assert len(supra) == n * L

    def test_same_mode(self):
        """Same mode should produce labels for all nodes."""
        n = 20
        kc = 2
        L = 2
        k = np.array([3, 3])
        kpl = np.array([1, 1])

        rng = np.random.default_rng(42)
        A1 = np.ones((n, n)); np.fill_diagonal(A1, 0)
        A2 = np.ones((n, n)); np.fill_diagonal(A2, 0)
        H = rng.random((n, kc))
        Hl = [rng.random((n, kpl[l])) for l in range(L)]

        Il, supra = assign_communities(
            [A1, A2], H, Hl, kc, k, kpl, L, mode='same', perc=0.8
        )
        assert len(Il) == L
        assert len(supra) == n * L

    def test_all_nodes_assigned(self):
        """No node should remain with label 0 after assignment."""
        n = 20
        kc = 2
        L = 2
        k = np.array([3, 3])
        kpl = np.array([1, 1])

        rng = np.random.default_rng(42)
        A1 = np.ones((n, n)); np.fill_diagonal(A1, 0)
        A2 = np.ones((n, n)); np.fill_diagonal(A2, 0)
        H = rng.random((n, kc))
        Hl = [rng.random((n, kpl[l])) for l in range(L)]

        Il, _ = assign_communities(
            [A1, A2], H, Hl, kc, k, kpl, L, mode='same', perc=0.5
        )
        for l in range(L):
            assert np.all(Il[l] > 0), f"Layer {l} has unassigned nodes"

    def test_invalid_mode(self):
        """Should raise ValueError for invalid mode."""
        with pytest.raises(ValueError, match="mode"):
            assign_communities(
                [np.eye(5)], np.ones((5, 2)), [np.ones((5, 1))],
                kc=2, k=np.array([3]), kpl=np.array([1]), L=1, mode='invalid'
            )
