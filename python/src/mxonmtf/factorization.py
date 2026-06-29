"""Core MX-ONMTF multiplex factorization algorithm.

Implements the multiplicative update rules for decomposing multiplex network
layers into common (H) and private (Hl) community factor matrices, with
diagonal scaling matrices Sl and Gl.

Each layer's adjacency is approximated as:
    Al ≈ H × Sl × H' + Hl × Gl × Hl'

Reference:
    M. Ortiz-Bouza and S. Aviyente, "Community Detection in Multiplex Networks
    Based on Orthogonal Nonnegative Matrix Tri-Factorization," IEEE Access,
    vol. 12, pp. 6423-6436, 2024.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from mxonmtf.metrics import normalized_mutual_information

_EPS = 1e-16  # small constant to avoid division by zero


@dataclass
class MXONMTFResult:
    """Container for MX-ONMTF results.

    Attributes
    ----------
    H : np.ndarray or None
        Common community factor matrix (n x kc).
    Hl : list[np.ndarray] or None
        Private community factor matrices, one per layer (each n x kpl[l]).
    labels_per_layer : list[np.ndarray]
        Community labels for each layer (each length n, 1-indexed).
    clusters_supra : np.ndarray or None
        Concatenated community labels across all layers (length n*L).
    averagedNMI : float or None
        Best NMI across runs (when ground_truth provided).
    """

    H: np.ndarray | None = None
    Hl: list[np.ndarray] | None = None
    labels_per_layer: list[np.ndarray] = field(default_factory=list)
    clusters_supra: np.ndarray | None = None
    averagedNMI: float | None = None

    def get_layer_labels(self, layer: int) -> np.ndarray:
        """Return community labels for a specific layer.

        Parameters
        ----------
        layer : int
            Layer index (0-indexed).

        Returns
        -------
        np.ndarray
            Community labels (length n, 1-indexed).
        """
        if not self.labels_per_layer:
            raise ValueError("No per-layer labels available.")
        return self.labels_per_layer[layer]

    def summary(self) -> str:
        """Return a human-readable summary of the results."""
        lines = ["MX-ONMTF Results"]
        lines.append("-" * 40)

        if self.labels_per_layer:
            L = len(self.labels_per_layer)
            lines.append(f"Layers: {L}")
            for l, lbl in enumerate(self.labels_per_layer):
                unique, counts = np.unique(lbl[lbl > 0], return_counts=True)
                sizes = ", ".join(f"{int(u)}:{int(c)}" for u, c in zip(unique, counts))
                lines.append(f"  Layer {l}: {len(unique)} communities ({sizes})")

        if self.averagedNMI is not None:
            lines.append(f"NMI: {self.averagedNMI:.4f}")

        if self.H is not None:
            lines.append(f"H shape: {self.H.shape}")

        return "\n".join(lines)


def mx_onmtf(
    Al: list[np.ndarray],
    k: np.ndarray,
    kc: int,
    kpl: np.ndarray,
    ground_truth: list[np.ndarray] | None = None,
    eta: float = 0.5,
    runs: int = 20,
    max_iter: int = 1000,
    conv_thresh: float = 1e-3,
    assign_mode: str = "same",
    assign_perc: float = 0.8,
    rng: np.random.Generator | None = None,
) -> MXONMTFResult:
    """Multiplex Orthogonal Nonnegative Matrix Tri-Factorization.

    When ground_truth is provided, uses NMI for solution selection (synthetic).
    When omitted, uses trace minimization (real).

    Update rules (Eq. from paper):
        Hl ← Hl × ((Al Hl Gl + Hl Hl' H Sl H' Hl Gl') /
                    (H Sl H' Hl Gl + Hl Hl' Al Hl Gl))^eta

        H  ← H  × ((Σ_l Al H Sl + H H' Hl Gl' Hl' H Sl) /
                    (Σ_l Hl Gl' Hl' H Sl + H H' Al H Sl))^eta

        Sl ← Sl × ((H' Al H) / (H' H Sl H' H + H' Hl Gl Hl' H))^eta

        Gl ← Gl × ((Hl' Al Hl) / (Hl' Hl Gl Hl' Hl + Hl' H Sl H' Hl))^eta

    Parameters
    ----------
    Al : list of np.ndarray
        List of L adjacency matrices (each n x n).
    k : array-like of int
        Total communities per layer (length L).
    kc : int
        Number of common communities.
    kpl : array-like of int
        Private communities per layer (length L).
    ground_truth : list of np.ndarray, optional
        Ground truth labels per layer (each length n). When provided, NMI is
        used for solution selection. When None, trace minimization is used.
    eta : float
        Learning rate in (0, 1]. Default 0.5.
    runs : int
        Number of independent runs. Default 20.
    max_iter : int
        Maximum iterations per run. Default 1000.
    conv_thresh : float
        Convergence threshold. Default 1e-3.
    assign_mode : str
        Community assignment mode ('flex' or 'same'). Default 'same'.
    assign_perc : float
        Threshold for common community assignment in [0, 1]. Default 0.8.
    rng : np.random.Generator, optional
        Random number generator for reproducibility.

    Returns
    -------
    MXONMTFResult
        Result object with H_best, Hl_best, Clusters, and NMI stats if applicable.
    """
    # Lazy import to avoid circular dependency
    from mxonmtf.community import assign_communities

    k = np.asarray(k, dtype=int)
    kpl = np.asarray(kpl, dtype=int)

    # Input validation
    if not (0 < eta <= 1):
        raise ValueError("eta must be in (0, 1].")
    if kc < 1:
        raise ValueError("kc must be a positive integer.")
    if len(k) != len(kpl):
        raise ValueError("k and kpl must have the same length.")
    if assign_mode not in ("flex", "same"):
        raise ValueError("assign_mode must be 'flex' or 'same'.")

    if rng is None:
        rng = np.random.default_rng()

    L = len(Al)
    n = Al[0].shape[0]
    use_nmi = ground_truth is not None

    result = MXONMTFResult()

    # Single realization
    NMI_sup = np.full(runs, -np.inf)
    trace_all = np.full(runs, np.inf)

    for j in range(runs):
        # Initialize Hl, H, Sl, Gl
        Hl = [rng.random((n, kpl[l])) for l in range(L)]
        Sl = [np.diag(rng.random(kc)) for _ in range(L)]
        Gl = [np.diag(rng.random(kpl[l])) for l in range(L)]
        H = rng.random((n, kc))

        # Iterative update
        for _ in range(max_iter):
            # Update Hl for each layer
            Hl_new = []
            for l in range(L):
                num = (Al[l] @ Hl[l] @ Gl[l]
                       + Hl[l] @ Hl[l].T @ H @ Sl[l] @ H.T @ Hl[l] @ Gl[l].T)
                den = (H @ Sl[l] @ H.T @ Hl[l] @ Gl[l]
                       + Hl[l] @ Hl[l].T @ Al[l] @ Hl[l] @ Gl[l])
                den = np.maximum(den, _EPS)
                Hl_new_l = Hl[l] * (num / den) ** eta
                if np.all(np.isnan(Hl_new_l)):
                    Hl_new_l = rng.random((n, kpl[l]))
                Hl_new.append(Hl_new_l)

            # Update H (summed across layers)
            numH = np.zeros((n, kc))
            denH = np.zeros((n, kc))
            for l in range(L):
                numH += Al[l] @ H @ Sl[l] + H @ H.T @ Hl[l] @ Gl[l].T @ Hl[l].T @ H @ Sl[l]
                denH += Hl[l] @ Gl[l].T @ Hl[l].T @ H @ Sl[l] + H @ H.T @ Al[l] @ H @ Sl[l]
            denH = np.maximum(denH, _EPS)
            H_new = H * (numH / denH) ** eta
            if np.all(np.isnan(H_new)):
                H_new = rng.random((n, kc))

            # Update Sl and Gl for each layer
            Sl_new = []
            Gl_new = []
            for l in range(L):
                # Sl update
                num_S = H.T @ Al[l] @ H
                den_S = H.T @ H @ Sl[l] @ (H.T @ H) + H.T @ Hl[l] @ Gl[l] @ Hl[l].T @ H
                den_S = np.maximum(den_S, _EPS)
                Sl_new_l = Sl[l] * (num_S / den_S) ** eta
                if np.all(np.isnan(Sl_new_l)):
                    Sl_new_l = np.diag(rng.random(kc))
                Sl_new.append(Sl_new_l)

                # Gl update
                num_G = Hl[l].T @ Al[l] @ Hl[l]
                den_G = (Hl[l].T @ Hl[l] @ Gl[l] @ (Hl[l].T @ Hl[l])
                         + Hl[l].T @ H @ Sl[l] @ H.T @ Hl[l])
                den_G = np.maximum(den_G, _EPS)
                Gl_new_l = Gl[l] * (num_G / den_G) ** eta
                if np.all(np.isnan(Gl_new_l)):
                    Gl_new_l = np.diag(rng.random(kpl[l]))
                Gl_new.append(Gl_new_l)

            # Convergence check (last layer + H)
            conv_Hl = max(np.linalg.norm(Hl[l] - Hl_new[l]) for l in range(L))
            conv_Sl = max(np.linalg.norm(Sl[l] - Sl_new[l]) for l in range(L))
            conv_H = np.linalg.norm(H - H_new)

            Hl = Hl_new
            H = H_new
            Sl = Sl_new
            Gl = Gl_new

            if conv_Hl < conv_thresh and conv_Sl < conv_thresh and conv_H < conv_thresh:
                break

        # Community assignment
        Il, clusters_supra = assign_communities(
            Al, H, Hl, kc, k, kpl, L, mode=assign_mode, perc=assign_perc
        )

        # Solution selection
        if use_nmi:
            GT_supra = np.concatenate(ground_truth)
            nmi_sup_j = normalized_mutual_information(GT_supra, clusters_supra)
            NMI_sup[j] = nmi_sup_j

            if nmi_sup_j >= np.max(NMI_sup[:j + 1]):
                result.H = H
                result.Hl = Hl
                result.labels_per_layer = Il
                result.clusters_supra = clusters_supra

            if abs(np.max(NMI_sup[:j + 1]) - 1.0) < 1e-3:
                break
        else:
            # Trace minimization
            tr_val = 0.0
            for l in range(L):
                recon = H @ Sl[l] @ H.T + Hl[l] @ Gl[l] @ Hl[l].T
                tr_val += np.linalg.norm(Al[l] - recon, 'fro')
                tr_val += np.trace(Hl[l].T @ Hl[l] - np.eye(kpl[l]))
            tr_val += np.trace(H.T @ H - np.eye(kc))
            trace_all[j] = tr_val

            if tr_val <= np.min(trace_all[:j + 1]):
                result.H = H
                result.Hl = Hl
                result.labels_per_layer = Il
                result.clusters_supra = clusters_supra

    if use_nmi:
        best_nmi = np.max(NMI_sup[NMI_sup > -np.inf])
        result.averagedNMI = best_nmi

    return result
