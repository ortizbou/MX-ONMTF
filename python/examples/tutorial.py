"""MX-ONMTF Tutorial: Python implementation.

This script demonstrates the full MX-ONMTF pipeline on a synthetic
multiplex network with known community structure.

Usage:
    cd python/
    pip install -e .
    python examples/tutorial.py
"""

import numpy as np
from mxonmtf import mx_onmtf, findingk, onmtf, normalized_mutual_information


def generate_planted_partition(labels, p_in=0.4, p_out=0.05, rng=None):
    """Generate a symmetric adjacency matrix from a planted partition model."""
    if rng is None:
        rng = np.random.default_rng()
    n = len(labels)
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            p = p_in if labels[i] == labels[j] else p_out
            if rng.random() < p:
                A[i, j] = 1.0
                A[j, i] = 1.0
    return A


def main():
    rng = np.random.default_rng(42)

    # -------------------------------------------------------------------------
    # Step 0: Generate synthetic multiplex network
    # -------------------------------------------------------------------------
    print("=" * 60)
    print("MX-ONMTF Tutorial")
    print("=" * 60)

    n = 60  # nodes
    # Layer 1: communities {1, 2, 3} — communities 1,2 are common, 3 is private
    labels1 = np.array([1] * 20 + [2] * 20 + [3] * 20)
    # Layer 2: communities {1, 2, 4} — communities 1,2 are common, 4 is private
    labels2 = np.array([1] * 20 + [2] * 20 + [4] * 20)

    A1 = generate_planted_partition(labels1, p_in=0.5, p_out=0.03, rng=rng)
    A2 = generate_planted_partition(labels2, p_in=0.5, p_out=0.03, rng=rng)
    Al = [A1, A2]

    print(f"\nGenerated 2-layer multiplex network with {n} nodes")
    print(f"  Layer 1: communities [1, 2, 3] (3 is private)")
    print(f"  Layer 2: communities [1, 2, 4] (4 is private)")
    print(f"  True common communities (kc): 2")

    # -------------------------------------------------------------------------
    # Step 1: Run MX-ONMTF with known parameters
    # -------------------------------------------------------------------------
    print("\n" + "-" * 60)
    print("Step 1: Running MX-ONMTF with known community structure")
    print("-" * 60)

    kc = 2  # common communities
    k = np.array([3, 3])  # total per layer
    kpl = np.array([1, 1])  # private per layer

    # Without ground truth (trace minimization)
    result = mx_onmtf(
        Al, k=k, kc=kc, kpl=kpl,
        runs=10, max_iter=500,
        rng=np.random.default_rng(42),
    )
    print(f"\n  Trace minimization mode (no ground truth):")
    print(f"  Clusters shape: {result.Clusters[0].shape}")

    # With ground truth (NMI selection)
    result_nmi = mx_onmtf(
        Al, k=k, kc=kc, kpl=kpl,
        ground_truth=[labels1, labels2],
        runs=10, max_iter=500,
        rng=np.random.default_rng(42),
    )
    print(f"\n  NMI selection mode (with ground truth):")
    print(f"  Best NMI: {result_nmi.averagedNMI:.4f}")

    # -------------------------------------------------------------------------
    # Step 2: Single-layer baseline comparison
    # -------------------------------------------------------------------------
    print("\n" + "-" * 60)
    print("Step 2: Single-layer ONMTF baseline comparison")
    print("-" * 60)

    for l, (A, gt) in enumerate(zip(Al, [labels1, labels2])):
        U1, nmi = onmtf(
            A, k=k[l], ground_truth=gt,
            runs=10, max_iter=500,
            rng=np.random.default_rng(42),
        )
        print(f"  Layer {l + 1} single-layer NMI: {nmi:.4f}")

    # -------------------------------------------------------------------------
    # Step 3: Automatic parameter selection
    # -------------------------------------------------------------------------
    print("\n" + "-" * 60)
    print("Step 3: Automatic parameter selection (findingk)")
    print("-" * 60)

    kc_auto, k_auto, kpl_auto = findingk(Al)
    print(f"  Detected kc:  {kc_auto} (true: 2)")
    print(f"  Detected k:   {k_auto} (true: [3, 3])")
    print(f"  Detected kpl: {kpl_auto} (true: [1, 1])")

    print("\n" + "=" * 60)
    print("Tutorial complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
