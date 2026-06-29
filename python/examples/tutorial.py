"""MX-ONMTF Tutorial: Python implementation.

This script demonstrates the full MX-ONMTF pipeline on a synthetic
multiplex network with known community structure.

Usage:
    cd python/
    pip install -e .
    python examples/tutorial.py
"""

import tempfile
from pathlib import Path

import numpy as np
from mxonmtf import (
    mx_onmtf, findingk, onmtf,
    load_layers, validate_multiplex, save_results,
)
from mxonmtf.evaluation import evaluate


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
    # Step 1: Save and reload data (demonstrating I/O)
    # -------------------------------------------------------------------------
    print("\n" + "-" * 60)
    print("Step 1: Data I/O — save and load multiplex network")
    print("-" * 60)

    with tempfile.TemporaryDirectory() as tmpdir:
        data_path = Path(tmpdir) / "layers.npz"
        np.savez(data_path, layer_0=A1, layer_1=A2)
        print(f"  Saved to {data_path}")

        Al_loaded = load_layers(data_path)
        validate_multiplex(Al_loaded)
        print(f"  Loaded {len(Al_loaded)} layers, {Al_loaded[0].shape[0]} nodes each")
        print(f"  Validation passed")

    # -------------------------------------------------------------------------
    # Step 2: Run MX-ONMTF with known parameters
    # -------------------------------------------------------------------------
    print("\n" + "-" * 60)
    print("Step 2: Running MX-ONMTF")
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
    print(f"  Supralayer clusters shape: {result.clusters_supra.shape}")
    print(f"  Per-layer labels: {len(result.labels_per_layer)} layers")

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
    # Step 3: Inspect results
    # -------------------------------------------------------------------------
    print("\n" + "-" * 60)
    print("Step 3: Result summary and per-layer access")
    print("-" * 60)

    print()
    print(result.summary())

    print(f"\n  Layer 0 labels (first 10): {result.get_layer_labels(0)[:10]}")
    print(f"  Layer 1 labels (first 10): {result.get_layer_labels(1)[:10]}")

    # -------------------------------------------------------------------------
    # Step 4: Evaluate without ground truth
    # -------------------------------------------------------------------------
    print("\n" + "-" * 60)
    print("Step 4: Quality evaluation (no ground truth needed)")
    print("-" * 60)

    metrics = evaluate(Al, result.labels_per_layer, k)
    print(f"\n  Mean Modularity Density: {metrics['mean_modularity_density']:.4f}")
    for l, pl in enumerate(metrics["per_layer"]):
        print(f"  Layer {l}: ModDen={pl['modularity_density']:.4f}, "
              f"ModDenNorm={pl['modularity_density_norm']:.4f}")
    print(f"\n  Cluster sizes: {metrics['cluster_sizes']}")

    # -------------------------------------------------------------------------
    # Step 5: Export results
    # -------------------------------------------------------------------------
    print("\n" + "-" * 60)
    print("Step 5: Save results")
    print("-" * 60)

    with tempfile.TemporaryDirectory() as tmpdir:
        csv_path = Path(tmpdir) / "communities.csv"
        save_results(result, csv_path, format="csv")
        print(f"  Saved CSV to {csv_path}")

        json_path = Path(tmpdir) / "communities.json"
        save_results(result, json_path, format="json")
        print(f"  Saved JSON to {json_path}")

    # -------------------------------------------------------------------------
    # Step 6: Single-layer baseline comparison
    # -------------------------------------------------------------------------
    print("\n" + "-" * 60)
    print("Step 6: Single-layer ONMTF baseline comparison")
    print("-" * 60)

    for l, (A, gt) in enumerate(zip(Al, [labels1, labels2])):
        U1, nmi = onmtf(
            A, k=k[l], ground_truth=gt,
            runs=10, max_iter=500,
            rng=np.random.default_rng(42),
        )
        print(f"  Layer {l + 1} single-layer NMI: {nmi:.4f}")

    # -------------------------------------------------------------------------
    # Step 7: Automatic parameter selection
    # -------------------------------------------------------------------------
    print("\n" + "-" * 60)
    print("Step 7: Automatic parameter selection (findingk)")
    print("-" * 60)

    kc_auto, k_auto, kpl_auto = findingk(Al)
    print(f"  Detected kc:  {kc_auto} (true: 2)")
    print(f"  Detected k:   {k_auto} (true: [3, 3])")
    print(f"  Detected kpl: {kpl_auto} (true: [1, 1])")

    # -------------------------------------------------------------------------
    # CLI reminder
    # -------------------------------------------------------------------------
    print("\n" + "-" * 60)
    print("CLI usage (alternative to scripting)")
    print("-" * 60)
    print("  python -m mxonmtf findingk --input layers.npz")
    print("  python -m mxonmtf run --input layers.npz --k 3,3 --kc 2 --kpl 1,1 --output results.csv")

    print("\n" + "=" * 60)
    print("Tutorial complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
