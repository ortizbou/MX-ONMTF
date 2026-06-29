"""CLI entry point for MX-ONMTF.

Usage:
    python -m mxonmtf run --input layers.npz --k 3,3 --kc 2 --kpl 1,1 --output results.csv
    python -m mxonmtf findingk --input layers.npz
"""

from __future__ import annotations

import argparse
import sys

import numpy as np


def _parse_int_list(s: str) -> np.ndarray:
    """Parse comma-separated integers."""
    return np.array([int(x.strip()) for x in s.split(",")])


def cmd_run(args):
    from mxonmtf import mx_onmtf
    from mxonmtf.evaluation import evaluate
    from mxonmtf.io import load_layers, save_results, validate_multiplex

    Al = load_layers(args.input)
    validate_multiplex(Al)

    k = _parse_int_list(args.k)
    kpl = _parse_int_list(args.kpl)

    result = mx_onmtf(
        Al, k=k, kc=args.kc, kpl=kpl,
        eta=args.eta, runs=args.runs, max_iter=args.max_iter,
        rng=np.random.default_rng(args.seed) if args.seed is not None else None,
    )

    print(result.summary())

    if result.labels_per_layer:
        metrics = evaluate(Al, result.labels_per_layer, k)
        print(f"\nMean Modularity Density: {metrics['mean_modularity_density']:.4f}")
        for l, pl in enumerate(metrics["per_layer"]):
            print(f"  Layer {l}: {pl['modularity_density']:.4f}")

    if args.output:
        save_results(result, args.output, format=args.format)
        print(f"\nResults saved to {args.output}")


def cmd_findingk(args):
    from mxonmtf import findingk
    from mxonmtf.io import load_layers, validate_multiplex

    Al = load_layers(args.input)
    validate_multiplex(Al)

    kc, k, kpl = findingk(Al)

    print(f"Common communities (kc):  {kc}")
    print(f"Total per layer (k):      {k}")
    print(f"Private per layer (kpl):  {kpl}")
    print(f"\nUse with: python -m mxonmtf run --input {args.input} "
          f"--k {','.join(map(str, k))} --kc {kc} --kpl {','.join(map(str, kpl))}")


def main():
    parser = argparse.ArgumentParser(
        prog="mxonmtf",
        description="MX-ONMTF: Community detection in multiplex networks",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # run command
    run_parser = subparsers.add_parser("run", help="Run MX-ONMTF on data")
    run_parser.add_argument("--input", required=True, help="Path to data file (.npz, .mat, .npy) or CSV directory")
    run_parser.add_argument("--k", required=True, help="Total communities per layer (comma-separated)")
    run_parser.add_argument("--kc", required=True, type=int, help="Number of common communities")
    run_parser.add_argument("--kpl", required=True, help="Private communities per layer (comma-separated)")
    run_parser.add_argument("--output", help="Output file path")
    run_parser.add_argument("--format", default="csv", choices=["csv", "json"], help="Output format (default: csv)")
    run_parser.add_argument("--eta", type=float, default=0.5, help="Learning rate (default: 0.5)")
    run_parser.add_argument("--runs", type=int, default=20, help="Number of runs (default: 20)")
    run_parser.add_argument("--max-iter", type=int, default=1000, help="Max iterations (default: 1000)")
    run_parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility")
    run_parser.set_defaults(func=cmd_run)

    # findingk command
    fk_parser = subparsers.add_parser("findingk", help="Auto-detect community counts")
    fk_parser.add_argument("--input", required=True, help="Path to data file")
    fk_parser.set_defaults(func=cmd_findingk)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
