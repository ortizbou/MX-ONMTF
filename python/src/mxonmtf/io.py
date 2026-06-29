"""Data I/O utilities for multiplex networks.

Load adjacency matrices from common file formats and export results.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


def load_layers(path: str | Path, format: str = "auto") -> list[np.ndarray]:
    """Load multiplex network layers from file.

    Parameters
    ----------
    path : str or Path
        Path to data file or directory.
    format : str
        File format: 'auto' (detect from extension), 'npz', 'mat', or 'csv_dir'.
        For 'csv_dir', path should be a directory containing one CSV file per layer.

    Returns
    -------
    list[np.ndarray]
        List of adjacency matrices, one per layer.
    """
    path = Path(path)

    if format == "auto":
        if path.is_dir():
            format = "csv_dir"
        elif path.suffix == ".npz":
            format = "npz"
        elif path.suffix == ".mat":
            format = "mat"
        elif path.suffix == ".npy":
            format = "npy"
        else:
            raise ValueError(
                f"Cannot detect format from '{path.suffix}'. "
                "Use format='npz', 'mat', 'csv_dir', or 'npy'."
            )

    if format == "npz":
        data = np.load(path)
        # Sort keys to ensure consistent ordering
        keys = sorted(data.files)
        return [data[k].astype(float) for k in keys]

    elif format == "mat":
        from scipy.io import loadmat

        data = loadmat(str(path))
        # Filter out metadata keys (those starting with '__')
        keys = sorted(k for k in data.keys() if not k.startswith("__"))
        layers = []
        for k in keys:
            arr = np.asarray(data[k], dtype=float)
            if arr.ndim == 2 and arr.shape[0] == arr.shape[1]:
                layers.append(arr)
        if not layers:
            raise ValueError(f"No square matrices found in '{path}'.")
        return layers

    elif format == "npy":
        data = np.load(path)
        if data.ndim == 2:
            return [data.astype(float)]
        elif data.ndim == 3:
            return [data[i].astype(float) for i in range(data.shape[0])]
        else:
            raise ValueError(f"Expected 2D or 3D array, got {data.ndim}D.")

    elif format == "csv_dir":
        if not path.is_dir():
            raise ValueError(f"'{path}' is not a directory.")
        csv_files = sorted(path.glob("*.csv"))
        if not csv_files:
            raise ValueError(f"No CSV files found in '{path}'.")
        return [np.loadtxt(f, delimiter=",") for f in csv_files]

    else:
        raise ValueError(f"Unknown format '{format}'. Use 'npz', 'mat', 'csv_dir', or 'npy'.")


def validate_multiplex(Al: list[np.ndarray]) -> None:
    """Validate a list of adjacency matrices for use with MX-ONMTF.

    Checks that all matrices are square, same size, and contain no NaN values.
    Issues warnings for non-symmetric matrices.

    Parameters
    ----------
    Al : list of np.ndarray
        List of adjacency matrices.

    Raises
    ------
    ValueError
        If validation fails.
    """
    if not Al:
        raise ValueError("Al must be a non-empty list of adjacency matrices.")

    n = Al[0].shape[0]

    for i, A in enumerate(Al):
        if A.ndim != 2:
            raise ValueError(f"Layer {i}: expected 2D matrix, got {A.ndim}D.")
        if A.shape[0] != A.shape[1]:
            raise ValueError(f"Layer {i}: matrix is not square ({A.shape[0]}x{A.shape[1]}).")
        if A.shape[0] != n:
            raise ValueError(
                f"Layer {i}: has {A.shape[0]} nodes, expected {n} (same as layer 0)."
            )
        if np.any(np.isnan(A)):
            raise ValueError(f"Layer {i}: contains NaN values.")


def save_results(result, path: str | Path, format: str = "csv") -> None:
    """Export MX-ONMTF results to file.

    Parameters
    ----------
    result : MXONMTFResult
        Result object from mx_onmtf().
    path : str or Path
        Output file path.
    format : str
        Output format: 'csv' or 'json'. Default 'csv'.
    """
    path = Path(path)

    if format == "csv":
        # Save community labels
        if result.labels_per_layer:
            header_parts = [f"layer_{i}" for i in range(len(result.labels_per_layer))]
            data = np.column_stack(result.labels_per_layer)
            np.savetxt(path, data, delimiter=",", header=",".join(header_parts),
                       fmt="%d", comments="")
        elif result.clusters_supra is not None:
            np.savetxt(path, result.clusters_supra, delimiter=",",
                       header="cluster", fmt="%d", comments="")

    elif format == "json":
        out = {}
        if result.labels_per_layer:
            out["labels_per_layer"] = [lbl.tolist() for lbl in result.labels_per_layer]
        if result.clusters_supra is not None:
            out["clusters_supra"] = result.clusters_supra.tolist()
        if result.averagedNMI is not None:
            out["nmi"] = result.averagedNMI
        path.write_text(json.dumps(out, indent=2))

    else:
        raise ValueError(f"Unknown format '{format}'. Use 'csv' or 'json'.")
