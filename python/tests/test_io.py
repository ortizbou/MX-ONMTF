"""Tests for mxonmtf.io."""

import json
import numpy as np
import pytest
from mxonmtf.io import load_layers, save_results, validate_multiplex
from mxonmtf.factorization import MXONMTFResult


class TestValidateMultiplex:
    def test_valid(self):
        """Should pass for valid multiplex."""
        Al = [np.eye(5), np.ones((5, 5))]
        validate_multiplex(Al)  # should not raise

    def test_empty(self):
        with pytest.raises(ValueError, match="non-empty"):
            validate_multiplex([])

    def test_non_square(self):
        with pytest.raises(ValueError, match="not square"):
            validate_multiplex([np.ones((3, 4))])

    def test_size_mismatch(self):
        with pytest.raises(ValueError, match="nodes"):
            validate_multiplex([np.eye(3), np.eye(5)])

    def test_nan(self):
        A = np.eye(3)
        A[0, 0] = np.nan
        with pytest.raises(ValueError, match="NaN"):
            validate_multiplex([A])

    def test_non_2d(self):
        with pytest.raises(ValueError, match="2D"):
            validate_multiplex([np.ones(5)])


class TestLoadLayers:
    def test_load_npz(self, tmp_path):
        """Should load from .npz file."""
        A1 = np.eye(5)
        A2 = np.ones((5, 5))
        np.savez(tmp_path / "layers.npz", layer_0=A1, layer_1=A2)
        result = load_layers(tmp_path / "layers.npz")
        assert len(result) == 2
        np.testing.assert_array_equal(result[0], A1)
        np.testing.assert_array_equal(result[1], A2)

    def test_load_npy_3d(self, tmp_path):
        """Should load 3D array from .npy."""
        data = np.random.rand(3, 5, 5)
        np.save(tmp_path / "layers.npy", data)
        result = load_layers(tmp_path / "layers.npy")
        assert len(result) == 3
        np.testing.assert_array_almost_equal(result[0], data[0])

    def test_load_csv_dir(self, tmp_path):
        """Should load from directory of CSV files."""
        A1 = np.eye(4)
        A2 = np.ones((4, 4))
        np.savetxt(tmp_path / "layer_0.csv", A1, delimiter=",")
        np.savetxt(tmp_path / "layer_1.csv", A2, delimiter=",")
        result = load_layers(tmp_path)
        assert len(result) == 2
        np.testing.assert_array_almost_equal(result[0], A1)

    def test_auto_detect_npz(self, tmp_path):
        """Auto format should detect .npz."""
        np.savez(tmp_path / "data.npz", a=np.eye(3))
        result = load_layers(tmp_path / "data.npz")
        assert len(result) == 1

    def test_auto_detect_dir(self, tmp_path):
        """Auto format should detect directory."""
        sub = tmp_path / "layers"
        sub.mkdir()
        np.savetxt(sub / "a.csv", np.eye(3), delimiter=",")
        result = load_layers(sub)
        assert len(result) == 1

    def test_unknown_format(self, tmp_path):
        """Should raise for unknown extension."""
        (tmp_path / "data.xyz").touch()
        with pytest.raises(ValueError, match="Cannot detect"):
            load_layers(tmp_path / "data.xyz")


class TestSaveResults:
    def _make_result(self):
        return MXONMTFResult(
            H=np.eye(4, 2),
            Hl=[np.ones((4, 1))],
            labels_per_layer=[np.array([1, 1, 2, 2]), np.array([1, 2, 2, 1])],
            clusters_supra=np.array([1, 1, 2, 2, 1, 2, 2, 1]),
            averagedNMI=0.85,
        )

    def test_save_csv(self, tmp_path):
        """Should save labels to CSV."""
        result = self._make_result()
        out = tmp_path / "results.csv"
        save_results(result, out, format="csv")
        data = np.loadtxt(out, delimiter=",", skiprows=1)
        assert data.shape == (4, 2)

    def test_save_json(self, tmp_path):
        """Should save to JSON with labels and NMI."""
        result = self._make_result()
        out = tmp_path / "results.json"
        save_results(result, out, format="json")
        with open(out) as f:
            data = json.load(f)
        assert "labels_per_layer" in data
        assert "nmi" in data
        assert data["nmi"] == pytest.approx(0.85)

    def test_unknown_format(self, tmp_path):
        with pytest.raises(ValueError, match="Unknown format"):
            save_results(self._make_result(), tmp_path / "out.xyz", format="xyz")
