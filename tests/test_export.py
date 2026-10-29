"""Executable release exports: real simulations, not mocked catalogue fixtures."""

import importlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from oceansignatures.passive import passive_case
from oceansignatures.metrics import welch_psd


def exporter():
    assert importlib.util.find_spec("oceansignatures.export") is not None, (
        "scientific exporter is missing"
    )
    return importlib.import_module("oceansignatures.export")


@pytest.mark.parametrize(
    "family, parameters",
    [("active", (150.5, -18)), ("active", (150, -18.5)), ("passive", ("whistle", 100.5, False))],
)
def test_fractional_case_downloads_preserve_the_complete_id(tmp_path, family, parameters):
    module = exporter()
    entry = getattr(module, f"export_{family}_case")(tmp_path, *parameters)
    for download in entry["downloads"]:
        path = tmp_path / download["href"]
        assert path.is_file(), f"Missing declared download: {path.name}"
        assert path.name.startswith(entry["id"])
    metadata = json.loads((tmp_path / "data" / (entry["id"] + ".json")).read_text())
    assert metadata["case_id"] == entry["id"]


@pytest.mark.parametrize(
    "family, index, parameters",
    [
        ("active", 0, (150.5, -18)),
        ("active", 1, (150, -18.5)),
        ("passive", 1, ("whistle", 100.5, False)),
    ],
)
def test_distinct_float_parameters_do_not_collide_or_replace_exports(
    tmp_path, family, index, parameters
):
    export = getattr(exporter(), f"export_{family}_case")
    first = export(tmp_path, *parameters)
    original = {d["href"]: (tmp_path / d["href"]).read_bytes() for d in first["downloads"]}
    adjacent = list(parameters)
    adjacent[index] = np.nextafter(parameters[index], np.inf)
    second = export(tmp_path, *adjacent)
    assert first["id"] != second["id"], "Distinct binary64 parameters collided in the case ID"
    assert set(original).isdisjoint(d["href"] for d in second["downloads"])
    assert all((tmp_path / name).read_bytes() == payload for name, payload in original.items())


def test_passive_export_preserves_exact_arrays_windows_and_spectra(tmp_path):
    module = exporter()
    assert hasattr(module, "export_passive_case"), "passive export feature is missing"
    entry = module.export_passive_case(tmp_path, "clicks", 100, True)
    case = passive_case("clicks", 100, True)
    assert entry["kind"] == "clicks" and entry["surface"] is True and entry["range_m"] == 100
    paths = {Path(d["href"]).suffix: tmp_path / d["href"] for d in entry["downloads"]}
    assert {".npz", ".csv", ".json"} <= paths.keys()
    arrays = np.load(paths[".npz"], allow_pickle=False)
    for key, value in case.items():
        if isinstance(value, np.ndarray):
            np.testing.assert_array_equal(arrays[key], value)
    metadata = json.loads(paths[".json"].read_text())
    assert metadata["model_metadata"] == case["metadata"]
    assert metadata["spectral_method"]["detrend"] == "constant"
    assert metadata["spectral_method"]["source_window_samples"] == [0, 96000]
    assert (
        metadata["spectral_method"]["receiver_window_samples"]
        == case["metadata"]["analysis_window_samples"]
    )
    start, stop = case["metadata"]["analysis_window_samples"]
    f, p = welch_psd(case["received_pa"][start:stop], case["fs"])
    table = np.genfromtxt(paths[".csv"], delimiter=",", names=True)
    np.testing.assert_array_equal(table["frequency_hz"], f)
    np.testing.assert_array_equal(table["received_pa2_hz"], p)
    assert all(np.isfinite(arrays[key]).all() for key in arrays.files)
    first = {p.name: p.read_bytes() for p in paths.values()}
    module.export_passive_case(tmp_path, "clicks", 100, True)
    assert first == {p.name: p.read_bytes() for p in paths.values()}


def test_case_figures_are_real_readable_exports(tmp_path):
    from PIL import Image

    module = exporter()
    assert hasattr(module, "render_case"), "case figure renderer is missing"
    for entry, expected in (
        (module.export_passive_case(tmp_path, "whistle", 100, True), 3),
        (module.export_active_case(tmp_path, 150, -18), 2),
    ):
        figures = module.render_case(tmp_path, entry)
        assert len(figures) == expected
        for fig in figures:
            assert set(("src", "title", "caption", "alt")) <= fig.keys()
            assert len(fig["caption"]) > 100
            image = Image.open(tmp_path / fig["src"])
            assert image.size[0] >= 1400 and image.size[1] >= 850
            assert np.asarray(image).std() > 20
    with np.load(tmp_path / "data/passive-whistle-100m-surface-true-spectrogram.npz") as arrays:
        assert {
            "source_psd_pa2_hz",
            "received_psd_pa2_hz",
            "source_frame_s",
            "received_frame_s",
        } <= set(arrays.files)
        assert all(np.isfinite(arrays[k]).all() for k in arrays.files)


def test_full_public_builder_coverage_schema_downloads_and_determinism(tmp_path):
    import hashlib
    import itertools
    import subprocess
    import sys
    import zipfile

    module = exporter()
    assert hasattr(module, "build_site"), "complete public site builder is missing"
    root = Path(__file__).resolve().parents[1]
    first = tmp_path / "first"
    first.mkdir()
    (first / "index.html").write_text("preserve existing frontend")
    catalog = module.build_site(first)
    assert (first / "index.html").read_text() == "preserve existing frontend"
    assert {(e["kind"], e["range_m"], e["surface"]) for e in catalog["passive_cases"]} == set(
        itertools.product(("whistle", "clicks", "vessel"), (25, 100, 400), (False, True))
    )
    assert {(e["range_m"], e["snr_db"]) for e in catalog["active_cases"]} == set(
        itertools.product((75, 150, 300), (-18, -6))
    )
    assert catalog["reference"] is None
    assert len(catalog["gallery"]) == 5
    assert len(list((first / "plots").glob("*.png"))) == 71
    for family, count in (("passive", 3), ("active", 2)):
        for entry in catalog[family + "_cases"]:
            assert len(entry["figures"]) == count
            for download in entry["downloads"]:
                assert (first / download["href"]).is_file()
            for figure in entry["figures"]:
                assert (first / figure["src"]).is_file()
    for path in (first / "data").glob("*.npz"):
        with np.load(path, allow_pickle=False) as arrays:
            assert all(np.isfinite(arrays[key]).all() for key in arrays.files), path
    for download in catalog["downloads"]:
        assert not download["href"].startswith(("/", "http", ".cache"))
        assert (first / download["href"]).is_file()
    robustness = np.genfromtxt(first / "data/active-robustness.csv", delimiter=",", names=True)
    assert robustness.shape == (72,)
    from oceansignatures.active import active_case

    for row in robustness:
        actual = active_case(150, row["snr_db"], seed=int(row["seed"]))
        assert row["estimated_range_m"] == actual["estimated_range_m"]
    assert np.max(abs(robustness["error_m"])) > 250
    assert all(c["status"] == "passed" for c in catalog["verification"]["checks"])
    subprocess.run(
        [
            "node",
            "--input-type=module",
            "-e",
            'import fs from "node:fs"; import {validateCatalog} from '
            + json.dumps((root / "site/app.mjs").as_uri())
            + '; validateCatalog(JSON.parse(fs.readFileSync(process.argv[1], "utf8")));',
            str(first / "data/catalog.json"),
        ],
        check=True,
        capture_output=True,
    )
    pack = first / "data/selected-plates.pdf"
    assert pack.read_bytes().startswith(b"%PDF-")
    manifest = json.loads((first / "data/plot-pack.json").read_text())
    assert len(manifest["pages"]) == 10
    with zipfile.ZipFile(first / "data/full-waveforms.zip") as archive:
        assert len(archive.namelist()) == 24
        assert all(info.date_time == (1980, 1, 1, 0, 0, 0) for info in archive.infolist())
    second = tmp_path / "second"
    subprocess.run(
        [sys.executable, str(root / "scripts/build_site.py"), "--output", str(second)],
        check=True,
        capture_output=True,
        cwd=root,
    )

    def hashes(folder):
        return {
            str(p.relative_to(folder)): hashlib.sha256(p.read_bytes()).hexdigest()
            for name in ("data", "plots")
            for p in (folder / name).rglob("*")
            if p.is_file()
        }

    assert hashes(first) == hashes(second)


def test_passive_listening_pair_uses_one_disclosed_gain(tmp_path):
    from scipy.io import wavfile

    entry = exporter().export_passive_case(tmp_path, "clicks", 400, True)
    wavs = [d for d in entry["downloads"] if d["href"].endswith(".wav")]
    assert len(wavs) == 2, "same-gain synthetic listening pair is missing"
    metadata = json.loads((tmp_path / "data" / (entry["id"] + ".json")).read_text())
    listening = metadata["listening"]
    assert listening["physical_playback_calibration"] is False
    gain = listening["digital_fraction_per_pa"]
    case = passive_case("clicks", 400, True)
    for name in ("source", "received"):
        fs, digital = wavfile.read(tmp_path / listening[name + "_href"])
        assert fs == case["fs"] and digital.dtype == np.int16
        np.testing.assert_allclose(
            digital.astype(float) / 32767, case[name + "_pa"] * gain, atol=0.5 / 32767
        )
    assert gain > 0 and max(abs(case["source_pa"] * gain)) <= 0.80000000001


def test_interference_neutral_colour_means_zero_db(tmp_path, monkeypatch):
    module = exporter()

    class StopAfterFirstFigure(Exception):
        pass

    def inspect(output, stem, fig, *args, **kwargs):
        field = fig.axes[0].collections[0]
        assert field.norm(0.0) == 0.5, (
            "white must represent 0 dB, not the midpoint of asymmetric limits"
        )
        import matplotlib.pyplot as plt

        plt.close(fig)
        raise StopAfterFirstFigure

    monkeypatch.setattr(module, "_save_figure", inspect)
    import pytest

    with pytest.raises(StopAfterFirstFigure):
        module._gallery(tmp_path, [], [], None)


def test_gallery_comparison_scales_preserve_small_errors(tmp_path, monkeypatch):
    module = exporter()
    passive = [
        module.export_passive_case(tmp_path, kind, distance, surface)
        for kind in ("whistle", "clicks", "vessel")
        for distance in (25, 100, 400)
        for surface in (False, True)
    ]
    active = [
        module.export_active_case(tmp_path, distance, snr)
        for distance in (75, 150, 300)
        for snr in (-18, -6)
    ]
    observed = []

    def inspect(output, stem, fig, title, caption, *args, **kwargs):
        if stem == "gallery-active-robustness":
            assert fig.axes[1].yaxis.get_transform().linthresh <= 0.02, (
                "default symlog hides carrier-lobe errors"
            )
            plotted = np.concatenate(
                [collection.get_offsets() for collection in fig.axes[0].collections]
            )
            assert plotted.shape == (72, 2)
            assert np.max(abs(plotted[:, 1])) > 500
        observed.append(stem)
        import matplotlib.pyplot as plt

        plt.close(fig)
        return {"src": stem, "title": title, "caption": caption, "alt": title}

    monkeypatch.setattr(module, "_save_figure", inspect)
    module._gallery(tmp_path, passive, active, None)
    assert len(observed) == 5


def test_active_export_keeps_signed_score_and_echo_support(tmp_path):
    from oceansignatures.active import active_case

    module = exporter()
    assert hasattr(module, "export_active_case"), "active export feature is missing"
    entry = module.export_active_case(tmp_path, 150, -18)
    case = active_case(150, -18)
    paths = {Path(d["href"]).suffix: tmp_path / d["href"] for d in entry["downloads"]}
    arrays = np.load(paths[".npz"], allow_pickle=False)
    for key, value in case.items():
        if isinstance(value, np.ndarray):
            np.testing.assert_array_equal(arrays[key], value)
    assert arrays["matched_filter"].min() < 0 < arrays["matched_filter"].max()
    metadata = json.loads(paths[".json"].read_text())
    assert metadata["model_metadata"] == case["metadata"]
    assert (
        metadata["spectral_method"]["receiver_window_samples"]
        == case["metadata"]["echo_support_samples"]
    )
    assert entry["snr_db"] == -18 and entry["range_m"] == 150
    assert all(np.isfinite(arrays[key]).all() for key in arrays.files)
