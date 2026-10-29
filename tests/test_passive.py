"""Numerical controls for the ideal point-source passive laboratory."""

import importlib
import importlib.util
import numpy as np
import pytest


def module():
    assert importlib.util.find_spec("oceansignatures.passive") is not None, "passive is missing"
    return importlib.import_module("oceansignatures.passive")


def test_whistle_reference_rms_taper_and_frequency_sweep():
    x = module().make_source("whistle", fs=48000, duration_s=0.5, rms_pa=2.0, seed=7)
    assert len(x) == 24000
    assert np.sqrt(np.mean(x * x)) == pytest.approx(2.0, rel=1e-14)
    assert abs(x[0]) < 1e-14 and abs(x[-1]) < 1e-14
    np.testing.assert_array_equal(x, module().make_source("whistle", 48000, 0.5, 2.0, 7))

    def peak(block):
        return np.argmax(np.abs(np.fft.rfft(block * np.hanning(len(block))))) * 48000 / len(block)

    assert 1800 < peak(x[3000:6000]) < 2400
    assert 3400 < peak(x[18000:21000]) < 4100
    for args in [
        ("unknown",),
        ("whistle", 8000),
        ("whistle", 48000, 0),
        ("whistle", 48000, 0.5, -1),
        ("whistle", 48000, 0.5, 1, -1),
    ]:
        with pytest.raises(ValueError):
            module().make_source(*args)


def test_click_train_is_seeded_sparse_and_reference_normalized():
    x = module().make_source("clicks", duration_s=0.5, seed=10)
    assert np.sqrt(np.mean(x * x)) == pytest.approx(1.0, rel=1e-14)
    assert np.max(np.abs(x)) > 10
    assert np.count_nonzero(np.abs(x) > 0.1) / len(x) < 0.03
    assert x[0] == 0 and x[-1] == 0
    np.testing.assert_array_equal(x, module().make_source("clicks", duration_s=0.5, seed=10))
    assert not np.array_equal(x, module().make_source("clicks", duration_s=0.5, seed=11))
    f = np.fft.rfftfreq(len(x), 1 / 48000)
    power = np.abs(np.fft.rfft(x)) ** 2
    assert power[(f > 5000) & (f < 11000)].sum() / power.sum() > 0.99


def test_vessel_has_harmonics_and_seeded_broadband():
    x = module().make_source("vessel", duration_s=1.0, seed=10)
    assert np.sqrt(np.mean(x * x)) == pytest.approx(1.0, rel=1e-14)
    assert x[0] == 0 and x[-1] == 0
    np.testing.assert_array_equal(x, module().make_source("vessel", duration_s=1.0, seed=10))
    assert not np.array_equal(x, module().make_source("vessel", duration_s=1.0, seed=11))
    power = np.abs(np.fft.rfft(x * np.hanning(len(x)))) ** 2
    for harmonic in (80, 160, 240, 320):
        assert power[harmonic] > 20 * np.median(power[harmonic - 30 : harmonic + 30])
    assert power[2000:6000].sum() / power.sum() > 0.001


def test_paths_use_slant_distance_signed_image_and_one_metre_gain():
    assert hasattr(module(), "path_parameters"), "path_parameters is missing"
    paths = module().path_parameters(25.0, source_depth_m=8.0, receiver_depth_m=12.0)
    assert [p["kind"] for p in paths] == ["direct", "surface"]
    for p, distance, sign in zip(paths, [np.sqrt(25**2 + 4**2), np.sqrt(25**2 + 20**2)], [1, -1]):
        assert p["distance_m"] == pytest.approx(distance)
        assert p["gain"] == pytest.approx(sign / distance)
        assert p["delay_s"] == pytest.approx(distance / 1500)
    assert module().path_parameters(25.0, surface=False) == paths[:1]
    for kwargs in (
        {"range_m": 0},
        {"range_m": -1},
        {"range_m": 1, "sound_speed": 0},
        {"range_m": 1, "source_depth_m": -1},
        {"range_m": 1, "surface": "false"},
    ):
        with pytest.raises(ValueError):
            module().path_parameters(**kwargs)


def test_propagation_retains_integer_delay_tail_and_signed_superposition():
    assert hasattr(module(), "propagate"), "propagate is missing"
    x = np.array([1.0, 0.0, -2.0])
    y = module().propagate(
        x, 1000.0, [{"delay_s": 0.1, "gain": 0.5}, {"delay_s": 0.2, "gain": -0.25}]
    )
    expected = np.zeros(len(y))
    expected[100:103] += 0.5 * x
    expected[200:203] -= 0.25 * x
    assert len(y) >= 203
    np.testing.assert_array_equal(y, expected)
    for paths in ([], [{"delay_s": -1, "gain": 1}], [{"delay_s": 0, "gain": np.nan}]):
        with pytest.raises(ValueError):
            module().propagate(x, 1000.0, paths)


@pytest.mark.parametrize("frequency", [300.0, 2000.0, 8000.0, 12000.0, 19200.0])
@pytest.mark.parametrize("fraction", [0.05, 0.25, 0.5, 0.9])
def test_fractional_delay_signed_phase_and_amplitude(frequency, fraction):
    fs, delay = 48000.0, 131 + fraction
    n = np.arange(4096)
    x = np.sin(2 * np.pi * frequency * n / fs + 0.3)
    y = module().propagate(x, fs, [{"delay_s": delay / fs, "gain": -0.7}])
    indices = np.arange(256, 4000)
    expected = -0.7 * np.sin(2 * np.pi * frequency * (indices - delay) / fs + 0.3)
    assert np.linalg.norm(y[indices] - expected) / np.linalg.norm(expected) < 5e-5
    assert abs(np.linalg.norm(y[indices]) / np.linalg.norm(expected) - 1) < 5e-5


def test_passive_case_applies_receiver_and_reports_actual_window_metrics():
    assert hasattr(module(), "passive_case"), "passive_case is missing"
    case = module().passive_case("clicks", 100.0, True, duration_s=0.3, seed=22)
    metadata = case["metadata"]
    n = len(case["time_s"])
    for key in ("source_pa", "received_clean_pa", "received_pa", "received_channel_pa", "noise_pa"):
        assert case[key].shape == (n,)
        assert np.all(np.isfinite(case[key]))
    np.testing.assert_array_equal(case["time_s"], np.arange(n) / 48000)
    source = module().make_source("clicks", duration_s=0.3, seed=22)
    channel = module().propagate(source, 48000, metadata["paths"])
    taps = np.asarray(metadata["receiver"]["fir_coefficients"])
    assert len(taps) == 257
    assert metadata["receiver"]["group_delay_samples"] == 128
    np.testing.assert_allclose(case["received_clean_pa"], np.convolve(channel, taps), atol=1e-17)
    np.testing.assert_allclose(
        case["received_pa"], case["received_clean_pa"] + case["noise_pa"], atol=1e-17
    )
    assert np.std(case["noise_pa"][:100]) > 0
    a, b = metadata["analysis_window_samples"]
    measured = metadata["measured"]
    clean_rms = np.sqrt(np.mean(case["received_clean_pa"][a:b] ** 2))
    noise_rms = np.sqrt(np.mean(case["noise_pa"][a:b] ** 2))
    assert measured["received_clean_rms_pa"] == pytest.approx(clean_rms)
    assert measured["noise_rms_pa"] == pytest.approx(noise_rms)
    assert measured["snr_db"] == pytest.approx(20 * np.log10(clean_rms / noise_rms))
    assert measured["source_rms_pa"] == pytest.approx(1.0)
    assert b == n and 0 < a < b
    assert len(case["source_pa"]) > len(source) + metadata["paths"][-1]["delay_s"] * 48000
    import json

    json.dumps(metadata, allow_nan=False)
    other = module().passive_case("clicks", 100.0, True, duration_s=0.3, seed=22)
    np.testing.assert_array_equal(case["received_pa"], other["received_pa"])


# Verification-only extensions: these audit implemented behavior, not new RED slices.
@pytest.mark.parametrize("kind", ["whistle", "clicks", "vessel"])
@pytest.mark.parametrize("range_m", [25.0, 100.0, 400.0])
@pytest.mark.parametrize("surface", [False, True])
def test_all_default_cases_are_finite_repeatable_and_complete(kind, range_m, surface):
    case = module().passive_case(kind, range_m, surface)
    repeated = module().passive_case(kind, range_m, surface)
    n = len(case["time_s"])
    assert (
        n
        == 96000
        + int(np.ceil(max(p["delay_s"] for p in case["metadata"]["paths"]) * 48000))
        + 64
        + 256
    )
    assert len(case["metadata"]["paths"]) == 1 + int(surface)
    for key in ("source_pa", "received_pa", "received_clean_pa", "received_channel_pa", "noise_pa"):
        assert case[key].shape == (n,)
        assert np.all(np.isfinite(case[key]))
        np.testing.assert_array_equal(case[key], repeated[key])
    assert case["metadata"]["measured"]["source_rms_pa"] == pytest.approx(1.0, abs=1e-14)


def test_fractional_impulse_keeps_tail_and_does_not_wrap():
    x = np.zeros(300)
    x[-1] = 1
    y = module().propagate(x, 48000, [{"delay_s": 100.5 / 48000, "gain": 1.0}])
    assert np.count_nonzero(y[:335]) == 0
    assert np.count_nonzero(y[400:]) > 50
    assert np.sum(y) == pytest.approx(1.0, rel=1e-14)
    assert len(y) == 465


def test_surface_pressure_release_cancels_at_zero_source_depth():
    x = np.random.default_rng(2).normal(size=1000)
    paths = module().path_parameters(100, source_depth_m=0.0)
    np.testing.assert_array_equal(module().propagate(x, 48000, paths), 0.0)


def test_surface_two_path_tone_matches_independent_signed_phasors():
    fs, f = 48000.0, 1373.0
    x = np.sin(2 * np.pi * f * np.arange(10000) / fs)
    paths = module().path_parameters(25.0)
    y = module().propagate(x, fs, paths)
    indices = np.arange(2000, 9000)
    distances = [np.sqrt(25**2 + 4**2), np.sqrt(25**2 + 20**2)]
    expected = sum(
        sign / d * np.sin(2 * np.pi * f * (indices / fs - d / 1500))
        for sign, d in zip([1, -1], distances)
    )
    np.testing.assert_allclose(y[indices], expected, rtol=0, atol=2e-6)


def test_sources_do_not_modify_global_rng_and_zero_rms_is_silence():
    np.random.seed(123)
    before = np.random.get_state()
    for kind in ("whistle", "clicks", "vessel"):
        np.testing.assert_array_equal(module().make_source(kind, duration_s=0.1, rms_pa=0), 0.0)
    after = np.random.get_state()
    assert before[0] == after[0] and before[2:] == after[2:]
    np.testing.assert_array_equal(before[1], after[1])
