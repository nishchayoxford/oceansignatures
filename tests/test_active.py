"""Independent analytic checks for the stationary active-echo benchmark."""

import importlib
import importlib.util
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]


def active_module():
    spec = importlib.util.find_spec("oceansignatures.active")
    assert spec is not None, "active benchmark module is missing"
    module = importlib.import_module("oceansignatures.active")
    assert module.__file__ is not None
    assert Path(module.__file__).resolve() == ROOT / "src/oceansignatures/active.py"
    return module


def test_stationary_echo_end_to_end_recovers_analytic_two_way_range():
    case = active_module().active_case(range_m=150.0, snr_db=60.0)
    required = {
        "time_s",
        "transmitted_pa",
        "received_pa",
        "lag_s",
        "matched_filter",
        "estimated_range_m",
        "true_range_m",
        "fs",
        "metadata",
    }
    assert required <= case.keys()
    assert case["fs"] == 48000
    assert case["true_range_m"] == 150.0
    np.testing.assert_allclose(case["time_s"], np.arange(len(case["received_pa"])) / 48000)
    assert case["transmitted_pa"].shape == case["received_pa"].shape
    peak = np.argmax(np.abs(case["matched_filter"]))
    assert case["lag_s"][peak] == pytest.approx(2 * 150 / 1500, abs=1 / 48000)
    assert case["estimated_range_m"] == case["lag_s"][peak] * 1500 / 2
    assert abs(case["estimated_range_m"] - 150) <= 1500 / (2 * 48000)


def test_seeded_white_noise_uses_echo_support_mean_square_snr():
    case = active_module().active_case(range_m=150, snr_db=-18, seed=2026)
    fs = case["fs"]
    delay_samples = round(2 * 150 / 1500 * fs)
    ping_samples = round(0.05 * fs)
    expected_echo = np.zeros_like(case["received_pa"])
    expected_echo[delay_samples : delay_samples + ping_samples] = (
        case["transmitted_pa"][:ping_samples] / 150**2
    )
    np.testing.assert_allclose(case["received_clean_pa"], expected_echo, atol=3e-16)
    echo_power = np.mean(expected_echo[delay_samples : delay_samples + ping_samples] ** 2)
    sigma = np.sqrt(echo_power / 10 ** (-18 / 10))
    expected_noise = np.random.default_rng(2026).normal(0, sigma, len(expected_echo))
    np.testing.assert_allclose(case["noise_pa"], expected_noise, rtol=1e-12, atol=1e-18)
    np.testing.assert_array_equal(case["received_pa"], case["received_clean_pa"] + case["noise_pa"])
    metadata = case["metadata"]
    assert metadata["noise_std_pa"] == pytest.approx(sigma, rel=1e-12)
    assert metadata["input_snr_db"] == -18
    assert (
        metadata["input_snr_definition"]
        == "10*log10(mean(clean_echo_pa[echo_support]**2)/noise_variance_pa2)"
    )
    assert metadata["seed"] == 2026


def test_matched_filter_is_signed_valid_correlation_in_noise_sigma_units():
    case = active_module().active_case(range_m=75, snr_db=-6)
    fs = case["fs"]
    ping = case["transmitted_pa"][: round(0.05 * fs)]
    count = len(case["received_pa"]) - len(ping) + 1
    np.testing.assert_array_equal(case["lag_s"], np.arange(count) / fs)
    assert case["lag_s"][-1] == 1.0  # Fixed search through 750 m, not a truth-centred gate.
    divisor = case["metadata"]["noise_std_pa"] * np.sqrt(np.sum(ping**2))
    for index in (0, 101, 4799, 4800, 4801, count - 1):
        direct = np.dot(case["received_pa"][index : index + len(ping)], ping) / divisor
        assert case["matched_filter"][index] == pytest.approx(direct, abs=1e-12)
    assert (
        case["metadata"]["matched_filter_units"]
        == "dimensionless; noise-only standard deviation = 1 at each lag"
    )
    peak = int(np.argmax(np.abs(case["matched_filter"])))
    assert case["metadata"]["peak_index"] == peak
    assert case["metadata"]["peak_abs_score"] == abs(case["matched_filter"][peak])


@pytest.mark.parametrize("fs", [48000, 96000])
def test_fractional_delay_matches_continuous_lfm_and_reports_measured_metadata(fs):
    distance = 123.4567
    case = active_module().active_case(range_m=distance, snr_db=80, fs=fs)
    metadata = case["metadata"]
    duration, f0, f1 = 0.05, 3000.0, 9000.0
    tau = 2 * distance / 1500
    template_time = np.arange(round(duration * fs)) / fs
    raw_template = np.sin(np.pi * template_time / duration) ** 2 * np.cos(
        2 * np.pi * (f0 * template_time + (f1 - f0) * template_time**2 / (2 * duration))
    )
    normalization = 1 / np.sqrt(np.mean(raw_template**2))
    np.testing.assert_allclose(
        case["transmitted_pa"][: len(template_time)], raw_template * normalization, atol=8e-13
    )
    assert np.count_nonzero(case["transmitted_pa"][len(template_time) :]) == 0
    local_time = case["time_s"] - tau
    support = (local_time >= 0) & (local_time < duration)
    expected = np.zeros_like(local_time)
    u = local_time[support]
    expected[support] = (
        normalization
        * np.sin(np.pi * u / duration) ** 2
        * np.cos(2 * np.pi * (f0 * u + (f1 - f0) * u**2 / (2 * duration)))
        / distance**2
    )
    np.testing.assert_allclose(case["received_clean_pa"], expected, atol=5e-16)
    assert abs(case["estimated_range_m"] - distance) <= 1500 / (4 * fs)
    assert metadata["expected_delay_s"] == tau
    assert metadata["sound_speed_m_s"] == 1500
    assert metadata["target_scalar_gain"] == 1.0
    assert metadata["two_way_pressure_gain"] == pytest.approx(1 / distance**2)
    assert metadata["reference_distance_m"] == 1.0
    assert metadata["ping_duration_s"] == duration
    assert metadata["ping_f_start_hz"] == f0
    assert metadata["ping_f_end_hz"] == f1
    assert metadata["ping_taper"] == "continuous Hann"
    assert metadata["ping_samples"] == len(template_time)
    assert metadata["ping_rms_pa"] == pytest.approx(1.0)
    assert metadata["record_samples"] == len(case["time_s"])
    assert metadata["range_bin_m"] == 1500 / (2 * fs)
    assert metadata["search_lag_s"] == [0.0, 1.0]
    assert metadata["range_error_m"] == case["estimated_range_m"] - distance
    start, stop = metadata["echo_support_samples"]
    np.testing.assert_array_equal(np.arange(start, stop), np.flatnonzero(support))
    clean_power = np.mean(case["received_clean_pa"][start:stop] ** 2)
    noise_power = np.mean(case["noise_pa"][start:stop] ** 2)
    assert metadata["echo_mean_square_pa2"] == clean_power
    assert metadata["realized_noise_mean_square_pa2"] == noise_power
    assert metadata["realized_input_snr_db"] == pytest.approx(
        10 * np.log10(clean_power / noise_power)
    )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"range_m": 0},
        {"range_m": -1},
        {"range_m": 0.5},
        {"range_m": 750.1},
        {"range_m": np.nan},
        {"range_m": np.inf},
        {"range_m": True},
        {"snr_db": np.nan},
        {"snr_db": np.inf},
        {"snr_db": -121},
        {"snr_db": 121},
        {"snr_db": "-6"},
        {"fs": 0},
        {"fs": 18000},
        {"fs": 23999},
        {"fs": 48000.5},
        {"fs": True},
        {"fs": 384001},
        {"seed": -1},
        {"seed": 1.5},
        {"seed": True},
    ],
)
def test_invalid_benchmark_parameters_raise_clear_value_error(kwargs):
    with pytest.raises(ValueError, match=next(iter(kwargs))):
        active_module().active_case(**kwargs)


# Verification-only coverage of the now-implemented analytic model: these are
# not represented as new RED/GREEN feature slices in the retained logs.
@pytest.mark.parametrize("distance", [75.0, 150.0, 300.0])
@pytest.mark.parametrize("snr", [-18.0, -6.0])
def test_prescribed_low_snr_examples_recover_delay_without_truth_gating(distance, snr):
    case = active_module().active_case(distance, snr, seed=2026)
    tau = 2 * distance / 1500
    assert case["metadata"]["expected_delay_s"] == tau
    assert abs(case["estimated_range_m"] - distance) <= 0.125
    assert case["lag_s"][np.argmax(abs(case["matched_filter"]))] == pytest.approx(
        tau, abs=0.125 * 2 / 1500
    )
    assert 10 * np.log10(
        case["metadata"]["echo_mean_square_pa2"] / case["metadata"]["noise_std_pa"] ** 2
    ) == pytest.approx(snr)
    assert case["metadata"]["search_lag_s"] == [0.0, 1.0]


@pytest.mark.parametrize("distance,fs", [(1, 24000), (750, 48001), (150, 384000)])
def test_valid_domain_boundaries_retain_complete_echo(distance, fs):
    case = active_module().active_case(distance, 80, fs=fs)
    assert case["lag_s"][-1] == 1.0
    start, stop = case["metadata"]["echo_support_samples"]
    assert 0 <= start < stop <= len(case["received_pa"])
    assert np.all(np.isfinite(case["matched_filter"]))
    assert abs(case["estimated_range_m"] - distance) <= 0.125


def test_default_case_is_reproducible_without_touching_global_random_state():
    state_before = np.random.get_state()
    assert isinstance(state_before, tuple)
    first = active_module().active_case()
    second = active_module().active_case()
    different = active_module().active_case(seed=2027)
    state_after = np.random.get_state()
    assert isinstance(state_after, tuple)
    assert state_before[0] == state_after[0]
    np.testing.assert_array_equal(state_before[1], state_after[1])
    assert state_before[2:] == state_after[2:]
    for key in ("received_pa", "noise_pa", "matched_filter", "transmitted_pa"):
        np.testing.assert_array_equal(first[key], second[key])
    assert first["metadata"] == second["metadata"]
    assert not np.array_equal(first["noise_pa"], different["noise_pa"])
    np.testing.assert_array_equal(first["received_clean_pa"], different["received_clean_pa"])


def test_fractional_delay_sweep_has_only_grid_quantization_at_high_snr():
    fs = 48000
    for fractional_samples in np.linspace(0.01, 0.99, 31):
        distance = (9600 + fractional_samples) * 1500 / (2 * fs)
        case = active_module().active_case(distance, 100, fs=fs)
        assert abs(case["estimated_range_m"] - distance) <= 1500 / (4 * fs) + 1e-12


def test_small_seeded_estimation_study_retains_low_snr_failures():
    # Conditional on a target being present; this is not a detection-rate test.
    errors = {}
    for snr in (-30, -18, -6):
        errors[snr] = np.array(
            [
                active_module().active_case(150, snr, seed=seed)["estimated_range_m"] - 150
                for seed in range(2026, 2050)
            ]
        )
    assert np.all(np.abs(errors[-6]) <= 0.125)
    assert np.count_nonzero(np.abs(errors[-30]) > 0.125) >= 12
    assert np.count_nonzero(np.abs(errors[-18]) <= 0.125) > np.count_nonzero(
        np.abs(errors[-30]) <= 0.125
    )
    assert np.count_nonzero(np.abs(errors[-18]) > 0.125) >= 1
    assert np.sqrt(np.mean(errors[-18] ** 2)) > 1.0  # Keep the catastrophic outlier in RMSE.


def test_metadata_carries_model_and_estimation_assumptions():
    metadata = active_module().active_case()["metadata"]
    assert (
        metadata["model"] == "stationary monostatic scalar point target; homogeneous 3-D free space"
    )
    assert (
        metadata["range_definition"]
        == "horizontal separation; target and transceiver at equal depth"
    )
    assert (
        metadata["receiver_model"]
        == "ideal pressure sampler; no receiver filter, clipping or transmit leakage"
    )
    assert (
        metadata["noise_model"]
        == "iid zero-mean real Gaussian; full sampled Nyquist band; not realization-normalized"
    )
    assert (
        metadata["fractional_delay_method"]
        == "evaluate continuous finite LFM at (sample_index - delay_s*fs)/fs"
    )
    assert (
        metadata["lag_definition"]
        == "template start lag k/fs; valid correlation; no filter-length offset"
    )
    assert (
        metadata["range_estimator"]
        == "global maximum absolute signed correlation; integer lag; no detection threshold"
    )
    assert metadata["omitted_physics"] == [
        "surface and seabed multipath",
        "refraction",
        "absorption",
        "Doppler",
        "frequency-dependent or extended-target scattering",
        "reverberation",
    ]
