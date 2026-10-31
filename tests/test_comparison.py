"""Independent comparison controls; fixtures are not observations or CFD."""

import importlib
import importlib.util

import numpy as np
import pytest


def comparison_module():
    assert importlib.util.find_spec("oceansignatures.comparison") is not None, (
        "prediction/observation comparison feature is missing"
    )
    return importlib.import_module("oceansignatures.comparison")


def spectrum_pair():
    prediction = {
        "frequency_hz": np.arange(5, dtype=float) * 2,
        "window_start_s": np.array([0.0, 2.0, 4.0]),
        "window_stop_s": np.array([2.0, 4.0, 6.0]),
        "psd": np.tile(np.arange(1.0, 6.0), (3, 1)),
        "psd_unit": "Pa2/Hz",
        "evidence_kind": "predicted_acoustic",
        "time_origin_utc": "2026-01-01T00:00:00Z",
        "valid_window": np.ones(3, dtype=bool),
    }
    observation = {
        key: value.copy() if isinstance(value, np.ndarray) else value
        for key, value in prediction.items()
    }
    observation["evidence_kind"] = "synthetic_observation"
    observation["psd"] *= 4
    return prediction, observation


def test_known_power_ratio_has_signed_db_residual_and_exact_support():
    prediction, observation = spectrum_pair()
    result = comparison_module().compare_prediction_observation(prediction, observation)
    expected = 20 * np.log10(2.0)
    np.testing.assert_allclose(result["observed_minus_predicted_db"], expected, rtol=1e-14)
    np.testing.assert_array_equal(result["linear_residual_pa2_hz"], 3 * prediction["psd"])
    assert result["summary"]["mean_bias_db"] == pytest.approx(expected)
    assert result["summary"]["rmse_db"] == pytest.approx(expected)
    assert result["summary"]["compared_windows"] == 3
    np.testing.assert_array_equal(result["window_start_s"], prediction["window_start_s"])
    np.testing.assert_array_equal(result["frequency_hz"], prediction["frequency_hz"])


@pytest.mark.parametrize(
    "key, value",
    [
        ("frequency_hz", np.arange(5, dtype=float) * 2 + 0.1),
        ("window_start_s", np.array([0.1, 2.0, 4.0])),
        ("window_stop_s", np.array([2.1, 4.0, 6.0])),
        ("psd_unit", "digital2/Hz"),
        ("evidence_kind", "hydrodynamic_pressure"),
    ],
)
def test_comparison_rejects_incompatible_support_units_or_evidence(key, value):
    prediction, observation = spectrum_pair()
    observation[key] = value
    with pytest.raises(ValueError, match="aligned|Pa2/Hz|acoustic observation"):
        comparison_module().compare_prediction_observation(prediction, observation)


def test_flagged_window_is_retained_but_excluded_from_errors():
    prediction, observation = spectrum_pair()
    observation["valid_window"][1] = False
    observation["psd"][1] = np.nan
    result = comparison_module().compare_prediction_observation(prediction, observation)
    np.testing.assert_array_equal(result["compared_window"], [True, False, True])
    assert np.isnan(result["observed_minus_predicted_db"][1]).all()
    assert np.isnan(result["linear_residual_pa2_hz"][1]).all()
    assert result["summary"]["compared_windows"] == 2
    assert result["summary"]["excluded_windows"] == 1
    assert result["summary"]["mean_bias_db"] == pytest.approx(20 * np.log10(2.0))


@pytest.mark.parametrize(
    "key, value",
    [
        ("frequency_hz", np.array([0.0, 2.0, 1.0, 6.0, 8.0])),
        ("window_start_s", np.array([0.0, 2.0, np.inf])),
        ("window_stop_s", np.array([0.0, 4.0, 6.0])),
        ("psd", np.ones((3, 4))),
        ("valid_window", np.array([1, 0, 1])),
    ],
)
def test_malformed_batches_rejected_before_arithmetic(key, value):
    prediction, observation = spectrum_pair()
    prediction[key] = value.copy()
    observation[key] = value.copy()
    with pytest.raises(ValueError):
        comparison_module().compare_prediction_observation(prediction, observation)


def test_physical_zero_is_not_replaced_by_a_logarithmic_floor():
    prediction, observation = spectrum_pair()
    prediction["psd"][0, [0, 3]] = 0
    observation["psd"][0, [0, 1]] = 0
    result = comparison_module().compare_prediction_observation(prediction, observation)
    delta = result["observed_minus_predicted_db"]
    assert np.isnan(delta[0, 0])
    assert np.isneginf(delta[0, 1])
    assert np.isposinf(delta[0, 3])
    finite_pairs = (prediction["psd"] > 0) & (observation["psd"] > 0)
    np.testing.assert_array_equal(result["finite_db_comparison"], finite_pairs)
    assert result["summary"]["finite_db_pairs"] == int(np.count_nonzero(finite_pairs))
    assert result["summary"]["mean_bias_db"] == pytest.approx(20 * np.log10(2.0))


def test_no_valid_windows_returns_unavailable_errors_not_success():
    prediction, observation = spectrum_pair()
    observation["valid_window"][:] = False
    observation["psd"][:] = np.nan
    result = comparison_module().compare_prediction_observation(prediction, observation)
    assert result["summary"]["mean_bias_db"] is None
    assert result["summary"]["rmse_db"] is None
    assert result["summary"]["compared_windows"] == 0


def test_equal_relative_windows_with_different_utc_origins_do_not_align():
    prediction, observation = spectrum_pair()
    observation["time_origin_utc"] = "2026-01-02T00:00:00Z"
    with pytest.raises(ValueError, match="UTC origin"):
        comparison_module().compare_prediction_observation(prediction, observation)


@pytest.mark.parametrize("origin", [None, "", "2026-01-01T00:00:00", "2026-01-01T00:00:00+01:00"])
def test_shared_unknown_or_non_utc_origin_does_not_prove_alignment(origin):
    prediction, observation = spectrum_pair()
    prediction["time_origin_utc"] = observation["time_origin_utc"] = origin
    with pytest.raises(ValueError, match="UTC origin"):
        comparison_module().compare_prediction_observation(prediction, observation)


def test_band_comparison_integrates_linear_density_before_decibels():
    prediction, observation = spectrum_pair()
    module = comparison_module()
    assert hasattr(module, "compare_band_levels"), "band comparison feature is missing"
    result = module.compare_band_levels(prediction, observation, {"low": (0, 4), "high": (4, 8)})
    expected = np.tile(np.sum([[1.0, 2.0], [3.0, 4.0]], axis=1) * 2, (3, 1))
    np.testing.assert_array_equal(result["predicted_power_pa2"], expected)
    np.testing.assert_array_equal(result["observed_power_pa2"], expected * 4)
    np.testing.assert_allclose(result["observed_minus_predicted_db"], 20 * np.log10(2.0))
    np.testing.assert_allclose(
        result["predicted_level_db_re_1upa"], 10 * np.log10(expected / 1e-12)
    )
    assert result["band_names"] == ["low", "high"]


@pytest.mark.parametrize("bands", [{}, {"bad": (-1, 4)}, {"bad": (0, 10)}, {"": (0, 4)}])
def test_band_comparison_rejects_undefined_or_outside_support(bands):
    prediction, observation = spectrum_pair()
    with pytest.raises(ValueError, match="band"):
        comparison_module().compare_band_levels(prediction, observation, bands)


def test_band_comparison_rejects_nonuniform_frequency_spacing():
    prediction, observation = spectrum_pair()
    prediction["frequency_hz"][3] = observation["frequency_hz"][3] = 7
    with pytest.raises(ValueError, match="uniform"):
        comparison_module().compare_band_levels(prediction, observation, {"low": (0, 4)})


def test_band_zero_and_missing_windows_keep_distinct_meanings():
    prediction, observation = spectrum_pair()
    prediction["psd"][0, :2] = observation["psd"][0, :2] = 0
    observation["valid_window"][1] = False
    observation["psd"][1] = np.nan
    result = comparison_module().compare_band_levels(prediction, observation, {"low": (0, 4)})
    assert result["predicted_power_pa2"][0, 0] == 0
    assert np.isneginf(result["predicted_level_db_re_1upa"][0, 0])
    assert np.isnan(result["observed_minus_predicted_db"][:2]).all()
    assert np.isfinite(result["predicted_power_pa2"][1, 0])
    assert np.isnan(result["observed_power_pa2"][1, 0])


def test_band_power_avoids_intermediate_sum_overflow():
    prediction, observation = spectrum_pair()
    for batch in (prediction, observation):
        batch["frequency_hz"] = np.arange(5) * 1e-308
        batch["psd"][:] = 1e308
    result = comparison_module().compare_band_levels(prediction, observation, {"low": (0, 2e-308)})
    np.testing.assert_allclose(result["predicted_power_pa2"], 2.0, rtol=1e-14)


def test_unrepresentable_band_power_is_rejected_cleanly():
    prediction, observation = spectrum_pair()
    prediction["psd"][:] = observation["psd"][:] = 1e308
    with pytest.raises(ValueError, match="band power"):
        comparison_module().compare_band_levels(prediction, observation, {"low": (0, 4)})
