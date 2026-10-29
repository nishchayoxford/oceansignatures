"""Tests for faithful plot reduction and analytic channel-reference curves."""

import importlib
import importlib.util

import numpy as np
import pytest


def reporting():
    assert importlib.util.find_spec("oceansignatures.reporting") is not None, "missing plot helpers"
    return importlib.import_module("oceansignatures.reporting")


def test_envelope_preserves_each_block_extremum_including_a_single_sample_click():
    module = reporting()
    x = np.zeros(1000)
    x[501] = 12
    x[502] = -7
    t, low, high = module.envelope(np.arange(1000) / 1000, x, max_bins=100)
    assert len(t) == 100
    assert np.max(high) == 12
    assert np.min(low) == -7
    assert np.all(low <= high)
    for index in range(100):
        assert high[index] == np.max(x[index * 10 : (index + 1) * 10])
        assert low[index] == np.min(x[index * 10 : (index + 1) * 10])


def test_channel_reference_has_correct_direct_phase_and_surface_dc_cancellation():
    module = reporting()
    assert callable(getattr(module, "channel_transfer", None)), "missing channel reference"
    f = np.array([0.0, 75.0, 150.0])
    direct = module.channel_transfer(f, 3.0, source_depth_m=1, receiver_depth_m=5, surface=False)
    np.testing.assert_allclose(direct, [0.2, -0.2j, -0.2], atol=1e-15)
    reflected = module.channel_transfer(np.array([0.0]), 3.0, source_depth_m=1, receiver_depth_m=5)
    np.testing.assert_allclose(reflected.real, [1 / 5 - 1 / np.sqrt(45)])


def test_envelope_rejects_mismatched_or_empty_data():
    with pytest.raises(ValueError):
        reporting().envelope([], [])
    with pytest.raises(ValueError):
        reporting().envelope([0, 1], [0])
