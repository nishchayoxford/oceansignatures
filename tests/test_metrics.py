"""Independent normalization tests for real pressure histories."""

import importlib
import importlib.util
import numpy as np
import pytest


def module():
    assert importlib.util.find_spec("oceansignatures.metrics") is not None, "metrics is missing"
    return importlib.import_module("oceansignatures.metrics")


def test_rms_preserves_dc_and_pressure_units():
    assert module().rms([3.0, 4.0]) == pytest.approx(np.sqrt(12.5))
    assert module().rms([2.0, 2.0]) == 2.0
    assert module().rms([0.0, 0.0]) == 0.0
    assert module().rms([1e300, -1e300]) == pytest.approx(1e300)
    for bad in ([], [[1.0, 2.0]], [np.nan], [np.inf], [1j]):
        with pytest.raises(ValueError):
            module().rms(bad)


def test_welch_matches_independent_hann_periodograms():
    assert hasattr(module(), "welch_psd"), "welch_psd is missing"
    fs, n = 8192.0, 256
    x = np.random.default_rng(33).normal(size=1024) + 4
    f, actual = module().welch_psd(x, fs, n)
    w = 0.5 - 0.5 * np.cos(2 * np.pi * np.arange(n) / n)
    expected = []
    for start in range(0, len(x) - n + 1, n // 2):
        block = x[start : start + n] - np.mean(x[start : start + n])
        p = np.abs(np.fft.rfft(block * w)) ** 2 / (fs * np.sum(w * w))
        p[1:-1] *= 2
        expected.append(p)
    np.testing.assert_allclose(f, np.arange(n // 2 + 1) * fs / n)
    np.testing.assert_allclose(actual, np.mean(expected, axis=0), rtol=1e-13)
    for args in ((x, 0, n), (x, fs, 1), (x, fs, 1.5), (x[:8], fs, 16)):
        with pytest.raises(ValueError):
            module().welch_psd(*args)


@pytest.mark.parametrize("n", [255, 256])
def test_spectrogram_parseval_and_half_open_frame_support(n):
    assert hasattr(module(), "spectrogram_psd"), "spectrogram_psd is missing"
    fs, overlap = 4096.0, 64
    x = np.random.default_rng(99).normal(size=1003) + 3
    t, f, p = module().spectrogram_psd(x, fs, n, overlap)
    starts = np.arange(0, len(x) - n + 1, n - overlap)
    np.testing.assert_allclose(t, (starts + n / 2) / fs)
    assert p.shape == (n // 2 + 1, len(starts))
    w = 0.5 - 0.5 * np.cos(2 * np.pi * np.arange(n) / n)
    for i, start in enumerate(starts):
        block = x[start : start + n] - np.mean(x[start : start + n])
        expected = np.abs(np.fft.rfft(block * w)) ** 2 / (fs * np.sum(w * w))
        expected[1 : -1 if n % 2 == 0 else None] *= 2
        np.testing.assert_allclose(p[:, i], expected, rtol=1e-12, atol=1e-17)
        assert np.sum(p[:, i]) * fs / n == pytest.approx(np.sum((block * w) ** 2) / np.sum(w * w))
    for bad in (-1, n, 1.5, True):
        with pytest.raises(ValueError):
            module().spectrogram_psd(x, fs, n, bad)
