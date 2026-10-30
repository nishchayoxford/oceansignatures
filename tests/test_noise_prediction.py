"""Independent numerical contracts for the environmental prediction extension."""

import importlib
import math

import numpy as np
import pytest


def model():
    return importlib.import_module("oceansignatures.noise_prediction")


def test_wind_receiver_psd_matches_published_equation():
    # Stojanovic Eq. 6, 1 kHz, 4 m/s: 65 - 40 log10(1.4).
    expected = 1e-12 * 10 ** ((65 - 40 * math.log10(1.4)) / 10)
    assert model().wind_psd([1000.0], wind_speed_m_s=4)[0] == pytest.approx(expected)


def test_wind_domain_and_calm_semantics():
    m = model()
    for bad in (-1, np.nan, np.inf, True, [4]):
        with pytest.raises(ValueError):
            m.wind_psd([1000], wind_speed_m_s=bad)
    for bad in ([0], [99], [100001], [np.nan], [], [[1000]], [1000, 100]):
        with pytest.raises(ValueError):
            m.wind_psd(bad, wind_speed_m_s=0)
    assert m.wind_psd([100], wind_speed_m_s=0)[0] > 0
    assert m.wind_psd(
        [20, 100], wind_speed_m_s=0, allow_low_frequency_extrapolation=True
    ).shape == (2,)


def test_jomopans_reference_spectra_and_scaling():
    m = model()
    # Independent algebra at f1=480/8=60 Hz, government baseline (Eq.10).
    base = m.vessel_source_psd([60, 1000], speed_kn=8, length_m=91.44, vessel_class="government")
    assert base[0] == pytest.approx(1e-12 * 10**19.1 / (60**2 * 9))
    faster = m.vessel_source_psd(
        [60, 1000], speed_kn=16, length_m=182.88, vessel_class="government"
    )
    np.testing.assert_allclose(faster / base, 256, rtol=1e-12)
    # Cargo low-frequency hump at its own f1=600/18.
    peak = m.vessel_source_psd(
        [600 / 18], speed_kn=18, length_m=91.44, vessel_class="containership"
    )[0]
    assert peak == pytest.approx(1e-12 * 10**20.8 / ((600 / 18) ** 3 * 0.8**2))
    for kwargs in (
        {"speed_kn": 0},
        {"speed_kn": np.inf},
        {"length_m": -1},
        {"vessel_class": "dredger"},
    ):
        args = dict(speed_kn=8, length_m=100, vessel_class="government") | kwargs
        with pytest.raises(ValueError):
            m.vessel_source_psd([1000], **args)
    for f in ([19], [20001], [np.nan]):
        with pytest.raises(ValueError):
            m.vessel_source_psd(f, speed_kn=8, length_m=100)


def test_spectral_reference_and_bandwidth_units():
    m = model()
    np.testing.assert_allclose(m.db_to_psd([0, 120]), [1e-12, 1])
    np.testing.assert_allclose(m.psd_to_db([1e-12, 1]), [0, 120])
    assert np.isneginf(m.psd_to_db([0])[0])
    np.testing.assert_allclose(m.psd_to_band_power([2, 3], [10, 100]), [20, 300])
    np.testing.assert_allclose(m.band_power_to_psd([20, 300], [10, 100]), [2, 3])
    assert m.fractional_octave_bandwidth([1000], fraction=3)[0] == pytest.approx(
        1000 * (10**0.05 - 10**-0.05)
    )
    for call in (
        lambda: m.db_to_psd([np.inf]),
        lambda: m.psd_to_db([-1]),
        lambda: m.psd_to_band_power([1], [0]),
        lambda: m.band_power_to_psd([np.nan], [10]),
        lambda: m.fractional_octave_bandwidth([1000], fraction=0),
    ):
        with pytest.raises(ValueError):
            call()


def test_absorption_and_spherical_range_benchmark():
    m = model()
    # Independent three-term evaluation at T=0, S=35, pH=8, depth=0.
    expected = 0.106 * 0.78 / (0.78**2 + 1) + 0.52 * 42 / (42**2 + 1) + 0.00049
    water = m.SeaWater(temperature_c=0, salinity_ppt=35, depth_m=0, ph=8)
    alpha = m.absorption_db_per_km([1000], water=water)[0]
    assert alpha == pytest.approx(expected, rel=1e-13)
    loss1 = m.propagation_loss([1000], distance_m=1000, water=water)[0]
    loss2 = m.propagation_loss([1000], distance_m=2000, water=water)[0]
    assert loss1 == pytest.approx(60 + expected)
    assert loss2 - loss1 == pytest.approx(20 * math.log10(2) + expected)
    for args in (
        {"temperature_c": 35},
        {"salinity_ppt": 5},
        {"ph": 8.3},
        {"depth_m": 7000},
        {"depth_m": -1},
        {"ph": np.nan},
    ):
        with pytest.raises(ValueError):
            m.SeaWater(**args)
    with pytest.raises(ValueError):
        m.absorption_db_per_km([20])
    assert m.absorption_db_per_km([20], allow_low_frequency_extrapolation=True)[0] > 0
    for distance in (0, -1, np.inf, np.nan):
        with pytest.raises(ValueError):
            m.propagation_loss([1000], distance_m=distance)


def test_passage_receiver_clock_slant_distance_and_linear_components():
    m = model()
    t = np.array([-20.0, 0, 20])
    kw = dict(
        speed_kn=10,
        length_m=100,
        closest_approach_m=300,
        source_depth_m=10,
        receiver_depth_m=410,
        sound_speed_m_s=1500,
        wind_speed_m_s=4,
    )
    result = m.predict_passage([100, 1000], t, **kw)
    np.testing.assert_array_equal(result["receiver_time_s"], t)
    emission = result["emission_time_s"]
    expected_distance = np.sqrt((10 * 1852 / 3600 * emission) ** 2 + 300**2 + 400**2)
    np.testing.assert_allclose(result["distance_m"], expected_distance, rtol=1e-13)
    np.testing.assert_allclose(emission + expected_distance / 1500, t, atol=1e-12)
    assert result["total_psd_pa2_hz"].shape == (3, 2)
    np.testing.assert_allclose(
        result["total_psd_pa2_hz"], result["wind_psd_pa2_hz"] + result["vessel_psd_pa2_hz"]
    )
    distant = m.predict_passage([100, 1000], t, **(kw | {"closest_approach_m": 600}))
    np.testing.assert_array_equal(distant["wind_psd_pa2_hz"], result["wind_psd_pa2_hz"])
    assert np.all(distant["vessel_psd_pa2_hz"] < result["vessel_psd_pa2_hz"])
    silence = m.predict_passage([100, 1000], t, **kw, include_wind=False, include_vessel=False)
    np.testing.assert_array_equal(silence["total_psd_pa2_hz"], 0)
    assert result["metadata"]["doppler"] == "omitted"
    assert result["metadata"]["clock"] == "receiver arrival time; retarded emission geometry"
    for update in (
        {"closest_approach_m": 0},
        {"source_depth_m": -1},
        {"sound_speed_m_s": np.nan},
        {"include_wind": 1},
    ):
        with pytest.raises(ValueError):
            m.predict_passage([1000], t, **(kw | update))
    for clock in ([1, 1], [np.inf], [], [[1]], [1, 0]):
        with pytest.raises(ValueError):
            m.predict_passage([1000], clock, **kw)


@pytest.mark.parametrize(
    "vessel_class,vc,d",
    [
        ("fishing", 6.4, 3),
        ("tug", 3.7, 3),
        ("naval", 11.1, 3),
        ("recreational", 10.6, 3),
        ("government", 8, 3),
        ("cruise", 17.1, 4),
        ("passenger", 9.7, 3),
        ("bulker", 13.9, 3),
        ("containership", 18, 3),
        ("vehicle_carrier", 15.8, 3),
        ("tanker", 12.4, 3),
        ("other", 7.4, 3),
    ],
)
def test_published_class_table_verification(vessel_class, vc, d):
    # Verification-only additional coverage after the class-law tracer bullet.
    expected = 1e-12 * 10**19.1 / ((480 / vc) ** 2 * ((1 - 100 / (480 / vc)) ** 2 + d * d))
    actual = model().vessel_source_psd(
        [100], speed_kn=vc, length_m=91.44, vessel_class=vessel_class
    )[0]
    assert actual == pytest.approx(expected, rel=2e-14)


def test_independent_environmental_absorption_and_clock_verification():
    m = model()
    # Scalar math oracle, nonzero pressure depth and all environmental factors.
    f, t, s, d, ph = 4.0, 18.0, 32.0, 1.5, 7.9
    b = 0.78 * math.sqrt(s / 35) * math.exp(t / 26)
    mg = 42 * math.exp(t / 17)
    expected = (
        0.106 * math.exp((ph - 8) / 0.56) * b / (1 + (b / f) ** 2)
        + 0.52 * (1 + t / 43) * (s / 35) * math.exp(-d / 6) * mg / (1 + (mg / f) ** 2)
        + 0.00049 * math.exp(-t / 27 - d / 17) * f * f
    )
    water = m.SeaWater(temperature_c=t, salinity_ppt=s, depth_m=1500, ph=ph)
    assert m.absorption_db_per_km([4000], water=water)[0] == pytest.approx(expected, rel=1e-14)
    # Exact forward construction: emission CPA at 17 s arrives 500/1500 s later.
    result = m.predict_passage(
        [1000],
        [17 + 500 / 1500],
        speed_kn=12,
        length_m=100,
        closest_approach_m=300,
        source_depth_m=10,
        receiver_depth_m=410,
        wind_speed_m_s=4,
        closest_approach_emission_time_s=17,
    )
    assert result["emission_time_s"][0] == pytest.approx(17, abs=1e-12)
    assert result["distance_m"][0] == pytest.approx(500, abs=1e-12)
    assert result["travel_time_s"][0] == pytest.approx(1 / 3, abs=1e-12)
    slow = m.vessel_source_psd([1000], speed_kn=6, length_m=100)
    fast = m.vessel_source_psd([1000], speed_kn=12, length_m=100)
    np.testing.assert_allclose(fast / slow, 64, rtol=1e-13)
