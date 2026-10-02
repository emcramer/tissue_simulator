"""Tests for density-aware layouts in ``tissue_simulator.density``."""

import json

import numpy as np
from scipy import ndimage as ndi
import pytest
from scipy.stats import ks_2samp

from tissue_simulator import Cell, TissueSection
from tissue_simulator.density import (
    DensityModel,
    _assign_compartments,
    _band_index,
    _calibrate_patches,
    _largest_remainder,
    _noise_shape,
    _same_label_profile,
    _sample_labels,
    pooled_patch_prior,
)
from tissue_simulator.packing import SpherePacker

SIZE = 300.0
NESTS = np.array([[70.0, 70.0], [220.0, 80.0], [80.0, 225.0], [215.0, 220.0]])
NEST_RADIUS = 45.0


def _rsa_points(seed, max_attempts=300):
    cells = SpherePacker((SIZE, SIZE, 0.0), {'c': (3.5, 5.0)}, min_spacing=0.3,
                         seed=seed).pack(max_attempts=max_attempts)
    return np.array([c.center[:2] for c in cells]), np.array([c.radius for c in cells])


def _uniform_region(seed=0):
    xy, r = _rsa_points(seed, 150)
    types = np.random.default_rng(seed).choice(['A', 'B', 'C'], len(r))
    return xy[:, 0], xy[:, 1], r, types


def _nest_region(seed=0):
    """Full-density tumor nests in stroma thinned to 25%."""
    rng = np.random.default_rng(seed)
    xy, r = _rsa_points(seed)
    inside = (np.linalg.norm(xy[:, None] - NESTS[None], axis=2) < NEST_RADIUS).any(axis=1)
    keep = inside | (rng.random(len(r)) < 0.25)
    types = np.where(inside, np.where(rng.random(len(r)) < 0.9, 'Tumor', 'CD8'),
                     np.where(rng.random(len(r)) < 0.8, 'Stroma', 'CD8'))
    return xy[keep, 0], xy[keep, 1], r[keep], types[keep]


@pytest.fixture(scope="module")
def nest_model():
    x, y, r, t = _nest_region()
    model = DensityModel.fit(x, y, r, t, bounds=(0, 0, SIZE, SIZE), n_compartments=2, seed=0)
    return model, len(x)


def test_fit_validates_inputs():
    with pytest.raises(ValueError):
        DensityModel.fit([1.0] * 20, [1.0] * 19, [1.0] * 20, ['A'] * 20)
    with pytest.raises(ValueError):
        DensityModel.fit([1.0] * 5, [1.0] * 5, [1.0] * 5, ['A'] * 5)


def test_uniform_region_is_homogeneous():
    x, y, r, t = _uniform_region()
    model = DensityModel.fit(x, y, r, t, bounds=(0, 0, SIZE, SIZE), seed=0)
    assert model.homogeneous
    layout = model.sample_layout(rng=1)
    assert layout.mode == "uniform"
    assert np.ptp(layout.intensity) == pytest.approx(0.0)
    assert layout.n_target == len(x)


def test_nest_region_compartments(nest_model):
    model, n = nest_model
    assert not model.homogeneous
    tumor = model.cell_types.index('Tumor')
    dense = int(np.argmax(model.compartment_composition[:, tumor]))
    true_fraction = len(NESTS) * np.pi * NEST_RADIUS ** 2 / SIZE ** 2
    assert model.compartment_fractions[dense] == pytest.approx(true_fraction, abs=0.05)
    assert model.compartment_composition[dense, tumor] > 0.8
    assert model.compartment_composition[1 - dense, tumor] < 0.35
    assert model.region_intensity.sum() * model.grid_step ** 2 == pytest.approx(n, rel=1e-6)


def test_copy_layout_reproduces_region_maps(nest_model):
    model, n = nest_model
    layout = model.sample_layout(rng=0, layout="copy")
    assert layout.mode == "copy" and layout.n_target == n
    np.testing.assert_array_equal(layout.compartment, model.region_compartments)
    np.testing.assert_allclose(layout.intensity, model.region_intensity.sum(axis=0))
    with pytest.raises(ValueError):
        model.sample_layout(layout="copy", width=SIZE + 50)
    with pytest.raises(ValueError):
        model.sample_layout(layout="tile")


def test_resampled_layout_matches_region_statistics(nest_model):
    model, _ = nest_model
    layout = model.sample_layout(rng=3)
    assert layout.mode == "resample"
    counts = np.bincount(layout.compartment.ravel(), minlength=model.n_compartments)
    np.testing.assert_array_equal(
        counts, _largest_remainder(model.compartment_fractions, layout.compartment.size))
    assert layout.intensity.sum() * layout.grid_step ** 2 == pytest.approx(layout.n_target)
    region = model.region_intensity.sum(axis=0).ravel()
    assert ks_2samp(layout.intensity.ravel(), region).statistic < 0.1
    assert not np.array_equal(layout.compartment, model.region_compartments)


def test_resampled_layout_keeps_boundary_margins():
    rng = np.random.default_rng(2)
    xy, r = _rsa_points(2)
    d = (np.linalg.norm(xy[:, None] - NESTS[None], axis=2) - NEST_RADIUS).min(axis=1)
    ring = (d >= 0) & (d < 12)
    keep = (d < 0) | (ring & (rng.random(len(r)) < 0.8)) | (~ring & (rng.random(len(r)) < 0.25))
    types = np.where(d < 0, 'Tumor', np.where(ring, 'CD8', 'Stroma'))
    model = DensityModel.fit(xy[keep, 0], xy[keep, 1], r[keep], types[keep],
                             bounds=(0, 0, SIZE, SIZE), n_null=0, n_compartments=2, seed=0)
    cd8 = model.cell_types.index('CD8')
    sparse = int(np.argmin(model.compartment_composition[:, model.cell_types.index('Tumor')]))
    assert model.band_composition[sparse, 0, cd8] > 1.5 * model.band_composition[sparse, -1, cd8]

    layout = model.sample_layout(rng=4)
    band = _band_index(layout.compartment, layout.grid_step)
    edge = (layout.compartment == sparse) & (band == 0)
    interior = (layout.compartment == sparse) & (band == band.max())
    assert layout.composition[cd8][edge].mean() > 1.5 * layout.composition[cd8][interior].mean()


def test_layout_sampling_is_seeded(nest_model):
    model, _ = nest_model
    a, b, c = (model.sample_layout(rng=s) for s in (5, 5, 6))
    np.testing.assert_array_equal(a.intensity, b.intensity)
    assert not np.array_equal(a.intensity, c.intensity)


def test_resampled_layout_supports_other_window(nest_model):
    model, _ = nest_model
    layout = model.sample_layout(rng=0, width=500.0, height=250.0)
    assert layout.intensity.shape == (50, 100)
    assert layout.n_target == round(model.density * 500 * 250)


def test_json_round_trip(nest_model):
    model, _ = nest_model
    restored = DensityModel.from_dict(json.loads(json.dumps(model.to_dict())))
    np.testing.assert_array_equal(restored.sample_layout(rng=2).intensity,
                                  model.sample_layout(rng=2).intensity)
    assert restored.cell_types == model.cell_types and restored.flags == model.flags


def test_mask_excludes_tissue_free_area():
    x, y, r, t = _nest_region()
    mask = np.ones((60, 60), dtype=bool)
    mask[:, 45:] = False
    model = DensityModel.fit(x, y, r, t, bounds=(0, 0, SIZE, SIZE), mask=mask,
                             n_null=0, n_compartments=2, seed=0)
    layout = model.sample_layout(layout="copy")
    assert np.all(layout.intensity[:, 45:] == 0)
    assert np.all(layout.compartment[:, 45:] == -1)
    assert model.n_cells == int(np.sum(x < 225.0))


def test_gradient_region_is_flagged_as_trend():
    rng = np.random.default_rng(4)
    xy, r = _rsa_points(4)
    keep = rng.random(len(r)) < 0.1 + 0.9 * xy[:, 0] / SIZE
    types = rng.choice(['A', 'B'], keep.sum())
    with pytest.warns(UserWarning, match="stationary"):
        model = DensityModel.fit(xy[keep, 0], xy[keep, 1], r[keep], types,
                                 bounds=(0, 0, SIZE, SIZE), seed=0)
    assert "trend" in model.flags and not model.stationary


def test_assign_compartments_hits_exact_counts():
    rng = np.random.default_rng(0)
    points = rng.normal(size=(1000, 2))
    centers = np.array([[-1.0, 0.0], [1.0, 0.5], [0.0, -1.5]])
    labels = _assign_compartments(points, centers, np.array([0.6, 0.3, 0.1]))
    np.testing.assert_array_equal(np.bincount(labels, minlength=3), [600, 300, 100])


def test_patch_calibration_recovers_patch_structure():
    shape, step, bandwidth, n_lags, pad = (80, 80), 5.0, 15.0, 27, 60
    centers, fractions = np.array([[-0.8, 0.0], [0.8, 0.0]]), np.array([0.7, 0.3])

    def labels_for(length, nu, seed):
        rng = np.random.default_rng(seed)
        noise = [rng.standard_normal(_noise_shape(shape, pad)) for _ in range(2)]
        return _sample_labels(shape, step, length, nu, bandwidth, 0.0, centers, fractions,
                              noise[0], noise[1], pad)[0].reshape(shape)

    truth = labels_for(40.0, 1.0, 0)
    length, nu, _, at_upper = _calibrate_patches(truth, 2, centers, fractions, 0.0, step,
                                                 bandwidth, step, 400.0 / 3,
                                                 np.random.default_rng(1))
    assert not at_upper

    def mean_profile(length_, nu_):
        return np.mean([_same_label_profile(labels_for(length_, nu_, s), 2, n_lags)
                        for s in (10, 11, 12)], axis=0)

    np.testing.assert_allclose(mean_profile(length, nu)[1:], mean_profile(40.0, 1.0)[1:],
                               atol=0.06)


def test_pooled_patch_prior_uses_unflagged_fits(nest_model):
    model, _ = nest_model
    prior = pooled_patch_prior([model, model])
    if prior:
        assert prior["length"] == pytest.approx(model.patch_length)
        assert prior["smoothness"] == model.patch_smoothness
    fixed = DensityModel.fit(*_nest_region(), bounds=(0, 0, SIZE, SIZE), n_null=0,
                             n_compartments=2, patch_prior={"length": 30.0, "smoothness": 2.0})
    assert (fixed.patch_length, fixed.patch_smoothness) == (30.0, 2.0)


def test_from_tissue_matches_fit():
    x, y, r, t = _nest_region()
    tissue = TissueSection(SIZE, SIZE, 1.0, {'c': (3.5, 5.0)})
    tissue.cells = [Cell((xi, yi, 0.5), ri, ti) for xi, yi, ri, ti in zip(x, y, r, t)]
    a = DensityModel.from_tissue(tissue, n_null=0, seed=0)
    b = DensityModel.fit(x, y, r, t, bounds=(0, 0, SIZE, SIZE), n_null=0, seed=0)
    np.testing.assert_array_equal(a.region_compartments, b.region_compartments)


# ---------------------------------------------------------------------------
# Adaptive strategy: estimation, organization, metadata, serialization
# ---------------------------------------------------------------------------

def _dict_equal(a, b):
    return json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)


def _gradient_region(seed=1):
    """Constant-density RSA cells whose P(type A) rises 0.1 -> 0.9 along x."""
    xy, r = _rsa_points(seed)
    rng = np.random.default_rng(seed)
    types = np.where(rng.random(len(r)) < 0.1 + 0.8 * xy[:, 0] / SIZE, 'A', 'B')
    return xy[:, 0], xy[:, 1], r, types


@pytest.fixture(scope="module")
def gradient_model():
    x, y, r, t = _gradient_region()
    return DensityModel.fit(x, y, r, t, bounds=(0, 0, SIZE, SIZE), n_null=0, seed=0,
                            strategy="adaptive")


def test_legacy_fit_defaults_are_unchanged():
    x, y, r, t = _nest_region()
    kw = dict(bounds=(0, 0, SIZE, SIZE), n_compartments=2, seed=0)
    default = DensityModel.fit(x, y, r, t, **kw)
    explicit = DensityModel.fit(x, y, r, t, strategy="legacy", bandwidth_range=(15.0, 80.0), **kw)
    strip = lambda m: {k: v for k, v in m.to_dict().items() if k != "estimation"}
    assert _dict_equal(strip(default), strip(explicit))
    assert default.strategy == "legacy"
    assert default.band_edges == (10.0, 20.0, 40.0)
    assert default.organization == {"model": "none"}
    assert default.estimation["bandwidth_range_source"] == "legacy_default"
    np.testing.assert_array_equal(default.type_bandwidths, default.bandwidth)
    with pytest.raises(ValueError):
        DensityModel.fit(x, y, r, t, strategy="fancy", **kw)


def test_v1_dict_loads_with_legacy_defaults(nest_model):
    model, _ = nest_model
    data = model.to_dict()
    assert data["format_version"] == 2
    for key in ("strategy", "type_bandwidths", "band_edges", "size_log_mu", "size_log_sigma",
                "estimation", "organization"):
        del data[key]
    data["format_version"] = 1
    old = DensityModel.from_dict(json.loads(json.dumps(data)))
    assert old.strategy == "legacy" and old.band_edges == (10.0, 20.0, 40.0)
    assert old.organization == {"model": "none"} and old.estimation == {}
    assert old.type_bandwidths.size == 0
    layout = old.sample_layout(rng=3)
    ref = model.sample_layout(rng=3)
    np.testing.assert_array_equal(layout.intensity, ref.intensity)


def test_adaptive_model_round_trips(gradient_model):
    clone = DensityModel.from_dict(json.loads(json.dumps(gradient_model.to_dict())))
    assert clone.strategy == "adaptive"
    assert clone.band_edges == gradient_model.band_edges
    assert _dict_equal(clone.organization, gradient_model.organization)
    np.testing.assert_array_equal(clone.size_log_mu, gradient_model.size_log_mu)
    np.testing.assert_array_equal(clone.type_bandwidths, gradient_model.type_bandwidths)


def test_homogeneous_control_has_no_organization():
    x, y, r, t = _uniform_region()
    model = DensityModel.fit(x, y, r, t, bounds=(0, 0, SIZE, SIZE), seed=0, strategy="adaptive")
    assert model.homogeneous and model.organization["model"] == "none"
    # Even with the homogeneity test skipped, no trend reaches the BIC threshold.
    free = DensityModel.fit(x, y, r, t, bounds=(0, 0, SIZE, SIZE), seed=0, n_null=0,
                            strategy="adaptive")
    assert free.organization["model"] == "none"
    assert free.organization["fallback"] == "no_candidate_reached_bic_threshold"
    assert not any(f.startswith("organization:") for f in free.flags)


def test_composition_gradient_selects_planar_composition(gradient_model):
    org = gradient_model.organization
    assert org["model"] in ("planar_composition", "planar_both")
    assert abs(org["direction"][0]) > 0.9
    assert org["delta_bic"] >= 10.0 and org["bic_threshold"] == 10.0
    assert org["gradient_strength"]["density_log_span"] < np.log(1.5)
    assert org["gradient_strength"]["composition_span"] > 0.5
    assert f"organization:{org['model']}" in gradient_model.flags
    # Profiles run monotonically along the (signed) direction for type A.
    comp = np.array(org["composition_profile"])
    knots = np.array(org["s_knots"])
    a = gradient_model.cell_types.index('A')
    slope = np.polyfit(knots, comp[:, a], 1)[0]
    assert slope * org["direction"][0] > 0
    np.testing.assert_allclose(comp.sum(axis=1), 1.0, atol=1e-9)
    # The trend was detrended away: the model is not flagged as a plain trend.
    assert "trend" not in gradient_model.flags
    assert {c["model"] for c in org["candidates"]} >= {"planar_both", "radial_both"}


def test_organization_selection_is_deterministic(gradient_model):
    x, y, r, t = _gradient_region()
    again = DensityModel.fit(x, y, r, t, bounds=(0, 0, SIZE, SIZE), n_null=0, seed=5,
                             strategy="adaptive")
    assert _dict_equal(again.organization, gradient_model.organization)


def test_planar_density_gradient_is_detected():
    xy, r = _rsa_points(3)
    rng = np.random.default_rng(3)
    keep = rng.random(len(r)) < 0.15 + 0.85 * xy[:, 1] / SIZE
    types = rng.choice(['A', 'B'], int(keep.sum()))
    model = DensityModel.fit(xy[keep, 0], xy[keep, 1], r[keep], types,
                             bounds=(0, 0, SIZE, SIZE), n_null=0, seed=0, strategy="adaptive")
    org = model.organization
    assert org["model"] == "planar_density"
    assert abs(org["direction"][1]) > 0.9
    assert org["gradient_strength"]["density_log_span"] > np.log(2.0)
    assert org["transition_widths"] == [None, None]


def test_rare_type_gets_pooled_bandwidth_and_valid_composition():
    xy, r = _rsa_points(2)
    rng = np.random.default_rng(2)
    types = np.where(rng.random(len(r)) < 0.5, 'A', 'B').astype(object)
    types[:6] = 'R'
    model = DensityModel.fit(xy[:, 0], xy[:, 1], r, types.astype(str),
                             bounds=(0, 0, SIZE, SIZE), n_null=0, seed=0, strategy="adaptive")
    rare = model.cell_types.index('R')
    assert model.estimation["type_bandwidth_evidence"]['R'] == "pooled_fallback"
    assert model.type_bandwidths[rare] == pytest.approx(model.bandwidth)
    assert np.isfinite(model.type_bandwidths).all() and (model.type_bandwidths > 0).all()
    np.testing.assert_allclose(model.compartment_composition.sum(axis=1), 1.0, atol=1e-9)
    np.testing.assert_allclose(model.band_composition.sum(axis=-1), 1.0, atol=1e-9)
    assert model.estimation["size_fallbacks"] == ['R']
    assert np.isfinite(model.size_log_mu).all() and (model.size_log_sigma > 0).all()


def test_bandwidth_at_bound_is_recorded():
    x, y, r, t = _nest_region()
    model = DensityModel.fit(x, y, r, t, bounds=(0, 0, 100, 100), bandwidth_range=(60, 80),
                             n_null=0, seed=0)
    assert model.estimation["bandwidth_at_bound"] is True
    assert model.estimation["bandwidth_range_source"] == "user"
    assert model.estimation["bandwidth_range"] == [60.0, 80.0]
    assert 60.0 <= model.bandwidth <= 80.0


def test_auto_bandwidth_range_and_learned_band_edges():
    x, y, r, t = _nest_region()
    kw = dict(bounds=(0, 0, SIZE, SIZE), n_compartments=2, seed=0)
    model = DensityModel.fit(x, y, r, t, strategy="adaptive", **kw)
    lo, hi = model.estimation["bandwidth_range"]
    assert model.estimation["bandwidth_range_source"] == "auto"
    assert lo >= 2 * np.median(r) and hi >= 3 * lo
    edges = model.band_edges
    assert list(edges) == sorted(set(edges)) and edges[0] >= model.grid_step
    assert model.band_density_quantiles.shape[1] == len(edges) + 1
    assert model.estimation["learned_band_edges"] is True
    layout = model.sample_layout(rng=1)
    assert np.isfinite(layout.intensity).all() and layout.mode in ("resample", "uniform")
    # Explicit "auto" works for the legacy strategy too.
    legacy = DensityModel.fit(x, y, r, t, bandwidth_range="auto", **kw)
    assert legacy.band_edges == (10.0, 20.0, 40.0)
    assert legacy.estimation["bandwidth_range_source"] == "auto"


def test_size_model_recovers_type_specific_radii():
    x, y, r, t = _nest_region()
    r = np.where(t == 'Tumor', 6.0, 3.0) * np.exp(0.05 * np.random.default_rng(0).standard_normal(len(r)))
    model = DensityModel.fit(x, y, r, t, bounds=(0, 0, SIZE, SIZE), n_compartments=2, seed=0)
    assert model.size_log_mu.shape == (len(model.cell_types), 5)
    tumor, stroma = model.cell_types.index('Tumor'), model.cell_types.index('Stroma')
    assert np.exp(model.size_log_mu[tumor]).mean() == pytest.approx(6.0, rel=0.1)
    assert np.exp(model.size_log_mu[stroma]).mean() == pytest.approx(3.0, rel=0.1)
    assert (model.size_log_sigma < 0.3).all()


def test_band_index_accepts_learned_edges():
    labels = np.zeros((20, 20), dtype=int)
    labels[:, 10:] = 1
    legacy = _band_index(labels, 5.0)
    custom = _band_index(labels, 5.0, (5.0,))
    assert legacy.max() == 3 and custom.max() == 1
    np.testing.assert_array_equal(_band_index(labels, 5.0, (10.0, 20.0, 40.0)), legacy)


# -- organized (directional) layouts -----------------------------------------

def _layered_region(seed=0, w=240.0, h=120.0):
    """Bands A | B | C along x (thirds), ~600 RSA cells."""
    cells = SpherePacker((h, w, 0.0), {'c': (2.3, 3.0)}, min_spacing=0.3,
                         seed=seed).pack(max_attempts=400)
    xy = np.array([c.center[:2] for c in cells])
    r = np.array([c.radius for c in cells])
    rng = np.random.default_rng(seed)
    third = np.minimum((xy[:, 0] / (w / 3)).astype(int), 2)
    pure = rng.random(len(r)) < 0.9
    types = np.where(pure, np.array(['A', 'B', 'C'])[third], rng.choice(['A', 'B', 'C'], len(r)))
    return xy[:, 0], xy[:, 1], r, types


@pytest.fixture(scope="module")
def layered_model():
    x, y, r, t = _layered_region(2)
    return DensityModel.fit(x, y, r, t, bounds=(0, 0, 240, 120), strategy="adaptive", seed=0,
                            units="none")  # tests the planar trend path, not units


def test_organized_layouts(layered_model):
    model = layered_model
    assert model.organization["model"] != "none"
    copy = model.sample_layout(layout="copy")
    layouts = [model.sample_layout(rng=s) for s in range(4)]
    thetas = []
    for lay in layouts:
        org = lay.organization
        assert org["accepted"] is True and org["fallback"] is None
        assert "organization_unsatisfied" not in lay.flags
        assert lay.quota_scale is not None and lay.strategy == "adaptive"
        thetas.append(org["theta"])
        d = np.array(org["direction"])
        ny, nx = lay.intensity.shape
        yy, xx = np.meshgrid((np.arange(ny) + .5) * lay.grid_step,
                             (np.arange(nx) + .5) * lay.grid_step, indexing="ij")
        s = xx * d[0] + yy * d[1]
        cent = []
        for t in range(3):
            w = lay.intensity * lay.composition[t]
            cent.append((w * s).sum() / w.sum())
        assert cent[0] < cent[1] < cent[2]
        w = (lay.composition * lay.intensity).sum(axis=(1, 2))
        np.testing.assert_allclose(w / w.sum(), model.proportions, rtol=0.15)
        assert np.all(np.isfinite(lay.composition))
        np.testing.assert_allclose(lay.composition.sum(axis=0), 1.0)
        assert np.corrcoef(lay.intensity.ravel(), copy.intensity.ravel())[0, 1] < 0.95
    for i in range(4):
        for j in range(i + 1, 4):
            diff = abs(thetas[i] - thetas[j]) % (2 * np.pi)
            assert np.degrees(min(diff, 2 * np.pi - diff)) > 10.0


def test_adaptive_without_organization_matches_legacy_resample():
    x, y, r, t = _nest_region()
    kw = dict(bounds=(0, 0, SIZE, SIZE), n_compartments=2, seed=0)
    legacy = DensityModel.fit(x, y, r, t, **kw)
    adaptive = DensityModel.fit(x, y, r, t, strategy="adaptive", organization=False,
                                bandwidth_range=(15.0, 80.0), per_type_bandwidth=False, **kw)
    a = legacy.sample_layout(rng=5)
    assert a.strategy == "legacy" and a.quota_scale is None and a.organization == {}
    # Legacy goes through the untouched resampler; adaptive with no trend model
    # must use exactly the same path (no extra draws).
    for model in (legacy, adaptive):
        lay = model.sample_layout(rng=5)
        ref = model._resampled_layout(np.random.default_rng(5), model.width, model.height)
        np.testing.assert_array_equal(lay.intensity, ref.intensity)
        np.testing.assert_array_equal(lay.composition, ref.composition)
        assert lay.organization == {}


def test_organized_layout_unsatisfied_on_narrow_window(layered_model):
    with pytest.warns(UserWarning, match="best-coverage"):
        lay = layered_model.sample_layout(rng=0, width=60.0, height=60.0, max_proposals=1)
    assert lay.organization["accepted"] is False
    assert lay.organization["fallback"] == "best_of_proposals"
    assert "organization_unsatisfied" in lay.flags
    assert np.all(np.isfinite(lay.intensity))


# -- radial organized layouts --------------------------------------------------

def _concentric_region(seed=0, size=200.0):
    """Core A (r<30), ring B (30-60), periphery C around the window centre."""
    cells = SpherePacker((size, size, 0.0), {'c': (2.3, 3.0)}, min_spacing=0.3,
                         seed=seed).pack(max_attempts=400)
    xy = np.array([c.center[:2] for c in cells])
    r = np.array([c.radius for c in cells])
    rng = np.random.default_rng(seed)
    d = np.linalg.norm(xy - size / 2, axis=1)
    ring = np.where(d < 30, 0, np.where(d < 60, 1, 2))
    pure = rng.random(len(r)) < 0.9
    types = np.where(pure, np.array(['A', 'B', 'C'])[ring], rng.choice(['A', 'B', 'C'], len(r)))
    return xy[:, 0], xy[:, 1], r, types


@pytest.fixture(scope="module")
def concentric_model():
    x, y, r, t = _concentric_region(2)
    return DensityModel.fit(x, y, r, t, bounds=(0, 0, 200, 200), strategy="adaptive", seed=0,
                            units="none")  # tests the radial path; units take precedence otherwise


def test_radial_organized_layouts(concentric_model):
    model = concentric_model
    assert model.organization["model"].startswith("radial_")
    centers = []
    for seed in range(3):
        lay = model.sample_layout(rng=seed)
        org = lay.organization
        assert org["accepted"] is True and org["fallback"] is None
        assert org["geometry"] == "radial"
        c = np.array(org["center"])
        centers.append(c)
        ny, nx = lay.intensity.shape
        yy, xx = np.meshgrid((np.arange(ny) + .5) * lay.grid_step,
                             (np.arange(nx) + .5) * lay.grid_step, indexing="ij")
        dist = np.hypot(xx - c[0], yy - c[1])
        mean_d = []
        for t in range(3):
            w = lay.intensity * lay.composition[t]
            mean_d.append((w * dist).sum() / w.sum())
        assert mean_d[0] < mean_d[1] < mean_d[2]
        w = (lay.composition * lay.intensity).sum(axis=(1, 2))
        np.testing.assert_allclose(w / w.sum(), model.proportions, rtol=0.15, atol=0.02)  # acceptance rule: max(0.02, 15 %)
    for i in range(3):
        for j in range(i + 1, 3):
            assert np.linalg.norm(centers[i] - centers[j]) > 15.0


def test_ambiguous_control_not_radial():
    rng = np.random.default_rng(3)
    xy, r = _rsa_points(1, 150)
    d = np.linalg.norm(xy - np.array([150.0, 150.0]), axis=1)
    keep = rng.random(len(r)) < (0.65 + 0.35 * np.exp(-(d / 60.0) ** 2))
    types = rng.choice(['A', 'B', 'C'], keep.sum())
    model = DensityModel.fit(xy[keep, 0], xy[keep, 1], r[keep], types,
                             bounds=(0, 0, SIZE, SIZE), strategy="adaptive", seed=0)
    org = model.organization
    assert not org["model"].startswith("radial_")
    assert org.get("weak_candidate") is None or isinstance(org.get("weak_candidate"), (dict, str))
    model.sample_layout(rng=0)


# -- trend null ------------------------------------------------------------

def _blob_region(nests, seed=0, size=300.0, radius=60.0):
    """Nest-type enrichment in thinned stroma, like ``_nest_region``, on a custom layout."""
    cells = SpherePacker((size, size, 0.0), {'c': (3.5, 5.0)}, min_spacing=0.3,
                         seed=seed).pack(max_attempts=300)
    xy = np.array([c.center[:2] for c in cells])
    r = np.array([c.radius for c in cells])
    rng = np.random.default_rng(seed)
    inside = (np.linalg.norm(xy[:, None] - np.array(nests)[None], axis=2) < radius).any(axis=1)
    keep = inside | (rng.random(len(r)) < 0.25)
    types = np.where(inside, np.where(rng.random(len(r)) < 0.9, 'Tumor', 'CD8'),
                     np.where(rng.random(len(r)) < 0.8, 'Stroma', 'CD8'))
    return xy[keep, 0], xy[keep, 1], r[keep], types[keep]


def _fit_blob(nests):
    x, y, r, t = _blob_region(nests)
    return DensityModel.fit(x, y, r, t, bounds=(0, 0, 300, 300), strategy="adaptive", seed=0)


def test_single_blob_is_accepted_as_radial():
    """One central blob is indistinguishable from a concentric trend in a single
    window, so it is accepted as radial (the capped-patch null cannot make it)."""
    org = _fit_blob([(150.0, 150.0)]).organization
    assert org["model"].startswith("radial_") and org["null"]["p_value"] <= 0.05


def test_two_blobs_are_not_radial():
    org = _fit_blob([(80.0, 80.0), (220.0, 220.0)]).organization
    assert org["model"] == "none" and org["fallback"] == "stationary_null"
    assert org["selected_model"] and org["null"]["p_value"] > 0.05


def test_null_skipped_when_n_null_is_zero():
    x, y, r, t = _gradient_region()
    m = DensityModel.fit(x, y, r, t, bounds=(0, 0, SIZE, SIZE), n_null=0, seed=0,
                         strategy="adaptive")
    assert m.organization["null"] is None and m.organization["model"].startswith("planar_")


def test_gradient_survives_trend_null():
    x, y, r, t = _gradient_region()
    org = DensityModel.fit(x, y, r, t, bounds=(0, 0, SIZE, SIZE), seed=0,
                           strategy="adaptive").organization
    assert org["model"].startswith("planar_")
    n = org["null"]
    assert {"n_null", "p_value", "observed_delta_bic", "null_delta_bic_quantiles",
            "ambiguous", "seconds"} <= set(n)
    assert n["p_value"] <= 0.05 and n["n_null"] == 19


# -- voids and lumens ----------------------------------------------------------

def _ring_lumen_region(seed=2, size=200.0):
    """No cells within 25 um of the centre (lumen); B 25-55; C outside."""
    x, y, r, _ = _concentric_region(seed, size)
    d = np.hypot(x - size / 2, y - size / 2)
    keep = d >= 25
    rng = np.random.default_rng(seed)
    ring = np.where(d < 55, 1, 2)
    pure = rng.random(len(r)) < 0.9
    types = np.where(pure, np.array(['A', 'B', 'C'])[ring], rng.choice(['A', 'B', 'C'], len(r)))
    return x[keep], y[keep], r[keep], types[keep]


def _two_lumen_region(seed=0, size=200.0):
    x, y, r, _ = _concentric_region(seed, size)
    centers = np.array([[60.0, 60.0], [140.0, 140.0]])
    keep = (np.linalg.norm(np.column_stack([x, y])[:, None] - centers[None], axis=2) >= 20).all(axis=1)
    types = np.random.default_rng(seed).choice(['A', 'B', 'C'], len(r))
    return x[keep], y[keep], r[keep], types[keep]


def _components(zero):
    from scipy import ndimage
    return ndimage.label(zero)


@pytest.fixture(scope="module")
def ring_lumen_model():
    x, y, r, t = _ring_lumen_region()
    return DensityModel.fit(x, y, r, t, bounds=(0, 0, 200, 200), strategy="adaptive",
                            seed=0, n_null=0, units="none")  # tests the voids path


@pytest.fixture(scope="module")
def two_lumen_model():
    x, y, r, t = _two_lumen_region()
    return DensityModel.fit(x, y, r, t, bounds=(0, 0, 200, 200), strategy="adaptive",
                            seed=0, n_null=0)


def test_ring_lumen_voids_are_inferred_and_placed(ring_lumen_model):
    import warnings
    model = ring_lumen_model
    assert model.voids["n_holes"] == 1
    assert model.estimation["voids"]["n_holes"] == 1
    diam = model.voids["equivalent_diameters"][0]
    assert abs(diam - 50.0) <= 0.2 * 50.0
    source_area = model.voids["areas"][0]
    off_center = 0
    for seed in range(3):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            lay = model.sample_layout(rng=seed)
        assert lay.voids["n_placed"] == 1 and "voids_unsatisfied" not in lay.flags
        labels, n = _components(lay.intensity == 0)
        assert n == 1
        area = float((labels > 0).sum()) * lay.grid_step ** 2
        assert abs(area - source_area) <= 0.3 * source_area
        c = np.array(lay.voids["centers"][0])
        off_center += int(np.linalg.norm(c - 100.0) > 5.0)
        # composition just outside the hole edge is dominated by B
        ny, nx = lay.intensity.shape
        yy, xx = np.meshgrid((np.arange(ny) + .5) * lay.grid_step,
                             (np.arange(nx) + .5) * lay.grid_step, indexing="ij")
        dist = np.hypot(xx - c[0], yy - c[1]) - lay.voids["radii"][0]
        ring = (dist > 0) & (dist <= 15.0) & (lay.intensity > 0)
        b = lay.composition[1][ring]
        assert (b * lay.intensity[ring]).sum() / lay.intensity[ring].sum() > 0.6
        if lay.organization and lay.organization.get("geometry") == "radial":
            assert np.linalg.norm(c - np.array(lay.organization["center"])) <= 10.0
    assert off_center >= 2


def test_two_lumens_are_placed_apart(two_lumen_model):
    model = two_lumen_model
    assert model.voids["n_holes"] == 2
    for seed in range(3):
        lay = model.sample_layout(rng=seed)
        assert lay.voids["n_placed"] == 2
        c = np.array(lay.voids["centers"])
        assert np.linalg.norm(c[0] - c[1]) >= 60.0


def test_uniform_control_has_no_voids_and_matches_none():
    x, y, r, t = _uniform_region(0)
    kw = dict(strategy="adaptive", seed=0, n_null=0, bounds=(0, 0, SIZE, SIZE))
    auto = DensityModel.fit(x, y, r, t, **kw)
    none = DensityModel.fit(x, y, r, t, voids="none", **kw)
    assert auto.voids["n_holes"] == 0 and none.voids == {}
    a, b = auto.sample_layout(rng=1), none.sample_layout(rng=1)
    assert np.array_equal(a.intensity, b.intensity)
    assert np.array_equal(a.composition, b.composition)
    assert a.voids == {}


def test_voids_serialization_and_defaults(ring_lumen_model):
    clone = DensityModel.from_dict(json.loads(json.dumps(ring_lumen_model.to_dict())))
    assert clone.voids["n_holes"] == 1
    data = ring_lumen_model.to_dict()
    del data["voids"]
    assert DensityModel.from_dict(data).voids == {}
    with pytest.raises(ValueError):
        DensityModel.fit(*_uniform_region(0), voids="bogus")


def test_legacy_default_draws_no_void_state():
    x, y, r, t = _nest_region(0)
    m = DensityModel.fit(x, y, r, t, n_null=0, seed=0)
    assert m.voids == {} and "voids" not in m.estimation
    assert m.sample_layout(rng=0).voids == {}


# ---------------------------------------------------------------------------
# Adaptive composition-only nests and homogeneous holes (WP-D)
# ---------------------------------------------------------------------------

_COMP_NESTS = np.array([[60.0, 70.0], [200.0, 60.0], [90.0, 220.0], [230.0, 210.0]])


def _composition_nests_region(seed=0):
    """Uniform density; four r=38 discs 90 % B inside, stroma 80/10/10 A/B/C."""
    rng = np.random.default_rng(seed)
    xy, r = _rsa_points(seed)
    inside = (np.linalg.norm(xy[:, None] - _COMP_NESTS[None], axis=2) < 38.0).any(axis=1)
    u = rng.random(len(r))
    out = np.where(u < 0.8, 'A', np.where(u < 0.9, 'B', 'C'))
    types = np.where(inside, np.where(rng.random(len(r)) < 0.9, 'B', 'A'), out)
    return xy[:, 0], xy[:, 1], r, types


@pytest.fixture(scope="module")
def nests_composition_only():
    x, y, r, t = _composition_nests_region()
    kw = dict(bounds=(0, 0, SIZE, SIZE), seed=0)
    return (DensityModel.fit(x, y, r, t, strategy="adaptive", **kw),
            DensityModel.fit(x, y, r, t, **kw))


def test_adaptive_detects_composition_only_nests(nests_composition_only):
    model, legacy = nests_composition_only
    het = model.heterogeneity
    assert "composition_bandwidths" in het
    assert not model.homogeneous
    assert het["p_composition"] <= 0.05
    assert len(set(model.region_compartments[model.mask].tolist())) >= 2
    lay = model.sample_layout(rng=1)
    b = list(lay.cell_types).index('B')
    comp_b = lay.composition[b][lay.compartment >= 0]
    assert (comp_b > 0.5).mean() >= 0.05
    assert (comp_b < 0.3).mean() >= 0.5
    # Legacy tests at the pooled bandwidth; its verdict may stay homogeneous.
    assert "composition_bandwidths" not in legacy.heterogeneity


def _glomeruli_region(seed=0):
    xy, r = _rsa_points(seed)
    centers = np.array([[70.0, 80.0], [200.0, 90.0], [140.0, 220.0]])
    keep = (np.linalg.norm(xy[:, None] - centers[None], axis=2) >= 18.0).all(axis=1)
    types = np.random.default_rng(seed).choice(['A', 'B'], len(r))
    return xy[keep, 0], xy[keep, 1], r[keep], types[keep]


def test_homogeneous_adaptive_model_keeps_holes():
    x, y, r, t = _glomeruli_region()
    model = DensityModel.fit(x, y, r, t, bounds=(0, 0, SIZE, SIZE), seed=0,
                             strategy="adaptive")
    assert model.voids["n_holes"] >= 1
    lay = model.sample_layout(rng=3)
    assert lay.mode == "uniform" or not model.homogeneous
    if model.homogeneous:
        assert lay.mode == "uniform"
        assert lay.voids["n_placed"] == model.voids["n_holes"] >= 1
        assert (lay.intensity == 0).any() and (lay.compartment == -1).any()
        assert lay.n_target == int(round(model.density * SIZE * SIZE
                                         * float((lay.compartment >= 0).mean())))


def test_edge_touching_voids_do_not_bias_cell_counts():
    x, y, r, t = _nest_region()
    kw = dict(bounds=(0, 0, SIZE, SIZE), seed=0, strategy="adaptive", n_null=0)
    auto = DensityModel.fit(x, y, r, t, voids="auto", **kw)
    none = DensityModel.fit(x, y, r, t, voids="none", **kw)
    assert np.isclose(auto.density, none.density)
    assert auto.sample_layout(rng=1).n_target == none.sample_layout(rng=1).n_target


def test_degenerate_trend_null_accepts_on_bic_and_flags(monkeypatch, gradient_model):
    from tissue_simulator import _organization
    assert gradient_model.organization["model"] != "none"
    x, y, r, t = _gradient_region()
    calls = []

    def fake(*a, **k):
        calls.append(1)
        return -np.inf
    monkeypatch.setattr(_organization, "trend_statistic", fake)
    model = DensityModel.fit(x, y, r, t, bounds=(0, 0, SIZE, SIZE), seed=0,
                             strategy="adaptive")
    assert calls
    assert model.organization["model"] != "none"
    assert model.organization["null"]["degenerate"] is True
    assert "trend_null_degenerate" in model.flags


# -- germ-grain unit layouts ------------------------------------------------

from .test_units import (  # noqa: E402  (module-scoped fixtures reused here)
    FOLLICLES, GLOMERULI, BLOBS, _fit as _unit_fit, _points as _unit_points,
    follicles, glomeruli, artery, blobs,
)


def _grain_s(layout, k):
    """Normalized distance map s = r / R_outer of placed grain k."""
    ny, nx = layout.intensity.shape
    gy, gx = np.mgrid[0:ny, 0:nx]
    c = layout.units["centers"][k]
    r = np.hypot((gx + 0.5) * layout.grid_step - c[0], (gy + 0.5) * layout.grid_step - c[1])
    return r / layout.units["outer_radii"][k]


def _sample_units(model, seeds=(0, 1, 2)):
    return [model.sample_layout(rng=s) for s in seeds]


def _basic_checks(lay):
    assert np.isfinite(lay.composition).all() and np.isfinite(lay.intensity).all()
    np.testing.assert_allclose(lay.composition.sum(axis=0), 1.0, atol=1e-9)
    assert (lay.intensity >= 0).all()


def test_follicle_layouts_place_new_c_core_b_ring_grains(follicles):
    assert follicles.units["n_units"] == 3 and "profiles" in follicles.units
    moved = 0
    for lay in _sample_units(follicles):
        _basic_checks(lay)
        assert lay.units["n_placed"] == 3 and "units_unsatisfied" not in lay.flags
        centers = np.array(lay.units["centers"])
        src = np.array(FOLLICLES, float)
        far = np.linalg.norm(centers[:, None] - src[None], axis=2).min(axis=1).max()
        moved += far > 10.0
        ci, bi = follicles.cell_types.index('C'), follicles.cell_types.index('B')
        for k in range(3):
            s = _grain_s(lay, k)
            core, ring = s < 0.3, (s > 0.5) & (s < 1.0)
            assert lay.composition[ci][core].mean() > 0.5
            assert lay.composition[bi][ring].mean() > 0.5
            assert lay.intensity[core].min() > 0  # a core is not a lumen
    assert moved >= 2


def test_glomerulus_lumens_are_unit_holes_not_duplicated(glomeruli):
    for lay in _sample_units(glomeruli):
        _basic_checks(lay)
        assert lay.units["n_placed"] == 3
        zero = lay.intensity == 0
        labels, n = ndi.label(zero)
        assert n == 3
        found = (np.array(ndi.center_of_mass(zero, labels, range(1, n + 1)))[:, ::-1] + 0.5) * lay.grid_step
        centers = np.array(lay.units["centers"])
        assert np.linalg.norm(found[:, None] - centers[None], axis=2).min(axis=1).max() < 5.0
        assert not lay.voids or lay.voids.get("n_placed", 0) == 0


def test_artery_unit_takes_precedence_over_radial_path(artery):
    for lay in _sample_units(artery):
        _basic_checks(lay)
        assert lay.units["n_placed"] == 1
        _, n = ndi.label(lay.intensity == 0)
        assert n == 1
        if artery.organization.get("model", "none") != "none":
            assert lay.organization == {"model": "none", "fallback": "units"}


def test_blob_layouts_are_b_dominant_inside(blobs):
    bi = blobs.cell_types.index('B')
    assert blobs.units["n_units"] == 4
    assert max(l.units["n_placed"] for l in _sample_units(blobs, range(6))) == 4
    for lay in _sample_units(blobs):
        _basic_checks(lay)
        # The source's min center separation (140 um in a 300 um window) can be
        # infeasible for four grains: a shortfall is reported, never silent.
        assert lay.units["n_placed"] + lay.units["shortfall"] == 4
        assert ("units_unsatisfied" in lay.flags) == (lay.units["shortfall"] > 0)
        for k in range(lay.units["n_placed"]):
            assert lay.composition[bi][_grain_s(lay, k) < 1.0].mean() > 0.5


def test_uniform_control_layout_matches_units_none():
    xy, r = _unit_points(5)
    t = np.random.default_rng(5).choice(['A', 'B', 'C'], len(r))
    a = _unit_fit(xy[:, 0], xy[:, 1], r, t)
    b = _unit_fit(xy[:, 0], xy[:, 1], r, t, units="none")
    la, lb = a.sample_layout(rng=3), b.sample_layout(rng=3)
    assert a.units == {} and la.units == {}
    assert np.array_equal(la.intensity, lb.intensity)
    assert np.array_equal(la.composition, lb.composition)


def test_units_profiles_round_trip(follicles):
    again = DensityModel.from_dict(json.loads(json.dumps(follicles.to_dict())))
    assert again.units == follicles.units and again.units["profiles"]
    a, b = follicles.sample_layout(rng=4), again.sample_layout(rng=4)
    assert np.array_equal(a.intensity, b.intensity)


def test_unit_free_voids_is_a_copy(glomeruli):
    before = json.dumps(glomeruli.voids, sort_keys=True, default=float)
    free = glomeruli._unit_free_voids()
    assert json.dumps(glomeruli.voids, sort_keys=True, default=float) == before
    assert free.get("n_holes", 0) == 0 or free is not glomeruli.voids
