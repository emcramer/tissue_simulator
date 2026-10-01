"""Tests for density-aware layouts in ``tissue_simulator.density``."""

import json

import numpy as np
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
    return DensityModel.fit(x, y, r, t, bounds=(0, 0, 240, 120), strategy="adaptive",
                            n_compartments=1, seed=0)


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
