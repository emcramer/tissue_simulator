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
