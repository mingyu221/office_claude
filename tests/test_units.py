import numpy as np
import pandas as pd
import pytest

from codhqsar import config as C
from codhqsar import evaluate as E
from codhqsar import labels as L
from codhqsar import splits as S
from codhqsar.domain import fit_domain
from codhqsar.features.composition import composition
from codhqsar.features.grid import center_dev, make_grid
from codhqsar.features.site import residue_properties


class TestLabels:
    def test_roundtrip(self):
        assert L.to_raw(L.to_score(-4.38)) == pytest.approx(-4.38)

    def test_score_is_positive(self):
        assert L.to_score(-4.38) > 0

    def test_missing_inversion_is_caught(self):
        with pytest.raises(AssertionError):
            L.assert_score_sign([-4.38, -5.0])

    def test_double_inversion_is_caught(self):
        twice = L.to_score(L.to_score(-4.38))
        with pytest.raises(AssertionError):
            L.assert_score_sign([twice])

    def test_noise_floor_gate(self):
        assert not L.resolvable(0.183)      # variant IQR, below the floor
        assert L.resolvable(1.020)          # homolog IQR, above it

    def test_rank_puts_strongest_first(self):
        s = pd.Series([4.0, 6.0, 5.0])
        assert L.rank_better(s).tolist() == [3.0, 1.0, 2.0]


class TestGrid:
    def test_geometry_matches_docking_box(self):
        g = make_grid()
        assert g.spacing == C.SPACING and g.npts == C.NPTS
        assert C.HALF_EDGE == pytest.approx(15.0)

    def test_mask_respects_radius(self):
        g = make_grid(radius=6.0)
        assert np.linalg.norm(g.offsets, axis=1).max() <= 6.0 + 1e-9

    def test_columns_are_stable_across_calls(self):
        a, b = make_grid(), make_grid()
        assert np.array_equal(a.index, b.index)

    def test_ref_origin_has_zero_deviation(self):
        assert center_dev(C.REF_CENTER) == pytest.approx(0.0)

    def test_stride_must_divide_npts(self):
        with pytest.raises(ValueError):
            make_grid(npts=80, stride=3)


class TestSplits:
    def _frame(self):
        return pd.DataFrame({
            "accession": [f"a{i}" for i in range(20)],
            "cluster": np.repeat(np.arange(10), 2),
            "cband": ["40-50"] * 20,
            "tier": ["core", "expansion"] * 10,
            "split": ["train"] * 16 + ["test"] * 4,
        })

    def test_leakage_detected(self):
        df = self._frame()
        with pytest.raises(AssertionError):
            S.assert_no_cluster_leakage(df.iloc[:12], df.iloc[8:])

    def test_expansion_in_test_rejected(self):
        df = self._frame()
        df.loc[df.index[-1], "tier"] = "expansion"
        df.loc[df.index[-1], "split"] = "test"
        with pytest.raises(AssertionError):
            S.split_from_targets(df)

    def test_folds_never_split_a_cluster(self):
        df = self._frame()
        for tr, te in S.cluster_folds(df, 5):
            assert not set(df.iloc[tr]["cluster"]) & set(df.iloc[te]["cluster"])

    def test_nested_check(self):
        S.assert_nested([{1}, {1, 2}, {1, 2, 3}])
        with pytest.raises(AssertionError):
            S.assert_nested([{1, 9}, {1, 2}])


class TestEvaluate:
    def test_ols_recovers_planted_effect(self):
        rng = np.random.default_rng(0)
        n = 200
        df = pd.DataFrame({"RMSD": rng.uniform(0.5, 2.5, n),
                           "center_dev": rng.uniform(0, 1.5, n)})
        df["score"] = 5.0 - 0.5 * df["RMSD"] + rng.normal(0, 0.1, n)
        fit = E.ols(df, "score", ["RMSD", "center_dev"])
        beta = dict(zip(fit.names, fit.beta))
        assert beta["RMSD"] == pytest.approx(-0.5, abs=0.05)

    def test_categorical_covariate_is_dummy_coded(self):
        rng = np.random.default_rng(1)
        n = 120
        df = pd.DataFrame({
            "RMSD": rng.uniform(0.5, 2.5, n),
            "center_dev": rng.uniform(0, 1.5, n),
            "mode": rng.choice(["global", "local"], n),
            "struct_source": rng.choice(["AFDB", "ColabFold"], n),
            "crit_clash": rng.poisson(1.0, n),
            "pLDDT": rng.uniform(93, 97, n),
            "cband": rng.choice(["40-50", "30-40"], n),
        })
        df["score"] = rng.normal(5, 1, n)
        fit = E.covariate_regression(df)
        assert any(t.startswith("mode[") for t in fit.names)
        assert np.isfinite(fit.r2)

    def test_covariate_regression_demands_confounders(self):
        df = pd.DataFrame({"score": [1.0, 2.0, 3.0], "cband": ["a", "b", "a"]})
        with pytest.raises(KeyError):
            E.covariate_regression(df)

    def test_thin_bands_are_pooled(self):
        rng = np.random.default_rng(2)
        n = 60
        df = pd.DataFrame({"cband": ["40-50"] * 25 + ["30-40"] * 25
                                    + ["≥60"] * 5 + ["<30"] * 5})
        df["score"] = rng.normal(size=n)
        df["pred"] = df["score"] + rng.normal(0, 0.3, n)
        bands = set(E.report_by_band(df, "score", "pred")["band"])
        assert "tails (pooled)" in bands
        assert "≥60" not in bands

    def test_noise_floor_and_snr(self):
        rng = np.random.default_rng(3)
        base = rng.normal(4.8, 1.0, 40)
        df = pd.DataFrame({
            "cluster": np.repeat(np.arange(40), 2),
            "tier": ["core", "expansion"] * 40,
            "score": np.concatenate([[b, b + rng.normal(0, 0.1)] for b in base]),
        })
        nf = E.noise_floor_from_pairs(df)
        assert nf.n_pairs == 40
        assert nf.median_within < 0.3
        assert nf.snr_canonical > 5


class TestDomain:
    def test_far_points_are_outside(self):
        rng = np.random.default_rng(0)
        tr = rng.normal(size=(80, 3))
        dom = fit_domain(tr)
        assert dom.inside(rng.normal(size=(20, 3))).mean() > 0.5
        assert not dom.inside(rng.normal(size=(20, 3)) + 20).any()


class TestFeatures:
    def test_composition_fractions_sum_to_one(self):
        v = composition("ACDEFGHIKLMNPQRSTVWY")
        assert v[:20].sum() == pytest.approx(1.0)

    def test_nonstandard_residues_count_in_length_only(self):
        v = composition("AAXX")
        assert v[20] == 4.0
        assert v[:20].sum() == pytest.approx(0.5)

    def test_residue_properties(self):
        assert residue_properties("ASP")[0] == -1.0
        assert residue_properties("PHE")[2] == 1.0
        assert residue_properties("TRP")[1] > residue_properties("GLY")[1]
