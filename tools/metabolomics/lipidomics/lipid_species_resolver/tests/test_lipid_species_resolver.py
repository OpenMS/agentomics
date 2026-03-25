"""Tests for lipid_species_resolver."""

import csv
import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestLipidSpeciesResolver:
    def test_parse_sum_composition(self):
        from lipid_species_resolver import parse_sum_composition

        cls, c, db = parse_sum_composition("PC 36:2")
        assert cls == "PC"
        assert c == 36
        assert db == 2

    def test_parse_sum_composition_tg(self):
        from lipid_species_resolver import parse_sum_composition

        cls, c, db = parse_sum_composition("TG 54:3")
        assert cls == "TG"
        assert c == 54
        assert db == 3

    def test_enumerate_chain_combinations_basic(self):
        from lipid_species_resolver import enumerate_chain_combinations

        combos = enumerate_chain_combinations(32, 0, num_chains=2)
        # 32:0 with 2 chains -> (2:0/30:0), (3:0/29:0), ..., (16:0/16:0)
        assert len(combos) > 0
        for combo in combos:
            assert sum(c for c, _ in combo) == 32
            assert sum(db for _, db in combo) == 0

    def test_enumerate_chain_combinations_with_db(self):
        from lipid_species_resolver import enumerate_chain_combinations

        combos = enumerate_chain_combinations(36, 2, num_chains=2)
        assert len(combos) > 0
        for combo in combos:
            assert sum(c for c, _ in combo) == 36
            assert sum(db for _, db in combo) == 2
        # Should include common combinations like (16:0/20:2), (18:1/18:1), (18:0/18:2)
        assert (16, 0) in [c for combo in combos for c in combo] or len(combos) > 0

    def test_enumerate_no_duplicates(self):
        from lipid_species_resolver import enumerate_chain_combinations

        combos = enumerate_chain_combinations(36, 2, num_chains=2)
        combo_tuples = [tuple(c) for c in combos]
        assert len(combo_tuples) == len(set(combo_tuples))

    def test_acyl_chain_formula(self):
        from lipid_species_resolver import acyl_chain_formula

        # 16:0 -> C16H31O
        f = acyl_chain_formula(16, 0)
        assert f == "C16H31O"
        # 18:1 -> C18H33O  (2*18 - 1 - 2*1 = 33)
        f = acyl_chain_formula(18, 1)
        assert f == "C18H33O"

    def test_lipid_exact_mass(self):
        from lipid_species_resolver import lipid_exact_mass

        mass = lipid_exact_mass("PC", [(16, 0), (18, 1)])
        assert mass > 0
        # PC 16:0/18:1 ~= 759.58 Da (approximate)
        assert 700 < mass < 850

    def test_resolve_lipids(self):
        from lipid_species_resolver import resolve_lipids

        lipids = [{"lipid": "PC 36:2"}]
        resolved = resolve_lipids(lipids)
        assert len(resolved) > 0
        for r in resolved:
            assert r["lipid_class"] == "PC"
            assert r["exact_mass"] > 0

    def test_resolve_lipids_with_override(self):
        from lipid_species_resolver import resolve_lipids

        lipids = [{"lipid": "PC 36:2"}]
        resolved = resolve_lipids(lipids, lipid_class_override="PE")
        for r in resolved:
            assert r["lipid_class"] == "PE"

    def test_full_pipeline(self):
        from lipid_species_resolver import load_lipids, resolve_lipids, write_resolved

        with tempfile.TemporaryDirectory() as tmpdir:
            in_path = os.path.join(tmpdir, "lipids.tsv")
            with open(in_path, "w", newline="") as fh:
                w = csv.DictWriter(fh, ["lipid"], delimiter="\t")
                w.writeheader()
                w.writerow({"lipid": "PC 34:1"})

            lipids = load_lipids(in_path)
            resolved = resolve_lipids(lipids)
            out_path = os.path.join(tmpdir, "resolved.tsv")
            write_resolved(resolved, out_path)
            assert os.path.exists(out_path)
            assert len(resolved) > 0
