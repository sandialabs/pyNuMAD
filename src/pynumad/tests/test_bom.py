import os
import unittest

import numpy as np
import pandas as pd

from pynumad.objects.bom import BillOfMaterials, find_layer_extents
from pynumad.objects.blade import Blade

test_data_dir = os.path.join(os.path.dirname(__file__), "test_data")


class TestBillOfMaterialsInit(unittest.TestCase):
    def test_init_creates_empty_dataframes(self):
        bom = BillOfMaterials()
        self.assertIsInstance(bom.hp, pd.DataFrame)
        self.assertIsInstance(bom.lp, pd.DataFrame)
        self.assertIsInstance(bom.sw, pd.DataFrame)
        self.assertTrue(bom.hp.empty)
        self.assertTrue(bom.lp.empty)
        self.assertTrue(bom.sw.empty)

    def test_init_bond_and_weight_attributes_none(self):
        bom = BillOfMaterials()
        self.assertIsNone(bom.lebond)
        self.assertIsNone(bom.tebond)
        self.assertIsNone(bom.total_dryweight)
        self.assertEqual(bom.swbonds, [])
        self.assertEqual(bom.mold, {})


class TestBillOfMaterialsGenerate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        yaml_path = os.path.join(test_data_dir, "blades", "blade.yaml")
        cls.blade = Blade(yaml_path)
        cls.blade.update_blade()

    def test_generate_returns_self(self):
        bom = BillOfMaterials()
        result = bom.generate(
            self.blade.definition.ispan,
            self.blade.definition.components,
            self.blade.definition.materials,
            self.blade.keypoints,
        )
        self.assertIs(result, bom)

    def test_generate_populates_hp_lp_with_expected_columns(self):
        bom = self.blade.bill_of_materials
        expected = [
            "layernum", "materialid", "name", "beginsta", "endsta",
            "maxwidth", "avgwidth", "area", "thickness", "dryweight", "angle",
            "sta_begin_idx", "sta_end_idx", "seg_start", "seg_end",
        ]
        for col in expected:
            self.assertIn(col, bom.hp.columns, msg=f"hp missing {col}")
            self.assertIn(col, bom.lp.columns, msg=f"lp missing {col}")

    def test_generate_sw_has_web_id(self):
        bom = self.blade.bill_of_materials
        self.assertIn("web_id", bom.sw.columns)

    def test_generate_total_dryweight_positive(self):
        bom = self.blade.bill_of_materials
        self.assertIsNotNone(bom.total_dryweight)
        self.assertGreater(bom.total_dryweight, 0)

    def test_generate_lebond_tebond_non_negative(self):
        bom = self.blade.bill_of_materials
        self.assertIsNotNone(bom.lebond)
        self.assertIsNotNone(bom.tebond)
        self.assertGreaterEqual(bom.lebond, 0)
        self.assertGreaterEqual(bom.tebond, 0)

    def test_generate_swbonds_length_matches_webs(self):
        bom = self.blade.bill_of_materials
        n_webs = len(self.blade.keypoints.web_areas) if self.blade.keypoints.web_areas else 0
        self.assertEqual(len(bom.swbonds), n_webs)


class TestBillOfMaterialsGenerateSegmented(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        yaml_path = os.path.join(test_data_dir, "blades", "blade.yaml")
        cls.blade = Blade(yaml_path)
        cls.blade.update_blade()
        cls.segments = [
            {"span_nd": [0.0, 0.5], "hp_extents": ["te", "le"], "lp_extents": ["le", "te"]},
            {"span_nd": [0.5, 1.0], "hp_extents": ["te", "le"], "lp_extents": ["le", "te"]},
        ]

    def test_generate_segmented_returns_self(self):
        bom = BillOfMaterials()
        result = bom.generate_segmented(
            self.segments,
            self.blade.definition.ispan,
            self.blade.definition.components,
            self.blade.definition.materials,
            self.blade.keypoints,
        )
        self.assertIs(result, bom)

    def test_generate_segmented_adds_extra_columns(self):
        bom = BillOfMaterials()
        bom.generate_segmented(
            self.segments,
            self.blade.definition.ispan,
            self.blade.definition.components,
            self.blade.definition.materials,
            self.blade.keypoints,
        )
        self.assertIn("segment_id", bom.hp.columns)
        self.assertIn("arc_length", bom.hp.columns)
        self.assertIn("perimeter", bom.hp.columns)

    def test_generate_segmented_lebond_tebond_are_lists(self):
        bom = BillOfMaterials()
        bom.generate_segmented(
            self.segments,
            self.blade.definition.ispan,
            self.blade.definition.components,
            self.blade.definition.materials,
            self.blade.keypoints,
        )
        self.assertIsInstance(bom.lebond, list)
        self.assertIsInstance(bom.tebond, list)
        self.assertEqual(len(bom.lebond), len(self.segments))
        self.assertEqual(len(bom.tebond), len(self.segments))

    def test_generate_segmented_populates_mold(self):
        bom = BillOfMaterials()
        bom.generate_segmented(
            self.segments,
            self.blade.definition.ispan,
            self.blade.definition.components,
            self.blade.definition.materials,
            self.blade.keypoints,
        )
        self.assertIsInstance(bom.mold, dict)
        for key in ["shell", "spar", "LEreinf", "TEreinf", "root", "web"]:
            if key in bom.mold:
                self.assertIsInstance(bom.mold[key], pd.DataFrame)


class TestFindLayerExtents(unittest.TestCase):
    def test_single_contiguous_region(self):
        layer_dist = np.array([0, 1, 1, 1, 0])
        begin, end = find_layer_extents(layer_dist, 1)
        self.assertEqual(begin, [1])
        self.assertEqual(end, [4])

    def test_layer_above_max_returns_empty(self):
        layer_dist = np.array([0, 1, 1, 0])
        begin, end = find_layer_extents(layer_dist, 2)
        self.assertEqual(begin, [])
        self.assertEqual(end, [])

    def test_full_span_at_layer(self):
        layer_dist = np.array([1, 1, 1])
        begin, end = find_layer_extents(layer_dist, 1)
        self.assertEqual(begin, [0])
        self.assertEqual(end, [2])


if __name__ == "__main__":
    unittest.main()
