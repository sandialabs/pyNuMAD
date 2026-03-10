# for type hints
from numpy import ndarray

import math
import numpy as np
import pandas as pd

from pynumad.io.yaml_to_blade import yaml_to_blade
from pynumad.io.excel_to_blade import excel_to_blade
from pynumad.utils.interpolation import interpolator_wrap
from pynumad.objects.geometry import Geometry
from pynumad.objects.keypoints import KeyPoints
from pynumad.objects.definition import Definition
from pynumad.objects.bom import BillOfMaterials
from pynumad.objects.materialdb import MaterialDatabase
from pynumad.objects.stackdb import StackDatabase

_TEA_COLS = [
    'BOM', 'Level', 'Component', 'Material', 'Type', 'Qty',
    'Thickness', 'Area', 'Length', 'Start_station', 'End_station', 'Perimeter',
]


class Blade:
    """Blade class
    
    Parameters
    ----------
    filename : str, optional
        Directory and filename of blade input file to load into the
        Blade object.
        
    Attributes
    ----------
    name : str
        Name of the blade
    definition : Definition
        Object containing the definition of the blade.
    geometry : Geometry
        Object containing the interpolated geometry of the blade
    keypoints : KeyPoints
        Object containing information about keypoint locations and areas
    bill_of_materials : BillOfMaterials
    stackdb : StackDatabase
    materialdb : MaterialDatabase
    mesh_size : float
        Target element size used by mesh generation utilities (default 0.45).

    Example
    -------
    
    blade = Blade("path/to/blade.yaml")
    """

    def __init__(self, filename: str = None):
        self.name: str = None
        self.definition: Definition = Definition()
        self.geometry: Geometry = Geometry()
        self.keypoints: KeyPoints = KeyPoints()
        self.bill_of_materials: BillOfMaterials = BillOfMaterials()
        self.stackdb: StackDatabase = StackDatabase()
        self.materialdb: MaterialDatabase = MaterialDatabase()
        self.mesh_size: float = 0.45  # target element size for mesh generation

        if filename:
            if "yaml" in filename or "yml" in filename:
                self.read_yaml(filename)
            elif "xls" in filename or "xlsx" in filename:
                self.read_excel(filename)
            else:
                raise Exception(
                    "Unknown filetype. Currently supported inputs are excel and yaml files."
                )
            self.name = filename.split(".")[0]
        else:
            self.name = "blade"

    @property
    def ispan(self) -> ndarray:
        """Interpolated span stations. Delegates to ``definition.ispan``."""
        return self.definition.ispan

    @ispan.setter
    def ispan(self, value: ndarray):
        self.definition.ispan = value

    def __str__(self):
        attributes = ""
        for attr_name, attr_value in vars(self).items():
            if isinstance(attr_value, list):
                attributes += f"{attr_name}={len(attr_value)}, "
            elif isinstance(attr_value, np.ndarray):
                attributes += f"{attr_name}={attr_value.shape}, "
            else:
                attributes += f"{attr_name}={attr_value}, "
        return f"Blade with {attributes[:-2]}"

    def read_yaml(self, filename: str):
        """Populate blade attributes with yaml file data

        Parameters
        ----------
        filename: str
            name of yaml file to be read

        Returns
        -------
        self

        """
        yaml_to_blade(self, filename)
        return self

    def read_excel(self, filename: str):
        """Populate blade attributes with excel file data

        Parameters
        ----------
        filename: str
            name of excel file to be read

        Returns
        -------
        self

        """
        excel_to_blade(self, filename)
        return self

    def update_blade(self):
        """Generates geometry, keypoints, bill of materials, 
        stack database, and material database based on the
        blade definition. 
        """
        self.geometry.generate(self.definition)
        self.keypoints.generate(self.definition, self.geometry)
        self.bill_of_materials.generate(
            self.definition.ispan,
            self.definition.components,
            self.definition.materials,
            self.keypoints,
        )
        self.stackdb.generate(self.keypoints, self.bill_of_materials)
        self.materialdb.generate(self.definition.materials, self.stackdb)
        return self

    def expand_blade_geometry_te(self, min_edge_lengths):
        """
        TODO: docstring
        """
        self.geometry.expand_blade_geometry_te(min_edge_lengths)
        self.keypoints.generate(self.definition, self.geometry)
        return

    def export_bom_tea(self, file, sheet_name, segments):
        """Generate a segmented BOM and export to Excel for TEA cost analysis.

        Parameters
        ----------
        file : str
            Output Excel file path.
        sheet_name : str
            Sheet name in the Excel file.
        segments : list of dict
            Blade segments. Each dict has:
              'span_nd'    : [start, end] normalized spanwise extents
              'hp_extents' : [key1, key2] HP chordwise bounds (e.g. ['le','te'])
              'lp_extents' : [key1, key2] LP chordwise bounds

        Returns
        -------
        pd.DataFrame
            The assembled bomTEA table.
        """
        bom = self.bill_of_materials
        definition = self.definition

        bom.generate_segmented(
            segments,
            self.ispan,
            definition.components,
            definition.materials,
            self.keypoints,
        )

        rows = []

        # Assembly-level rows (finishing steps, root hardware, bonds, LPS, etc.)
        rows.extend(self._build_assembly_rows(segments))

        # Ordered composite + mold rows per segment
        mold_types = ['shell', 'spar', 'LEreinf', 'TEreinf', 'root']
        mold_suffixes = [
            '_shell_mold', '_spar_prefab', '_LEreinf_preform',
            '_TEreinf_prefab', '_root_prefab',
        ]
        mold_levels = [5, 6, 6, 6, 6]

        for ks in range(len(segments)):
            for side in ('hp', 'lp'):
                region_str = f'{side.upper()}{ks + 1}'
                side_df = bom.hp if side == 'hp' else bom.lp
                if side_df.empty:
                    continue
                seg_rows = side_df[side_df['segment_id'] == region_str]

                for mm, (mold_field, mold_suffix) in enumerate(
                    zip(mold_types, mold_suffixes)
                ):
                    if mm == 0:
                        is_prefab = (
                            seg_rows['name'].str.contains('spar', case=False, na=False)
                            | seg_rows['name'].str.contains('reinf', case=False, na=False)
                            | seg_rows['name'].str.contains('root', case=False, na=False)
                        )
                        layer_rows = seg_rows[~is_prefab]
                    elif mm == 1:
                        layer_rows = seg_rows[
                            seg_rows['name'].str.contains('spar', case=False, na=False)
                        ]
                    elif mm == 2:
                        layer_rows = seg_rows[
                            seg_rows['name'].str.contains('le', case=False, na=False)
                            & seg_rows['name'].str.contains('reinf', case=False, na=False)
                        ]
                    elif mm == 3:
                        layer_rows = seg_rows[
                            seg_rows['name'].str.contains('te', case=False, na=False)
                            & seg_rows['name'].str.contains('reinf', case=False, na=False)
                        ]
                    else:
                        layer_rows = seg_rows[
                            seg_rows['name'].str.contains('root', case=False, na=False)
                        ]

                    if layer_rows.empty:
                        continue

                    level = mold_levels[mm]
                    mold_str = region_str + mold_suffix
                    mold_field_df = bom.mold.get(mold_field, pd.DataFrame())
                    if not mold_field_df.empty:
                        for _, m_row in mold_field_df[
                            mold_field_df['Component'].str.contains(
                                mold_str, case=False, na=False
                            )
                        ].iterrows():
                            rows.append(_mold_row(m_row, level))

                    for _, l_row in layer_rows.iterrows():
                        mat_name = definition.materials[l_row['materialid']].name
                        rows.append(_layer_row(l_row, mat_name, level + 1))

            # Shear web rows for this segment
            if not bom.sw.empty:
                for kw in sorted(bom.sw['web_id'].unique()):
                    kw_rows = bom.sw[bom.sw['web_id'] == kw]
                    seg_ids = [
                        sid for sid in kw_rows['segment_id'].unique()
                        if sid.endswith(f'_SW{kw + 1}')
                        and sid.split('_SW')[0][-1] == str(ks + 1)
                    ]
                    if not seg_ids:
                        continue
                    seg_id = seg_ids[0]
                    web_seg_rows = kw_rows[kw_rows['segment_id'] == seg_id]

                    web_mold_df = bom.mold.get('web', pd.DataFrame())
                    if not web_mold_df.empty:
                        for _, m_row in web_mold_df[
                            web_mold_df['Component'].str.contains(
                                seg_id, case=False, na=False
                            )
                        ].iterrows():
                            rows.append(_mold_row(m_row, 5))

                    for _, l_row in web_seg_rows.iterrows():
                        mat_name = definition.materials[l_row['materialid']].name
                        rows.append(_layer_row(l_row, mat_name, 6))

        bom_tea = pd.DataFrame(rows)
        bom_tea.insert(0, 'BOM', range(1, len(bom_tea) + 1))
        bom_tea = bom_tea.reindex(columns=_TEA_COLS)
        bom_tea.to_excel(file, sheet_name=sheet_name, index=False, engine='openpyxl')
        bom.bom_tea = bom_tea
        return bom_tea

    def _build_assembly_rows(self, segments, **kwargs):
        """Build assembly-level BOM rows (hardware, adhesive, LPS, etc.).

        Mirrors MATLAB's calc_bom_assembly. All dimensional parameters can be
        overridden via keyword arguments.

        Returns
        -------
        list of dict
        """
        bom = self.bill_of_materials
        span = self.ispan[-1]
        root_diameter = float(self.geometry.ichord[0])

        p = {
            'overlaminate_LE': 1 / 3,
            'overlaminate_TE': 1 / 3,
            'overlaminate_ply': 3,
            'overlaminate_width_mm': 500,
            'LPS_type': 'cable',
            'balance_box_span': 35,
            'meter_per_root_bolt': np.pi * 3.386 / 90,
            'le_bond_thick_mm': 5,
            'le_bond_width_mm': 100,
            'te_bond_thick_mm': 5,
            'te_bond_width_mm': 100,
            'sw_bond_thick_mm': 5,
            'sw_bond_width_mm': 100,
            'sw_clips_material': 'Saertex(DB)',
            'sw_meter_per_clip': 0.5,
        }
        p.update(kwargs)

        asm = 'MP0_'
        rows = []

        def _asm_row(**fields):
            r = {k: float('nan') for k in _TEA_COLS if k not in ('BOM', 'Component')}
            r['Component'] = ''
            r['Material'] = ''
            r['Type'] = ''
            r.update(fields)
            return r

        # Blade finishing hierarchy
        for level, component in [
            (1, f'{asm}Finished Blade'),
            (2, f'{asm}Machined Rotor Blade'),
            (3, f'{asm}Wet Finished Laminated Blade'),
            (4, f'{asm}Bonded Assembly'),
        ]:
            rows.append(_asm_row(Level=level, Component=component, Type='Assembly', Qty=1))

        # Root hardware
        num_root_bolts = math.floor(np.pi * root_diameter / p['meter_per_root_bolt'])
        for level, component, material, qty_type, qty in [
            (2, f'{asm}Root Close Out System', '', 'Assembly', 1),
            (3, f'{asm}RCO Platform', 'RCO PLATFORM', 'Diameter', root_diameter),
            (3, f'{asm}RCO Attachment', 'RCO FLANGE', 'Diameter', root_diameter),
            (3, f'{asm}RCO Hatch', 'HATCH', 'Each', 1),
            (2, f'{asm}Root Attachment System', '', 'Assembly', 1),
            (3, f'{asm}Barrel Nuts', 'BARREL NUT', 'Each', num_root_bolts),
            (3, f'{asm}Root Bolts', 'ROOT BOLT', 'Each', num_root_bolts),
            (3, f'{asm}O-Ring', 'O-RING', 'Each', num_root_bolts),
            (3, f'{asm}Sealant', 'SEALANT', 'Cubic Meters', float('nan')),
            (2, f'{asm}Zero Degree Kit', 'ZERO DEGREE KIT', 'Each', 1),
            (2, f'{asm}Labels', 'LABELS', 'Each', 15),
            (2, f'{asm}Nameplates', 'NAMEPLATE', 'Each', 5),
        ]:
            rows.append(_asm_row(
                Level=level, Component=component, Material=material,
                Type=qty_type, Qty=qty,
            ))

        # Overlaminate placeholders
        for level, component in [
            (4, f'{asm}External Overlaminates'),
            (99, f'{asm}Missing: Shear Web Bonding'),
            (99, f'{asm}Missing: Blade Closing Materials'),
        ]:
            rows.append(_asm_row(Level=level, Component=component, Type='Kit', Qty=1))

        # LE/TE adhesive per segment
        for ks in range(len(segments)):
            asm_ks = f'MP{ks + 1}_'
            le_m = bom.lebond[ks] / 1000 if ks < len(bom.lebond) else 0.0
            te_m = bom.tebond[ks] / 1000 if ks < len(bom.tebond) else 0.0
            if le_m > 0:
                le_area = p['le_bond_width_mm'] / 1000 * le_m
                rows.append(_asm_row(
                    Level=5,
                    Component=f'{asm_ks}Leading Edge Bond Paste',
                    Material='Paste',
                    Type='Cubic Meters',
                    Qty=le_area * p['le_bond_thick_mm'] / 1000,
                    Thickness=p['le_bond_thick_mm'],
                    Area=le_area,
                    Length=le_m,
                ))
            if te_m > 0:
                te_area = p['te_bond_width_mm'] / 1000 * te_m
                rows.append(_asm_row(
                    Level=5,
                    Component=f'{asm_ks}Trailing Edge Bond Paste',
                    Material='Paste',
                    Type='Cubic Meters',
                    Qty=te_area * p['te_bond_thick_mm'] / 1000,
                    Thickness=p['te_bond_thick_mm'],
                    Area=te_area,
                    Length=te_m,
                ))

        # Balance boxes
        n_balance = math.floor(span / p['balance_box_span'])
        rows.append(_asm_row(
            Level=5, Component=f'{asm}Balance Box',
            Type='Assembly', Qty=n_balance,
        ))

        # LPS
        for level, component, material, qty_type, qty in [
            (2, f'{asm}Lightning Protection System', '', 'Assembly', 1),
            (3, f'{asm}LPS Tip Receptor', 'LPS Tip Receptor', 'Each', 1),
            (3, f'{asm}LPS Cable', 'Electrical Cable', 'Meter', span),
            (3, f'{asm}LPS Root Hardware', '', 'Kit', 1),
        ]:
            rows.append(_asm_row(
                Level=level, Component=component, Material=material,
                Type=qty_type, Qty=qty,
            ))

        # Shear web clips and adhesive per web per segment
        for kw, sw_bond_arr in enumerate(bom.swbonds):
            for ks in range(len(segments)):
                asm_ks = f'MP{ks + 1}_'
                hp_len = sw_bond_arr[0, ks] / 1000
                lp_len = sw_bond_arr[1, ks] / 1000
                for side_len, side_label in [(hp_len, 'Pressure'), (lp_len, 'Suction')]:
                    if side_len > 0:
                        bond_area = p['sw_bond_width_mm'] / 1000 * side_len
                        rows.append(_asm_row(
                            Level=5,
                            Component=f'{asm_ks}Shear Web #{kw + 1} Bond Paste, {side_label} Side',
                            Material='Paste',
                            Type='Cubic Meters',
                            Qty=bond_area * p['sw_bond_thick_mm'] / 1000,
                            Thickness=p['sw_bond_thick_mm'],
                            Area=bond_area,
                            Length=side_len,
                        ))
                        rows.append(_asm_row(
                            Level=5,
                            Component=f'{asm_ks}Shear Web #{kw + 1} Clips, {side_label} Side',
                            Material=p['sw_clips_material'],
                            Type='Piece',
                            Qty=math.ceil(side_len / p['sw_meter_per_clip']),
                            Length=side_len,
                        ))

        return rows

    def add_interpolated_station(self, span_location: float):
        """Adds an interpolated station to blade geometry

        Parameters
        ----------
        span_location : float
            location along span between 0 and 1.

        Returns
        -------
        int
            integer index where the new span was inserted
        """
        x0 = self.definition.ispan.copy()

        if span_location < self.definition.ispan[-1] and span_location > 0:
            for i_span, spanLocation in enumerate(self.definition.ispan[1:]):
                if span_location < spanLocation:
                    insertIndex = i_span + 1
                    break
        else:
            raise ValueError(
                f"A new span location with value {span_location} is not possible."
            )

        new_ispan = np.insert(x0, insertIndex, span_location)

        # Interpolate all ispan-grid columns in one loop, then rebuild the DataFrame.
        new_data = {}
        for col in self.definition.ispan_data.columns:
            old_values = self.definition.ispan_data[col].to_numpy(dtype=float)
            new_data[col] = interpolator_wrap(x0, old_values, new_ispan)
        self.definition.ispan_data = pd.DataFrame(new_data, index=new_ispan)

        self.update_blade()
        return insertIndex


def _mold_row(m_row, level):
    """Format a mold/prefab DataFrame row for the bomTEA output."""
    return {
        'Level': level,
        'Component': m_row['Component'],
        'Material': '',
        'Type': 'Assembly',
        'Qty': 1,
        'Thickness': float('nan'),
        'Area': m_row['Area'],
        'Length': m_row['Length'],
        'Start_station': m_row['Start_station'],
        'End_station': m_row['End_station'],
        'Perimeter': m_row['Perimeter'],
    }


def _layer_row(l_row, mat_name, level):
    """Format a composite layer BOM row for the bomTEA output."""
    return {
        'Level': level,
        'Component': f"{l_row['segment_id']}_{l_row['name']}",
        'Material': mat_name,
        'Type': 'Piece',
        'Qty': 1,
        'Thickness': l_row['thickness'],
        'Area': l_row['area'],
        'Length': l_row['arc_length'],
        'Start_station': l_row['beginsta'],
        'End_station': l_row['endsta'],
        'Perimeter': l_row['perimeter'],
    }
