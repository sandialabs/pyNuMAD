# for type hints
from numpy import ndarray

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
        self.bill_of_materials.generate(self.definition, self.keypoints)
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
