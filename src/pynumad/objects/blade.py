import re, warnings
import numpy as np
from copy import deepcopy

from pynumad.io.yaml_to_blade import yaml_to_blade
from pynumad.io.excel_to_blade import excel_to_blade
from pynumad.utils.interpolation import interpolator_wrap
from pynumad.objects.geometry import Geometry
from pynumad.objects.settings import BladeSettings
from pynumad.objects.keypoints import KeyPoints
from pynumad.objects.definition import Definition
from pynumad.objects.bom import BillOfMaterials
from pynumad.objects.materialdb import MaterialDatabase
from pynumad.objects.stackdb import StackDatabase

# for type hints
from numpy import ndarray


class Blade:
    """BladeDef A class definition for wind & water turbine blades.

    Attributes
    ----------
    name : str
        Name of the blade
    definition : Definition
        Object containing the definition of the blade.
    ispan : ndarray
        Array defining the interpolated stations from 0 to 1
    geometry : Geometry
        Object containing the interpolated geometry of the blade
    keypoints : KeyPoints
        Object containing information about keypoint locations and areas
    bill_of_materials : BillOfMaterials
    stackdb : StackDatabase
    materialdb : MaterialDatabase
    settings : BladeSettings

    Example
    -------
    blade = Blade("path/to/blade.yaml")
    """

    def __init__(self, filename: str = None):
        """
        Parameters
        ----------
        filename : str, optional
            Directory and filename of blade input file to load into the
            Blade object.
        """
        self.name: str = None
        self.definition: Definition = Definition()
        self.ispan: ndarray = None
        self.geometry: Geometry = Geometry()
        self.keypoints: KeyPoints = KeyPoints()
        self.bill_of_materials: BillOfMaterials = BillOfMaterials()
        self.stackdb: StackDatabase = StackDatabase()
        self.materialdb: MaterialDatabase = MaterialDatabase()
        self.settings: BladeSettings = BladeSettings()

        # read input file
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
        return

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
        stacks database, and materials database based on the
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
        x0 = self.ispan

        if span_location < self.ispan[-1] and span_location > 0:
            for iSpan, spanLocation in enumerate(self.ispan[1:]):
                if span_location < spanLocation:
                    insertIndex = iSpan + 1
                    break
        else:
            raise ValueError(
                f"A new span location with value {span_location} is not possible."
            )

        self.ispan = np.insert(self.ispan, insertIndex, np.array([span_location]))

        self.leband = interpolator_wrap(x0, self.leband, self.ispan)
        self.teband = interpolator_wrap(x0, self.teband, self.ispan)
        self.sparcapwidth_hp = interpolator_wrap(x0, self.sparcapwidth_hp, self.ispan)
        self.sparcapwidth_lp = interpolator_wrap(x0, self.sparcapwidth_lp, self.ispan)
        self.sparcapoffset_hp = interpolator_wrap(x0, self.sparcapoffset_hp, self.ispan)
        self.sparcapoffset_lp = interpolator_wrap(x0, self.sparcapoffset_lp, self.ispan)

        self.update_blade()
        return insertIndex
