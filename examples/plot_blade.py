"""Demonstrate PyVista blade geometry and mesh plots.

Requires: uv pip install pyvista

  - plot_blade_geometry_pv : interactive 3D OML surface
  - plot_shell_mesh_pv     : interactive shell FE mesh with edge wireframe
"""
import os

import pynumad as pynu
from pynumad.mesh_gen.mesh_gen import get_shell_mesh
from pynumad.graphics.graphics import plot_blade_geometry_pv, plot_shell_mesh_pv

yaml_file = os.path.join(os.path.dirname(__file__), "example_data", "blade.yaml")

blade = pynu.Blade()
blade.read_yaml(yaml_file)
blade.update_blade()

plot_blade_geometry_pv(blade)

mesh_data = get_shell_mesh(blade, includeAdhesive=False, elementSize=0.45)
print(f"Shell mesh: {len(mesh_data['nodes'])} nodes, {len(mesh_data['elements'])} elements")

plot_shell_mesh_pv(mesh_data)
