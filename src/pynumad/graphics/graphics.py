import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from pynumad.utils.interpolation import interpolator_wrap


def plot_airfoil(airfoil):
    """Plot airfoil"""
    fig, ax = plt.subplots()
    # ax[0].plot(self.x,self.y,'.-')
    ax.plot(airfoil.coordinates[:, 0], airfoil.coordinates[:, 1], ".-", color="black")
    ax.plot(airfoil.c, airfoil.camber, color="red")
    # mtx = self.maxthick * np.array([1,1])
    # kn = find(self.c >= self.maxthick,1)
    # mty = self.camber(kn) + self.thickness(kn) * np.array([0.5,- 0.5])
    # line(mtx,mty,'LineStyle',':','Color','k')
    # else:
    fig.show()
    return ax


def plot_regions(blade):
    fig, ax = plt.subplots()

    n = blade.keypoints.shape[0]
    for kn in range(n):
        blade.hgKeypoints[kn] = ax.plot(
            blade.keypoints[kn, 2, :],
            blade.keypoints[kn, 0, :],
            blade.keypoints[kn, 1, :],
        )
    fig.show()
    return ax


def plot_geometry(self):
    fig, ax = plt.subplots()
    n = self.geometry.shape[2]
    for k in range(n):
        self.hgGeometry[k] = ax.plot(
            self.geometry[:, 2, k], self.geometry[:, 0, k], self.geometry[:, 1, k]
        )
    fig.show()
    return ax


def plot_profile(blade, k):
    fig, ax = plt.subplots()
    ax.plot(blade.profiles[:, 0, k], blade.profiles[:, 1, k], ".-")
    fig.show()
    return ax


def plot_blade_geometry(blade):
    """Plot blade cross-section profiles in 3D.

    Each interpolated spanwise station is drawn as a closed profile line.
    Axes are (x, y, z) matching the blade coordinate system stored in
    ``geometry.coordinates`` with shape ``(n_pts, 3, n_stations)``.
    """
    coords = blade.geometry.coordinates
    n_stations = coords.shape[2]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    for k in range(n_stations):
        ax.plot(coords[:, 0, k], coords[:, 1, k], coords[:, 2, k], color="steelblue", lw=0.6)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z (span)")
    ax.set_title("Blade geometry")

    x_range = np.ptp(coords[:, 0, :])
    y_range = np.ptp(coords[:, 1, :])
    z_range = np.ptp(coords[:, 2, :])
    ax.set_box_aspect([x_range, y_range, z_range])
    return ax



def plot_segment_mass(bom, segments):
    """Plot per-segment mass breakdown as a stacked bar chart.

    Parameters
    ----------
    bom : BillOfMaterials
        Must have been populated by ``generate_segmented`` so that
        ``bom.hp``, ``bom.lp``, and ``bom.sw`` contain a ``segment_id``
        column.
    segments : list of dict
        Segment definitions (used only for ordering / labelling).
    """
    g_to_kg = 0.001
    n_segs = len(segments)
    seg_labels = [f"Seg {i + 1}\n{s['span_nd']}" for i, s in enumerate(segments)]

    hp_mass = np.zeros(n_segs)
    lp_mass = np.zeros(n_segs)
    sw_mass = np.zeros(n_segs)

    for i in range(n_segs):
        seg_num = str(i + 1)
        if not bom.hp.empty:
            mask = bom.hp["segment_id"].str.endswith(seg_num)
            hp_mass[i] = bom.hp.loc[mask, "dryweight"].sum()
        if not bom.lp.empty:
            mask = bom.lp["segment_id"].str.endswith(seg_num)
            lp_mass[i] = bom.lp.loc[mask, "dryweight"].sum()
        if not bom.sw.empty:
            mask = bom.sw["segment_id"].str.contains(f"_{seg_num}_SW|^{seg_num}_SW", regex=True)
            sw_mass[i] = bom.sw.loc[mask, "dryweight"].sum()

    x = np.arange(n_segs)
    fig, ax = plt.subplots()
    ax.bar(x, hp_mass, label="HP shell")
    ax.bar(x, lp_mass, bottom=hp_mass, label="LP shell")
    ax.bar(x, sw_mass, bottom=hp_mass + lp_mass, label="Shear web")
    ax.set_xticks(x)
    ax.set_xticklabels(seg_labels, fontsize=8)
    ax.set_ylabel("Mass (kg)")
    ax.set_title("Segment mass breakdown")
    ax.legend()
    return ax


def plot_blade_geometry_pv(blade, **plotter_kwargs):
    """Render the blade OML surface in PyVista as a structured grid.

    ``geometry.coordinates`` has shape ``(n_pts, 3, n_stations)``.  The points
    are arranged as a structured surface: n_pts along the profile direction and
    n_stations along the span, so a ``pv.StructuredGrid`` can be built directly.

    Parameters
    ----------
    blade : Blade
    **plotter_kwargs
        Forwarded to ``pv.Plotter`` (e.g. ``off_screen=True``).

    Returns
    -------
    pv.Plotter
    """
    import pyvista as pv

    coords = blade.geometry.coordinates          # (n_pts, 3, n_stations)
    n_pts, _, n_sta = coords.shape
    # StructuredGrid expects (n_sta, n_pts) meshgrid ordering with separate x/y/z
    x = coords[:, 0, :].T   # (n_sta, n_pts)
    y = coords[:, 1, :].T
    z = coords[:, 2, :].T
    grid = pv.StructuredGrid(x, y, z)

    plotter = pv.Plotter(**plotter_kwargs)
    plotter.add_mesh(grid, color="steelblue", show_edges=False, opacity=0.85)
    plotter.add_axes()
    plotter.show_grid()
    plotter.show()
    return plotter


def plot_blade_mesh_pv(mesh3d, scalars=None, scalar_name="", **plotter_kwargs):
    """Render a Mesh3D object in PyVista as an unstructured grid.

    Supports hex, wedge, and tet elements.  An optional per-element scalar
    array (e.g. laminate thickness, material ID) can be passed for coloring.

    Parameters
    ----------
    mesh3d : Mesh3D
    scalars : array-like, optional
        Per-element scalar values for color mapping.
    scalar_name : str
        Label shown on the scalar bar.
    **plotter_kwargs
        Forwarded to ``pv.Plotter``.

    Returns
    -------
    pv.Plotter
    """
    import pyvista as pv

    nodes = np.array(mesh3d.nodes)
    cells = []
    cell_types = []

    if mesh3d.numHexEls > 0:
        for el in mesh3d.hexElements:
            cells.extend([8] + list(el))
            cell_types.append(pv.CellType.HEXAHEDRON)

    if mesh3d.numWedgeEls > 0:
        for el in mesh3d.wedgeElements:
            cells.extend([6] + list(el))
            cell_types.append(pv.CellType.WEDGE)

    if mesh3d.numTetEls > 0:
        for el in mesh3d.tetElements:
            cells.extend([4] + list(el))
            cell_types.append(pv.CellType.TETRA)

    grid = pv.UnstructuredGrid(np.array(cells), np.array(cell_types), nodes)
    if scalars is not None:
        grid[scalar_name or "scalars"] = np.asarray(scalars)

    plotter = pv.Plotter(**plotter_kwargs)
    plotter.add_mesh(
        grid,
        scalars=scalar_name or ("scalars" if scalars is not None else None),
        show_edges=True,
        cmap="viridis",
    )
    plotter.add_axes()
    plotter.show_grid()
    plotter.show()
    return plotter


def plot_shell_mesh_pv(mesh_data, **plotter_kwargs):
    """Render a shell mesh dict (from ``get_shell_mesh``) in PyVista.

    Elements may be quads (4 unique node indices) or triangles (el[3] == -1).
    Node indices in ``mesh_data['elements']`` are 0-based.

    Parameters
    ----------
    mesh_data : dict
        Dict with keys ``'nodes'`` (N×3) and ``'elements'`` (M×4).
    **plotter_kwargs
        Forwarded to ``pv.Plotter``.

    Returns
    -------
    pv.Plotter
    """
    import pyvista as pv

    nodes = np.array(mesh_data["nodes"], dtype=float)
    faces = []
    for el in mesh_data["elements"]:
        if el[3] == -1:
            faces.extend([3, el[0], el[1], el[2]])
        else:
            faces.extend([4, el[0], el[1], el[2], el[3]])

    mesh = pv.PolyData(nodes, np.array(faces))
    plotter = pv.Plotter(**plotter_kwargs)
    plotter.add_mesh(mesh, color="steelblue", show_edges=True, edge_color="white", line_width=0.5)
    plotter.add_axes()
    plotter.show_grid()
    plotter.show()
    return plotter


def plot_component(component):
    """
    TODO docstring
    """
    fig, ax = plt.subplots()
    cpx, cpy = component.getcp()
    ax.plot(cpx, cpy)
    x = np.linspace(0, 1, 100)
    y = np.round(interpolator_wrap(cpx, cpy, x, "pchip", 0))
    ax.plot(x, y)
    plt.title(component.name)
    fig.show()
    return ax
