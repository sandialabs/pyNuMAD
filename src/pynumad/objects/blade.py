import re, warnings
import numpy as np
from copy import deepcopy

from pynumad.io.yaml_to_blade import yaml_to_blade
from pynumad.io.excel_to_blade import excel_to_blade
from pynumad.utils.interpolation import interpolator_wrap
from pynumad.utils.affinetrans import rotation, translation
from pynumad.objects.elements import MatDBentry, Layer, Shearweb, BOM, Ply, Stack
from pynumad.objects.geometry import Geometry
from pynumad.objects.settings import BladeSettings
from pynumad.objects.keypoints import KeyPoints
from pynumad.objects.definition import Definition

# for type hints
from numpy import ndarray


class Blade:
    """BladeDef A class definition for wind & water turbine blades.

    Parameters
    ----------
    filename : string

    Attributes
    ----------
    aerocenter : array
        Aerodynamic center of airfoil (used only by NuMAD->FAST)
    chord : array
        Chord distribution [m]
    chordoffset : array
        Chordwise offset (in addition to natural offset)
    components : list
        Blade components such as spar, panels, etc., refer to ``ComponentDef``
    degreestwist : array
        Twist distribution [degrees]
    ispan : array
        Spanwise locations of interpolated output
    leband : float
        Location of keypoint a
    materials : list
        Material properties, refer to ``MaterialDef``
    mesh : float
        Approximate element edge size for FE model [m]
    percentthick : array
        Percent thickness of airfoil [%]
    prebend : array
        Blade prebend, reference axis location along x2 [m]
    span : array
        Spanwise location of distributed properties [m]
    sparcapoffset : array
        (Does Nothing)
    sparcapwidth : array
        Locations of keypoints b & c, defines distance
        between keypoints b & c [mm]. First entry is the HP spar cap.
        Second entry is the LP spar cap
    stations : list
        Blade Stations, define the camber and thickness along the blade,
        refer to ``StationDef``
    sweep : array
        Blade Sweep, Reference axis location along x1 [m]
    teband : float
    geometry : array
        actual x,y,z geometry
    arclength : array
        surface distance from L.E.
    cpos : array
        chordwise position
    bom : dict
    bomIndices : dict
    stacks : array
        array of StackDef
    swstacks : list
        contains StackDefs for shearweb
    matdb : dict
        Composite definition for each region at each station
    te_type : list
        trailing edge type; assigned in update_keypoints
    shearweb : list
    bomPlot : dict
    hgGeometry : list
    hgKeypoints : list
    job_name : string
    paths : dict
    ansys : dict
        generate ANSYS settings
    write_airfoils : bool

    Example
    -------
    blade = BladeDef()
    """

    def __init__(self, filename: str = None):
        self.geometry: Geometry = None
        self.keypoints: KeyPoints = None
        self.definition: Definition = None
        self.ispan: ndarray = None
        self.bom: dict = None
        self.bomIndices: dict = None
        self.matdb: dict = None
        self.settings = BladeSettings()

        # read input file
        try:
            if "yaml" in filename or "yml" in filename:
                self.read_yaml(filename)
            elif "xls" in filename or "xlsx" in filename:
                self.read_excel(filename)
            else:
                raise Exception(
                    "Unknown filetype. Currently supported inputs are excel and yaml files."
                )

        # To handle when filename == None
        except TypeError:
            pass

        return

    ### Magic methods

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

    ### IO

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

    ### Update methods

    def update_blade(self):
        """
        TODO docstring
        """
        self.generate_geometry()
        self.generate_keypoints()
        self.generate_bom()
        return self

    def generate_geometry(self):
        """This method updates the interpolated blade parameters"""

        # update the interpolated station profiles
        definition = self.definition
        stations = definition.stations

        n_stations = len(stations)
        if n_stations > 0:
            n_points = len(stations[0].airfoil.c)
        else:
            raise Exception(
                "BladeDef must have at least one station before updating geometry."
            )

        # first station must be at blade root to prevent extrapolation
        assert stations[0].spanlocation == 0, "first station must be at the blade root"

        # initialize geometry object
        geometry = Geometry(shape=(n_points, n_stations))
        self.geometry = geometry
        geometry.ispan = self.ispan
        geometry._natural_offset = definition.natural_offset
        geometry._rotorspin = definition.rotorspin
        geometry._swtwisted = definition.swtwisted

        # Collect parameter tables from the stations.
        spanlocation = np.array(
            [stations[i].spanlocation for i in range(len(stations))]
        )
        tetype = ["_"] * n_stations
        for k in range(n_stations):
            station = stations[k]
            assert (
                len(station.airfoil.c) == n_points
            ), "Station airfoils must have same number of samples."

            geometry.c[:, k] = station.airfoil.c
            geometry.camber[:, k] = station.airfoil.camber
            geometry.thickness[:, k] = station.airfoil.thickness
            tetype[k] = station.airfoil.te_type

        # fix numerical issue due to precision on camber calculation
        # camber should start and end at y-values of zero
        geometry.camber[0, :] = np.zeros((1, n_stations))
        geometry.camber[-1, :] = np.zeros((1, n_stations))

        # Interpolate the station parameter tables.
        # Each column corresponds to an interpolated station.

        ## ic
        geometry.ic = interpolator_wrap(
            spanlocation, geometry.c, self.ispan, "pchip", axis=1
        )

        ## cpos
        geometry.cpos = np.concatenate(
            (
                -geometry.ic[-1, :].reshape(1, -1),
                -np.flipud(geometry.ic),
                geometry.ic[1:, :],
                geometry.ic[-1, :].reshape(1, -1),
            ),
            axis=0,
        )

        ## icamber
        geometry.icamber = interpolator_wrap(
            spanlocation, geometry.camber, self.ispan, "pchip", axis=1
        )

        ## ithickness
        geometry.ithickness = interpolator_wrap(
            spanlocation, geometry.thickness, self.ispan, "pchip", axis=1
        )
        # Adjust the thickness profiles based on te_type of stations.
        # This is mainly for transitions to flatbacks were the
        # interpolated airfoil needs to look like a round.
        for k in range(len(self.ispan)):
            try:
                ind = np.argwhere(self.ispan[k] < spanlocation)[0][0]
                # maybe better: ind = np.flatnonzero(self.ispan[k] < spanlocation)[0]
            except:
                continue
            else:
                if ind == 1:
                    continue
            if tetype[ind] == "flat" and tetype[ind - 1] == "round":
                geometry.ithickness[-1, k] = 0

        # Interpolate the blade parameter curves.
        ## idegreestwist
        geometry.idegreestwist = interpolator_wrap(
            definition.span, definition.degreestwist, self.ispan, "pchip"
        )

        ## ichord
        geometry.ichord = interpolator_wrap(
            definition.span, definition.chord, self.ispan, "pchip"
        )

        ## ipercentthick
        absolutethick = np.multiply(definition.percentthick, definition.chord) / 100
        iabsolutethick = interpolator_wrap(
            definition.span, absolutethick, self.ispan, "pchip"
        )
        geometry.ipercentthick = iabsolutethick / geometry.ichord * 100
        # ensure that the interpolation doesn't reduce the percent
        # thickness beneath the thinnest airfoil
        geometry.ipercentthick[
            geometry.ipercentthick < np.amin(definition.percentthick)
        ] = np.amin(definition.percentthick)

        ## ichordoffset
        geometry.ichordoffset = interpolator_wrap(
            definition.span, definition.chordoffset, self.ispan, "pchip"
        )

        ## iaerocenter
        geometry.iaerocenter = interpolator_wrap(
            definition.span, definition.aerocenter, self.ispan, "pchip"
        )

        ## isweep
        if len(definition.sweep) == 0:
            definition.sweep = np.zeros((self.span.shape, self.span.shape))
        if len(definition.prebend) == 0:
            definition.prebend = np.zeros((self.span.shape, self.span.shape))

        geometry.isweep = interpolator_wrap(
            definition.span, definition.sweep, self.ispan, "pchip"
        )

        ## iprebend
        geometry.iprebend = interpolator_wrap(
            definition.span, definition.prebend, self.ispan, "pchip"
        )

        # Generate the blade surface geometry.
        N = np.asarray(self.ispan).size
        M = n_points * 2 + 1
        geometry.profiles = np.zeros((M, 2, N))
        geometry.coordinates = np.zeros((M, 3, N))
        geometry.xoffset = np.zeros((1, N))
        geometry.LEindex = n_points

        for k in range(N):
            self.geometry.update_airfoil_profile(k)
            mtindex = np.argmax(geometry.ithickness[:, k])
            geometry.xoffset[0, k] = geometry.ic[mtindex, k]
            geometry.update_oml_geometry(k)

        # Calculate the arc length of each curve
        geometry.arclength = np.zeros((M, N))
        geometry.HParcx0 = np.zeros((1, N))
        geometry.LParcx0 = np.zeros((1, N))

        LE = geometry.LEindex
        for k in range(N):
            xx = geometry.coordinates[:, 0, k]
            yy = geometry.coordinates[:, 1, k]
            zz = geometry.coordinates[:, 2, k]
            arclen = np.sqrt(np.diff(xx) ** 2 + np.diff(yy) ** 2 + np.diff(zz) ** 2)
            arclen = np.concatenate((np.array([0]), np.cumsum(arclen)), axis=0)
            geometry.arclength[:, k] = arclen
            LEarcsum = geometry.arclength[LE, k]
            geometry.arclength[:, k] = geometry.arclength[:, k] - LEarcsum

            # find where x=0 intersects the surface
            geometry.HParcx0[0, k] = (
                interpolator_wrap(xx[1 : LE + 1], arclen[1 : LE + 1], 0) - LEarcsum
            )
            geometry.LParcx0[0, k] = (
                interpolator_wrap(xx[-2 : LE - 1 : -1], arclen[-2 : LE - 1 : -1], 0)
                - LEarcsum
            )

        return self

    def generate_keypoints(self):
        """This method updates the keypoints (a,b,c,...) which define the blade
        regions.

        Returns
        -------

        self

        Example:
          ``blade.update_keypoints``

        find the curves which bound each blade region
        """
        geometry = self.geometry
        definition = self.definition
        mm_to_m = 0.001

        # number of interpolated span stations
        num_istations = self.ispan.size

        # number of areas around airfoil profile; must be even (see calc of web areas)
        num_areas = 12

        # initialize keypoints
        keypoints = KeyPoints(shape=(num_areas, num_istations))
        self.keypoints = keypoints

        # start and finish indices in geometry/arcs
        ns = 1
        nf = geometry.coordinates.shape[0] - 2

        # keypoints, keyarcs, keycpos
        self.te_type = []  # reset te_type
        for k in range(num_istations):
            # allow for separate definitions of HP and LP spar cap
            # width and offset [HP LP]
            n1 = mm_to_m * definition.leband[k]  # no foam width
            n2 = mm_to_m * definition.teband[k]  # no foam width

            scwidth_hp = mm_to_m * definition.sparcapwidth_hp[k]  # type: float
            scwidth_lp = mm_to_m * definition.sparcapwidth_lp[k]  # type: float

            scoffset_hp = mm_to_m * definition.sparcapoffset_hp[k]  # type: float
            scoffset_lp = mm_to_m * definition.sparcapoffset_lp[k]  # type: float

            tempTE = geometry.get_profile_te_type(k)
            if definition.te_type:
                definition.te_type.append(tempTE)
            else:
                definition.te_type = []
                definition.te_type.append(tempTE)
            if definition.swtwisted:
                # get angle of each xy pair w.r.t. pitch axis (0,0)
                xyangle = np.zeros(geometry.coordinates.shape[0])
                for j in range(len(xyangle)):
                    xy = geometry.coordinates[j, 0:2, k]
                    xyangle[j] = np.arctan2(definition.rotorspin * xy[1], xy[0])
                # unwrap and center around 0
                xyangle = np.unwrap(xyangle)
                xyangle = xyangle - np.pi * np.round(xyangle[self.LEindex] / np.pi)

            k_arclen = self.geometry.arclength[ns : nf + 1, k]
            k_geom = self.geometry.coordinates[ns : nf + 1, :, k]
            k_cpos = self.geometry.cpos[ns : nf + 1, k]

            # ==================== HP surface ====================
            if self.definition.swtwisted:
                # find arclength where xyangle equals normal to chord
                # angle normal to chord line
                twistnorm = np.pi / 180 * (-self.idegreestwist[k] - 90)
                z = interpolator_wrap(xyangle[ns : nf + 1], k_arclen, twistnorm)
            else:
                z = geometry.HParcx0[0, k]
            z0 = z
            z = z - scoffset_hp
            a = np.amax(((0 - n1), 0.1 * geometry.arclength[ns, k]))  # type: float
            a = np.amin((a, 0.01 * geometry.arclength[ns, k]))
            b = np.amin(((z + 0.5 * scwidth_hp), 0.15 * geometry.arclength[ns, k]))
            c = np.amax(((z - 0.5 * scwidth_hp), 0.8 * geometry.arclength[ns, k]))
            d = np.amin(
                ((geometry.arclength[0, k] + n2), 0.85 * geometry.arclength[ns, k])
            )
            d = np.amax((d, 0.98 * geometry.arclength[ns, k]))
            if str(definition.te_type[k]) == "flat":
                e = geometry.arclength[ns, k]
                keypoints.key_points[0, :, k] = geometry.coordinates[ns, :, k]
                keypoints.key_cpos[1, k] = -1
            else:
                # e = 0.5 * (d + geometry.arclength(ns,k));
                e = 0.99 * geometry.arclength[ns, k]
                keypoints.key_points[0, :, k] = interpolator_wrap(k_arclen, k_geom, e)
                keypoints.key_cpos[1, k] = interpolator_wrap(k_arclen, k_cpos, e)

            # 1 -> e
            keypoints.key_points[1, :, k] = interpolator_wrap(k_arclen, k_geom, d)
            keypoints.key_points[2, :, k] = interpolator_wrap(k_arclen, k_geom, c)
            # keypoints.key_points(  ,:,k) = interpolator_wrap(geometry.arclength(ns:nf,k),self.geometry(ns:nf,:,k),z);
            keypoints.key_points[3, :, k] = interpolator_wrap(k_arclen, k_geom, b)
            keypoints.key_points[4, :, k] = interpolator_wrap(k_arclen, k_geom, a)
            keypoints.key_arcs[0, k] = geometry.arclength[ns, k]
            keypoints.key_arcs[1, k] = e
            keypoints.key_arcs[2, k] = d
            keypoints.key_arcs[3, k] = c
            # keypoints.key_arcs(  ,k)   = z;
            keypoints.key_arcs[4, k] = b
            keypoints.key_arcs[5, k] = a
            keypoints.key_arcs[6, k] = 0  # le
            keypoints.key_cpos[0, k] = geometry.cpos[ns, k]  # te, hp surface
            #            2   -> e
            keypoints.key_cpos[2, k] = interpolator_wrap(k_arclen, k_cpos, d)
            keypoints.key_cpos[3, k] = interpolator_wrap(k_arclen, k_cpos, c)
            #                 keypoints.key_cpos(  ,k) = interpolator_wrap(geometry.arclength(ns:nf,k),self.cpos(ns:nf,k),z);
            keypoints.key_cpos[4, k] = interpolator_wrap(k_arclen, k_cpos, b)
            keypoints.key_cpos[5, k] = interpolator_wrap(k_arclen, k_cpos, a)
            keypoints.key_cpos[6, k] = interpolator_wrap(k_arclen, k_cpos, 0)

            # ==================== LP surface ====================
            if self.definition.swtwisted:
                # angle normal to chord line
                twistnorm = np.pi / 180 * (-self.idegreestwist[k] + 90)
                z = interpolator_wrap(xyangle[ns : nf + 1], k_arclen, twistnorm)
            else:
                z = geometry.LParcx0[0, k]
            z0 = z  # ble: location where airfoil surface crosses Xglobal=0
            z = z + scoffset_lp  # positive scoffset moves z toward t.e.
            a = np.amin(((0 + n1), 0.1 * geometry.arclength[nf, k]))
            a = np.amax((a, 0.01 * geometry.arclength[nf, k]))
            b = np.amax(((z - 0.5 * scwidth_lp), 0.15 * geometry.arclength[nf, k]))
            c = np.amin((z + 0.5 * scwidth_lp, 0.8 * geometry.arclength[nf, k]))
            d = np.amax(
                (geometry.arclength[-1, k] - n2, 0.85 * geometry.arclength[nf, k])
            )
            d = np.amin((d, 0.96 * geometry.arclength[nf, k]))
            if str(definition.te_type[k]) == str("flat"):
                e = geometry.arclength[nf, k]
                keypoints.key_points[9, :, k] = geometry.coordinates[nf, :, k]
                keypoints.key_cpos[11, k] = 1
            else:
                # e = 0.5 * (d + geometry.arclength(nf,k));
                e = 0.98 * geometry.arclength[nf, k]
                keypoints.key_points[9, :, k] = interpolator_wrap(k_arclen, k_geom, e)
                keypoints.key_cpos[11, k] = interpolator_wrap(k_arclen, k_cpos, e)
            keypoints.key_points[5, :, k] = interpolator_wrap(k_arclen, k_geom, a)
            keypoints.key_points[6, :, k] = interpolator_wrap(k_arclen, k_geom, b)
            # keypoints.key_points(  ,:,k) = interpolator_wrap(geometry.arclength(ns:nf,k),self.geometry(ns:nf,:,k),z);
            keypoints.key_points[7, :, k] = interpolator_wrap(k_arclen, k_geom, c)
            keypoints.key_points[8, :, k] = interpolator_wrap(k_arclen, k_geom, d)
            # 10   -> e
            keypoints.key_arcs[7, k] = a
            keypoints.key_arcs[8, k] = b
            # keypoints.key_arcs( ,k)   = z;
            keypoints.key_arcs[9, k] = c
            keypoints.key_arcs[10, k] = d
            keypoints.key_arcs[11, k] = e
            keypoints.key_arcs[12, k] = geometry.arclength[nf, k]
            keypoints.key_cpos[7, k] = interpolator_wrap(k_arclen, k_cpos, a)
            keypoints.key_cpos[8, k] = interpolator_wrap(k_arclen, k_cpos, b)
            # keypoints.key_cpos(  ,k) = interpolator_wrap(geometry.arclength(ns:nf,k),self.cpos(ns:nf,k),z);
            keypoints.key_cpos[9, k] = interpolator_wrap(k_arclen, k_cpos, c)
            keypoints.key_cpos[10, k] = interpolator_wrap(k_arclen, k_cpos, d)
            # 12   -> e
            keypoints.key_cpos[12, k] = geometry.cpos[nf, k]  # te, lp surface

        # find the points used by each shear web
        component_groups = [
            definition.components[name].group for name in definition.components
        ]
        keypoints.web_indices = []
        keypoints.web_arcs = []
        keypoints.web_cpos = []
        keypoints.web_points = []
        keypoints.web_areas = []
        keypoints.web_width = []
        keypoints.web_bonds = []
        for ksw in range(max(component_groups)):  # for each shear web
            # pre-allocating arrays
            keypoints.web_indices.append([])
            keypoints.web_arcs.append(np.ndarray((2, num_istations)))
            keypoints.web_cpos.append(np.ndarray((2, num_istations)))
            keypoints.web_points.append(np.ndarray((2, 3, num_istations)))
            keypoints.web_areas.append(np.ndarray((num_istations - 1)))
            keypoints.web_width.append(np.ndarray(num_istations))
            keypoints.web_bonds.append(np.ndarray((2, num_istations - 1)))

            # find the components that are part of the shear web
            ksw_cmpts = [
                definition.components[comp]
                for comp in definition.components
                if definition.components[comp].group == ksw + 1
            ]

            # get hp extents
            hpextents = np.unique([comp.hpextents for comp in ksw_cmpts]).tolist()

            # get the lp extents
            lpextents = np.unique([comp.lpextents for comp in ksw_cmpts]).tolist()
            assert (
                len(hpextents) == 1
            ), f"HP Extents for components in group {ksw} must be identical and contain no spaces or commas"
            assert (
                len(lpextents) == 1
            ), f"LP Extents for components in group {ksw} must be identical and contain no spaces or commas"
            # match extents that have form of either '0.5b-c' or
            # 'b+/-100' or 'b' or 'z+/-100'
            # pat = '(?<fraction>\d*[\.]?\d*)(?<pt1>[a-zA-Z]+)-(?<pt2>[a-zA-Z]+)|(?<pt3>[a-zA-Z]+)(?<mm_offset>[+-]\d+)|(?<pt>[a-zA-Z])'
            pat = "(?P<fraction>\d*[\.]?\d*)(?P<pt1>[a-zA-Z]+)-(?P<pt2>[a-zA-Z]+)|(?P<pt3>[a-zA-Z]+)(?P<mm_offset>[+-]\d+)|(?P<pt>[a-zA-Z])"

            hp = re.search(pat, hpextents[0]).groupdict()
            lp = re.search(pat, lpextents[0]).groupdict()
            try:
                le = keypoints.key_labels.index("le")
            except ValueError:
                print(f"HP extent label \"{hp['pt']}\" not defined.")
            # get shear web placement on HP side
            if hp["pt"]:
                try:
                    n = keypoints.key_labels[0 : le + 1].index(hp["pt"])  ## EMA
                except ValueError:
                    print(f"HP extent label \"{hp['pt']}\" not defined.")
                keypoints.web_indices[ksw].append(n)
                keypoints.web_arcs[ksw][0, :] = keypoints.key_arcs[n, :]
                keypoints.web_cpos[ksw][0, :] = keypoints.key_cpos[n, :]
                n = n - 1
                keypoints.web_points[ksw][0, :, :] = keypoints.key_points[n, :, :]
            elif hp["pt1"]:
                f = float(hp["fraction"])
                if f <= 0 or f >= 1:
                    raise Exception(
                        f"Component group {ksw}: HP extent fraction={f}, which is outside range (0..1)"
                    )
                try:
                    n1 = keypoints.key_labels[0 : le + 1].index(hp["pt1"])
                except:
                    print(f"HP extent label \"{hp['pt1']}\" not defined.")
                try:
                    n2 = keypoints.key_labels[0 : le + 1].index(hp["pt2"])
                except:
                    print(f"HP extent label \"{hp['pt2']}\" not defined.")
                self.web_indices[ksw].append(np.nan)
                p1 = keypoints.key_arcs[n1, :]
                p2 = keypoints.key_arcs[n2, :]
                p = (1 - f) * p1 + f * p2
                keypoints.web_arcs[ksw][0, :] = p
                for k in range(num_istations):
                    keypoints.web_cpos[ksw][0, k] = interpolator_wrap(
                        k_arclen, k_cpos, p[k]
                    )
                    keypoints.web_points[ksw][0, :, k] = interpolator_wrap(
                        k_arclen, k_geom, p[k]
                    )
            elif hp["pt3"]:
                try:
                    n3 = keypoints.key_labels[0 : le + 1].index(hp["pt3"])
                except:
                    print(f"HP extent label \"{hp['pt3']}\" not defined.")
                self.web_indices[ksw].append(np.nan)
                p3 = keypoints.key_cpos[n3, :]
                p = p3 - float(hp["mm_offset"]) / 1000
                iMax = keypoints.key_labels[0, le + 1].index("d")
                # NOTE potential for error here - array shapes TBD -kb
                pMax = np.multiply(
                    keypoints.key_cpos[iMax, :], np.transpose(self.ichord)
                )
                p[np.abs(p) > np.abs(pMax)] = pMax[np.abs(p) > np.abs(pMax)]
                iMin = keypoints.key_labels[0 : le + 1].index("a")
                # NOTE same issue here -kb
                pMin = np.multiply(
                    keypoints.key_cpos[iMin, :], np.transpose(self.ichord)
                )
                p[np.abs(p) < np.abs(pMin)] = pMin[np.abs(p) < np.abs(pMin)]
                keypoints.web_cpos[ksw][0, :] = p
                for k in range(num_istations):
                    keypoints.web_arcs[ksw][0, k] = interpolator_wrap(
                        self.cpos[ns : nf + 1, :, k],
                        geometry.arclength[ns : nf + 1, :, k],
                        p[k],
                    )
                    keypoints.web_points[ksw][0, :, k] = interpolator_wrap(
                        k_cpos, k_geom, p[k]
                    )
            else:
                raise Exception(
                    "Shear web geometry HP extents not defined correctly (e.g., 0.5b-c, b, b+200)"
                )
            # get shear web placement on LP side
            if lp["pt"]:
                try:
                    n = keypoints.key_labels[le:].index(lp["pt"]) + le
                    keypoints.web_indices[ksw].append(n)
                    keypoints.web_arcs[ksw][1, :] = keypoints.key_arcs[n, :]
                    keypoints.web_cpos[ksw][1, :] = keypoints.key_cpos[n, :]
                    keypoints.web_points[ksw][1, :, :] = keypoints.key_points[n, :, :]
                except:
                    print(f"LP extent label \"{lp['pt']}\" not defined.")

            elif lp["pt1"]:
                f = float(lp["fraction"])
                if f < 0 or f > 1:
                    raise Exception(
                        f"Component group {ksw}: LP extent fraction={f}, which is outside range [0..1]"
                    )
                try:
                    n1 = keypoints.key_labels[le:].index(lp["pt1"]) + le
                except:
                    print(f"LP extent label \"{lp['pt1']}\" not defined.")
                try:
                    n2 = keypoints.key_labels[le:].index(lp["pt2"]) + le
                except:
                    print(f"LP extent label \"{lp['pt2']}\" not defined.")
                self.web_indices[ksw].append(np.nan)
                p1 = keypoints.key_arcs[n1, :]
                p2 = keypoints.key_arcs[n2, :]
                p = (1 - f) * p1 + f * p2
                keypoints.web_arcs[ksw][1, :] = p
                for k in range(num_istations):
                    keypoints.web_cpos[ksw][1, k] = interpolator_wrap(
                        k_arclen, k_cpos, p[k]
                    )
                    keypoints.web_points[ksw][1, :, k] = interpolator_wrap(
                        k_arclen, k_geom, p[k]
                    )
            elif lp["pt3"]:
                try:
                    n3 = keypoints.key_labels[le:].index(lp["pt3"]) + le
                except:
                    print(f"LP extent label \"{lp['pt3']}\" not defined.")
                self.web_indices[ksw].append(np.nan)
                p3 = keypoints.key_cpos[n3, :]
                p = p3 + float(lp["mm_offset"]) / 1000
                iMax = keypoints.key_labels[le:].index("d") + le
                pMax = np.multiply(
                    keypoints.key_cpos[iMax, :], np.transpose(self.ichord)
                )
                p[np.abs(p) > np.abs(pMax)] = pMax[np.abs(p) > np.abs(pMax)]
                iMin = keypoints.key_labels[le:].index("a") + le
                pMin = np.multiply(
                    keypoints.key_cpos[iMin, :], np.transpose(self.ichord)
                )
                p[np.abs(p) < np.abs(pMin)] = pMin[np.abs(p) < np.abs(pMin)]
                keypoints.web_cpos[ksw][1, :] = p
                for k in range(num_istations):
                    keypoints.web_arcs[ksw][1, k] = interpolator_wrap(
                        k_cpos, k_arclen, p[k]
                    )
                    keypoints.web_points[ksw][1, :, k] = interpolator_wrap(
                        k_cpos, k_geom, p[k]
                    )
            else:
                raise Exception(
                    "Shear web geometry LP extents not defined correctly (e.g., 0.5b-c, b, b+200)"
                )

        # calculate shell areas
        for kc in range(num_istations - 1):
            for kr in range(num_areas):
                # choose number of points to use in area calculation
                # jcb: I decided to base this on the number of points
                # in the interpolated station profile found within the region
                #  of interest.
                npts = sum(
                    np.logical_and(
                        geometry.arclength[:, kc] >= keypoints.key_arcs[kr, kc],
                        geometry.arclength[:, kc] <= keypoints.key_arcs[kr + 1, kc],
                    )
                )

                # need at least two points
                npts = np.amax((npts, 2))

                # inboard curve arclengths
                ibarc = np.linspace(
                    keypoints.key_arcs[kr, kc], keypoints.key_arcs[kr + 1, kc], npts
                )

                # outboard curve arclengths
                obarc = np.linspace(
                    keypoints.key_arcs[kr, kc + 1],
                    keypoints.key_arcs[kr + 1, kc + 1],
                    npts,
                )

                # inboard xyz
                ib = interpolator_wrap(
                    geometry.arclength[ns : nf + 1, kc],
                    geometry.coordinates[ns : nf + 1, :, kc],
                    ibarc,
                )

                # outboard xyz
                ob = interpolator_wrap(
                    geometry.arclength[ns : nf + 1, kc + 1],
                    geometry.coordinates[ns : nf + 1, :, kc + 1],
                    obarc,
                )

                # "ds" in the span direction
                dspan = np.sqrt(np.sum((ob - ib) ** 2, 1))

                # treat each "rectangular" area as two triangles
                t1 = 0.5 * np.dot(
                    np.sqrt(np.sum(np.diff(ib, 1, axis=0) ** 2, 1)), dspan[0:-1]
                )
                t2 = 0.5 * np.dot(
                    np.sqrt(np.sum(np.diff(ob, 1, axis=0) ** 2, 1)), dspan[1:]
                )
                keypoints.key_areas[kr, kc] = t1 + t2
                if kr == 0:
                    keypoints.te_bond[kc] = dspan[0]
                if (num_areas / 2 + 1) == (kr + 1):
                    keypoints.le_bond[kc] = dspan[0]

        # calculate areas used by shear webs
        # jcb: note that these areas come purely from the geometry and
        # do not take into account the thickness of the shell or
        # sparcap layup.
        for ksw in range(len(keypoints.web_points)):
            for kc in range(num_istations - 1):
                ib = keypoints.web_points[ksw][:, :, kc]
                ob = keypoints.web_points[ksw][:, :, kc + 1]
                # treat each "rectangular" area as two triangles
                b1 = np.diff(ib, axis=0)
                b2 = np.diff(ob, axis=0)
                base1 = np.sqrt(np.sum(b1**2, 1))[0]
                base2 = np.sqrt(np.sum(b2**2, 1))[0]
                b1 = b1 / base1
                b2 = b2 / base2
                h1 = np.abs(np.dot((ob[0, :] - ib[0, :]), (1 - np.transpose(b1))))
                h2 = np.abs(np.dot((ib[1, :] - ob[1, :]), (1 - np.transpose(b2))))
                keypoints.web_areas[ksw][kc] = 0.5 * (base1 * h1 + base2 * h2)
                keypoints.web_width[ksw][kc] = base1
                # calculate edge (bond-line) lengths
                keypoints.web_bonds[ksw][0:2, kc] = np.sqrt(np.sum((ob - ib) ** 2, 1))
            keypoints.web_width[ksw][num_istations - 1] = base2

        return self

    def generate_bom(self):
        """This method generates the Bill-of-Materials
        See datatypes.BOM

        Returns
        -------

        None

        """
        # raise DeprecationWarning("update_bom currently deprecated. Please do not use.")
        # set conversion constants
        G_TO_KG = 0.001
        M_TO_MM = 1000.0
        MM_TO_M = 0.001

        keypoints = self.keypoints
        geometry = self.geometry
        definition = self.definition
        materials = definition.materials
        components = definition.components

        # initialize structures
        self.bom = {
            "hp": [],
            "lp": [],
            "sw": [],
            "lebond": [],
            "tebond": [],
            "swbonds": [],
            "dryweight": [],
        }
        self.bomIndices = {"hp": [], "lp": [], "sw": []}

        # calculate non-dimensional span
        ndspan = (self.ispan - self.ispan[0]) / (self.ispan[-1] - self.ispan[0])

        hprow = 0
        lprow = 0
        outer_shape_comps = [name for name in components if components[name].group == 0]
        for comp_name in outer_shape_comps:
            comp = components[comp_name]
            mat = materials[comp.materialid]
            hpRegion, lpRegion = comp.find_region_extents()
            num_layers = comp.get_num_layers(ndspan)
            num_layers = np.round(num_layers)

            for k_layer in range(1, int(np.max(num_layers)) + 1):
                begin_station, end_station = self.find_layer_extents(
                    num_layers, k_layer
                )
                ks_max = np.amin((len(begin_station), len(end_station)))
                # situation that beginSta/endSta is longer than 1
                for ks in range(ks_max):
                    ## END
                    if hpRegion:
                        areas = keypoints.key_areas[
                            hpRegion[0] : hpRegion[1],
                            begin_station[ks] : end_station[ks],
                        ]
                        regionarea = sum(areas.flatten())
                        arcs = (
                            keypoints.key_arcs[
                                hpRegion[1], begin_station[ks] : end_station[ks] + 1
                            ]
                            - keypoints.key_arcs[
                                hpRegion[0], begin_station[ks] : end_station[ks] + 1
                            ]
                        )
                        cur_bom = BOM()
                        cur_bom.layernum = hprow
                        cur_bom.materialid = comp.materialid
                        cur_bom.name = comp.name
                        cur_bom.beginsta = self.ispan[begin_station[ks]]
                        cur_bom.endsta = self.ispan[end_station[ks]]
                        cur_bom.maxwidth = np.amax(arcs)
                        cur_bom.avgwidth = np.mean(arcs)
                        cur_bom.area = regionarea
                        cur_bom.thickness = mat.layerthickness
                        cur_bom.weight = mat.drydensity * regionarea
                        self.bomIndices["hp"].append(
                            [begin_station[ks], end_station[ks], *hpRegion]
                        )
                        self.bom["hp"].append(cur_bom)
                        hprow = hprow + 1

                    if lpRegion:
                        areas = keypoints.key_areas[
                            lpRegion[0] : lpRegion[1],
                            begin_station[ks] : end_station[ks],
                        ]
                        regionarea = sum(areas.flatten())
                        arcs = (
                            keypoints.key_arcs[
                                lpRegion[1], begin_station[ks] : end_station[ks] + 1
                            ]
                            - keypoints.key_arcs[
                                lpRegion[0], begin_station[ks] : end_station[ks] + 1
                            ]
                        )
                        cur_bom = BOM()
                        cur_bom.layernum = lprow
                        cur_bom.materialid = comp.materialid
                        cur_bom.name = comp.name
                        cur_bom.beginsta = self.ispan[begin_station[ks]]
                        cur_bom.endsta = self.ispan[end_station[ks]]
                        cur_bom.maxwidth = np.amax(arcs)
                        cur_bom.avgwidth = np.mean(arcs)
                        cur_bom.area = regionarea
                        cur_bom.thickness = mat.layerthickness
                        cur_bom.weight = mat.drydensity * regionarea
                        self.bomIndices["lp"].append(
                            [begin_station[ks], end_station[ks], *lpRegion]
                        )
                        self.bom["lp"].append(cur_bom)
                        lprow = lprow + 1

        # shearwebs
        swnum = None
        swrow = 0
        sw_begin_station = []
        sw_end_station = []
        sw_comps = [comp for comp in components.values() if comp.group > 0]

        def sorter(e):
            return e.group

        # enforce an ordering on components based on group number
        # to ensure below loop functions correctly
        # probably should re-write the loop at some point...
        sw_comps.sort(key=sorter)
        for comp in sw_comps:
            mat = materials[comp.materialid]
            num_layers = comp.get_num_layers(ndspan)
            num_layers = np.round(num_layers)

            for k_layer in range(1, int(np.max(num_layers)) + 1):
                begin_station, end_station = self.find_layer_extents(
                    num_layers, k_layer
                )
                ks_max = np.amin((len(begin_station), len(end_station)))
                # situation that beginSta/endSta is longer than 1
                for ks in range(ks_max):
                    if swnum != comp.group - 1:
                        swnum = comp.group - 1
                        swrow = 0
                        sw_begin_station.append(begin_station[0])
                        sw_end_station.append(end_station[0])
                        self.bom["sw"].append([])
                        self.bomIndices["sw"].append([])
                    sw_begin_station[swnum] = np.amin(
                        [*begin_station, sw_begin_station[swnum]]
                    )
                    sw_end_station[swnum] = np.amax(
                        [*end_station, sw_end_station[swnum]]
                    )
                    areas = keypoints.web_areas[swnum][
                        begin_station[ks] : end_station[ks]
                    ]
                    regionarea = sum(areas.flatten())
                    cur_bom = BOM()
                    cur_bom.layernum = swrow
                    cur_bom.materialid = comp.materialid
                    cur_bom.name = comp.name
                    cur_bom.beginsta = self.ispan[begin_station[ks]]
                    cur_bom.endsta = self.ispan[end_station[ks]]
                    cur_bom.maxwidth = np.amax(keypoints.web_width[swnum])
                    cur_bom.avgwidth = np.mean(keypoints.web_width[swnum])
                    cur_bom.area = regionarea
                    cur_bom.thickness = mat.layerthickness
                    cur_bom.weight = mat.drydensity * regionarea
                    self.bom["sw"][swnum].append(cur_bom)
                    self.bomIndices["sw"][swnum].append(
                        [begin_station[ks], end_station[ks]]
                    )
                    swrow = swrow + 1

        # compute lebond, tebond, and dryweight
        self.bom["lebond"] = sum(keypoints.le_bond) * M_TO_MM
        self.bom["tebond"] = sum(keypoints.te_bond) * M_TO_MM
        hp_dw = sum([L.weight for L in self.bom["hp"]])
        lp_dw = sum([L.weight for L in self.bom["lp"]])
        self.bom["dryweight"] = G_TO_KG * (hp_dw + lp_dw)

        nsw = len(self.bom["sw"])
        self.bom["swbonds"] = [None] * nsw
        for k in range(nsw):
            sw_dw = sum([L.weight for L in self.bom["sw"][k]])
            self.bom["dryweight"] = self.bom["dryweight"] + sw_dw
            C = keypoints.web_bonds[k][:, sw_begin_station[k] : sw_end_station[k]]
            self.bom["swbonds"][k] = M_TO_MM * np.sum(C, 1)

        self._generate_stacks()
        self._generate_matdb()
        return self

    def _generate_stacks(self):
        """_summary_"""
        keypoints = self.keypoints
        # build the material stack for each area
        n_segments = keypoints.key_areas.shape[0]
        n_stations = keypoints.key_areas.shape[1]
        n_webs = len(self.bomIndices["sw"])
        segment_labels = [
            "HP_TE_FLAT",
            "HP_TE_REINF",
            "HP_TE_PANEL",
            "HP_SPAR",
            "HP_LE_PANEL",
            "HP_LE",
            "LP_LE",
            "LP_LE_PANEL",
            "LP_SPAR",
            "LP_TE_PANEL",
            "LP_TE_REINF",
            "LP_TE_FLAT",
        ]

        # initialize stack array
        self.stacks = np.empty(shape=(n_segments, n_stations), dtype=object)
        for k_seg in range(n_segments):
            for k_stat in range(n_stations):
                self.stacks[k_seg, k_stat] = Stack()
                self.stacks[k_seg, k_stat].name = "{:02d}_{:02d}_{}".format(
                    k_seg, k_stat, segment_labels[k_seg]
                )
                self.stacks[k_seg, k_stat].indices = [
                    k_stat,
                    k_stat + 1,
                    k_seg,
                    k_seg + 1,
                ]

        for k in range(len(self.bom["hp"])):
            # for each row in the BOM, get the ply definition ...
            cur_ply = Ply()
            cur_ply.component = self.bom["hp"][k].name
            cur_ply.materialid = self.bom["hp"][k].materialid
            cur_ply.thickness = self.bom["hp"][k].thickness
            cur_ply.angle = 0  # TODO, set to 0 for now, self.bom['lp'](k, );
            cur_ply.nPlies = 1  # default to 1, modified in addply() if necessary

            # ... and add the ply to every area that is part of the region
            ind = self.bomIndices["hp"][k]
            for k_seg in range(ind[2], ind[3]):
                for k_stat in range(ind[0], ind[1]):
                    # deepcopy is important to keep make ply object unique in each stack
                    self.stacks[k_seg, k_stat].addply(deepcopy(cur_ply))  

        for k in range(len(self.bom["lp"])):
            # for each row in the BOM, get the ply definition ...
            cur_ply = Ply()
            cur_ply.component = self.bom["lp"][k].name
            cur_ply.materialid = self.bom["lp"][k].materialid
            cur_ply.thickness = self.bom["lp"][k].thickness
            cur_ply.angle = 0  # TODO, set to 0 for now, self.bom['lp'](k, );
            cur_ply.nPlies = 1  # default to 1, modified in addply() if necessary

            # ... and add the ply to every area that is part of the region
            ind = self.bomIndices["lp"][k]
            for k_seg in range(ind[2], ind[3]):
                for k_stat in range(ind[0], ind[1]):
                    self.stacks[k_seg, k_stat].addply(deepcopy(cur_ply))

        self.swstacks = np.empty(shape=(n_webs, n_stations), dtype=object)
        for k_web in range(n_webs):
            for k_stat in range(n_stations):
                self.swstacks[k_web, k_stat] = Stack()
                # name the stacks <webnumber+1>_<stationnumber+1>_SW
                self.swstacks[k_web, k_stat].name = "{:02d}_{:02d}_SW".format(
                    k_web, k_stat
                )
                # currently, the shearweb indices do not change down the span
                ind = keypoints.web_indices[k_web]
                self.swstacks[k_web][k_stat].indices = [
                    k_stat,
                    k_stat + 1,
                    ind[0],
                    ind[1],
                ]
        for k_web in range(n_webs):
            for k in range(len(self.bom["sw"][k_web])):
                # for each row in the BOM, get the ply definition ...
                cur_ply = Ply()
                cur_ply.component = self.bom["sw"][k_web][k].name
                cur_ply.materialid = self.bom["sw"][k_web][k].materialid
                cur_ply.thickness = self.bom["sw"][k_web][k].thickness
                cur_ply.angle = 0  # TODO, set to 0 for now, self.bom['lp'](k, );
                cur_ply.nPlies = 1  # default to 1, modified in addply() if necessary
                
                # ... and add the ply to every area that is part of the region
                ind = self.bomIndices["sw"][k_web][k]
                for k_stat in range(ind[0], ind[1]):
                    self.swstacks[k_web, k_stat].addply(deepcopy(cur_ply))
        

    def _generate_matdb(self):
        """Adds material and composites information to MatDB
        """
        MM_TO_M = 0.001
        materials = self.definition.materials
        n_segments = self.keypoints.key_areas.shape[0]
        n_stations = self.keypoints.key_areas.shape[1]
        n_webs = len(self.bomIndices["sw"])

        # prepare material database ==========================================
        self.matdb = dict()
        
        # add base materials
        for mat_name in materials:
            cur_entry = MatDBentry()
            cur_material = materials[mat_name]
            cur_entry.name = cur_material.name
            cur_entry.type = cur_material.type
            cur_entry.ex = cur_material.ex
            cur_entry.ey = cur_material.ey
            cur_entry.ez = cur_material.ez
            cur_entry.gxy = cur_material.gxy
            cur_entry.gyz = cur_material.gyz
            cur_entry.gxz = cur_material.gxz
            if cur_entry.type == "isotropic":
                cur_entry.nuxy = cur_material.prxy
            else:
                cur_entry.prxy = cur_material.prxy
                cur_entry.pryz = cur_material.pryz
                cur_entry.prxz = cur_material.prxz
            cur_entry.dens = cur_material.density
            cur_entry.reference = cur_material.reference
            self.matdb[mat_name] = cur_entry

        # add component stacks
        flat_stacks = self.stacks.flatten("F")
        for k in range(self.stacks.size):
            cur_entry = MatDBentry()
            cur_entry.name = flat_stacks[k].name
            cur_entry.type = "composite"
            cur_entry.reference = "Reference text"
            cur_entry.thicknessType = "Constant"
            cur_entry.uniqueLayers = len(flat_stacks[k].plygroups)
            cur_entry.symmetryType = "none"
            cur_entry.layer = [None] * cur_entry.uniqueLayers
            for j in range(cur_entry.uniqueLayers):
                cur_layer = Layer()
                matid = flat_stacks[k].plygroups[j].materialid
                cur_layer.layerName = self.matdb[matid].name
                cur_layer.thicknessA = MM_TO_M * flat_stacks[k].plygroups[j].thickness
                cur_layer.thicknessB = cur_layer.thicknessA
                cur_layer.quantity = flat_stacks[k].plygroups[j].nPlies
                cur_layer.theta = flat_stacks[k].plygroups[j].angle
                cur_entry.layer[j] = cur_layer
            self.matdb[cur_entry.name] = cur_entry

        # add shearweb stacks
        for k_web in range(n_webs):
            for k_stat in range(n_stations):
                cur_entry = MatDBentry()
                cur_entry.name = self.swstacks[k_web, k_stat].name
                cur_entry.type = "composite"
                cur_entry.reference = "Reference text"
                cur_entry.thicknessType = "Constant"
                try:
                    cur_entry.uniqueLayers = len(self.swstacks[k_web, k_stat].plygroups)
                except TypeError:
                    cur_entry.uniqueLayers = 0
                cur_entry.symmetryType = "none"
                cur_entry.layer = [None] * cur_entry.uniqueLayers
                for j in range(cur_entry.uniqueLayers):
                    cur_layer = Layer()
                    matid = self.swstacks[k_web, k_stat].plygroups[j].materialid
                    cur_layer.layerName = self.matdb[matid].name
                    cur_layer.thicknessA = (
                        MM_TO_M * self.swstacks[k_web, k_stat].plygroups[j].thickness
                    )
                    cur_layer.thicknessB = cur_layer.thicknessA
                    cur_layer.quantity = self.swstacks[k_web, k_stat].plygroups[j].nPlies
                    cur_layer.theta = self.swstacks[k_web, k_stat].plygroups[j].angle
                    cur_entry.layer[j] = cur_layer
                self.matdb[cur_entry.name] = cur_entry
        
        # shearweb information from NuMAD v1 is formatted in a specific
        # way, recreating that here
        # recreating data.shearweb ====================================
        # NOTE: do we care about this anymore? -kb
        ctr = 0
        self.shearweb = []
        for k_web in range(n_webs):
            ind = self.keypoints.web_indices[k_web]
            for k_stat in range(n_stations):
                if self.swstacks[k_web, k_stat].plygroups:
                    cur_sw = Shearweb()
                    cur_sw.Material = self.swstacks[k_web, k_stat].name
                    cur_sw.BeginStation = self.swstacks[k_web, k_stat].indices[0]  # =k
                    cur_sw.EndStation = self.swstacks[k_web, k_stat].indices[1]  # =k+1
                    cur_sw.Corner = [
                        ind[1] - 1,
                        ind[0] - 1,
                        ind[0] - 1,
                        ind[1] - 1,
                    ]  # dp number is offset by 1 in NuMAD v1
                    self.shearweb.append(cur_sw)
                    ctr += 1
        return self
        
    # Supporting function for update_bom
    def find_layer_extents(self, layer_dist, layer_n):
        """
        TODO docstring
        """
        assert np.isscalar(layer_n), 'second argument "layer_n" must be a scalar'
        sta_logical = layer_dist >= layer_n
        prev = 0
        begin_station = []
        end_station = []
        for k in range(len(sta_logical)):
            if sta_logical[k] == 1 and prev == 0:
                begin_station.append(k)
            if sta_logical[k] == 0 and prev == 1:
                end_station.append(k)
            elif k == len(sta_logical) - 1 and prev == 1:
                end_station.append(k)
            prev = sta_logical[k]

        return begin_station, end_station

    def expand_blade_geometry_te(self, min_edge_lengths):
        """
        TODO: docstring
        """
        self.geometry.expand_blade_geometry_te(min_edge_lengths)
        self.generate_keypoints()
        return

    ### Shell

    def edit_stacks_for_solid_mesh(self):
        """_summary_

        Returns
        -------
        _type_
            _description_
        """
        numSec, numStat = self.stacks.shape
        for i in range(numSec):
            for j in range(numStat):
                pg = self.stacks[i, j].plygroups
                if len(pg) == 4:
                    ply1 = deepcopy(pg[1])
                    ply2 = deepcopy(pg[2])
                    ply3 = deepcopy(pg[3])
                    newPg = np.array([ply1, ply2, ply3])
                else:
                    if len(pg) == 3:
                        # newPg = np.array([pg[1],pg[1],pg[2]])
                        ply1 = deepcopy(pg[1])
                        ply2 = deepcopy(pg[1])
                        ply3 = deepcopy(pg[2])
                        t2 = ply1.thickness
                        t3 = ply3.thickness
                        ply2.thickness = 0.3333333 * (t2 + t3)
                        ply1.thickness = 0.6666666 * t2
                        ply3.thickness = 0.6666666 * t3
                        newPg = np.array([ply1, ply2, ply3])
                    else:
                        if len(pg) == 2:
                            ply1 = deepcopy(pg[0])
                            ply2 = deepcopy(pg[0])
                            ply3 = deepcopy(pg[1])
                            # newPg = np.array([pg[0],pg[0],pg[1]])
                            t1 = ply1.thickness
                            t2 = ply3.thickness
                            ply2.thickness = 0.3333333 * (t1 + t2)
                            ply1.thickness = 0.6666666 * t1
                            ply3.thickness = 0.6666666 * t2
                            newPg = np.array([ply1, ply2, ply3])
                        else:
                            ply1 = deepcopy(pg[0])
                            ply2 = deepcopy(pg[0])
                            ply3 = deepcopy(pg[0])
                            # newPg = np.array([pg[0],pg[0],pg[0]])
                            t1 = ply1.thickness
                            ply2.thickness = 0.3333333 * t1
                            ply1.thickness = 0.3333333 * t1
                            ply3.thickness = 0.3333333 * t1
                            newPg = np.array([ply1, ply2, ply3])
                self.stacks[i, j].plygroups = newPg

        for i in range(2):
            stackLst = self.swstacks[i]
            for j in range(len(stackLst)):
                pg = stackLst[j].plygroups
                if len(pg) == 2:
                    ply1 = deepcopy(pg[0])
                    ply2 = deepcopy(pg[0])
                    ply3 = deepcopy(pg[1])
                    # newPg = np.array([pg[0],pg[0],pg[1]])
                    t1 = ply1.thickness
                    t2 = ply3.thickness
                    ply2.thickness = 0.3333333 * (t1 + t2)
                    ply1.thickness = 0.6666666 * t1
                    ply3.thickness = 0.6666666 * t2
                    newPg = np.array([ply1, ply2, ply3])
                    self.swstacks[i][j].plygroups = newPg
                elif len(pg) == 1:
                    ply1 = deepcopy(pg[0])
                    ply2 = deepcopy(pg[0])
                    ply3 = deepcopy(pg[0])
                    # newPg = np.array([pg[0],pg[0],pg[0]])
                    t1 = ply1.thickness
                    ply2.thickness = 0.3333333 * t1
                    ply1.thickness = 0.3333333 * t1
                    ply3.thickness = 0.3333333 * t1
                    newPg = np.array([ply1, ply2, ply3])
                    self.swstacks[i][j].plygroups = newPg
        return self

    ### Blade modification API

    def add_interpolated_station(self, span_location):
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
