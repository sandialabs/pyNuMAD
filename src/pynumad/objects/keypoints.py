import re

import numpy as np
from numpy import ndarray

from pynumad.utils.interpolation import interpolator_wrap
from pynumad.objects.definition import Definition
from pynumad.objects.geometry import Geometry



class KeyPoints:
    """Keypoints class

    Attributes
    ----------
    key_labels: list
        list of key label names
    key_points: array
        keypoints in xyz geometry
    key_arcs: array
        surface arclength distance of keypoints from LE
    key_chordpos: array
        chordwise position of keypoints
    key_areas: array
        surface area of regions created by keypoints
    le_bond: array
    te_bond: array

    web_indices: list
    web_points: list
    web_arcs: list
    web_cpos: list
    web_areas: list
    web_width: list
    web_bonds: list

    """

    def __init__(self):
        self.key_labels = [
            "te",
            "e",
            "d",
            "c",
            "b",
            "a",
            "le",
            "a",
            "b",
            "c",
            "d",
            "e",
            "te",
        ]
        self.key_points: ndarray = None
        self.key_arcs: ndarray = None
        self.key_cpos: ndarray = None
        self.key_areas: ndarray = None

        self.le_bond: ndarray = None
        self.te_bond: ndarray = None

        self.web_indices: list = None
        self.web_points: list = None
        self.web_arcs: list = None
        self.web_cpos: list = None
        self.web_areas: list = None
        self.web_width: list = None
        self.web_bonds: list = None
        
        
    def __eq__(self, other):
        attrs = [
            a
            for a in dir(self)
            if not a.startswith("__") and not callable(getattr(self, a))
        ]
        for attr in attrs:
            self_attr = getattr(self, attr)
            other_attr = getattr(other, attr)
            if isinstance(self_attr, (int,)):
                if self_attr != other_attr:
                    return False
            elif isinstance(self_attr, ndarray):
                if (self_attr != other_attr).any():
                    return False
        return True
    
        
    def initialize(self,num_areas,num_stations):
        self.key_points = np.zeros((num_areas - 2, 3, num_stations))
        self.key_arcs = np.zeros((num_areas + 1, num_stations))
        self.key_cpos = np.zeros((num_areas + 1, num_stations))
        self.key_areas = np.zeros((num_areas, num_stations - 1))
        self.le_bond = np.zeros((num_stations - 1))
        self.te_bond = np.zeros((num_stations - 1))
        return self
    
    def generate(self, definition: Definition, geometry: Geometry):
        """This method generates the keypoints (a,b,c,...) which define the blade
        regions given by a definition and geometry object.
        
        Parameters
        ----------
        defintion : Definition
        geometry : Geometry

        Returns
        -------
        Self
        """
        mm_to_m = 0.001

        # number of interpolated span stations
        num_istations = geometry.ispan.size

        # number of areas around airfoil profile; must be even (see calc of web areas)
        num_areas = 12
        self.initialize(num_areas, num_istations)
        # initialize keypoints

        # start and finish indices in geometry/arcs
        ns = 1
        nf = geometry.coordinates.shape[0] - 2

        # keypoints, keyarcs, keycpos
        te_types = []  # reset te_type
        i_leband_start=np.min(np.nonzero(definition.leband))
        i_teband_start=np.min(np.nonzero(definition.teband))

        i_leband_end=np.max(np.nonzero(definition.leband))
        i_teband_end=np.max(np.nonzero(definition.teband))

        for k in range(num_istations):
            # allow for separate definitions of HP and LP spar cap
            # width and offset [HP LP]
            n1 = mm_to_m * definition.leband[k]  # no foam width
            n2 = mm_to_m * definition.teband[k]  # no foam width

            ### 
            #In order to avoid abrupt changes in geometry when the le/te bands
            #begin, set the le/te band width equal to the first nonzero value.
            #This algorithm will not work as well if small numbers exist in the
            #le/te band widths since it is based on nonzero values.
            
            
            if k < i_leband_start :
                n1= mm_to_m * definition.leband[i_leband_start]

            if k < i_teband_start :
                n2= mm_to_m * definition.teband[i_teband_start]

            #In order to avoid abrupt changes in geometry when the le/te bands
            #end, set the le/te band width is tapered starting with the last nonzero.
            #value. This algorithm will not work as well if small numbers exist in the
            #le/te band widths since it is based on nonzero values.
            
            if k > i_leband_end :
                n1= mm_to_m * definition.leband[k-1]*0.75
                definition.leband[k]=n1/mm_to_m

            if k > i_teband_end :
                n2= mm_to_m * definition.teband[k-1]*0.75
                definition.teband[k]=n2/mm_to_m
            ###

            scwidth_hp = mm_to_m * definition.sparcapwidth_hp[k]  # type: float
            scwidth_lp = mm_to_m * definition.sparcapwidth_lp[k]  # type: float

            scoffset_hp = mm_to_m * definition.sparcapoffset_hp[k]  # type: float
            scoffset_lp = mm_to_m * definition.sparcapoffset_lp[k]  # type: float

            tempTE = geometry.get_profile_te_type(k)
            if te_types:
                te_types.append(tempTE)
            else:
                te_types = []
                te_types.append(tempTE)
            if definition.swtwisted:
                # get angle of each xy pair w.r.t. pitch axis (0,0)
                xyangle = np.zeros(geometry.coordinates.shape[0])
                for j in range(len(xyangle)):
                    xy = geometry.coordinates[j, 0:2, k]
                    xyangle[j] = np.arctan2(definition.rotorspin * xy[1], xy[0])
                # unwrap and center around 0
                xyangle = np.unwrap(xyangle)
                xyangle = xyangle - np.pi * np.round(xyangle[self.LEindex] / np.pi)

            k_arclen = geometry.arclength[ns : nf + 1, k]
            k_geom = geometry.coordinates[ns : nf + 1, :, k]
            k_cpos = geometry.cpos[ns : nf + 1, k]

            # ==================== HP surface ====================
            if definition.swtwisted:
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
            if str(te_types[k]) == "flat":
                e = geometry.arclength[ns, k]
                self.key_points[0, :, k] = geometry.coordinates[ns, :, k]
                self.key_cpos[1, k] = -1
            else:
                # e = 0.5 * (d + geometry.arclength(ns,k));
                e = 0.99 * geometry.arclength[ns, k]
                self.key_points[0, :, k] = interpolator_wrap(k_arclen, k_geom, e)
                self.key_cpos[1, k] = interpolator_wrap(k_arclen, k_cpos, e)

            # 1 -> e
            self.key_points[1, :, k] = interpolator_wrap(k_arclen, k_geom, d)
            self.key_points[2, :, k] = interpolator_wrap(k_arclen, k_geom, c)
            # self.key_points(  ,:,k) = interpolator_wrap(geometry.arclength(ns:nf,k),self.geometry(ns:nf,:,k),z);
            self.key_points[3, :, k] = interpolator_wrap(k_arclen, k_geom, b)
            self.key_points[4, :, k] = interpolator_wrap(k_arclen, k_geom, a)
            self.key_arcs[0, k] = geometry.arclength[ns, k]
            self.key_arcs[1, k] = e
            self.key_arcs[2, k] = d
            self.key_arcs[3, k] = c
            # self.key_arcs(  ,k)   = z;
            self.key_arcs[4, k] = b
            self.key_arcs[5, k] = a
            self.key_arcs[6, k] = 0  # le
            self.key_cpos[0, k] = geometry.cpos[ns, k]  # te, hp surface
            #            2   -> e
            self.key_cpos[2, k] = interpolator_wrap(k_arclen, k_cpos, d)
            self.key_cpos[3, k] = interpolator_wrap(k_arclen, k_cpos, c)
            #                 self.key_cpos(  ,k) = interpolator_wrap(geometry.arclength(ns:nf,k),self.cpos(ns:nf,k),z);
            self.key_cpos[4, k] = interpolator_wrap(k_arclen, k_cpos, b)
            self.key_cpos[5, k] = interpolator_wrap(k_arclen, k_cpos, a)
            self.key_cpos[6, k] = interpolator_wrap(k_arclen, k_cpos, 0)

            # ==================== LP surface ====================
            if definition.swtwisted:
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
            if str(te_types[k]) == str("flat"):
                e = geometry.arclength[nf, k]
                self.key_points[9, :, k] = geometry.coordinates[nf, :, k]
                self.key_cpos[11, k] = 1
            else:
                # e = 0.5 * (d + geometry.arclength(nf,k));
                e = 0.98 * geometry.arclength[nf, k]
                self.key_points[9, :, k] = interpolator_wrap(k_arclen, k_geom, e)
                self.key_cpos[11, k] = interpolator_wrap(k_arclen, k_cpos, e)
            self.key_points[5, :, k] = interpolator_wrap(k_arclen, k_geom, a)
            self.key_points[6, :, k] = interpolator_wrap(k_arclen, k_geom, b)
            # self.key_points(  ,:,k) = interpolator_wrap(geometry.arclength(ns:nf,k),self.geometry(ns:nf,:,k),z);
            self.key_points[7, :, k] = interpolator_wrap(k_arclen, k_geom, c)
            self.key_points[8, :, k] = interpolator_wrap(k_arclen, k_geom, d)
            # 10   -> e
            self.key_arcs[7, k] = a
            self.key_arcs[8, k] = b
            # self.key_arcs( ,k)   = z;
            self.key_arcs[9, k] = c
            self.key_arcs[10, k] = d
            self.key_arcs[11, k] = e
            self.key_arcs[12, k] = geometry.arclength[nf, k]
            self.key_cpos[7, k] = interpolator_wrap(k_arclen, k_cpos, a)
            self.key_cpos[8, k] = interpolator_wrap(k_arclen, k_cpos, b)
            # self.key_cpos(  ,k) = interpolator_wrap(geometry.arclength(ns:nf,k),self.cpos(ns:nf,k),z);
            self.key_cpos[9, k] = interpolator_wrap(k_arclen, k_cpos, c)
            self.key_cpos[10, k] = interpolator_wrap(k_arclen, k_cpos, d)
            # 12   -> e
            self.key_cpos[12, k] = geometry.cpos[nf, k]  # te, lp surface

        # find the points used by each shear web
        component_groups = [
            definition.components[name].group for name in definition.components
        ]
        self.web_indices = []
        self.web_arcs = []
        self.web_cpos = []
        self.web_points = []
        self.web_areas = []
        self.web_width = []
        self.web_bonds = []
        for ksw in range(max(component_groups)):  # for each shear web
            # pre-allocating arrays
            self.web_indices.append([])
            self.web_arcs.append(np.ndarray((2, num_istations)))
            self.web_cpos.append(np.ndarray((2, num_istations)))
            self.web_points.append(np.ndarray((2, 3, num_istations)))
            self.web_areas.append(np.ndarray((num_istations - 1)))
            self.web_width.append(np.ndarray(num_istations))
            self.web_bonds.append(np.ndarray((2, num_istations - 1)))

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
                le = self.key_labels.index("le")
            except ValueError:
                print(f"HP extent label \"{hp['pt']}\" not defined.")
            # get shear web placement on HP side
            if hp["pt"]:
                try:
                    n = self.key_labels[0 : le + 1].index(hp["pt"])  ## EMA
                except ValueError:
                    print(f"HP extent label \"{hp['pt']}\" not defined.")
                self.web_indices[ksw].append(n)
                self.web_arcs[ksw][0, :] = self.key_arcs[n, :]
                self.web_cpos[ksw][0, :] = self.key_cpos[n, :]
                n = n - 1
                self.web_points[ksw][0, :, :] = self.key_points[n, :, :]
            elif hp["pt1"]:
                f = float(hp["fraction"])
                if f <= 0 or f >= 1:
                    raise Exception(
                        f"Component group {ksw}: HP extent fraction={f}, which is outside range (0..1)"
                    )
                try:
                    n1 = self.key_labels[0 : le + 1].index(hp["pt1"])
                except:
                    print(f"HP extent label \"{hp['pt1']}\" not defined.")
                try:
                    n2 = self.key_labels[0 : le + 1].index(hp["pt2"])
                except:
                    print(f"HP extent label \"{hp['pt2']}\" not defined.")
                self.web_indices[ksw].append(np.nan)
                p1 = self.key_arcs[n1, :]
                p2 = self.key_arcs[n2, :]
                p = (1 - f) * p1 + f * p2
                self.web_arcs[ksw][0, :] = p
                for k in range(num_istations):
                    self.web_cpos[ksw][0, k] = interpolator_wrap(
                        k_arclen, k_cpos, p[k]
                    )
                    self.web_points[ksw][0, :, k] = interpolator_wrap(
                        k_arclen, k_geom, p[k]
                    )
            elif hp["pt3"]:
                try:
                    n3 = self.key_labels[0 : le + 1].index(hp["pt3"])
                except:
                    print(f"HP extent label \"{hp['pt3']}\" not defined.")
                self.web_indices[ksw].append(np.nan)
                p3 = self.key_cpos[n3, :]
                p = p3 - float(hp["mm_offset"]) / 1000
                iMax = self.key_labels[0, le + 1].index("d")
                # NOTE potential for error here - array shapes TBD -kb
                pMax = np.multiply(
                    self.key_cpos[iMax, :], np.transpose(self.ichord)
                )
                p[np.abs(p) > np.abs(pMax)] = pMax[np.abs(p) > np.abs(pMax)]
                iMin = self.key_labels[0 : le + 1].index("a")
                # NOTE same issue here -kb
                pMin = np.multiply(
                    self.key_cpos[iMin, :], np.transpose(self.ichord)
                )
                p[np.abs(p) < np.abs(pMin)] = pMin[np.abs(p) < np.abs(pMin)]
                self.web_cpos[ksw][0, :] = p
                for k in range(num_istations):
                    self.web_arcs[ksw][0, k] = interpolator_wrap(
                        self.cpos[ns : nf + 1, :, k],
                        geometry.arclength[ns : nf + 1, :, k],
                        p[k],
                    )
                    self.web_points[ksw][0, :, k] = interpolator_wrap(
                        k_cpos, k_geom, p[k]
                    )
            else:
                raise Exception(
                    "Shear web geometry HP extents not defined correctly (e.g., 0.5b-c, b, b+200)"
                )
            # get shear web placement on LP side
            if lp["pt"]:
                try:
                    n = self.key_labels[le:].index(lp["pt"]) + le
                    self.web_indices[ksw].append(n)
                    self.web_arcs[ksw][1, :] = self.key_arcs[n, :]
                    self.web_cpos[ksw][1, :] = self.key_cpos[n, :]
                    self.web_points[ksw][1, :, :] = self.key_points[n, :, :]
                except:
                    print(f"LP extent label \"{lp['pt']}\" not defined.")

            elif lp["pt1"]:
                f = float(lp["fraction"])
                if f < 0 or f > 1:
                    raise Exception(
                        f"Component group {ksw}: LP extent fraction={f}, which is outside range [0..1]"
                    )
                try:
                    n1 = self.key_labels[le:].index(lp["pt1"]) + le
                except:
                    print(f"LP extent label \"{lp['pt1']}\" not defined.")
                try:
                    n2 = self.key_labels[le:].index(lp["pt2"]) + le
                except:
                    print(f"LP extent label \"{lp['pt2']}\" not defined.")
                self.web_indices[ksw].append(np.nan)
                p1 = self.key_arcs[n1, :]
                p2 = self.key_arcs[n2, :]
                p = (1 - f) * p1 + f * p2
                self.web_arcs[ksw][1, :] = p
                for k in range(num_istations):
                    self.web_cpos[ksw][1, k] = interpolator_wrap(
                        k_arclen, k_cpos, p[k]
                    )
                    self.web_points[ksw][1, :, k] = interpolator_wrap(
                        k_arclen, k_geom, p[k]
                    )
            elif lp["pt3"]:
                try:
                    n3 = self.key_labels[le:].index(lp["pt3"]) + le
                except:
                    print(f"LP extent label \"{lp['pt3']}\" not defined.")
                self.web_indices[ksw].append(np.nan)
                p3 = self.key_cpos[n3, :]
                p = p3 + float(lp["mm_offset"]) / 1000
                iMax = self.key_labels[le:].index("d") + le
                pMax = np.multiply(
                    self.key_cpos[iMax, :], np.transpose(self.ichord)
                )
                p[np.abs(p) > np.abs(pMax)] = pMax[np.abs(p) > np.abs(pMax)]
                iMin = self.key_labels[le:].index("a") + le
                pMin = np.multiply(
                    self.key_cpos[iMin, :], np.transpose(self.ichord)
                )
                p[np.abs(p) < np.abs(pMin)] = pMin[np.abs(p) < np.abs(pMin)]
                self.web_cpos[ksw][1, :] = p
                for k in range(num_istations):
                    self.web_arcs[ksw][1, k] = interpolator_wrap(
                        k_cpos, k_arclen, p[k]
                    )
                    self.web_points[ksw][1, :, k] = interpolator_wrap(
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
                        geometry.arclength[:, kc] >= self.key_arcs[kr, kc],
                        geometry.arclength[:, kc] <= self.key_arcs[kr + 1, kc],
                    )
                )

                # need at least two points
                npts = np.amax((npts, 2))

                # inboard curve arclengths
                ibarc = np.linspace(
                    self.key_arcs[kr, kc], self.key_arcs[kr + 1, kc], npts
                )

                # outboard curve arclengths
                obarc = np.linspace(
                    self.key_arcs[kr, kc + 1],
                    self.key_arcs[kr + 1, kc + 1],
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
                self.key_areas[kr, kc] = t1 + t2
                if kr == 0:
                    self.te_bond[kc] = dspan[0]
                if (num_areas / 2 + 1) == (kr + 1):
                    self.le_bond[kc] = dspan[0]

        # calculate areas used by shear webs
        # jcb: note that these areas come purely from the geometry and
        # do not take into account the thickness of the shell or
        # sparcap layup.
        for ksw in range(len(self.web_points)):
            for kc in range(num_istations - 1):
                ib = self.web_points[ksw][:, :, kc]
                ob = self.web_points[ksw][:, :, kc + 1]
                # treat each "rectangular" area as two triangles
                b1 = np.diff(ib, axis=0)
                b2 = np.diff(ob, axis=0)
                base1 = np.sqrt(np.sum(b1**2, 1))[0]
                base2 = np.sqrt(np.sum(b2**2, 1))[0]
                b1 = b1 / base1
                b2 = b2 / base2
                h1 = np.abs(np.dot((ob[0, :] - ib[0, :]), (1 - np.transpose(b1))))
                h2 = np.abs(np.dot((ib[1, :] - ob[1, :]), (1 - np.transpose(b2))))
                self.web_areas[ksw][kc] = 0.5 * (base1 * h1 + base2 * h2)
                self.web_width[ksw][kc] = base1
                # calculate edge (bond-line) lengths
                self.web_bonds[ksw][0:2, kc] = np.sqrt(np.sum((ob - ib) ** 2, 1))
            self.web_width[ksw][num_istations - 1] = base2

        return self
