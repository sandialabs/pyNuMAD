import numpy as np
from numpy import ndarray

from pynumad.utils.interpolation import interpolator_wrap
from pynumad.objects.airfoil import (
    get_airfoil_normals,
    get_airfoil_normals_angle_change,
)
from pynumad.utils.affinetrans import rotation, translation
from pynumad.objects.settings import BladeSettings


class Geometry:
    """Contains the geometry of a blade object

    Parameters
    ----------
    shape: string

    Attributes
    ----------
    aerocenter: array
        Aerodynamic center of airfoil (used only by NuMAD->FAST)
    chord : array
        Chord distribution [m]
    chordoffset : array
        Chordwise offset (in addition to natural offset)
    degreestwist : array
        Twist distribution [degrees]
    ispan : array
        Spanwise locations of interpolated output
    leband : float
        Location of keypoint a
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
    idegreestwist : array
        interpolated twist
    ichord : array
        interpolated chord
    ipercentthick : array
        interpolated thickness
    ic : array
    icamber : array
    ithickness : array
    ichordoffset : array
        interpolated offset
    iaerocenter : array
        interpolated aerocenter
    isweep : array
        interpolated sweep
    iprebend : array
        interpolated prebend
    xoffset : array
        natural offset
    profiles : array
        normalized airfoil profiles
    geometry : array
        actual x,y,z geometry
    arclength : array
        surface distance from L.E.
    cpos : array
        chordwise position
    hgGeometry : list
    hgKeypoints : list

    Example
    -------
    blade = BladeDef()
    """

    def __init__(self, shape, settings=None):
        n_points, n_stations = shape
        self.c: ndarray = np.zeros((n_points, n_stations))
        # nondimensional (0-1) value of x axis #
        self.camber: ndarray = np.zeros((n_points, n_stations))
        self.thickness: ndarray = np.zeros((n_points, n_stations))
        self.ic: ndarray = np.zeros((n_points, n_stations))
        self.icamber: ndarray = np.zeros((n_points, n_stations))
        self.ithickness: ndarray = np.zeros((n_points, n_stations))
        self.cpos = None
        # cpos x-axis parametrization for airfoil curve
        self.idegreestwist = None
        self.ichord = None
        self.ichordoffset: ndarray = None
        self.iaerocenter: ndarray = None
        self.idegreestwist: ndarray = None
        self.ipercentthick = None
        self.profiles = None
        self.coordinates = None
        self.xoffset = None
        self.LEindex = None
        self.iprebend: ndarray = None
        self.isweep: ndarray = None

        # init properties
        self._natural_offset: int = 1
        self._rotorspin: int = 1
        self._swtwisted: int = 0

    def update_airfoil_profile(self, k):
        """_summary_

        Parameters
        ----------
        k : _type_
            _description_

        Returns
        -------
        _type_
            _description_
        """
        thickness = self.ithickness[:, k]
        percentthick = self.ipercentthick[k]
        camber = self.icamber[:, k]
        c = self.ic[:, k]
        # jcb: note that I'm using max thickness about camber
        # instead of overall thickness of airfoil. We may need to
        # change this definition.
        maxthick = np.amax(thickness)
        tratio = percentthick / (maxthick * 100)
        thick = thickness * tratio
        hp = camber - 0.5 * thick
        lp = camber + 0.5 * thick
        profile1 = np.concatenate(([c[-1]], np.flipud(c), c[1:], [c[-1]]))
        profile2 = np.concatenate(([0], np.flipud(hp), lp[1:], [0]))
        profile = np.stack((profile1, profile2), axis=1)
        self.profiles[:, :, k] = profile
        return self

    def expand_blade_geometry_te(self, min_edge_lengths):
        """Adjusts trailing edges for solid models.

        Opens trailing edge
        TODO: docstring
        """
        nStations = self.coordinates.shape[2]

        for i_station in range(nStations):
            first_point = self.ichord[i_station] * self.profiles[-2, :, i_station]
            second_point = self.ichord[i_station] * self.profiles[1, :, i_station]
            edge_length = np.linalg.norm(second_point - first_point)
            # fprintf('station #i, edgeLength: #f\n',iStation,edgeLength*1000)

            maxthick = np.amax(self.ithickness[:, i_station])
            mtindex = np.argmax(self.ithickness[:, i_station])
            tratio = self.ipercentthick[i_station] / (maxthick * 100)
            airfoil_thickness = self.ithickness[:, i_station] * tratio
            onset = self.ic[mtindex, i_station]
            if edge_length < min_edge_lengths[i_station]:
                print(
                    f"Updating station: {i_station} TE separation from {edge_length} to {min_edge_lengths[i_station]}"
                )
                tet = (min_edge_lengths[i_station] - edge_length) / self.ichord[
                    i_station
                ]
                tes = 5 / 3 * tet  # slope of TE adjustment; 5/3*tet is "natural"
                # continuous first & second derivatives at 'onset'
                # maintain second & third derivative at mc==1 (TE)
                # adjust slope at mc==1 (TE) by tes
                A = np.array(
                    [[1, 1, 1, 1], [3, 4, 5, 6], [6, 12, 20, 30], [6, 24, 60, 120]]
                )
                d = np.array([[tet], [tes], [0], [0]])
                p = np.linalg.solve(A, d)
                # onset = self(k).maxthick;  # start of TE modification, measured from LE
                vec = (self.ic[:, i_station] - onset) / (1 - onset)
                mc = np.maximum(vec, np.zeros(vec.shape))
                temod = np.vstack([mc**3, mc**4, mc**5, mc**6]).T @ p
                temod = temod.reshape(-1)
                airfoil_thickness = airfoil_thickness + temod

                self.ithickness[:, i_station] = airfoil_thickness / tratio
                self.update_airfoil_profile(i_station)

                mtindex = np.argmax(self.ithickness[:, i_station])
                self.xoffset[0, i_station] = self.ic[mtindex, i_station]
                self.update_oml_geometry(i_station)
                # firstPoint=self.ichord(iStation)*self.profiles(end-1,:,iStation);
                # secondPont=self.ichord(iStation)*self.profiles(2,:,iStation);
                # edgeLength2=norm(secondPont-firstPoint);
                # fprintf('station #i, edgeLength: #f, New edgeLength=#f, percent diff: #f\n',iStation,edgeLength*1000,edgeLength2*1000,(edgeLength2-edgeLength)/edgeLength2*100)
        return self

    def update_oml_geometry(self, k):
        """ """
        x = self.profiles[:, 0, k]
        y = self.profiles[:, 1, k]

        # self.xoffset[0,k] = c[mtindex]
        if self._natural_offset:
            x = x - self.xoffset[0, k]
        x = x - self.ichordoffset[k]  # apply chordwise offset
        x = x * self.ichord[k] * -1 * self._rotorspin  # scale by chord
        y = y * self.ichord[k]  # scale by chord
        twist = -1 * self._rotorspin * self.idegreestwist[k]

        # prepare for hgtransform rotate & translate
        coords = np.zeros((len(x), 4))
        coords[:, 0] = np.cos(np.deg2rad(twist)) * x - np.sin(np.deg2rad(twist)) * y
        coords[:, 1] = np.sin(np.deg2rad(twist)) * x + np.cos(np.deg2rad(twist)) * y
        coords[:, 2] = np.zeros(len(x))
        coords[:, 3] = np.ones(len(x))

        # use the generating line to translate and rotate the coordinates
        # NOTE currently, rotation is not assigned from blade properties
        # and defaults to 0
        prebend_rot = 0
        sweep_rot = 0
        """
        jcb: This code, copied from NuMAD 2.0, causes each section to rotate out
        of plane so that its normal follows the generating line direction. Need
        to replace 'twistFlag' with '-1*self.rotorspin' and calculate the slopes
        based on the available data. For now, default to parallel sections.
        if isequal(blade.PresweepRef.method,'normal')
            sweep_slope = ppval(blade.PresweepRef.dpp,sta.LocationZ);
            sweep_rot = atan(sweep_slope*twistFlag);
        end
        if isequal(blade.PrecurveRef.method,'normal')
            prebend_slope = ppval(blade.PrecurveRef.dpp,sta.LocationZ);
            prebend_rot = atan(-prebend_slope);
        endc
        """
        transX = -1 * self._rotorspin * self.isweep[k]
        transY = self.iprebend[k]
        transZ = self.ispan[k]
        Ry = rotation("y", sweep_rot)
        Rx = rotation("x", prebend_rot)
        R = Ry @ Rx
        T = translation(transX, transY, transZ)
        coords = coords @ np.transpose(R) @ np.transpose(T)
        # save the transformed coordinates
        self.coordinates[:, :, k] = coords[:, 0:3]
        return self

    def get_profile_te_type(self, k: int):
        """Get trailing edge type for the profile at the kth span location

        Parameters
        ----------
        k

        Return
        ------
        tetype : str
        """
        xy = self.profiles[:, :, k]
        unitNormals = get_airfoil_normals(xy)
        angleChange = get_airfoil_normals_angle_change(unitNormals)
        disconts = np.flatnonzero(angleChange > 30)

        if np.std(angleChange) < 2:
            tetype = "round"
        elif len(disconts) > 1:
            tetype = "flat"
        else:
            tetype = "sharp"
        return tetype
