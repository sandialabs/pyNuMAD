import numpy as np
from numpy import ndarray

from pynumad.utils.interpolation import interpolator_wrap
from pynumad.objects.airfoil import (
    get_airfoil_normals,
    get_airfoil_normals_angle_change,
)
from pynumad.utils.affinetrans import rotation, translation


class Geometry:
    """Geometry class
    
    Contains the interpolated geometry of a blade

    Attributes
    ----------
    c : ndarray
        nondimensional (0-1) value of x axis
    camber : ndarray
    thickness : ndarray
    ic : ndarray
    icamber : ndarray
    ithickness : array
        interpolated thickness
    cpos : ndarray
        x-axis parametrization for airfoil curve
    ichord : ndarray
        interpolated chord
    ichordoffset : ndarray
        interpolated offset
    iaerocenter : ndarray
        interpolated aerocenter
    idegreestwist : ndarray
        interpolated twist
    ipercentthick : ndarray
        interpolated percent thickness
    profiles : array
        normalized airfoil profiles
    coordinates : ndarray
        actual x,y,z geometry
    xoffset : ndarray
        natural offset
    LEindex : ndarray
    iprebend : ndarray
        interpolated prebend
    isweep : ndarray
        interpolated sweep
    """

    def __init__(self, settings=None):
        self.c: ndarray = None
        self.camber: ndarray = None
        self.thickness: ndarray = None
        self.ic: ndarray = None
        self.icamber: ndarray = None
        self.ithickness: ndarray = None
        self.cpos: ndarray = None
        self.ichord: ndarray = None
        self.ichordoffset: ndarray = None
        self.iaerocenter: ndarray = None
        self.idegreestwist: ndarray = None
        self.ipercentthick: ndarray = None
        self.profiles: ndarray = None
        self.coordinates: ndarray = None
        self.xoffset: ndarray = None
        self.LEindex: ndarray = None
        self.iprebend: ndarray = None
        self.isweep: ndarray = None

        # init properties
        self._natural_offset: int = 1
        self._rotorspin: int = 1
        self._swtwisted: int = 0
        
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
    
        
    def initialize_arrays(self, num_points: int, num_stations: int, num_istations: int):
        """_summary_

        Parameters
        ----------
        size : tuple[int]
            num_points, num_stations = size
        """
        size = (num_points, num_stations)
        isize = (num_points, num_istations)        
        self.c = np.zeros(size)
        self.camber = np.zeros(size)
        self.thickness = np.zeros(size)
        self.ic = np.zeros(isize)
        self.icamber = np.zeros(isize)
        self.ithickness = np.zeros(isize)
        self.idegreestwist = np.zeros(isize)
        self.ichord = np.zeros(isize)
        self.ichordoffset = np.zeros(isize)
        self.iaerocenter = np.zeros(isize)
        self.idegreestwist = np.zeros(isize)
        self.ipercentthick = np.zeros(isize)
        self.iprebend = np.zeros(isize)
        self.isweep = np.zeros(isize)
        
        
    def generate(self, definition):
        """Populates geometry attributes based on a given blade defintion

        Parameters
        ----------
        definition : Definition

        Returns
        -------
        self
        """
        self.ispan = definition.ispan
        num_istations = self.ispan.size
        stations = definition.stations
        num_stations = len(stations)
        if num_stations > 0:
            num_points = len(stations[0].airfoil.c)
        else:
            raise Exception(
                "BladeDef must have at least one station before updating geometry."
            )

        # first station must be at blade root to prevent extrapolation
        assert stations[0].spanlocation == 0, "first station must be at the blade root"

        self.initialize_arrays(num_points, num_stations, num_istations)

        self._natural_offset = definition.natural_offset
        self._rotorspin = definition.rotorspin
        self._swtwisted = definition.swtwisted

        # Collect parameter tables from the stations.
        spanlocation = np.array(
            [stations[i].spanlocation for i in range(len(stations))]
        )
        tetype = ["_"] * num_stations
        for k in range(num_stations):
            station = stations[k]
            assert (
                len(station.airfoil.c) == num_points
            ), "Station airfoils must have same number of samples."

            self.c[:,k] = station.airfoil.c
            self.camber[:, k] = station.airfoil.camber
            self.thickness[:, k] = station.airfoil.thickness
            tetype[k] = station.airfoil.te_type

        # fix numerical issue due to precision on camber calculation
        # camber should start and end at y-values of zero
        self.camber[0, :] = np.zeros((1, num_stations))
        self.camber[-1, :] = np.zeros((1, num_stations))

        # Interpolate the station parameter tables.
        # Each column corresponds to an interpolated station.

        ## ic
        self.ic = interpolator_wrap(
            spanlocation, self.c, self.ispan, "pchip", axis=1
        )

        ## cpos
        self.cpos = np.concatenate(
            (
                -self.ic[-1, :].reshape(1, -1),
                -np.flipud(self.ic),
                self.ic[1:, :],
                self.ic[-1, :].reshape(1, -1),
            ),
            axis=0,
        )

        ## icamber
        self.icamber = interpolator_wrap(
            spanlocation, self.camber, self.ispan, "pchip", axis=1
        )

        ## ithickness
        self.ithickness = interpolator_wrap(
            spanlocation, self.thickness, self.ispan, "pchip", axis=1
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
                self.ithickness[-1, k] = 0

        # Interpolate the blade parameter curves.
        ## idegreestwist
        self.idegreestwist = interpolator_wrap(
            definition.span, definition.degreestwist, self.ispan, "pchip"
        )

        ## ichord
        self.ichord = interpolator_wrap(
            definition.span, definition.chord, self.ispan, "pchip"
        )

        ## ipercentthick
        absolutethick = np.multiply(definition.percentthick, definition.chord) / 100
        iabsolutethick = interpolator_wrap(
            definition.span, absolutethick, self.ispan, "pchip"
        )
        self.ipercentthick = iabsolutethick / self.ichord * 100
        # ensure that the interpolation doesn't reduce the percent
        # thickness beneath the thinnest airfoil
        self.ipercentthick[
            self.ipercentthick < np.amin(definition.percentthick)
        ] = np.amin(definition.percentthick)

        ## ichordoffset
        self.ichordoffset = interpolator_wrap(
            definition.span, definition.chordoffset, self.ispan, "pchip"
        )

        ## iaerocenter
        self.iaerocenter = interpolator_wrap(
            definition.span, definition.aerocenter, self.ispan, "pchip"
        )

        ## isweep
        if len(definition.sweep) == 0:
            definition.sweep = np.zeros((self.span.shape, self.span.shape))
        if len(definition.prebend) == 0:
            definition.prebend = np.zeros((self.span.shape, self.span.shape))

        self.isweep = interpolator_wrap(
            definition.span, definition.sweep, self.ispan, "pchip"
        )

        ## iprebend
        self.iprebend = interpolator_wrap(
            definition.span, definition.prebend, self.ispan, "pchip"
        )

        # Generate the blade surface self.
        n_istations = np.asarray(self.ispan).size
        n_areas = num_points * 2 + 1
        self.profiles = np.zeros((n_areas, 2, n_istations))
        self.coordinates = np.zeros((n_areas, 3, n_istations))
        self.xoffset = np.zeros((1, n_istations))
        self.LEindex = num_points

        for k in range(n_istations):
            self.update_airfoil_profile(k)
            mtindex = np.argmax(self.ithickness[:, k])
            self.xoffset[0, k] = self.ic[mtindex, k]
            self.update_oml_geometry(k)

        # Calculate the arc length of each curve
        self.arclength = np.zeros((n_areas, n_istations))
        self.HParcx0 = np.zeros((1, n_istations))
        self.LParcx0 = np.zeros((1, n_istations))

        le_index = self.LEindex
        for k in range(n_istations):
            xx = self.coordinates[:, 0, k]
            yy = self.coordinates[:, 1, k]
            zz = self.coordinates[:, 2, k]
            arclen = np.sqrt(np.diff(xx) ** 2 + np.diff(yy) ** 2 + np.diff(zz) ** 2)
            arclen = np.concatenate((np.array([0]), np.cumsum(arclen)), axis=0)
            self.arclength[:, k] = arclen
            LEarcsum = self.arclength[le_index, k]
            self.arclength[:, k] = self.arclength[:, k] - LEarcsum

            # find where x=0 intersects the surface
            self.HParcx0[0, k] = (
                interpolator_wrap(xx[1 : le_index + 1], arclen[1 : le_index + 1], 0) - LEarcsum
            )
            self.LParcx0[0, k] = (
                interpolator_wrap(xx[-2 : le_index - 1 : -1], arclen[-2 : le_index - 1 : -1], 0)
                - LEarcsum
            )

        return self

    def update_airfoil_profile(self, k):
        """Updates kth airfoil profile

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

    def expand_blade_geometry_te(self, min_edge_lengths: ndarray):
        """Adjusts trailing edges for solid models.

        Opens trailing edge

        Parameters
        ----------
        min_edge_lengths : ndarray

        Returns
        -------
        Self
        """
        nStations = self.coordinates.shape[2]

        for i_station in range(nStations):
            first_point = self.ichord[i_station] * self.profiles[-2, :, i_station]
            second_point = self.ichord[i_station] * self.profiles[1, :, i_station]
            edge_length = np.linalg.norm(second_point - first_point)
            # fprintf('station #i, edgeLength: #f\n',i_station,edgeLength*1000)

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
                # first_point=self.ichord(i_station)*self.profiles(end-1,:,i_station);
                # second_point=self.ichord(i_station)*self.profiles(2,:,i_station);
                # edgeLength2=norm(second_point-first_point);
                # fprintf('station #i, edgeLength: #f, New edgeLength=#f, percent diff: #f\n',i_station,edgeLength*1000,edgeLength2*1000,(edgeLength2-edgeLength)/edgeLength2*100)
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
