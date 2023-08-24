# for type hints
from numpy import ndarray
from pynumad.objects.airfoil import Airfoil
from pynumad.objects.station import Station


class Definition:
    def __init__(self):
        self.components: list = None
        self.shearweb: list = None
        self.materials: list = None
        self.stacks: ndarray = None
        self.swstacks: list = None
        self.ispan: ndarray = None
        self.aerocenter: ndarray = None
        self.chord: ndarray = None
        self.chordoffset: ndarray = None
        self.degreestwist: ndarray = None
        self.leband: float = None
        self.percentthick: ndarray = None
        self.prebend: ndarray = None
        self.span: ndarray = None
        self.sparcapoffset: ndarray = None
        self.sparcapwidth: ndarray = None
        self.stations: list = []
        self.sweep: ndarray = None
        self.teband: float = None
        self.te_type: list = None

        # init properties
        self._natural_offset: int = 1
        self._rotorspin: int = 1
        self._swtwisted: int = 0

    @property
    def natural_offset(self):
        """
        1 = offset by max thickness location,
        0 = do not offset to max thickness
        """
        return self._natural_offset

    @natural_offset.setter
    def natural_offset(self, new_natural_offset):
        if not (new_natural_offset == 0 or new_natural_offset == 1):
            raise Exception("natural_offset must be 0 or 1")
        else:
            self._natural_offset = new_natural_offset

    @property
    def rotorspin(self):
        """
        Rotor Spin, 1= CW rotation looking downwind, -1= CCW rotation
        """
        return self._rotorspin

    @rotorspin.setter
    def rotorspin(self, new_rotorspin):
        if not (new_rotorspin == 1 or new_rotorspin == -1):
            raise Exception("rotorspin must be 1 (cw) or -1 (ccw)")
        else:
            self._rotorspin = new_rotorspin

    @property
    def swtwisted(self):
        """
        Shear Web,
        0 = planar shear webs,
        1 = shear webs twisted by blade twist
        """
        return self._swtwisted

    @swtwisted.setter
    def swtwisted(self, new_swtwisted):
        if not (new_swtwisted == 0 or new_swtwisted == 1):
            raise Exception("swtwisted must be 0 or 1")
        else:
            self._swtwisted = new_swtwisted

    def _compare(self, other):
        """
        Parameters
        ----------
        other : Definition

        Returns
        -------
        bool
        """

        # compare components

        # compare shearweb

        # compare materials

    def add_station(self, af: Airfoil, spanlocation: float):
        """This method adds a station

        Specifically, the station object is created
        and appended to self.stations.

        Parameters
        ----------
        af : airfoil
        spanlocation : float

        Returns
        -------
        None

        Example
        -------
        ``blade.add_station(af,spanlocation)`` where  ``af`` = airfoil filename
        or ``AirfoilDef`` object
        """
        new_station = Station(af)
        new_station.spanlocation = spanlocation
        new_station.parent = self
        self.stations.append(new_station)

        return self
