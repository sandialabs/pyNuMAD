# for type hints
from numpy import ndarray

from pynumad.objects.airfoil import Airfoil
from pynumad.objects.station import Station
from pynumad.objects.component import Component
from pynumad.objects.material import Material


class Definition:
    """Definition Class
    
    A Definition object is designed to contain the
    basic defining features of a blade, which can
    be populated manually or by reading a blade file.

    Attributes
    ----------
    components : dict
        Dictionary of components indexed by component name
    shearweb : list
        List of shearwebs. len(shearweb) == number of shearwebs.
        shearweb[i] gives an array defining the ith shearweb.
    materials : dict
        Dictionary of material objects indexed by material name.
    stacks : ndarray
    swstacks : array
    ispan : ndarray
        Spanwise locations of interpolated output
    aerocenter : ndarray
        Aerodynamic center of airfoil (used only by NuMAD->FAST)
    chord : ndarray
        Chord distribution [m]
    chordoffset : ndarray
        Chordwise offset (in addition to natural offset)
    degreestwist : ndarray
        Twist distribution [degrees]
    percentthick : ndarray
        Percent thickness of airfoil [%]
    prebend : ndarray
        Blade prebend, reference axis location along x2 [m]
    span : ndarray
        Spanwise location of distributed properties [m]
    sparcapoffset : ndarray
    sparcapwidth : ndarray
        Locations of keypoints b & c, defines distance
        between keypoints b & c [mm]. First entry is the HP spar cap.
        Second entry is the LP spar cap
    stations : list
        List of station objects
    sweep : ndarray
        Blade sweep, reference axis location along x1 [m]
    teband : float
        Location of keypoint e
    leband : float
        Location of keypoint a
    te_type : list
    """
    def __init__(self):
        self.components: dict[str,Component] = None
        self.shearweb: list = None
        self.materials: dict[str,Material] = None
        self.stations: list[Station] = None
        self.te_type: list[str] = None
        self.stacks: ndarray = None
        self.swstacks: ndarray = None
        self.ispan: ndarray = None
        self.aerocenter: ndarray = None
        self.chord: ndarray = None
        self.chordoffset: ndarray = None
        self.degreestwist: ndarray = None
        self.percentthick: ndarray = None
        self.prebend: ndarray = None
        self.span: ndarray = None
        self.sparcapoffset: ndarray = None
        self.sparcapwidth: ndarray = None
        self.sweep: ndarray = None
        self.teband: float = None
        self.leband: float = None

        # init properties
        self._natural_offset: int = 1
        self._rotorspin: int = 1
        self._swtwisted: int = 0
        
    def __eq__(self, other):
        attrs = vars(self).keys()
        for attr in attrs:
            self_attr = getattr(self, attr)
            other_attr = getattr(other, attr)
            if isinstance(self_attr, (int, float, str, list, dict)):
                if self_attr != other_attr:
                    return False
            elif isinstance(self_attr, ndarray):
                if (self_attr != other_attr).any():
                    return False
        return True

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
        Rotor Spin, 
        1 = CW rotation looking downwind, 
        -1 = CCW rotation
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

    def add_station(self, af: Airfoil, spanlocation: float):
        """This method adds a station

        Specifically, a station object is created
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
        or ``Airfoil`` object
        """
        new_station = Station(af)
        new_station.spanlocation = spanlocation
        self.stations.append(new_station)

        return self
