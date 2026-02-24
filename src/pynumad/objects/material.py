from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class Material:
    """Material class

    Attributes
    ----------
    name : str
        User selected name of the material
    type : str
        Two options: isotropic or orthotropic
    layerthickness : float
        Layer thickness [mm]
    ex : float
        Longitudinal elastic modulus [Pa]
    ey : float
        Transverse elastic modulus [Pa]
    ez : float
        Through-the-thickness elastic modulus in the
        principal material coordinates [Pa]
    gxy : float
        In-plane shear modulus [Pa]
    gyz : float
        Transverse shear modulus [Pa]
    gxz : float
        Transverse shear modulus [Pa]
    prxy : float
        In-plane Poisson ratio [ ]
    pryz : float
        Transverse Poisson ratio [ ]
    prxz : float
        Transverse Poisson ratio [ ]
    density : float
        Cured mass density [kg/m2]
    drydensity : float
        Density of fabric
    uts : float
        1 x 3 array of ultimate tensile strength design values.
        Sequence: SL , ST, Sz, 1 x 1 for isotropic.
    ucs : float
        1 x 3 array of ultimate compressive strength design values.
        Sequence: SL , ST, Sz, 1 x 1 for isotropic.
    uss : float
        1 x 3 array of ultimate shear strength design values.
        Sequence: SLT , STz, SLz, 1 x 1 for isotropic.
    xzit : float
        Lz tensile inclination parameter for Puck failure index
    xzic : float
        Lz compressive inclination parameter for Puck failure index
    yzit : float
        Tz tensile inclination parameter for Puck failure index
    yzic : float
        Tz compressive inclination parameter for Puck failure index
    g1g2 : float
        Fracture toughness ratio between GI (mode I) and GII (mode II) [ ]
    alp0 : float
        Fracture angle under pure transverse compression [degrees]
    etat : float
        Transverse friction coefficient for Larc [ ]
    etal : float
        Longitudinal friction coefficient for Larc [ ]
    m : list
        Fatigue slope exponent [ ]
    gamma_mf : list
        From DNL-GL standard, fatigue strength reduction factor
    gamma_ms : list
        From DNV-GL standard, short term strength reduction factor
    reference : str
    """

    name: Optional[str] = None
    type: Optional[str] = None
    reference: Optional[str] = None
    layerthickness: Optional[float] = None
    ex: Optional[float] = None
    ey: Optional[float] = None
    ez: Optional[float] = None
    gxy: Optional[float] = None
    gyz: Optional[float] = None
    gxz: Optional[float] = None
    prxy: Optional[float] = None
    pryz: Optional[float] = None
    prxz: Optional[float] = None
    density: Optional[float] = None
    drydensity: Optional[float] = None
    uts: Optional[float] = None
    ucs: Optional[float] = None
    uss: Optional[float] = None
    xzit: Optional[float] = None
    xzic: Optional[float] = None
    yzit: Optional[float] = None
    yzic: Optional[float] = None
    g1g2: Optional[float] = None
    alp0: Optional[float] = None
    etat: Optional[float] = None
    etal: Optional[float] = None
    m: Optional[list] = None
    gamma_mf: Optional[list] = None
    gamma_ms: Optional[list] = None


@dataclass
class Layer:
    """A single layer entry within a :class:`~pynumad.objects.materialdb.MaterialDatabaseEntry`.

    Attributes
    ----------
    layer_name : str
        Name of the base material used in this layer.
    thicknessA : float
        Nominal ply thickness [m].
    thicknessB : float
        Secondary thickness (e.g. for tapered plies) [m].
    quantity : int
        Number of plies.
    theta : float
        Fibre orientation angle [degrees].
    """

    layer_name: Optional[str] = None
    thicknessA: Optional[float] = None
    thicknessB: Optional[float] = None
    quantity: Optional[int] = None
    theta: Optional[float] = None
