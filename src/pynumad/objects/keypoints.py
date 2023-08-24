import numpy as np


class KeyPoints:
    """This class holds keypoint information of a blade

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

    def __init__(self, shape=(0, 0)):
        M, N = shape
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
        self.key_points = np.zeros((M - 2, 3, N))  #
        self.key_arcs = np.zeros((M + 1, N))
        self.key_cpos = np.zeros((M + 1, N))
        self.key_areas = np.zeros((M, N - 1))

        self.le_bond = np.zeros((N - 1))
        self.te_bond = np.zeros((N - 1))

        self.web_indices: list = None
        self.web_points: list = None
        self.web_arcs: list = None
        self.web_cpos: list = None
        self.web_areas: list = None
        self.web_width: list = None
        self.web_bonds: list = None
