from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import ArrayLike


class StiffnessMatrix(ABC):

    _base_shape = (6, 6)

    def __init__(self, n_states: int):
        self.n_states = n_states
        self._actual_shape = self._base_shape + (n_states,)
        self.stiffness = np.zeros(self._actual_shape)
        self.fill_stiffness()

    @abstractmethod
    def fill_stiffness(self):
        """fill in the stiffness matrix"""


class TransverseIsotropicStiffness(StiffnessMatrix):
    """

    Following sayers 1992 notation (eqs 24-29):
    Sayers, Colin M. "Elastic anisotropy of short-fibre reinforced
        composites." International journal of solids and structures 29,
        no. 23 (1992): 2933-2944. https://doi.org/10.1016/0020-7683(92)90150-R


    Parameters
    ----------
    E_longitudinal :  longitudinal Youngs
    E_transverse: transverse Youngs
    G_inplane :  in-plane shear modulus
    G_outplane: out of plane shear modulus
    nu_major: major poisson ratio
    K_plain: plane-strain bulk modulus

    """

    def __init__(
        self,
        E_longitudinal: ArrayLike,
        E_transverse: ArrayLike,
        G_inplane: ArrayLike,
        G_outplane: ArrayLike,
        nu_major: ArrayLike,
        K_plain: ArrayLike,
    ):

        self.E_longitudinal = np.asarray(E_longitudinal)
        self.E_transverse = np.asarray(E_transverse)
        self.G_inplane = np.asarray(G_inplane)
        self.G_outplane = np.asarray(G_outplane)
        self.nu_major = np.asarray(nu_major)
        self.K_plain = np.asarray(K_plain)

        in_shape = self.E_longitudinal.shape
        for mat in (
            "E_longitudinal",
            "G_inplane",
            "E_transverse",
            "G_outplane",
            "nu_major",
            "K_plain",
        ):
            ndim = getattr(self, mat).ndim
            if ndim != 1:
                raise ValueError(f"All input arrays must be 1D. {mat} is {ndim}.")

            matshape = getattr(self, mat).shape
            if matshape != in_shape:
                raise ValueError(
                    f"All input arrays must have the same shape: "
                    f"{mat}, {matshape}, does not match {in_shape}."
                )

        super().__init__(in_shape[0])

    def fill_stiffness(self):
        # fill the stiffness matrix

        # orig:
        mu12 = self.G_outplane
        mu31 = self.G_inplane

        C11 = mu12 + self.K_plain
        nusq = self.nu_major**2
        C33 = self.E_longitudinal + 4 * nusq * self.K_plain
        C12 = self.K_plain - mu12
        C13 = 2 * self.nu_major * self.K_plain
        C44 = mu31
        C66 = mu12

        # now fill in the stiffness matrix (note index offset for python)

        # diagonals
        self.stiffness[0, 0, :] = C11
        self.stiffness[1, 1, :] = C11
        self.stiffness[2, 2, :] = C33
        self.stiffness[3, 3, :] = C44
        self.stiffness[4, 4, :] = C44
        self.stiffness[5, 5, :] = C66

        # off diagonals (upper-right)
        self.stiffness[0, 1, :] = C12
        self.stiffness[0, 2, :] = C13
        self.stiffness[1, 2, :] = C13

        # symmetric lower-left
        self.stiffness[1, 0, :] = C12
        self.stiffness[2, 0, :] = C13
        self.stiffness[2, 1, :] = C13

    @property
    def a_1(self):
        # anisotropy factor a_1
        C = self.stiffness
        t1 = C[0, 0] + C[2, 2]
        t2 = -2.0 * C[0, 2]
        t3 = -4.0 * C[3, 3]
        return t1 + t2 + t3

    @property
    def a_2(self):
        # anisotropy factor a_2
        C = self.stiffness
        return C[0, 0] - 3 * C[0, 1] + 2 * C[0, 2] - 2 * C[3, 3]

    @property
    def a_3(self):
        # anisotropy factor a_3
        C = self.stiffness
        return 4 * C[0, 0] - 3 * C[2, 2] - C[0, 2] - 2 * C[3, 3]
