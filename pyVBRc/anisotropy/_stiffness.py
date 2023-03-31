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
        pass


class TransverseIsotropicStiffness(StiffnessMatrix):
    """

    Following sayers 1992 notation (eqs 24-29):

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
        # see, e.g., Nejati 2019

        # # pop out some vars for convenience
        # E0, E_transverse = self.E_longitudinal, self.E_transverse
        # nu_major = self.nu_major
        # K_plain = self.K_plain

        # self.E_longitudinal = np.asarray(E_longitudinal)
        # self.E_transverse = np.asarray(E_transverse)
        # self.G_inplane = np.asarray(G_inplane)
        # self.G_outplane = np.asarray(G_outplane)
        # self.nu_major = np.asarray(nu_major)
        # self.K_plain = np.asarray(K_plain)

        C11 = self.G_outplane + self.K_plain
        C33 = self.E_longitudinal + 4 * self.nu_major * self.K_plain
        C12 = self.K_plain - self.G_outplane
        C13 = 2 * self.nu_major * self.K_plain
        C44 = self.G_inplane
        C66 = self.G_outplane

        # denom_c11_c12 = (nu0*nu0-1)*E_transverse + 2 * E_longitudinal * (1 + nu0) * nu1 * nu1
        # C11 = E_longitudinal * (E_longitudinal * nu1* nu1 - E_transverse) / denom_c11_c12
        # C12 = - E_longitudinal * (E_longitudinal * nu1 * nu1 + nu0 * E_transverse) / denom_c11_c12
        #
        # C13 = E_longitudinal * E_transverse * nu1 / ((nu0 * nu0 - 1)* E_transverse + 2 * E_longitudinal * nu1 * nu1)
        # C33 = (nu0 - 1) * E_transverse * E_transverse / ((nu0 - 1) * E_transverse + 2 * E_longitudinal * nu1 * nu1)
        # C44 = self.G_outplane
        # C66 = self.G_inplane

        # now fill in the stiffness matrix (note index offset for python)

        # diagonals
        self.stiffness[0, 0, :] = C11
        self.stiffness[1, 1, :] = C11
        self.stiffness[2, 2, :] = C33
        self.stiffness[3, 3, :] = C44
        self.stiffness[4, 4, :] = C44
        self.stiffness[5, 5, :] = C66

        # off diagonals (upper-right)
        self.stiffness[0, 1] = C12
        self.stiffness[0, 2] = C13
        self.stiffness[1, 2] = C13

        # symmetric lower-left
        self.stiffness[1, 0] = C12
        self.stiffness[2, 0] = C13
        self.stiffness[2, 1] = C13
