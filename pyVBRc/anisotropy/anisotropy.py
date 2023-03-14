from abc import ABC
from typing import Union

import numpy as np
from numpy.typing import ArrayLike
from scipy.optimize import fsolve


class IsotropicMedium:
    def __init__(self, poisson_ratio: float, modulus: float, modulus_type: str):
        self.poisson_ratio = poisson_ratio

        if modulus_type == "shear":
            self.shear_modulus = modulus
        elif modulus_type == "youngs":
            self._youngs_modulus = modulus
            self.shear_modulus = modulus / (2 * (1 + poisson_ratio))
        elif modulus_type == "bulk":
            self._bulk_modulus = modulus
            mu = 3 * modulus * (1 - 2 * poisson_ratio) / (2 * (1 + poisson_ratio))
            self.shear_modulus = mu

    _bulk_modulus = None

    @property
    def bulk_modulus(self):
        if self._bulk_modulus is None:
            nu = self.poisson_ratio
            K = 2 * self.shear_modulus * (1 + nu) / (3 * (1 - 2 * nu))
            self._bulk_modulus = K
        return self._bulk_modulus

    _youngs_modulus = None

    @property
    def youngs_modulus(self):
        if self._youngs_modulus is None:
            E = 2 * self.shear_modulus * (1 + self.poisson_ratio)
            self._youngs_modulus = E
        return self._youngs_modulus

    _lame_first = None

    @property
    def lame_first_parameter(self):
        # commonly lambda
        if self._lame_first is None:
            nu = self.poisson_ratio
            lam = 2 * self.shear_modulus * nu / (1 - 2 * nu)
            self._lame_first = lam
        return self._lame_first

    _pwave_M = None

    @property
    def pwave_effective_modulus(self):
        if self._pwave_M is None:
            nu = self.poisson_ratio
            M = 2 * self.shear_modulus * (1 - nu) / (1 - 2 * nu)
            self._pwave_M = M
        return self._pwave_M


class AnisotropicMedium(ABC):

    model_name: str = None
    model_reference: dict = None


class AlignedInclusions(AnisotropicMedium):

    model_name = "AlignedInclusions"
    model_reference = dict(
        title="The effect of aspect ratio of inclusions on the elastic properties "
        "of unidirectionally aligned composites",
        authors=("Tandon, Gp P", "Weng, Gj J"),
        journal="Polymer Composites",
        volume=5,
        number=4,
        pages="327--333",
        year=1984,
    )

    def __init__(self, aspect_ratio: Union[ArrayLike, float]):

        aspect_ratio = np.asarray(aspect_ratio)
        lt_1 = np.sum(aspect_ratio < 1)
        eq_1 = np.sum(aspect_ratio == 1)
        gt_1 = np.sum(aspect_ratio > 1)

        n_aspects = aspect_ratio.size
        if lt_1 == n_aspects:
            self.inclusion_type = "discs"
        elif eq_1 == n_aspects:
            self.inclusion_type = "spheres"
        elif gt_1 == n_aspects:
            self.inclusion_type = "fibers"
        else:
            raise ValueError(
                "Detected a mix of aspect ratio types. Ensure that "
                "all aspect ratio values are in the same range: <1, =1 or >1. "
                "Split calculation to match these bounds."
            )

        self.aspect_ratio = aspect_ratio

    def _disc_g(self):
        alpha = self.aspect_ratio
        costerm = np.arccos(alpha) - alpha * np.sqrt(1 - alpha**2)
        alphaterm = alpha / ((1 - alpha**2) ** (3 / 2))
        return alphaterm * costerm

    def _fiber_g(self):
        alpha = self.aspect_ratio
        al2 = alpha**2
        al2m1 = al2 - 1
        costerm = alpha * np.sqrt(al2m1) - np.arccosh(alpha)
        alphaterm = alpha / (al2m1 ** (3.0 / 2.0))
        gval = alphaterm * costerm
        return gval

    def _empty_eshelby(self):
        return {}  # defaultdict(lambda: 0.0)

    def _eshelby_spherical(self, poisson_0: float):
        S = self._empty_eshelby()
        p0fac = 15 * (1 - poisson_0)
        S[1111] = (7 - 5 * poisson_0) / p0fac
        S[1122] = (5 * poisson_0 - 1) / p0fac
        S[1212] = (4 - 5 * poisson_0) / p0fac
        S[2222] = S[1111]
        S[3333] = S[1111]

        S[2233] = S[1122]
        S[3311] = S[1122]
        S[2211] = S[3311]
        S[1133] = S[1122]

        S[2323] = S[1212]
        S[3131] = S[1212]
        S[1313] = S[1212]
        S[3232] = S[2323]

        return S

    def _eshelby_disc_fiber(self, poisson_0: float):
        if self.inclusion_type == "discs":
            print("getting disc g")
            g = self._disc_g()
        elif self.inclusion_type == "fibers":
            print("getting fiber g")
            g = self._fiber_g()
        else:
            raise RuntimeError(f"Unexpected inclusion_type: {self.inclusion_type}.")

        print("getting eshelby tensor entries")
        S = self._empty_eshelby()
        nu0 = poisson_0
        al = self.aspect_ratio
        al2 = al**2
        al2m1 = al2 - 1
        nu0m1 = 1 - nu0

        S[1111] = (
            1.0
            - 2.0 * nu0
            + (3 * al2 - 1) / al2m1
            - (1 - 2 * nu0 + 3 * al2 / al2m1) * g
        )
        S[1111] = S[1111] / (2 * nu0m1)

        S[2222] = 3 * al2 / (8 * nu0m1 * al2m1)
        S[2222] += (1 - 2 * nu0 - 9 / (4 * al2m1)) * g / (4 * nu0m1)
        S[3333] = S[2222]

        S[2233] = (al2 / (2 * al2m1) - (1 - 2 * nu0 + 3 / (4 * al2m1)) * g) / (
            4 * nu0m1
        )
        S[3322] = S[2233]

        S[2211] = -al2 / (2 * nu0m1 * al2m1) + (3 * al2 / al2m1 - (1 - 2 * nu0)) * g / (
            4 * nu0m1
        )
        S[3311] = S[2211]

        S[1122] = (1 - 2 * nu0 + 1 / al2m1) / (2 * nu0m1) + (
            1 - 2 * nu0 + 3 / (2 * al2m1)
        ) * g / (2 * nu0m1)
        S[1133] = S[1122]

        S[2323] = (al2 / (2 * al2m1) + (1 - 2 * nu0 - 3 / (4 * al2m1)) * g) / (
            4 * nu0m1
        )
        S[3232] = S[2323]

        S[1212] = (
            1
            - 2 * nu0
            - (al2 + 1) / al2m1
            - (1 - 2 * nu0 - 3 * (al2 + 1) / al2m1) * g / 2.0
        )
        S[1212] = S[1212] / (4 * nu0m1)
        S[1313] = S[1212]

        return S

    def _eshelby_tensor(self, poisson_0: float):
        # Sijkl
        if self.inclusion_type == "spheres":
            return self._eshelby_spherical(poisson_0)
        else:
            return self._eshelby_disc_fiber(poisson_0)

    def effective_shear_modulus(
        self,
        matrix_material: IsotropicMedium,
        inclusion_material: IsotropicMedium,
        vol_fraction: float,
    ):
        S = self._eshelby_tensor(matrix_material.poisson_ratio)

        # in-plane
        mu0 = matrix_material.shear_modulus
        mu1 = inclusion_material.shear_modulus
        c = vol_fraction
        mu_12 = (1 + c / (mu0 / (mu1 - mu0) + 2 * (1 - c) * S[1212])) * mu0

        # out-plane
        mu_23 = (1 + c / (mu0 / (mu1 - mu0) + 2 * (1 - c) * S[2323])) * mu0

        return mu_12, mu_23

    def _D_values(
        self,
        matrix_material: IsotropicMedium,
        inclusion_material: IsotropicMedium,
    ):
        lam0 = matrix_material.lame_first_parameter
        lam1 = inclusion_material.lame_first_parameter
        dlam = lam1 - lam0
        mu0 = matrix_material.shear_modulus
        mu1 = inclusion_material.shear_modulus

        D1 = 1.0 + 2.0 * (mu1 - mu0) / dlam
        D2 = (lam0 + 2.0 * mu0) / dlam
        D3 = lam0 / dlam
        return D1, D2, D3

    def _B_values(
        self,
        matrix_material: IsotropicMedium,
        D1,
        D2,
        D3,
        vol_fraction: Union[ArrayLike, float],
    ):

        S = self._eshelby_tensor(matrix_material.poisson_ratio)
        c = vol_fraction
        B1 = c * D1 + D2 + (1.0 - c) * (D1 * S[1111] + 2 * S[2211])
        B2 = c + D3 + (1.0 - c) * (D1 * S[1122] + S[2222] + S[2233])
        B3 = c + D3 + (1.0 - c) * (S[1111] + (1 + D1) * S[2211])
        B4 = c * D1 + D2 + (1.0 - c) * (S[1122] + D1 * S[2222] + S[2233])
        B5 = c + D3 + (1.0 - c) * (S[1122] + S[2222] + D1 * S[2233])

        # A_first = 2. * B2 * B3
        # A_sub = B1 * (B4 + B5)
        # A = A_first - A_sub = 0 at some point
        return B1, B2, B3, B4, B5

    def _A_coefficients(
        self,
        matrix_material: IsotropicMedium,
        inclusion_material: IsotropicMedium,
        vol_fraction: Union[ArrayLike, float],
    ):

        D1, D2, D3 = self._D_values(matrix_material, inclusion_material)
        B1, B2, B3, B4, B5 = self._B_values(matrix_material, D1, D2, D3, vol_fraction)

        A1 = D1 * (B4 + B5) - 2.0 * B2
        A2 = (1.0 + D1) * B2 - (B4 + B5)
        A3 = B1 - D1 * B3
        A4 = (1.0 + D1) * B1 - 2.0 * B3
        A5 = (1.0 - D1) / (B4 - B5)

        A_sub = B1 * (B4 + B5)
        A_first = 2.0 * B2 * B3
        A = A_first - A_sub

        return A, A1, A2, A3, A4, A5

    def effective_youngs_modulus(
        self,
        matrix_material: IsotropicMedium,
        inclusion_material: IsotropicMedium,
        vol_fraction: float,
    ):

        if np.any(self.aspect_ratio > 1):
            raise NotImplementedError(
                "effective_youngs_modulus is only valid for aspect ratios <= 1."
            )

        A_i = self._A_coefficients(matrix_material, inclusion_material, vol_fraction)
        c = vol_fraction
        E0 = matrix_material.youngs_modulus
        # longitudinal youngs
        nu0 = matrix_material.poisson_ratio
        e11f = A_i[1] + 2 * nu0 * A_i[2]
        e11f2 = c * e11f / A_i[0]
        E11 = E0 / (1 + e11f2)

        # transverse youngs
        c1 = -2.0 * nu0 * A_i[3]
        c2 = (1.0 - nu0) * A_i[4]
        c3 = (1.0 + nu0) * A_i[5] * A_i[0]
        e22f = c1 + c2 + c3
        c4 = e22f / (2.0 * A_i[0])
        cc4 = c * c4
        denom = 1.0 + cc4
        E22 = E0 / denom

        return E11, E22

    def effective_bulk_modulus(
        self,
        matrix_material: IsotropicMedium,
        inclusion_material: IsotropicMedium,
        vol_fraction: float,
    ):

        if np.any(self.aspect_ratio > 1):
            raise NotImplementedError(
                "effective_youngs_modulus is only valid for aspect ratios <= 1."
            )

        A_i = self._A_coefficients(matrix_material, inclusion_material, vol_fraction)
        c = np.asarray(vol_fraction)
        if c.ndim == 0:
            c = np.array(
                [
                    c,
                ]
            )

        _, mu_23 = self.effective_shear_modulus(
            matrix_material, inclusion_material, vol_fraction
        )
        E11, E22 = self.effective_youngs_modulus(
            matrix_material, inclusion_material, vol_fraction
        )
        nu0 = matrix_material.poisson_ratio
        # plane-strain bulk modulus:
        K0bar = matrix_material.lame_first_parameter + matrix_material.shear_modulus
        init_guess = mu_23

        K23 = np.empty(c.shape)
        for ic, cval in enumerate(c):
            all_the_args = (
                K0bar,
                nu0,
                A_i[3][ic],
                A_i[4][ic],
                A_i[0][ic],
                E11[ic],
                E22[ic],
                mu_23[ic],
                cval,
            )
            result = fsolve(_bulk_mod_solver_func, init_guess[ic], args=all_the_args)
            K23[ic] = result
        return K23


def _bulk_mod_solver_func(K23, K0bar, nu0, A3, A4, A, E11, E22, mu23, c):
    # eqs 36 and 37
    nu12 = E11 / E22 - E11 / 4 * (1 / mu23 + 1 / K23)
    rhsf1 = c * (2 * (nu12 - nu0) * A3 + (1 - nu0 * (1 + 2 * nu12)) * A4) / A
    rhs_denom = 1 - nu0 * (1 + 2 * nu12) + rhsf1
    rhs = (1 + nu0) * (1 - 2 * nu0) / rhs_denom
    lhs = K23 / K0bar
    return lhs - rhs  # will equal 0 at true K23
