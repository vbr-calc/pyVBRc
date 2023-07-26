from abc import ABC, abstractmethod
from typing import Optional, Tuple, Union

import numpy as np
from numpy.typing import ArrayLike
from unyt import unyt_array

from pyVBRc.anisotropy._stiffness import TransverseIsotropicStiffness
from pyVBRc.vbrc_structure import VBRCstruct


class IsotropicMedium:
    """
    An end-member isotropic medium

    Parameters
    ----------
    poisson_ratio
        the poisson ratio of the material
    modulus
        either the shear, Youngs or bulk modulus of the material
    modulus_type : str
        one of "shear", "youngs" or "bulk" depending on the previous parameter
    density : Optional
        the density of the material. only needed for calculating velocities, can
        omit if you only want moduli.

    """

    def __init__(
        self, poisson_ratio: float, modulus: float, modulus_type: str, density=None
    ):
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

        self.density = density

    _bulk_modulus = None

    @property
    def bulk_modulus(self):
        """the bulk modulus"""
        if self._bulk_modulus is None:
            nu = self.poisson_ratio
            K = 2 * self.shear_modulus * (1 + nu) / (3 * (1 - 2 * nu))
            self._bulk_modulus = K
        return self._bulk_modulus

    _youngs_modulus = None

    @property
    def youngs_modulus(self):
        """the Youngs modulus"""
        if self._youngs_modulus is None:
            E = 2 * self.shear_modulus * (1 + self.poisson_ratio)
            self._youngs_modulus = E
        return self._youngs_modulus

    _lame_first = None

    @property
    def lame_first_parameter(self):
        """Lame's first parameter (lambda)"""
        # commonly lambda
        if self._lame_first is None:
            nu = self.poisson_ratio
            lam = 2 * self.shear_modulus * nu / (1 - 2 * nu)
            self._lame_first = lam
        return self._lame_first

    _pwave_M = None

    @property
    def pwave_effective_modulus(self):
        """Effective p-wave Modulus, 2 * G (1 - nu)/ (1 - 2 * nu) = K + 4/3 G"""
        if self._pwave_M is None:
            nu = self.poisson_ratio
            self._pwave_M = 2 * self.shear_modulus * (1 - nu) / (1 - 2 * nu)
        return self._pwave_M

    @property
    def plain_strain_bulk(self):
        lame = self.lame_first_parameter
        shear = self.shear_modulus
        return lame + shear

    _v_p = None

    @property
    def v_p(self):
        """p-wave velocity"""
        if self._v_p is None:
            if self.density is None:
                raise ValueError("calculating v_p requires density to be set.")
            self._v_p = np.sqrt(self.pwave_effective_modulus / self.density)
        return self._v_p

    _v_s = None

    @property
    def v_s(self):
        """shear wave velocity"""
        if self._v_s is None:
            if self.density is None:
                raise ValueError("calculating v_s requires density to be set.")
            self._v_s = np.sqrt(self.shear_modulus / self.density)
        return self._v_s

    def set_density(self, density):
        """
        Set the density for the material

        Parameters
        ----------
        density: scalar or array
        """
        self.density = density


class _AnisotropicMedium(ABC):

    model_name: str = None
    model_reference: dict = None

    @abstractmethod
    def get_stiffness_matrix(self):
        """returns the stiffness matrix"""


class AlignedInclusions(_AnisotropicMedium):

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
        self.matrix_material: IsotropicMedium = None
        self.inclusion_material: IsotropicMedium = None
        self.volume_fraction: ArrayLike = None
        self.composite_density: ArrayLike = None
        self._unravel_shape: tuple = None

    def set_material(
        self,
        matrix_material: IsotropicMedium,
        inclusion_material: IsotropicMedium,
        vol_fraction: ArrayLike,
    ):

        self.matrix_material = matrix_material
        self.inclusion_material = inclusion_material
        self.volume_fraction = vol_fraction

        d_matrix = matrix_material.density
        d_includ = inclusion_material.density
        if d_matrix is not None and d_includ is not None:
            self.composite_density = (
                vol_fraction * d_includ + (1 - vol_fraction) * d_matrix
            )
        else:
            self.composite_density = None

    def set_materials_from_vbrc_structures(
        self,
        vbrc_matrix: VBRCstruct,
        vbrc_inclusion: VBRCstruct,
        shear_modulus_location: Tuple[str],
        vol_frac: ArrayLike,
        ifreq: Optional[int] = None,
    ):
        """
        set the material properties using VBRc structures

        Parameters
        ----------
        vbrc_matrix: VBRCstruct
            the VBRc structure for the matrix phase
        vbrc_inclusion: VBRCstruct
            the VBRc structure for the inclusion phase
        shear_modulus_location: Tuple[str]
            the output shear modulus of the VBRc structure to use. e.g,
            ["elastic", "anharmonic", "Gu"] or ["anelastic", "andrade_psp", "M"]
        vol_frac: ArrayLike
            the volume fraction of vbrc_inclusion
        ifreq: Optional[int]
            the frequency index to select if the selected shear modulus is
            frequency dependent.
        """

        # TODO: allow subselections

        vol_frac = np.asarray(vol_frac)
        if vol_frac.ndim == 0:
            vol_frac = _promote_to_1d_array(vol_frac)

        mat_inc = []
        for iv, v in enumerate((vbrc_matrix, vbrc_inclusion)):

            # always frequency independent:
            rho = v.input.SV.rho

            nu = v.input.elastic.anharmonic.nu  # constant of length 1

            # this one might be frequency dependent
            G = v._get_nested_output_field(shear_modulus_location)
            if G.ndim - rho.ndim == 1:
                # frequency dependent!
                if ifreq is None:
                    vbrloc = ".".join(shear_modulus_location)
                    raise ValueError(
                        f"{vbrloc} is frequency dependent, "
                        f"supply a frequency index using ifreq"
                    )
                G = G[..., ifreq]

            if isinstance(G, unyt_array):
                G = G.value
                rho = rho.value

            if iv == 0:
                matrix_shape = rho.shape
            else:
                if rho.shape != matrix_shape:
                    raise RuntimeError(
                        "The matrix and inclusion arrays do not have"
                        f"the same shape: {matrix_shape}, {rho.shape}"
                    )

            # must be 1D
            G = G.ravel()
            rho = rho.ravel()
            nu = np.full(G.shape, nu)

            m = IsotropicMedium(nu, G, modulus_type="shear", density=rho)
            mat_inc.append(m)

        self._unravel_shape = matrix_shape

        if vol_frac.shape == matrix_shape:
            vol_frac = vol_frac.ravel()
        elif len(vol_frac) == 1:
            pass
        else:
            raise RuntimeError(
                "When material properties are arrays, vol_frac must "
                "be constant or the same shape. Found shape(vol_frac)"
                f"={vol_frac.shape} with material shape of {matrix_shape}."
            )

        self.set_material(mat_inc[0], mat_inc[1], vol_frac)

    def _require_material(self, func: str):
        if any(
            [
                attr is None
                for attr in (
                    self.matrix_material,
                    self.inclusion_material,
                    self.volume_fraction,
                )
            ]
        ):
            raise RuntimeError(
                f"You must set the material properties with `.set_material` before calling {func}"
            )

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
            g = self._disc_g()
        elif self.inclusion_type == "fibers":
            g = self._fiber_g()

        # print("getting eshelby tensor entries")
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

        S[1122] = -(1 - 2 * nu0 + 1 / al2m1) / (2 * nu0m1) + (
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

    @property
    def _eshelby_tensor(self):
        # Sijkl
        poisson_0 = self.matrix_material.poisson_ratio
        if self.inclusion_type == "spheres":
            return self._eshelby_spherical(poisson_0)
        else:
            return self._eshelby_disc_fiber(poisson_0)

    def _shear_moduli(self):
        S = self._eshelby_tensor

        # in-plane
        mu0 = self.matrix_material.shear_modulus
        mu1 = self.inclusion_material.shear_modulus
        c = self.volume_fraction
        mu_12 = (1 + c / (mu0 / (mu1 - mu0) + 2 * (1 - c) * S[1212])) * mu0

        # out-plane
        mu_23 = (1 + c / (mu0 / (mu1 - mu0) + 2 * (1 - c) * S[2323])) * mu0

        if isinstance(mu_12, np.ndarray) is False:
            mu_12 = _promote_to_1d_array(mu_12)
            mu_23 = _promote_to_1d_array(mu_23)

        return mu_12, mu_23

    def _D_values(
        self,
    ):
        lam0 = self.matrix_material.lame_first_parameter
        lam1 = self.inclusion_material.lame_first_parameter
        dlam = lam1 - lam0
        mu0 = self.matrix_material.shear_modulus
        mu1 = self.inclusion_material.shear_modulus

        D1 = 1.0 + 2.0 * (mu1 - mu0) / dlam
        D2 = (lam0 + 2.0 * mu0) / dlam
        D3 = lam0 / dlam
        return D1, D2, D3

    def _B_values(
        self,
        D1,
        D2,
        D3,
    ):

        S = self._eshelby_tensor
        c = self.volume_fraction
        B1 = c * D1 + D2 + (1.0 - c) * (D1 * S[1111] + 2 * S[2211])
        B2 = c + D3 + (1.0 - c) * (D1 * S[1122] + S[2222] + S[2233])
        B3 = c + D3 + (1.0 - c) * (S[1111] + (1 + D1) * S[2211])
        B4 = c * D1 + D2 + (1.0 - c) * (S[1122] + D1 * S[2222] + S[2233])
        B5 = c + D3 + (1.0 - c) * (S[1122] + S[2222] + D1 * S[2233])

        return B1, B2, B3, B4, B5

    def _A_coefficients(
        self,
    ):

        D1, D2, D3 = self._D_values()
        B1, B2, B3, B4, B5 = self._B_values(D1, D2, D3)

        A1 = D1 * (B4 + B5) - 2.0 * B2
        A2 = (1.0 + D1) * B2 - (B4 + B5)
        A3 = B1 - D1 * B3
        A4 = (1.0 + D1) * B1 - 2.0 * B3
        A5 = (1.0 - D1) / (B4 - B5)

        A_sub = B1 * (B4 + B5)
        A_first = 2.0 * B2 * B3
        A = A_first - A_sub

        if isinstance(A1, np.ndarray) is False:
            A = _promote_to_1d_array(A)
            A1 = _promote_to_1d_array(A1)
            A2 = _promote_to_1d_array(A2)
            A3 = _promote_to_1d_array(A3)
            A4 = _promote_to_1d_array(A4)
            A5 = _promote_to_1d_array(A5)

        return A, A1, A2, A3, A4, A5

    def _youngs_moduli(
        self,
    ):
        A_i = self._A_coefficients()
        c = self.volume_fraction
        E0 = self.matrix_material.youngs_modulus
        # longitudinal youngs
        nu0 = self.matrix_material.poisson_ratio
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

        return E11, E22, A_i

    def _poisson_bulk_moduli(self, mu23, E11, E22, A_i):
        c = np.asarray(self.volume_fraction)
        if c.ndim == 0:
            c = _promote_to_1d_array(c)

        nu0 = self.matrix_material.poisson_ratio

        # plane-strain bulk modulus of the matrix:
        K0bar = (
            self.matrix_material.lame_first_parameter
            + self.matrix_material.shear_modulus
        )

        # Zhao et al 1989 have explicit form for nu_12 (see Sayers 1992 too)
        numer = nu0 * (A_i[1] + 2 * nu0 * A_i[2]) + (A_i[3] - nu0 * A_i[4])
        denom = A_i[0] + c * (A_i[1] + 2 * nu0 * A_i[2])
        nu12 = nu0 - c * numer / denom

        # and now we can calculate K23
        term1 = 2 * (nu12 - nu0) * A_i[3] + (1 - nu0 * (1 + 2 * nu12)) * A_i[4]
        term2 = c * term1 / A_i[0]
        rhs_denom = 1 - nu0 * (1 + 2 * nu12) + term2
        K23 = K0bar * (1 + nu0) * (1 - 2 * nu0) / rhs_denom
        return nu12, K23

    def get_moduli(
        self,
    ):

        mu12, mu23 = self._shear_moduli()
        E11, E22, A_i = self._youngs_moduli()
        nu_12, K23 = self._poisson_bulk_moduli(mu23, E11, E22, A_i)

        return E11, E22, mu12, mu23, nu_12, K23

    def get_stiffness_matrix(self):
        self._require_material("get_stiffness_matrix")
        E11, E22, mu12, mu23, nu12, K23 = self.get_moduli()
        stiff = TransverseIsotropicStiffness(E11, E22, mu12, mu23, nu12, K23)
        return stiff

    def velocities(self, theta: ArrayLike):
        # theta can either be size 1 or same size as density

        density = self.composite_density
        if density is None:
            raise ValueError("calculating velocities requires densities.")

        theta = np.asarray(theta)
        if theta.ndim == 0:
            theta = _promote_to_1d_array(theta)

        if self._unravel_shape is not None:
            if theta.shape == self._unravel_shape:
                theta = theta.ravel()
            elif len(theta) == 1:
                pass
            else:
                raise RuntimeError(
                    "When material properties are arrays, theta must "
                    "be constant or the same shape. Found shape(theta)"
                    f"={theta.shape} with material shape of {self._unravel_shape}."
                )

        stiffness = self.get_stiffness_matrix()
        tc = ThomsenCalculator(density, stiffness.stiffness)
        tc.set_theta(theta)

        v_p, v_sv, v_sh = tc.v_p, tc.v_sv, tc.v_sh
        if self._unravel_shape is not None:
            v_p = np.reshape(v_p, self._unravel_shape)
            v_sh = np.reshape(v_sh, self._unravel_shape)
            v_sv = np.reshape(v_sv, self._unravel_shape)

        return v_p, v_sv, v_sh


class ThomsenCalculator:
    def __init__(self, density: ArrayLike, stiffness: ArrayLike, approx: bool = False):

        # Thomsen, Leon. "Weak elastic anisotropy." Geophysics 51.10 (1986): 1954-1966.
        # following form in Kendall 2000
        self.density = density
        self.stiffness = stiffness

        C = self.stiffness
        # epsilon = (C11 - C33) / (2 * C33)
        epsilon = (C[0, 0] - C[2, 2]) / (2 * C[2, 2])
        # gamma = (C66 - C44) / (2 * C44)
        gamma = (C[5, 5] - C[3, 3]) / (2 * C[3, 3])

        dstar_denom = 2 * C[2, 2] * C[2, 2]
        t1 = 2 * (C[0, 2] + C[3, 3]) ** 2
        t2 = -(C[2, 2] - C[3, 3]) * (C[0, 0] + C[2, 2] - 2 * C[3, 3])
        dstar = (t1 + t2) / dstar_denom

        self.epsilon = epsilon
        self.gamma = gamma
        self.dstar = dstar

        dweak_denom = 2 * C[2, 2] * (C[2, 2] - C[3, 3])
        self.d_weak = (
            (C[0, 2] + C[3, 3]) ** 2 - (C[2, 2] - C[3, 3]) ** 2
        ) / dweak_denom
        self.alpha_o = np.sqrt(C[2, 2] / density)
        self.beta_o = np.sqrt(C[3, 3] / density)
        self.theta: np.ndarray = None
        self.Dstar_theta: np.ndarray = None
        self.approx = approx

    def set_theta(self, theta: ArrayLike):
        self.theta = np.asarray(theta)
        self.Dstar_theta = self._dstar_theta()
        if self.approx:
            self.v_p = self._v_p()
            self.v_sv = self._vsv()
            self.v_sh = self._vsh()
        else:
            self.v_p = self._v_p_full()
            self.v_sv = self._vsv_full()
            self.v_sh = self._vsh_full()

    def _dstar_theta(self):
        theta = self.theta
        b_a_2 = 1 - (self.beta_o / self.alpha_o) ** 2

        sin_cos_2 = np.sin(theta) ** 2 * np.cos(theta) ** 2
        sin_4 = np.sin(theta) ** 4

        epsi = self.epsilon
        term_1 = 4 * self.dstar / (b_a_2 * b_a_2) * sin_cos_2
        term_2 = 4 * (b_a_2 + epsi) * epsi / (b_a_2 * b_a_2) * sin_4
        sqrt_term = np.sqrt(1 + term_1 + term_2)
        D_star = 0.5 * b_a_2 * (sqrt_term - 1)
        return D_star

    def _v_p_full(self):
        theta = self.theta
        D_star = self.Dstar_theta

        v_p = self.alpha_o * np.sqrt(1 + self.epsilon * np.sin(theta) ** 2 + D_star)
        return v_p

    def _vsv_full(self):
        theta = self.theta
        b_a_2 = (self.alpha_o / self.beta_o) ** 2
        epsi = self.epsilon
        vsv = self.beta_o * np.sqrt(
            1 + b_a_2 * epsi * np.sin(theta) ** 2 - b_a_2 * self.Dstar_theta
        )
        return vsv

    def _vsh_full(self):
        theta = self.theta
        vsh = self.beta_o * np.sqrt(1 + 2 * self.gamma * np.sin(theta) ** 2)
        return vsh

    def _v_p(self):
        theta = self.theta
        d = self.d_weak
        sin_cos_2 = np.sin(theta) ** 2 * np.cos(theta) ** 2
        sin_4 = np.sin(theta) ** 4
        v_p = self.alpha_o * (1 + d * sin_cos_2 + self.epsilon * sin_4)
        return v_p

    def _vsv(self):
        theta = self.theta
        a_b_2 = (self.alpha_o / self.beta_o) ** 2
        epsi = self.epsilon
        sin_cos_2 = np.sin(theta) ** 2 * np.cos(theta) ** 2
        vsv = self.beta_o * (1 + a_b_2 * (epsi - self.d_weak) * sin_cos_2)
        return vsv

    def _vsh(self):
        theta = self.theta
        vsh = self.beta_o * (1 + self.gamma * np.sin(theta) ** 2)
        return vsh


def _promote_to_1d_array(input):
    return np.array(
        [
            input,
        ]
    )
