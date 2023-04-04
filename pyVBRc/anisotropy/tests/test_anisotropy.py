import numpy as np
import pytest

import pyVBRc.anisotropy.materials as pam
from pyVBRc.anisotropy._stiffness import TransverseIsotropicStiffness


def test_anisotropy():
    # initial test, just try some calls that hit a bunch
    matrix = pam.IsotropicMedium(0.25, 60 * 1e9, "shear", 3300)
    inclusions = pam.IsotropicMedium(0.25, 40 * 1e9, "shear", 2700)
    degs = np.linspace(0, 180, 10)
    theta = degs * np.pi / 180.0

    vol_frac = 0.05
    ars = (0.01, 0.1, 1.0, 10.0, 100.0)
    for case_i, ar in enumerate(ars):
        mixture = pam.AlignedInclusions(ar)
        mixture.set_material(matrix, inclusions, vol_frac)
        _ = mixture.velocities(theta)


def test_aligned_inclusion_errors():

    with pytest.raises(ValueError, match="Detected a mix"):
        _ = pam.AlignedInclusions((0.01, 10))

    ai = pam.AlignedInclusions(0.01)
    with pytest.raises(RuntimeError, match="You must set"):
        _ = ai.get_stiffness_matrix()

    with pytest.raises(ValueError, match="calculating velocities"):
        _ = ai.velocities((0.0,))

    matrix = pam.IsotropicMedium(0.25, 60 * 1e9, "shear")
    inclusions = pam.IsotropicMedium(0.25, 40 * 1e9, "shear")
    ai.set_material(matrix, inclusions, 0.05)
    assert ai.composite_density is None


def test_isotropic_medium():
    m1 = pam.IsotropicMedium(0.25, 60 * 1e9, "shear")
    m2 = pam.IsotropicMedium(0.25, 80 * 1e9, "shear")

    assert m1.bulk_modulus < m2.bulk_modulus
    assert m1.pwave_effective_modulus < m2.pwave_effective_modulus

    m3 = pam.IsotropicMedium(0.25, m1.bulk_modulus, "bulk")
    assert m3.shear_modulus == m1.shear_modulus

    m3 = pam.IsotropicMedium(0.25, m1.youngs_modulus, "youngs")
    assert m3.shear_modulus == m1.shear_modulus


def test_thomsen_calculator():
    matrix = pam.IsotropicMedium(0.25, 60 * 1e9, "shear", 3300)
    inclusions = pam.IsotropicMedium(0.25, 40 * 1e9, "shear", 2700)
    mixture = pam.AlignedInclusions(0.01)
    mixture.set_material(matrix, inclusions, 0.01)
    C = mixture.get_stiffness_matrix()
    dens = mixture.composite_density

    t_c = pam.ThomsenCalculator(dens, C.stiffness)

    t_c.set_theta(
        np.array(
            [
                0.0,
            ]
        )
    )

    # check that they run
    assert t_c.v_p > 0.0
    assert t_c.v_sh > 0.0
    assert t_c.v_sv > 0.0

    assert t_c._v_p_full() > 0.0
    assert t_c._vsh_full() > 0.0
    assert t_c._vsv_full() > 0.0


def test_transverse_isotropic_stiffness():

    e_l = np.ones((5,))
    e_t = np.ones((5,))
    g_in = np.ones((5,))
    g_out = np.ones((5,))
    nu = np.full((5,), 0.25)
    K = np.ones((5,))
    tis = TransverseIsotropicStiffness(e_l, e_t, g_in, g_out, nu, K)
    assert tis.stiffness.shape == (6, 6, 5)

    K = np.ones((2, 2))
    with pytest.raises(ValueError, match="All input arrays must be 1D"):
        _ = TransverseIsotropicStiffness(e_l, e_t, g_in, g_out, nu, K)

    K = np.ones((6,))
    with pytest.raises(ValueError, match="All input arrays must have the same"):
        _ = TransverseIsotropicStiffness(e_l, e_t, g_in, g_out, nu, K)
