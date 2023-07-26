import numpy as np
import pytest
from unyt import unyt_array, unyt_quantity

import pyVBRc.anisotropy.materials as pam
import pyVBRc.sample_data as sd
from pyVBRc.anisotropy._stiffness import TransverseIsotropicStiffness


def test_aligned_inclusions():
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

    with pytest.raises(ValueError, match="calculating v_p requires"):
        _ = m3.v_p
    with pytest.raises(ValueError, match="calculating v_s requires"):
        _ = m3.v_s

    m3.set_density(3000)
    assert m3.v_p > m3.v_s


def test_thomsen_calculator():
    matrix = pam.IsotropicMedium(0.25, 60 * 1e9, "shear", 3300)
    inclusions = pam.IsotropicMedium(0.25, 40 * 1e9, "shear", 2700)
    mixture = pam.AlignedInclusions(0.01)
    mixture.set_material(matrix, inclusions, 0.01)
    C = mixture.get_stiffness_matrix()
    dens = mixture.composite_density

    t_c = pam.ThomsenCalculator(dens, C.stiffness)

    theta = np.array(
        [
            0.0,
        ]
    )
    t_c.set_theta(theta)

    # check that they run
    assert t_c.v_p > 0.0
    assert t_c.v_sh > 0.0
    assert t_c.v_sv > 0.0

    assert t_c._v_p_full() > 0.0
    assert t_c._vsh_full() > 0.0
    assert t_c._vsv_full() > 0.0

    inclusions = pam.IsotropicMedium(0.25, 50 * 1e9, "shear", 2700)
    mixture = pam.AlignedInclusions(10)
    mixture.set_material(matrix, inclusions, 0.01)
    C = mixture.get_stiffness_matrix()
    dens = mixture.composite_density
    t_c2 = pam.ThomsenCalculator(dens, C.stiffness, approx=True)
    t_c2.set_theta(theta)
    diff = np.abs(t_c.v_p - t_c2.v_p) / t_c2.v_p
    assert diff < 0.01  # approximate agrees to within a percent when little aniso


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


def test_transverse_isotropic_stiffness_values():
    # Tando and Weng, Sayers test case
    # AlignedInclusions()
    # epoxy matrix:
    matrix = pam.IsotropicMedium(0.35, 2.76 * 1e9, "youngs")
    # glass fibres
    inclusions = pam.IsotropicMedium(0.2, 72.4 * 1e9, "youngs")
    ar = 2.0
    ai = pam.AlignedInclusions(ar)
    ai.set_material(matrix, inclusions, 0.0)

    # check that when vol fraction is 0, we get the matrix exactly
    E11, E22, mu12, mu23, nu12, K23 = ai.get_moduli()

    assert E11 == E22
    assert E11 == matrix.youngs_modulus
    assert mu12 == mu23
    assert mu12 == matrix.shear_modulus
    assert nu12 == matrix.poisson_ratio

    # nu = matrix.poisson_ratio
    # Kiso = matrix.plain_strain_bulk * ()
    assert np.allclose(K23, matrix.plain_strain_bulk)

    # check some stiffness matrix values
    stiffness = ai.get_stiffness_matrix()
    C = stiffness.stiffness
    assert C[3, 3] == C[4, 4]
    assert C[3, 3] == C[5, 5]
    assert C[0, 0] == C[2, 2]
    assert C[0, 0] == C[1, 1]
    # assert C[2, 2] == E11 + 2 * C[2, 0]
    assert C[0, 1] == C[1, 0]
    assert C[0, 2] == C[2, 0]
    assert C[0, 2] == C[2, 1]
    assert C[0, 2] == C[1, 2]

    # vol fracgtion of 0 should reduce to isotropic case
    # in which case
    assert C[0, 1] == C[0, 2]

    # and anisotropy factors should go to 0

    assert np.allclose(stiffness.a_1, 0.0, atol=1e-6)
    assert np.allclose(stiffness.a_2, 0.0, atol=1e-6)
    assert np.allclose(stiffness.a_3, 0.0, atol=1e-6)


def _get_offset_vbr_structs(file):
    vbr_1 = sd.load_sample_structure(file)
    vbr_2 = sd.load_sample_structure(file)

    if isinstance(vbr_2.output.elastic.anharmonic.Gu, unyt_array):
        offset = unyt_quantity(10, "GPa")
    else:
        offset = 10 * 1e9

    vbr_2.output.elastic.anharmonic.Gu += offset
    vbr_2.output.anelastic.andrade_psp.M += offset

    shear1 = vbr_2.output.elastic.anharmonic.Gu
    shear2 = vbr_1.output.elastic.anharmonic.Gu
    assert np.all(shear1 != shear2)

    return vbr_1, vbr_2


_VBRfiles = [
    "VBRc_sample_LUT.mat",
    "VBRc_sample_LUT_R2021a.mat",
    "VBRc_sample_LUT_v1pt0pt0.mat",
]


@pytest.mark.parametrize("file", _VBRfiles)
def test_load_from_vbrc(file):
    # set_materials_from_vbrc_structures
    vbr_1, vbr_2 = _get_offset_vbr_structs(file)
    input_shape = vbr_2.output.elastic.anharmonic.Gu.shape

    ar = 0.1
    mixture = pam.AlignedInclusions(ar)

    mixture.set_materials_from_vbrc_structures(
        vbr_1,
        vbr_2,
        ["elastic", "anharmonic", "Gu"],
        0.25,
    )

    shear1 = mixture.matrix_material.shear_modulus
    shear2 = mixture.inclusion_material.shear_modulus
    assert np.all(shear1 != shear2)

    _ = mixture.get_stiffness_matrix()
    vp, vsv, vsh = mixture.velocities(0.0)
    assert all([v.shape == input_shape for v in (vp, vsv, vsh)])

    # try other shapes, different combos
    mixture.set_materials_from_vbrc_structures(
        vbr_1,
        vbr_2,
        ["elastic", "anharmonic", "Gu"],
        np.full(input_shape, 0.25),
    )
    _ = mixture.velocities(0.0)
    _ = mixture.velocities(np.full(input_shape, 0.5))

    mixture.set_materials_from_vbrc_structures(
        vbr_1,
        vbr_2,
        ["elastic", "anharmonic", "Gu"],
        0.25,
    )
    _ = mixture.velocities(np.full(input_shape, 0.5))

    with pytest.raises(RuntimeError, match="When material properties"):
        mixture.set_materials_from_vbrc_structures(
            vbr_1,
            vbr_2,
            ["elastic", "anharmonic", "Gu"],
            np.full((3, 3), 0.25),
        )

    mixture.set_materials_from_vbrc_structures(
        vbr_1,
        vbr_2,
        ["elastic", "anharmonic", "Gu"],
        0.25,
    )
    with pytest.raises(RuntimeError, match="When material properties"):
        _ = mixture.velocities(np.full((4, 5), 0.5))

    mixture.set_materials_from_vbrc_structures(
        vbr_1,
        vbr_2,
        ["anelastic", "andrade_psp", "M"],
        0.25,
        ifreq=0,
    )
    _ = mixture.velocities(np.full(input_shape, 0.5))

    with pytest.raises(ValueError, match="anelastic.andrade_psp.M is freq"):
        mixture.set_materials_from_vbrc_structures(
            vbr_1,
            vbr_2,
            ["anelastic", "andrade_psp", "M"],
            0.25,
        )

    # make it fail (do this last)
    vbr_1.input.SV.rho = vbr_1.input.SV.rho[0:2, 0:3]
    with pytest.raises(RuntimeError, match="The matrix and inclusion"):
        mixture.set_materials_from_vbrc_structures(
            vbr_1, vbr_2, ["anelastic", "andrade_psp", "M"], 0.25, ifreq=0
        )
