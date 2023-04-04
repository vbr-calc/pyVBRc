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
