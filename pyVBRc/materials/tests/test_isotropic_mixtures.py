import numpy as np
import pytest
from unyt import unyt_array, unyt_quantity

import pyVBRc.sample_data as sd
from pyVBRc.materials import materials


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
def test_load_isotropic_from_vbrc(file):
    # set_materials_from_vbrc_structures
    vbr_1, vbr_2 = _get_offset_vbr_structs(file)
    sml = ("elastic", "anharmonic", "Gu")
    m1 = materials.load_isotropic_medium(vbr_1, shear_modulus_location=sml)
    m2 = materials.load_isotropic_medium(vbr_2, shear_modulus_location=sml)

    props = (0.5, 0.5)
    m = materials.IsotropicMixture([m1, m2], props)
    assert np.all(m.shear_velocity() < m.compressional_velocity())

    sml = ("anelastic", "andrade_psp", "M")
    m1 = materials.load_isotropic_medium(vbr_1, shear_modulus_location=sml, ifreq=0)
    m2 = materials.load_isotropic_medium(vbr_2, shear_modulus_location=sml, ifreq=0)

    props = (0.5, 0.5)
    m = materials.IsotropicMixture([m1, m2], props)
    assert np.all(m.shear_velocity() < m.compressional_velocity())

    with pytest.raises(ValueError, match=".".join(sml) + " is frequency"):
        _ = materials.load_isotropic_medium(vbr_2, shear_modulus_location=sml)


@pytest.mark.parametrize("method", ("voigt", "reuss", "voigt-reuss"))
def test_isotropic_mixture(method):
    m1 = materials.IsotropicMedium(0.25, 80 * 1e9, modulus_type="shear", density=3000)
    m2 = materials.IsotropicMedium(0.25, 50 * 1e9, modulus_type="shear", density=2000)

    props = (0.5, 0.5)
    m = materials.IsotropicMixture([m1, m2], props)

    assert m.density() == 2500

    vs = m.shear_velocity(method=method)
    vp = m.compressional_velocity(method=method)
    assert np.all(vs < vp)


def test_mixtures_from_arrays():
    shp = (3, 4)
    nu = np.full(shp, 0.25)
    G1 = np.full(shp, 80 * 1e9)
    G2 = np.full(shp, 40 * 1e9)
    d1 = np.full(shp, 3000)
    d2 = np.full(shp, 2000)
    m1 = materials.IsotropicMedium(nu, G1, modulus_type="shear", density=d1)
    m2 = materials.IsotropicMedium(nu, G2, modulus_type="shear", density=d2)

    props = (0.5, 0.5)
    m = materials.IsotropicMixture([m1, m2], props)
    assert np.all(m.density() == 2500)
    assert np.all(m.shear_modulus(method="voigt") == 60 * 1e9)


def test_some_errors():
    m1 = materials.IsotropicMedium(0.25, 80 * 1e9, modulus_type="shear", density=3000)
    m2 = materials.IsotropicMedium(0.25, 50 * 1e9, modulus_type="shear", density=2000)

    props = (0.5, 0.5, 0.1)
    with pytest.raises(
        RuntimeError, match="length of materials does not match length of proportions"
    ):
        _ = materials.IsotropicMixture([m1, m2], props)

    props = (0.5, 0.25)
    with pytest.raises(RuntimeError, match="proportions should sum to 1"):
        _ = materials.IsotropicMixture([m1, m2], props)

    props = (0.25, 0.75)
    m = materials.IsotropicMixture([m1, m2], props)

    with pytest.raises(RuntimeError, match="unexpected shapes"):
        _ = m._array_mult_and_sum(np.random.random((3,)), np.random.random((2, 3)))
