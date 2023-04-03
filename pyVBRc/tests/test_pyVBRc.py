import os

import numpy as np
import pytest
from unyt import unyt_array

from pyVBRc.vbrc_structure import VBRCstruct

_VBRfiles = [
    "VBRc_sample_LUT.mat",
    "VBRc_sample_LUT_R2021a.mat",
    "VBRc_sample_LUT_v1pt0pt0.mat",
]


@pytest.mark.parametrize("fname", _VBRfiles)
def test_mat_load(fname):
    sample_file = os.path.join("pyVBRc/sample_data/", fname)

    vbr = VBRCstruct(sample_file, lut_dimensions=["T_K", "phi", "dg_um"])

    assert vbr.input.SV.T_K.ndim == 3


@pytest.mark.parametrize("fname", _VBRfiles)
def test_interpolator(fname):

    sample_file = os.path.join("pyVBRc/sample_data/", fname)

    vbr = VBRCstruct(sample_file, lut_dimensions=["T_K", "phi", "dg_um"])

    interp = vbr.interpolator(
        ("anelastic", "andrade_psp", "V"), 0, log_vars=["phi", "dg_um"]
    )

    # evaluate at one arbitrary point
    # (T in K, log(phi), log(dg in micrometer))
    target = (1333 + 273.0, np.log10(0.0012), np.log10(1131))
    Vs_interp = interp(target)
    assert Vs_interp > 0

    # resample at a higher resolution
    T = vbr.input.SV.T_K[:, 0, 0]
    nT = len(T)
    T_targets = np.linspace(T.min(), T.max(), nT * 2)
    phival = vbr.input.SV.phi.min()
    dgval = vbr.input.SV.dg_um.min()
    phi_targets = np.full(T_targets.shape, np.log10(phival))
    dg_targets = np.full(T_targets.shape, np.log10(dgval))
    targets = np.column_stack((T_targets, phi_targets, dg_targets))
    Vs_interp = interp(targets)
    assert len(Vs_interp) == nT * 2


def test_units():
    sample_file = "pyVBRc/sample_data/VBRc_sample_LUT_v1pt0pt0.mat"

    vbr = VBRCstruct(sample_file)
    assert isinstance(vbr.input.SV.T_K, unyt_array)
    assert isinstance(vbr.output.elastic.anh_poro.Gu, unyt_array)
    assert isinstance(vbr.output.anelastic.eburgers_psp.Q, unyt_array)

    vbr = VBRCstruct(sample_file, attach_units=False)
    assert isinstance(vbr.input.SV.T_K, unyt_array) is False


@pytest.mark.parametrize("fname", _VBRfiles)
def test_interpolator_further(fname):
    sample_file = os.path.join("pyVBRc/sample_data/", fname)

    # vbr = VBRCstruct(sample_file, lut_dimensions=["T_K", "phi", "dg_um"])
    vbr = VBRCstruct(sample_file)
    with pytest.raises(RuntimeError, match="Please call set_lut_dimensions"):
        _ = vbr.interpolator(("anelastic", "andrade_psp", "V"), 0)

    vbr.set_lut_dimensions(["T_K", "phi", "dg_um"])
    with pytest.raises(RuntimeError, match="the data_field selection"):
        _ = vbr.interpolator(("anelastic", "andrade_psp"), 0)

    with pytest.raises(ValueError, match="must supply a frequency dimension "):
        _ = vbr.interpolator(("anelastic", "andrade_psp", "V"))

    _ = vbr.interpolator(
        ("anelastic", "andrade_psp", "V"),
        0,
        log_vars=[
            "V",
        ],
    )

    lut_dims = vbr._get_lut_dimensions()
    assert len(lut_dims) == 3

    with pytest.raises(ValueError, match="the number of fields in lut_dimensions"):
        vbr.set_lut_dimensions(
            [
                "T_K",
            ]
        )
