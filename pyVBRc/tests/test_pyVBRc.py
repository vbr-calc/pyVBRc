import os

import numpy as np
import pytest

from pyVBRc.vbrc_structure import VBRCstruct

_VBRfiles = [
    "VBRc_sample_LUT.mat",
    "VBRc_sample_LUT_R2021a.mat",
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
