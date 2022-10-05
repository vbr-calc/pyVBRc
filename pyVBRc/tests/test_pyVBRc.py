#!/usr/bin/env python

"""Tests for `pyVBRc` package."""

from pyVBRc.vbrc_structure import VBRCstruct


def test_mat_load():
    sample_file = "pyVBRc/sample_data/VBRc_sample_LUT.mat"

    vbr = VBRCstruct(sample_file, lut_dimensions=["T_K", "phi", "dg_um"])

    assert vbr.input.SV.T_K.ndim == 3
