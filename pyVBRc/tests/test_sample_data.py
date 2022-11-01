import os

import pytest
from packaging.version import Version

import pyVBRc.sample_data as sd

_VBRfiles = [
    "VBRc_sample_LUT.mat",
    "VBRc_sample_LUT_R2021a.mat",
    "VBRc_sample_LUT_v1pt0pt0.mat",
]


@pytest.mark.parametrize("file", _VBRfiles)
def test_load_sample_structure(file):
    vbr = sd.load_sample_structure(file)
    assert vbr.input.SV.T_K.ndim == 3


def test_list_sample_data():
    available_data = sd.list_sample_files()
    for f in _VBRfiles:
        assert f in available_data
        vbr = sd.load_sample_structure(f)
        assert vbr.input.SV.T_K.ndim == 3


def test_get_sample_filename():
    for f in _VBRfiles:
        fullfi = sd.get_sample_filename(f)
        assert os.path.isfile(fullfi)


def test_version():
    vbr = sd.load_sample_structure("VBRc_sample_LUT_v1pt0pt0.mat")
    assert vbr.vbrc_version == Version("1.0.0")


def test_missing_version():
    vbr = sd.load_sample_structure("VBRc_sample_LUT_R2021a.mat")
    assert vbr.vbrc_version is None
