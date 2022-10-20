from pyVBRc.sample_data import (
    load_sample_structure,
    list_sample_files,
    get_sample_filename,
)
import pytest
import os

_VBRfiles = [
    "VBRc_sample_LUT.mat",
    "VBRc_sample_LUT_R2021a.mat",
]


@pytest.mark.parametrize("file", _VBRfiles)
def test_load_sample_structure(file):
    vbr = load_sample_structure(file)
    assert vbr.input.SV.T_K.ndim == 3


def test_list_sample_data():
    available_data = list_sample_files()
    for f in _VBRfiles:
        assert f in available_data
        vbr = load_sample_structure(f)
        assert vbr.input.SV.T_K.ndim == 3


def test_get_sample_filename():
    for f in _VBRfiles:
        fullfi = get_sample_filename(f)
        assert os.path.isfile(fullfi)