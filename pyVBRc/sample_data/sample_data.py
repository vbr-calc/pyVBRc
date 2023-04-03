import os
from importlib import resources
from pathlib import PosixPath
from typing import List, Optional

from pyVBRc.logging import pyvbrc_log
from pyVBRc.vbrc_structure import VBRCstruct


def load_sample_structure(
    file: str, lut_dimensions: Optional[List[str]] = None
) -> VBRCstruct:
    fullfile = _sample_dir / file
    lut_dims = lut_dimensions or ["T_K", "phi", "dg_um"]
    pyvbrc_log.info(f"Loading sample file {file}.")
    vbr = VBRCstruct(fullfile, lut_dimensions=lut_dims)
    return vbr


_sample_dir = resources.files("pyVBRc") / "sample_data"


def list_sample_files() -> List[str]:
    return [f for f in os.listdir(_sample_dir) if f.endswith(".mat")]


def get_sample_filename(file: str) -> PosixPath:
    return _sample_dir / file
