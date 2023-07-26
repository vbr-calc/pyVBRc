from pathlib import PosixPath
from typing import List, Tuple, Union

import numpy as np
from packaging.version import Version
from scipy.interpolate import RegularGridInterpolator
from scipy.io import loadmat
from unyt import unyt_array
from unyt.unit_object import UnitParseError

from pyVBRc.logging import pyvbrc_log


class VBRCstruct:
    def __init__(
        self,
        filename: Union[str, PosixPath],
        vbr_name: str = "VBR",
        lut_dimensions: List[str] = None,
        attach_units: bool = True,
    ):
        """
        A data class to load a VBRc structure saved by the VBRc

        Parameters
        ----------
        filename: str | PosixPath
            the file to load VBRc structure from
        vbr_name: str = "VBR"
            the variable name of the top-level VBR structure in the file.
        lut_dimensions: list[str] = None
            the look-up table (lut) dimensions. The variables listed here
            are taken to be the primary state variable dimensions. Used
            when constructing interpolaters.
        attach_units: bool = True
            if True (default), will attempt a recursive walk through the VBRc
            structure and convert the arrays to unyt arrays.
        """

        self.filename = filename

        raw_data = loadmat(
            filename,
            uint16_codec="ascii",
            struct_as_record=False,
            squeeze_me=True,
        )
        self.input = getattr(raw_data[vbr_name], "in")
        self.output = getattr(raw_data[vbr_name], "out")
        vbrc_raw_version = getattr(raw_data[vbr_name], "version_used", None)
        self.vbrc_version = None
        if vbrc_raw_version:
            self.vbrc_version = Version(vbrc_raw_version.version)

        self.lut_dimensions = None
        if lut_dimensions is not None:
            self.set_lut_dimensions(lut_dimensions)

        if attach_units and self.vbrc_version and self.vbrc_version >= Version("1.0.0"):
            self._attach_units()

        pyvbrc_log.info(f"{filename} loaded.")

    def _attach_units(self):

        # first the inputs
        self.input = _recursive_unitfication(self.input, "input")

        # now the outputs
        self.output = _recursive_unitfication(self.output, "output")

    def set_lut_dimensions(self, lut_dimensions: List[str]):
        # check that number of dimensions equals number of dimensions
        # in the state variables before storing the lut_dimensions
        nstate_dim = self.input.SV.T_K.ndim
        ndims = len(lut_dimensions)
        if nstate_dim != ndims:
            raise ValueError(
                f"the number of fields in lut_dimensions ({ndims}) must match"
                f"the dimensionality of the state variables ({nstate_dim})"
            )

        self.lut_dimensions = lut_dimensions

    def interpolator(
        self,
        data_field: Tuple[str],
        i_frequency: int = None,
        log_vars: List[str] = None,
        method: Union[str, None] = "linear",
        bounds_error: Union[bool, None] = True,
        fill_value: Union[float, None] = np.nan,
    ) -> RegularGridInterpolator:
        """

        Parameters
        ----------
        data_field: Tuple[str]
            the field to interpolate. The tuple should point to nested attributes
            in the VBRCstruct.output. For example, ('anelastic', 'andrade_psp', 'V'),
            will extract VBRCstruct.output.anelastic.andrade_psp.V as the data
            field to hand off to the interopolator.
        i_frequency: int
            if the data_field is frequency-dependent, you must supply the frequency
            index to use.
        log_vars: list[str]
            any variable names appearing in this list will be passed through
            np.log10 before handing to the interpolator. Rmember that you will
            then need to supply the target poitns in log-space for these
            variables as well.

        Remaining arguments are passed to scipy.interpolate.RegularGridInterpolator:

        method: str | None = "linear"
        bounds_error: bool | None = True
        fill_value: float | None = np.nan

        Returns
        -------
        scipy.interpolate.RegularGridInterpolator instance

        """

        if log_vars is None:
            log_vars = []

        if self.lut_dimensions is None:
            raise RuntimeError(
                "Please call set_lut_dimensions prior to building an interpolator."
            )

        data = self._get_nested_output_field(data_field)

        if isinstance(data, np.ndarray) is False:
            raise RuntimeError("the data_field selection did not return an array.")

        if data.ndim - len(self.lut_dimensions) == 1:
            # frequency dependent
            if i_frequency is None:
                msg = (
                    f"must supply a frequency dimension for frequency "
                    f"dependent field {data_field}"
                )
                raise ValueError(msg)
            data = data[..., i_frequency]

        if data_field[-1] in log_vars:
            data = np.log10(data)

        pts = self._get_lut_dimensions(log_vars)
        return RegularGridInterpolator(
            pts, data, method=method, bounds_error=bounds_error, fill_value=fill_value
        )

    def _get_lut_dimensions(self, log_vars: list[str] = None) -> List[np.ndarray]:
        # extract the look up table dimension arrays, taking the log if necessary
        if log_vars is None:
            log_vars = []

        # build the point grid
        points = []
        for slice_dim, var in enumerate(self.lut_dimensions):
            # pull out the full n-d state variable array
            pts = getattr(self.input.SV, var)

            # fix the index to slice except along the current dimensions
            indxs = [0] * pts.ndim
            indxs[slice_dim] = slice(None)
            indxs = tuple(indxs)
            pts = pts[indxs]
            assert pts.ndim == 1

            if var in log_vars:
                pts = np.log10(pts)
            points.append(pts)
        return points

    def _get_nested_output_field(self, the_field: Tuple[str]):
        data = self.output
        for df in the_field:
            data = getattr(data, df)
        return data


def _recursive_unitfication(vbrc_sub_struct, struct_name: str):
    # recursively walks a VBRc structure, replacing all arrays with unyt arrays
    # if there is a units field at the same level as the arrays in the structure.
    if hasattr(vbrc_sub_struct, "_fieldnames"):
        for field in vbrc_sub_struct._fieldnames:
            field_value = getattr(vbrc_sub_struct, field)
            if isinstance(field_value, np.ndarray):
                # it is an array! check if units are known
                if hasattr(vbrc_sub_struct, "units"):
                    if hasattr(vbrc_sub_struct.units, field) is False:
                        pyvbrc_log.warning(f"{field} in {struct_name} has no units")
                    units = getattr(vbrc_sub_struct.units, field, "")
                    if isinstance(units, np.ndarray) and len(units) == 0:
                        units = ""
                    try:
                        new = unyt_array(field_value, units)
                    except UnitParseError:
                        pyvbrc_log.warning(
                            f"{field} in {struct_name} has "
                            f"unsupported units ({units}), using nondimensional"
                        )
                        new = unyt_array(field_value, "")
                    setattr(vbrc_sub_struct, field, new)
            elif hasattr(field_value, "_fieldnames"):
                # need to go deeper:
                new_value = _recursive_unitfication(field_value, field)
                setattr(vbrc_sub_struct, field, new_value)

    return vbrc_sub_struct
