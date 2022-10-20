from pathlib import PosixPath
from typing import List, Tuple, Union

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.io import loadmat


class VBRCstruct:
    def __init__(
        self,
        filename: Union[str, PosixPath],
        vbr_name: str = "VBR",
        lut_dimensions: List[str] = None,
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
        """

        self.filename = filename

        raw_data = loadmat(
            filename,
            uint16_codec="ascii",  # why this?
            struct_as_record=False,
            squeeze_me=True,
        )
        self.input = getattr(raw_data[vbr_name], "in")
        self.output = getattr(raw_data[vbr_name], "out")

        self.lut_dimensions = None
        if lut_dimensions is not None:
            self.set_lut_dimensions(lut_dimensions)

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

        data = self.output
        for df in data_field:
            data = getattr(data, df)

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
