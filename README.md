# pyVBRc

pyVBRc is a python package for working with output from the [Very Broadband Rheology Calculator](https://github.com/vbr-calc/vbr). pyVBRc is currently (very) experimental and is likely to change drastically between early versions, but feel free to use it.

## installation


## basic usage

The primary gateway into working with VBRc data in the `VBRstruct`. Note that the following examples rely on a tiny VBR output file, `VBRc_sample_LUT.mat`, included in pyVBRc to allow simple testing of functionality. For any real application, you should generate your own VBR structure (in MATLAB or Octave).

The most basic usage of the `VBRCstruct` object is to load a `.mat` file:

```python
from pyVBRc.vbrc_structure import VBRCstruct
from pyVBRc.sample_data import get_sample_filename
file = get_sample_filename("VBRc_sample_LUT.mat")
vbr = VBRCstruct(file)
```

The `vbr` object will mirror the standard VBR structure (see [here](https://vbr-calc.github.io/vbr/gettingstarted/)) with the exception that the `in` and `out` structures of the matlab VBR structure are replaced with `input` and `output` in the python class.

Once loaded, you can access the data just as in matlab. This `.mat` file contains a `VBR` structure that was run with varying temperature (`T_K`), melt fraction (`phi`) and grain size(`dg_um`) at two frequencies. So to pull out the temperature dependence of intrinsic attenuation at a given melt fraction, grain size and frequency:

```
import matplotlib.pyplot as plt
T_K = vbr.input.SV.T_K[:, 0, 0]
Qinv = vbr.output.anelastic.andrade_psp.Qinv[:, 0, 0, 0]

plt.semilogy(T_K, Qinv, '.k')
plt.xlabel("T [Kelvin]")
plt.ylabel("Q^{-1}, andrade pseudo-period")
plt.show()
```

will result in the following figure

![](./examples/andrade_pspa_T_dep.png)

