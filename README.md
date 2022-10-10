# pyVBRc

pyVBRc is a python package for working with output from the [Very Broadband Rheology Calculator](https://github.com/vbr-calc/vbr). pyVBRc is currently (very) experimental and is likely to change drastically between early versions, but feel free to use it.

## installation


(not working yet)
```
pip install pyVBRc
```

Currently working:

```
pip install git+https://github.com/vbr-calc/pyVBRc
```

## usage

As of now, you can use `pyVBRc` to:

* load a VBR structure
* interpolate between calculated values in the VBR structure

### loading a VBR structure

The primary gateway into working with VBRc data in the `VBRstruct`. Note that the following examples rely on a tiny VBR output file, `VBRc_sample_LUT.mat`, included in pyVBRc to allow simple testing of functionality. For any real application, you should generate your own VBR structure (in MATLAB or Octave).

The most basic usage of the `VBRCstruct` object is to load a `.mat` file:

```python
from pyVBRc.vbrc_structure import VBRCstruct
from pyVBRc.sample_data import get_sample_filename
file = get_sample_filename("VBRc_sample_LUT.mat")
vbr = VBRCstruct(file)
```

The `vbr` object will mirror the standard VBR structure (see [here](https://vbr-calc.github.io/vbr/gettingstarted/)) with the exception that the `in` and `out` structures of the matlab VBR structure are replaced with `input` and `output` in the python class.

Once loaded, you can access the data just as in matlab. This `.mat` file contains a `VBR` structure that was run with varying temperature (`T_K`), melt fraction (`phi`) and grain size(`dg_um`) at two frequencies. The following example plots the temperature dependence of attenuation for the Andrade pseuodo-period scaling at a fixed grain size, melt fraction and frequency:

```
import matplotlib.pyplot as plt
T_K = vbr.input.SV.T_K[:, 0, 0]
Qinv = vbr.output.anelastic.andrade_psp.Qinv[:, 0, 0, 0]

plt.semilogy(T_K, Qinv, '.k')
plt.xlabel("T [Kelvin]")
plt.ylabel("Q^{-1}, andrade pseudo-period")
plt.show()
```

![](https://raw.githubusercontent.com/vbr-calc/pyVBRc/main/examples/andrade_psp_T_dep.png)

## interpolating

In some situations, it may be useful to interpolate VBRc results at arbitrary points between calculated points. The sample file,`"VBRc_sample_LUT.mat"` was calculated for a very coarse grid of temperature (`T_K`), melt fraction (`phi`) and grain size(`dg_um`) at two frequencies to serve as a simple look-up table (LUT). You can provide these variable names to the `lut_dimensions` keyword argument when initializing the `VBRCstruct`:

```python
from pyVBRc.vbrc_structure import VBRCstruct
from pyVBRc.sample_data import get_sample_filename

file = get_sample_filename("VBRc_sample_LUT.mat")
vbr = VBRCstruct(file, lut_dimensions=["T_K", "phi", "dg_um"])
```

Now, you can build a grid interpolator for a single method-array and a single frequency index with the `interpolator` method. This method also accepts the `log_vars` keyword argument that you can use to specify LUT fields that you want to interpolate in log-space:

```python
import numpy as np
interp = vbr.interpolator(
    ("anelastic", "andrade_psp", "V"), 0, log_vars=["phi", "dg_um"]
)

# evaluate at one arbitrary point
# (T in K, log(phi), log(dg in micrometer))
target = (1333 + 273.0, np.log10(0.0012), np.log10(1131))

Vs_interp = interp(target)
print(Vs_interp)

```

Note that the array to interpolate is specified as a tuple of output fieldnames: `("anelastic", "andrade_psp", "V")` corresponds to accessing `VBR.output.anelastic.andrade_psp.V`. To interpolate at more than one point, you can provide an array of points:

```python
import matplotlib.pyplot as plt
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

# compare to actual Vs along curve
actual_Vs = vbr.output.anelastic.andrade_psp.V[:, 0, 0, 0]

plt.plot(T_targets, Vs_interp, label="interpolated output", marker=".")
plt.plot(T, actual_Vs, ".k", label="actual VBRc output", markersize=12)
plt.legend()
plt.xlabel("T [K]")
plt.ylabel("Vs [m/s]")
plt.show()
```

which yields

![](https://raw.githubusercontent.com/vbr-calc/pyVBRc/main/examples/interpolate_example.png)


**NOTE** that the validity of any interpolation will depend on the underlying VBR structure that you have built.

## Getting help

## Contributing
### installing from source
### running the test suite


