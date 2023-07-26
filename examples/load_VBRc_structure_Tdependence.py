import matplotlib.pyplot as plt

from pyVBRc.sample_data import get_sample_filename
from pyVBRc.vbrc_structure import VBRCstruct

file = get_sample_filename("VBRc_sample_LUT.mat")
vbr = VBRCstruct(file)
T_K = vbr.input.SV.T_K[:, 0, 0]
Qinv = vbr.output.anelastic.andrade_psp.Qinv[:, 0, 0, 0]

plt.semilogy(T_K, Qinv, ".k")
plt.xlabel("T [Kelvin]")
plt.ylabel("Q$^{-1}$, andrade pseudo-period")
plt.savefig("andrade_psp_T_dep.png")
