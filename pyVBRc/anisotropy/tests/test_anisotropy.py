import numpy as np

from pyVBRc.anisotropy.materials import AlignedInclusions, IsotropicMedium


def test_anisotropy():
    # initial test, just try some calls that hit a bunch
    matrix = IsotropicMedium(0.25, 60 * 1e9, "shear", 3300)
    inclusions = IsotropicMedium(0.25, 40 * 1e9, "shear", 2700)
    degs = np.linspace(0, 180, 10)
    theta = degs * np.pi / 180.0

    vol_frac = 0.05
    ars = np.linspace(0.0001, 0.2, 10)
    for case_i, ar in enumerate(ars):
        mixture = AlignedInclusions(ar)
        mixture.set_material(matrix, inclusions, vol_frac)
        v_p, v_sv, v_sh = mixture.velocities(theta)
