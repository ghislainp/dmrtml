

import numpy as np
import dmrtml


def test_dmrt():
    freq = 37.0e9                 # frequency : Hz
    height = np.array([100.0])       # height : m
    temp = np.array([250.0])         # temperature : K
    density = np.array([300.0])      # density : kg.m-3
    dist = False                  # if True => use RAYLEIGH size distribution of particles
    soilp = None

    radius = np.arange(100, 1000, 200)

    tbv = []
    tbh = []

    for r in radius:
        radius = np.array([r]) * 1e-6
        res = dmrtml.run(freq, 64, height, density, radius, temp,
                         tau=dmrtml.NONSTICKY, dist=dist, soilp=soilp)

        tbv.append(res.TbV(53))
        tbh.append(res.TbH(53))

    print(tbv)
    print(tbh)
    assert np.allclose(tbv, [249.52730150104676, 239.57642259920044, 213.63542567993633,
                             181.29329484342307, 151.30705941776813], atol=0.01)
    assert np.allclose(tbh, [236.903725052563, 224.8419978870332, 196.4295944758548,
                             164.71891382020473, 137.04967795060455], atol=0.01)
