---
layout: default
title: Tutorial
---

In this tutorial we will set-up a simulation for a 2D fault with heterogeneous fault properties, using *rate-and-state friction* (RSF) as the underlying fault rheology.

Create a wrapper file (e.g. `single_asperity.py`) with the following imports:

```Python
import os
import sys

# Modify the line below to point to the QDYN source directory
qdyn_src_dir = "/path/to/qdyn/src"

sys.path.append(qdyn_src_dir)
from pyqdyn import qdyn
```

Next, we create two dictionaries with simulation parameters and settings. The first dictionary (`set_dict`) contains general settings, the second (`set_dict_RSF`) contains base parameters defining the fault rheology. The base parameters save as default values when the mesh is generated, and we will overwrite them later on when heterogeneity is introduced.

```Python
# General simulation parameters
set_dict = {
    "FAULT_TYPE": 1,
    "ACC": 1e-10,
    "SOLVER": 2,
    "MU": 3e10,
    "TMAX": 20,
    "DTTRY": 1e-6,
    "MESHDIM": 1,
    "NTOUT": 10000,
    "VS": 3000,
    "SIGMA": 5e6,
    "L": 1,
    "W": 1,
    "FEAT_STRESS_COUPL": 0,
    "FEAT_TP": 0,
    "FEAT_LOCALISATION": 0,
    "D": 0,
    "HD": 0,
    "FRICTION_MODEL": "RSF",
}

# RSF parameters
set_dict_RSF = {
    "RNS_LAW": 0,
    "THETA_LAW": 1,
    "A": 0.001,
    "B": 0.0015,
    "DC": 1e-5,
    "V1": 1e-2,
    "V2": 1e-7,
    "MU_SS": 0.6,
    "V_SS": 1e-6,
    "CO": 0,
}
```

The rheological parameters are then merged into the main settings dictionary:
```Python
set_dict["SET_DICT_RSF"] = set_dict_RSF
```

We are now ready to initiate the QDYN API, upload our simulation settings, and render the fault mesh:
```Python
p = qdyn()
p.settings(set_dict)
p.render_mesh()
```
The class instance `p` now contains a dictionary called `mesh_dict` that stores the parameter values at each mesh node, initially set to the default values defined in `set_dict`. To introduce the heterogeneity, we modify the value of `B` in a central portion of the mesh:
```Python
start = 100
stop = 800
p.mesh_dict["B"][start:stop] = 0.02
```