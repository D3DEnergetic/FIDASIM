title: Advanced Settings

This page documents FIDASIMs advanced settings. All settings can be found in the Namelist File.

[TOC]

---

# Setting Random Number Generator Seed
```
seed = -1 ! -1 = Default (pick random seed)
```

# Verbose Output
```
verbose = 1
```

# Loading Pre-computed Neutrals File
```
load_neutrals = 1 ! Default = 0
neutral_file = 'neutrals.h5'
```

# Output Individual Stark Components
```
stark_components = 1 ! Default = 0
```

# Control Finite Larmor Radius Correction Order
```
flr = 2 !! 0 => No correction, 1 => First order correction, 2 => Second Order correction (Default)
```

# Output NPA Particles
```
calc_npa = 2
calc_pnpa = 2
```

# Spatially Averaged FIDA Velocity-space Weight Functions
```
calc_fida_wght = 2
```

# Output Extra NPA Velocity-space Weight Function Variables
```
calc_npa_wght = 2
```
