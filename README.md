# Corrugated Web Structural Optimization

This repository contains MATLAB scripts for topology optimization and buckling analysis of corrugated steel webs and flat plates, developed for the research project "Structural Optimization of Corrugated Plates for Improved Elastic Shear Buckling Capacity".

---

## External Code Dependencies
The code relies on the following external libraries and subroutines:
1.  **`fsparse` function**
    From the `stenglib` library by Stefan Engblom.
    Source: https://github.com/stefanengblom/stenglib

2.  **`ocUpdate` optimization subroutine**
    Available in Ref. [1].

3.  **`mmasub_new.m`**
    A modified version of the original Method of Moving Asymptotes (MMA) code written by Krister Svanberg [2-4].
    Original source: http://www.smoptit.se/

All third-party dependencies are included in the `3rd/` directory for convenience.

---

## Repository Structure
```text
root/
├── Case_trapezoidal_corrugated_web_op.m   # Main script: trapezoidal corrugated web optimization
├── Case_sinusoidal_corrugated_web_op.m    # Main script: sinusoidal corrugated web optimization
├── Case_plate_web_op.m                    # Main script: flat plate web optimization
├── initialData.mat                        # Precomputed initial data/parameters
├── CSWBuck.m                              # Core buckling analysis function
├── README.md                              # This file
├── PostProcessingScript/
│   ├── PostProcessingScript_3D_Geometry.m      # 3D geometry reconstruction & visualization
│   └── PostProcessingScript_2D_Waveform_GIF.m  # 2D waveform plotting & GIF generation
├── initialData/
│   ├── sinusoidal_corrugated_web_1600x220x34mm_1600x300pix.jpg
│   ├── FlatPlate_Vf_0.14_1600x300.bmp
│   └── BCSW1600_430x370x220x36mm.jpg
├── 3rd/
│   ├── fsparse.mexw64                     # Precompiled stenglib mex function
│   ├── ocUpdate.m                         # External optimization subroutine
│   ├── subsolv.m                          # Support routine for MMA solver
│   ├── fsparse.m                          # Source code of fsparse
│   └── mmasub_new.m                       # Modified MMA optimization solver
├── Mesh18x96x24_202606111173945/          # Generated mesh/result folder 1
└── Mesh18x96x24_20260612091023/           # Generated mesh/result folder 2
```

---

## References
[1] F. Ferrari, O. Sigmund, J.K. Guest, Topology optimization with linearized buckling criteria in 250 lines of Matlab, Struct Multidisc Optim 63 (2021) 3045–3066. https://doi.org/10.1007/s00158-021-02854-x.  

[2] Svanberg K. A Class of Globally Convergent Optimization Methods Based on Conservative Convex Separable Approximations. SIAM J Optim 2002;12:555–73. https://doi.org/10.1137/S1052623499362822.  

[3] Svanberg K. The method of moving asymptotes—a new method for structural optimization. International Journal for Numerical Methods in Engineering 1987;24:359–73. https://doi.org/10.1002/nme.1620240207.  

[4] K. Svanberg, MMA and GCMMA – two methods for nonlinear optimization.  
