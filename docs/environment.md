# Computational Environment and Dependencies

This project was developed and validated with the following software configuration:

## Operating System
- macOS Sonoma / Sequoia (Apple Silicon and Intel compatible)
- Linux (Ubuntu 20.04/22.04 LTS tested)
- Windows 10/11 compatible

## MATLAB Version
- MATLAB R2022b (Version 9.13) or higher.

## MATLAB Toolboxes
- Signal Processing Toolbox (spectral analysis, filtering)
- Statistics and Machine Learning Toolbox (distribution fitting, random distributions)

## External Toolboxes
1. **WAFO (Wave Analysis for Fatigue and Oceanography)**:
   - Version 2.6 / 2022b compatibility.
   - Available at: https://github.com/wafo-project/wafo
2. **jLab (Dr. Jonathan Lilly)**:
   - Used for `packfig` subplot arrangement. Bundled in `src/utils/packfig.m`.
   - Available at: http://www.jmlilly.net/software.html
3. **shadedErrorBar (Rob Campbell)**:
   - Bundled in `src/utils/shadedErrorBar.m`.
   - Available on MATLAB Central File Exchange.

## Estimator Algorithms Included
- **BDM-NR**: Hashimoto & Kobune (1988), Hashimoto et al. (2003) via Newton-Raphson.
- **BDM-NNLS**: Non-negative least squares formulation of BDM.
- **EMEM**: Extended Maximum Entropy Method (Kobune & Hashimoto 1986).
- **EWDM**: Extended Wavelet/Direct Method (Peláez et al. 2024).
- **IMLM**: Iterative Maximum Likelihood Method (Pawka 1983; DIWASP formulation).
