# Measuring Directional Wave Spectra with Random Spatial Arrays

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.placeholder.svg)](https://doi.org/10.5281/zenodo.placeholder)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2022b%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)

Official code and data repository for:

> **Estimating Directional Wave Spectra with Randomly Distributed Spatial Arrays**  
> Clarence O. Collins III, Alexei Skvortsov, Alexander Babanin, and Ian Young  
> *Submitted to Ocean Engineering* (Elsevier).

This repository provides everything needed to:
1. **Instantly reproduce Figures 1–13** from the manuscript using packaged, figure-ready datasets (<80 MB).
2. **Re-run the full Monte Carlo wave field simulations** and sensitivity experiments from scratch.
3. **Run standalone directional wave spectrum estimators** (BDM-NR, BDM-NNLS, EWDM, EMEM, and IMLM) on arbitrary random or fixed spatial arrays.

---

## Quickstart: Reproduce All Paper Figures (1 Command)

Open MATLAB in this directory and run:

```matlab
reproduce_all_figures
```

In approximately **30–40 seconds**, all 13 publication figures (in both high-resolution PNG and vector EPS formats) will be generated and saved into the `./figures/` folder:

| Figure | Description | Key Result |
| :--- | :--- | :--- |
| **Figure 1** | Directional distributions vs dimensionless scale $l^*$ (Tests 1–11) | Identifies the optimal array scale $l^* \approx 10^{-1.5}$ ($l_{rms} \approx 1.8$ m). |
| **Figure 2** | Precision $\sigma(\theta_p)$ and MIAE vs $l^*$ ($N=4, 10$) | Demonstrates the U-shaped error curve and wider operational zone for $N=10$. |
| **Figure 3** | Test 4 ($l^* = 10^{-2.5}$) Configurations & Co-arrays | Array collapse and resolution limitations of undersized arrays. |
| **Figure 4** | Test 8 ($l^* = 10^{-0.5}$) Configurations & Co-arrays | Spatial aliasing and elongated baseline geometry. |
| **Figure 5** | Test 7 ($l^* = 10^{-1.0}$) Configurations & Co-arrays | Robust baseline geometry near the optimal performance valley. |
| **Figure 6** | Error metrics vs sensor count ($N = 3 \dots 20$) | Major jump from $N=3 \to 4$; sweet spot at $N \approx 6\text{--}8$. |
| **Figure 7** | Directional distributions vs sensor count ($N = 3 \dots 20$) | Ensemble narrowing and reduction in realization variance. |
| **Figure 8** | 3-Element Arrays ($N=3$) Configurations & Co-arrays | Geometric instability and baseline clustering for minimal arrays. |
| **Figure 9** | BDM-NR Sensitivity to Uncorrelated Elevation Noise | Robustness from $0.1\%$ to $10\%$; breakdown at $>50\%$ noise. |
| **Figure 10** | BDM-NR Sensitivity to Temporal Desynchronization (Lag) | Phase error behavior across $0.1\% \dots 10.0\% \, T_p$. |
| **Figure 11** | BDM-NR Sensitivity to Sensor Spatial Offset | Position perturbation response across $0.1\% \dots 10.0\% \, \lambda_p$. |
| **Figure 12** | Estimator Comparison on Baseline Array ($N=7, l^* = 10^{-1.5}$) | Evaluation of BDM-NR, EMEM, EWDM, IMLM, and BDM-NNLS. |
| **Figure 13** | Estimator Metrics across Sensor Count ($N = 3, 7, 15$) | Scaling behavior: EWDM refinement, EMEM stability, IMLM irregular array sensitivity. |

---

## Repository Structure

```text
random-array-wave-spectra/
├── README.md                      # Documentation, quickstart, and citation
├── LICENSE                        # MIT Open Source License
├── reproduce_all_figures.m        # Master 1-click replication script
├── data/                          # Compact figure datasets (~76 MB total)
│   ├── targetL2/
│   │   ├── N4/                    # 11 scale sweep datasets (Tests 1-11, N=4)
│   │   └── N10/                   # 11 scale sweep datasets (Tests 1-11, N=10)
│   ├── N/                         # 10 sensor count datasets (N=3, 4, ..., 20)
│   ├── sensitivity/
│   │   ├── noise/                 # 6 noise level tests (0.1% to 100%)
│   │   ├── lag/                   # 5 time lag tests (0.1% to 10.0% Tp)
│   │   └── offset/                # 5 spatial offset tests (0.1% to 10.0% Lp)
│   └── methods/
│       ├── N3/                    # 5 estimators benchmarked on N=3 arrays
│       ├── N7/                    # 5 estimators benchmarked on N=7 arrays (baseline)
│       └── N15/                   # 5 estimators benchmarked on N=15 arrays
├── src/
│   ├── estimators/                # Core directional wave spectrum algorithms
│   │   ├── bdm.m                  # Bayesian Directional Method (Newton-Raphson)
│   │   ├── bdm_nnls_method.m      # Bayesian Directional Method (NNLS)
│   │   ├── emem.m                 # Extended Maximum Entropy Method
│   │   ├── ewdm_method.m          # Extended Wavelet/Direct Method
│   │   ├── imlm.m                 # Iterated Maximum Likelihood Method
│   │   ├── mlm.m                  # Maximum Likelihood Method
│   │   └── getcrossspectra.m      # Cross-spectral matrix calculation
│   ├── simulation/                # Sea surface & sensor synthesis
│   │   ├── simulateSensorTimeseries.m
│   │   ├── randomParticleLocationsFromL2.m
│   │   └── randomParticleLocations.m
│   ├── geometry/                  # Array and co-array analysis
│   │   ├── coarray.m
│   │   ├── coarraySubplot.m
│   │   ├── polar_from_Cartesian.m
│   │   └── circular_tukey.m
│   └── utils/                     # Plotting and formatting utilities
│       ├── shadedErrorBar.m
│       └── packfig.m
├── scripts/                       # Scripts to re-run simulations from scratch
│   ├── sim_scale_targetL2.m       # Re-runs scale sweep (Tests 1-11)
│   ├── sim_sensor_count.m         # Re-runs sensor count sweep (N=3..20)
│   ├── sim_noise_sensitivity.m    # Re-runs noise perturbation sweep
│   ├── sim_lag_sensitivity.m      # Re-runs clock desynchronization sweep
│   ├── sim_offset_sensitivity.m   # Re-runs spatial GPS error sweep
│   └── sim_estimator_benchmark.m  # Re-runs multi-estimator benchmark
├── figures/                       # Destination folder for generated figures
└── docs/
    └── environment.md             # Environment and dependency details
```

---

## Running Full Simulations From Scratch

If you wish to regenerate the underlying Monte Carlo ensemble datasets from scratch:

1. **Scale Sweep (Tests 1–11):**
   ```matlab
   run('scripts/sim_scale_targetL2.m')
   ```
2. **Sensor Count Sweep ($N=3\dots20$):**
   ```matlab
   run('scripts/sim_sensor_count.m')
   ```
3. **Sensitivity Experiments (Noise, Lag, Spatial Offset):**
   ```matlab
   run('scripts/sim_noise_sensitivity.m')
   run('scripts/sim_lag_sensitivity.m')
   run('scripts/sim_offset_sensitivity.m')
   ```
4. **Directional Estimator Benchmark:**
   ```matlab
   run('scripts/sim_estimator_benchmark.m')
   ```

*Note: The Monte Carlo simulations generate 100 random geometric configurations evaluated over 30 wavefield realizations per test case. Depending on your CPU and parallel worker pool, full simulations take between 15 minutes and several hours.*

---

## Software Requirements & Toolboxes

- **MATLAB:** R2022b or newer recommended.
- **Required MATLAB Toolboxes:**
  - Signal Processing Toolbox
  - Statistics and Machine Learning Toolbox
- **Third-Party Libraries:**
  - [WAFO (Wave Analysis for Fatigue and Oceanography)](https://github.com/wafo-project/wafo)
  - Standalone estimator implementations and utility functions (`shadedErrorBar`, `packfig`) are bundled directly inside `./src/` for self-contained execution.

---

## Citation

If you use this code or data in your research, please cite:

```bibtex
@article{collins2026estimating,
  author    = {Collins, III, Clarence O. and Skvortsov, Alexei and Babanin, Alexander and Young, Ian},
  title     = {Estimating Directional Wave Spectra with Randomly Distributed Spatial Arrays},
  journal   = {Ocean Engineering},
  year      = {2026},
  note      = {Submitted},
  doi       = {10.5281/zenodo.placeholder}
}
```

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
