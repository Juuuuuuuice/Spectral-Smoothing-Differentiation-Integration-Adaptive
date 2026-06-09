# Spectral-Smoothing-Differentiation-Integration-Adaptive

Code implementation for the paper: **"A fully adaptive method for spectral smoothing, differentiation, and integration"** DOI: [10.1016/j.measurement.2026.122204](https://doi.org/10.1016/j.measurement.2026.122204)

---

## Overview
This repository provides a fully adaptive framework for spectral data processing, covering:

* **Adaptive Smoothing**: Robust noise reduction that preserves signal features.
* **Differentiation**: High-precision numerical differentiation for spectral signals.
* **Integration**: Stable and accurate numerical integration techniques.

## Highlights
* A unified spectral framework based on time and frequency-domain modeling.
* The algorithm unifies smoothing, differentiation, and integration functions.
* A redefined metric ensures high signal fidelity and metrological rigor.
* Maintains linear computational complexity with physically interpretable results.

## Implementation Details

### Core Script: `ZXQuick.m`
This file contains the complete source code for the proposed **Zone eXtraction (ZX)** method. It includes:
* Hierarchical framework for manifold decomposition.
* Adaptive spectral smoothing algorithm.

The script is provided to facilitate the reproduction of the signal reconstruction results presented in the main text.

### Requirements & Dependencies
* **MATLAB**: R2023b or later recommended.
* **Toolbox**: Signal Processing Toolbox.

---

## Citation
If you find this work useful in your research, please consider citing our paper:

```bibtex
@article{YAO2026122204,
title = {A fully adaptive method for spectral smoothing, differentiation, and integration},
journal = {Measurement},
pages = {122204},
year = {2026},
issn = {0263-2241},
doi = {https://doi.org/10.1016/j.measurement.2026.122204},
url = {https://www.sciencedirect.com/science/article/pii/S0263224126019135},
author = {Zhixiang Yao and Ju Yao and Xiaowei Chen and Xiaocheng Huang and Hui Su},
keywords = {Spectral preprocessing, Adaptive denoising, Generalized Gaussian distribution, Fractional-order differentiation, Frequency-domain analysis, Metrological fidelity},
}
