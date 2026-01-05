# NL-LRD-2025-Data
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18147053.svg)](https://doi.org/10.5281/zenodo.18147053)

This repository contains the data products and analysis codes used in  
**Zhang et al. (2025), "[JWST Insights into Narrow-line Little Red Dots](https://arxiv.org/abs/2506.04350v1)"**.

The purpose of this repository is to ensure transparency and reproducibility of the
main results presented in the paper, including sample selection, spectral fitting,
SED modeling, and figure generation.

---

## 1. Repository Structure

The repository is organized as follows:

```
NL-LRD-2025-Data/
├── BL_simulation/           # Broad-line detectability simulation results (Section 4.2.1 and Figure 13)
├── Code/                    # Core analysis scripts and Jupyter notebooks
├── Figure/                  # figures in the paper
├── catalog/                 # catalogs used in this work (ASTRODEEP_full_combine.fits and all_field_final_cat.fits are large and can be download at https://disk.pku.edu.cn/link/AAC1D514B635784784B122F258A89D7242)
├── jades_galaxy/            # Line measurements for comparison galaxy samples from Curti+2024
├── nakajima_line_fit/       # Line measurements for comparison galaxy samples from Nakajima+2023
├── nircam_throughputs/      # JWST/NIRCam filter transmission curves
├── Narrowline_LRD_image/    # JWST image cutouts of narrow-line LRDs and GalfitM fitting result (Figure 4 and Figure 5)
├── Xray_stacking/           # X-ray stacking analysis, which is based on the code of Yue+2024 (https://github.com/cosmicdawn-mit/xray_lrd).
├── cigale_sed_fitting/      # CIGALE SED fitting data, configurations, and results
├── literature_data/         # Literature data used for comparison
├── narrow_Ha_result_lsf_z2/ # Line measurements for all LRD candidates
├── narrow_Ha_result_all_src/ # Line measurements for all sources
├── paper.mplstyle           # Matplotlib style file used in all figures
└── README.md
```

---

## 2. Code Directory

The `Code/` directory contains the main analysis notebooks used in this paper:

```
Code/
├── LRD_selection_plot.ipynb    # LRD selection plot (Figure 1 and Figure 2)
├── LRD_selection.ipynb         # code used to select LRD candidates (Figure 2)
├── property_compare.ipynb      # Comparison with star-forming galaxies (Figure 3, Figure 6, Figure 7, Figure 10, and Figure 12)
├── plot_RGB.ipynb              # RGB image generation (Figure 4 upper panels)
├── test_fit.ipynb              # Fit the emission lines of medium/high-resolution NIRSpec spectra (Figure 4 lower panels)
├── galfit_sample.ipynb         # Multi-band GALFITM morphology fitting script and fitting (Figure 5)
├── spectra_stacking.ipynb      # Prism spectrum stacking (Section 3.5 and Figure 8)
├── plot_SED.ipynb              # SED visualization (Figure 9)
├── Figure11_sfr.ipynb          # Comparison of UV- and Hα-based SFRs (Figure 11)
├── simulate_broad_Ha_revise.ipynb  # Broad Hα detectability simulations (Section 4.2.1)
├── M_BH-Lbol_revise.ipynb      # Broad Hα detectability simulations plot; Black hole mass – bolometric luminosity relation (Figure 13)
├── M_BH-M_stellar.ipynb        # Black hole mass – stellar mass relation (Figure 14 and Figure 16)
├── test_cal_NII_upperlimit.ipynb # estimate the rms around NII line position (for BPT plot)
├── select_narrowHa_allgal.ipynb # Narrow-line LRD selection (Figure 17)
└── basicfunc.py               # Common utility functions

```

Each notebook is self-documented and can be run independently, provided the required
data files are available.

---

## 3. Data Availability and Reproducibility

- All **derived data products** (line measurements, catalogs, SED fitting results)
  necessary to reproduce the figures and tables in the paper are included.
- Raw JWST imaging and spectroscopy data are **not redistributed**, but are publicly
  available from the JWST archive and the DAWN JWST Archive (DJA) or in our future data release.
- The source IDs used throughout this repository correspond to DJA UIDs.

---

## 4. Software Requirements

The analysis was performed using Python 3.10+.  
Key dependencies include:

- numpy
- scipy
- astropy
- matplotlib
- lmfit
- cigale
- bagpipes

---

## 5. License

This repository is released under the **CC0 1.0 Universal** license.

All data products and analysis codes are made available without restriction
to facilitate reproducibility and reuse.

---

## 6. Citation

If you use any code or data, please cite:

> Zhang, Z., Jiang, L., Liu, W., Ho, L. C., & Inayoshi, K. 2025, *ApJ*, acceptted

Zenodo DOI: [10.5281/zenodo.18147053](https://doi.org/10.5281/zenodo.18147053)

