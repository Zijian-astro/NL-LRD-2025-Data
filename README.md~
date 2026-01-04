# NL-LRD-2025-Data
<<<<<<< HEAD
Data and code of Zhang+2025 Narrow-line LRD paper:
https://arxiv.org/abs/2506.04350v1

All selected LRDs in this work are in catalog/all_nirspechighz_LRD_z2.ecsv.
=======
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
├── BL_simulation/           # Broad-line detectability simulations (Section 3.2, 4.2)
├── Code/                    # Core analysis scripts and Jupyter notebooks
├── Figure/                  # Scripts and intermediate files used to generate paper figures
├── catalog/                 # Source catalogs used in this work
├── jades_galaxy/            # Comparison samples from JADES
├── nakajima_line_fit/       # Line measurements for comparison galaxy samples
├── nircam_throughputs/      # JWST/NIRCam filter transmission curves
├── Narrowline_LRD_image/    # Image cutouts of narrow-line LRDs (Figure 4)
├── Xray_stacking/           # X-ray stacking analysis
├── cigale_sed_fitting/      # CIGALE SED fitting configurations and results
├── literature_data/         # Literature data used for comparison
├── narrow_Ha_result_lsf_z2/ # Hα fitting results with LSF correction
├── paper.mplstyle           # Matplotlib style file used in all figures
└── README.md
```

---

## 2. Code Directory

The `Code/` directory contains the main analysis notebooks used in this paper:

```
Code/
├── Figure11_sfr.ipynb          # Comparison of UV- and Hα-based SFRs (Figure 11)
├── M_BH-Lbol_revise.ipynb      # Black hole mass – bolometric luminosity relation
├── M_BH-M_stellar.ipynb       # Black hole mass – stellar mass relation
├── spectra_stacking.ipynb     # Prism spectrum stacking (Section 3.5)
├── simulate_broad_Ha_revise.ipynb  # Broad Hα detectability simulations
├── galfit_sample.ipynb        # Multi-band GALFITM morphology fitting
├── plot_SED.ipynb             # SED visualization
├── plot_RGB.ipynb             # RGB image generation
├── property_compare.ipynb     # Comparison with star-forming galaxies
├── test_cal_NII_upperlimit.ipynb
├── basicfunc.py               # Common utility functions
└── mysel_kor.fits / mysel_kov.fits  # LRD selection results from literature criteria
```

Each notebook is self-documented and can be run independently, provided the required
data files are available.

---

## 3. Data Availability and Reproducibility

- All **derived data products** (line measurements, catalogs, SED fitting results)
  necessary to reproduce the figures and tables in the paper are included.
- Raw JWST imaging and spectroscopy data are **not redistributed**, but are publicly
  available from the JWST archive and the DAWN JWST Archive (DJA).
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

A Zenodo DOI will be provided upon acceptance.
>>>>>>> 9fa247f (update README)
