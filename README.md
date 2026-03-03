# SpecAnalysisALMA

## Overview

**SpecAnalysisALMA** is a Python-based toolkit designed for analyzing ALMA FITS data cubes. It enables users to:
* Plot spectra for selected targets
* Create moment maps for target regions
* Perform Gaussian fitting to derive:
  * Spectroscopic redshift
  * Line flux (Gaussian & integrated)
  * Luminosity (Gaussian & integrated)
  * Line width (FWHM)
  * RMS noise
 
## How to Use GaussFitDict.py
Contains user-defined dictionaries for:
* Rest frequencies of common emission lines
* Associated FITS file paths of .pbcor for each field/line combination
* Associated FITS file paths of .pb for each field/line combination
* Associated file paths (FITS/CASA) of .image for each field/line combination

**Note**: modify or append entries as needed to match your dataset and directory structure. 

## GaussFit.py Version 3
**Note**: Use this if aperture selection based on moment 0 map is preferred. This version also calculates luminosity.  
### Workflow
The analysis flow is as follows:
1. **Initial Moment 0 Map**: creates a moment 0 map for every target source and line.
2. **Selection of Extraction Aperture Size (automated)**: for a given sigma, fits a circular aperture on the moment 0 map for every SPW. Then takes the maximum aperture size as common across all SPW for a given line.
3. **Region File Creation**: automatically generates region files for each target source and line with the selected extraction aperture size.
4. **RMS Noise Measurement**: uses the corresponding **.pb** file for each cube to define circular regions for 95% and 50% of the PB response. Places 100 random apertures (0.5" each) within this annulus on the .image files, and calculates mean RMS across apertures.
5. **FITS Conversion**: converts FITS files into CASA-friendly image format.
6. **Flux Extraction**: uses `imstat` to extract the flux spectrum for each region.
7. **Gaussian Fitting**: fits both single and double Gaussians to each spectral window (SPW). Selects the best fits based on reduced chi-squared statistics.
8. **Best Fit Selection**: among all SPWs, selects the best fit based on a scoring function that weights reduced chi-squared, SNR, FWHM (higher weight) and retains only fits above the user-defined SNR threshold.
9. **Moment Map Generation**: for the chosen SPW, generates moment 0, 1, 2 maps with `immoments`, over a range of ±2σ around the fitted line peak.
10. **Plots**: creates and saves (optional) [mJy] vs. frequency [GHz] plots for each line.
11. **Monte Carlo Analysis**: runs 1000 iterations of Gaussian fitting to estimate Gaussian flux, integrated flux, Gaussian luminosity, integrated luminosity, redshift, and linewidth. It also produces corner plots to verify parameter distributions.
12. **Output**: saves all extracted parameters (redshift, line flux, line luminosity, linewidth, RMS) into a CSV file.

### How to Use
This is the main script. It reads the dictionaries from `GaussFitDict.py` and runs the analysis. 

**Inputs**:
* Field of interest
* Lines of interest
  * If multiple lines from different FITS files: `["line1", "line2", ...]`
  * If multiple lines in the same FITS file: `"line1 & line2"`
  * For a single line: `"line1"`
* Plain text file listing your target coordinates in the format: `# id ra dec zphot`
* SNR threshold (useful if multiple SPW)
* Moment-0 contour level: by default `nsigma = 1` to fit an extraction aperture.
* Extraction aperture selection method: for automatic aperture selection use 'auto' else input numeric value in arcsec.
* Aperture size for moment maps (in arcsec)
* Base directory: required to store folders with additional data.
