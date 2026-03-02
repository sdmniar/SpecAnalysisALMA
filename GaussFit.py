from collections import defaultdict
from casatasks import importfits, exportfits, imstat, immoments 
from casatools import image
from astropy.cosmology import Planck18 as cosmo 
from astropy.wcs import WCS

from scipy.optimize import curve_fit
from scipy.special import erf
from scipy.ndimage import label
from scipy import stats

from astropy.io import fits 
from astropy.coordinates import SkyCoord
import astropy.units as u
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

import numpy as np
import re, time, os, csv, random, math
import corner

from GaussFitDict import freqs, fits_list, pb_list, noPBcor_list
import matplotlib as mpl
from pathlib import Path

timestamp = int(time.time())

# ----- INPUT PARAMETERS -----

field = " " 

line_list =  " " # if multiple entries use ["line1", "line2", ...] else "line1" or "line1 & line2" for same entry

coords_file = " " # path for .txt file with # id ra dec zphot

snr_threshold = 0 # if multiple SPW

nsigma = 1 # for aperture fitting around moment 0 map

rms_diagnose = True

aperture_selection = 0.8 # for automatic aperture selection use 'auto' else input numeric value in arcsec

mom_reg_size = "3.0000" # region size for moment maps

base_dir = Path(" ") # base directory to store data

# ----------------------------

# Create folders
dirs = {
    "inimommap": base_dir / "Data" / "ini_mommap" / f"{field}",
    "momregion": base_dir / "Data" / "region_mommap" / f"{field}",
    "exregion": base_dir / "Data" / "region_extract" / f"{field}",
    "finalmommap": base_dir / "Data" / "final_mommap" / f"{field}",
    "results": base_dir / "Results" / f"{field}",
    "fig": base_dir / "Figures" / f"{field}",
}

for path in dirs.values():
    path.mkdir(parents=True, exist_ok=True)

try:
    if isinstance(line_list, str):
        if "&" in line_list:
            line_name = [ln.strip() for ln in line_list.split("&")]
            combined_key = line_list.strip()
            rest_frequencies = {ln: freqs[ln] for ln in line_name}
            rest_frequencies_list = list(rest_frequencies.values())
            line_ratio = max(rest_frequencies_list) / min(rest_frequencies_list)          
            fits_files = {line_list: fits_list[field][line_list]}
            
        else:
            line_name = [line_list]
            combined_key = None
            rest_frequencies = {line_list: freqs[line_list]}
            line_ratio = None
            fits_files = {ln: fits_list[field][ln] for ln in line_name}

    elif isinstance(line_list, list):
        line_name = line_list
        combined_key = None
        rest_frequencies = {ln: freqs[ln] for ln in line_name}
        line_ratio = None
        fits_files = {ln: fits_list[field][ln] for ln in line_name}
    else:
        raise ValueError("Please use the recommended format.")

except RuntimeError:
    print("Please use the recommended format.")

# Region selection for RMS estimation
def pbres(pbfile, fractions=(0.95, 0.50)):
    pb_data = fits.getdata(pbfile)
    header = fits.getheader(pbfile)

    pb_data = np.nan_to_num(pb_data)

    if pb_data.ndim > 2:
        pb_2d = pb_data[0, 0, :, :]
    else:
        pb_2d = pb_data

    ny, nx = pb_2d.shape
    cy, cx = ny // 2, nx // 2

    peak = pb_2d[cy, cx]
    pb_norm = pb_2d / peak 

    y, x = np.indices(pb_2d.shape)

    r_pixels = np.sqrt((x - cx)**2 + (y - cy)**2)

    pix_arc_scale = abs(header['CDELT1']) * 3600

    results = {}

    for frac in fractions:
        mask = pb_norm <= frac
        if np.any(mask):
            r_pix = r_pixels[mask].min()
            results[frac] = r_pix * pix_arc_scale
        else:
            results[frac] = None

    return results
    

# Extract boundaries of each FITS file
def image_bounds(imagename):
    ia = image()
    try:
        ia.open(imagename)
        csys = ia.coordsys()
        shape = ia.shape()
        w1 = csys.toworld([0, 0, 0, 0])['numeric']
        w2 = csys.toworld([shape[0]-1, shape[1]-1, 0, 0])['numeric']
        ra1 = np.degrees(w1[0])
        dec1 = np.degrees(w1[1])
        ra2 = np.degrees(w2[0])
        dec2 = np.degrees(w2[1])
        ra_vals = [ra1, ra2]
        dec_vals = [dec1, dec2]
        ra_min, ra_max = min(ra_vals), max(ra_vals)
        dec_min, dec_max = min(dec_vals), max(dec_vals)
        ia.close()
        return ra_min, ra_max, dec_min, dec_max
    except Exception as e:
        try:
            ia.close()
        except: pass
        print(f"Error in obtaining image bounds for {imagename}: {e}")
        return None, None, None, None

# Check if coordinates are within image boundaries
def coord_bounds(ra, dec, ra_min, ra_max, dec_min, dec_max):
    if None in (ra_min, ra_max, dec_min, dec_max):
        return False
    # Normalise
    ra = ra % 360
    ra_min = ra_min % 360
    ra_max = ra_max % 360
    if ra_min <= ra_max:
        ra_ok = (ra_min <= ra <= ra_max)
    else:
        ra_ok = (ra >= ra_min) or (ra <= ra_max)
    dec_ok = (dec_min <= dec <= dec_max)
    return ra_ok and dec_ok

# Create circular string of regions around target to extract RMS noise
def rms_noise(imagename, pbfile,
              bkg_radius,
              num_regions=50,
              seed=42, diagnostic=None):

    import numpy as np
    np.random.seed(seed)

    radii = pbres(pbfile)
    inner_radius = radii[0.95]
    outer_radius = radii[0.50]

    ra_min, ra_max, dec_min, dec_max = image_bounds(imagename)

    ia = image()
    ia.open(imagename)

    shape = ia.shape()
    nchan = shape[2]

    csys = ia.coordsys()
    nx, ny = shape[0], shape[1]

    center_x = nx / 2
    center_y = ny / 2
    world = csys.toworld([[center_x, center_y]])

    center_ra_deg  = np.degrees(world['numeric'][0])
    center_dec_deg = np.degrees(world['numeric'][1])

    flux_spectra = []

    for i in range(num_regions):
        for j in range(10):
            angle = np.random.uniform(0, 2*np.pi)
            r = np.random.uniform(inner_radius, outer_radius)

            ra_bkg  = center_ra_deg  + (r * np.cos(angle)) / 3600.0
            dec_bkg = center_dec_deg + (r * np.sin(angle)) / 3600.0

            if not coord_bounds(ra_bkg, dec_bkg,
                                ra_min, ra_max, dec_min, dec_max):
                continue

            region = f"circle[[{ra_bkg}deg, {dec_bkg}deg], {bkg_radius}arcsec]"

            try:
                # imstat returns per-channel values when axes=[0,1]
                stats_result = imstat(imagename=imagename,
                                region=region,
                                axes=[0, 1])

                flux = np.array(stats_result["flux"])

                if flux.size != nchan:
                    continue

                flux_spectra.append(flux)
                break

            except Exception:
                continue

    ia.close()

    if len(flux_spectra) == 0:
        return np.full(nchan, np.nan)

    flux_spectra = np.array(flux_spectra)  # (N_ap, N_chan)

    # RMS across apertures per channel
    rms_per_channel = np.nanstd(flux_spectra, axis=0)

    if diagnostic:
        flat_flux_spectra = flux_spectra.flatten()
        flat_flux_spectra = flat_flux_spectra[~np.isnan(flat_flux_spectra)]
        plt.figure(figsize=(12, 5))
        plt.subplot(1,2,1)
        
        plt.hist(flat_flux_spectra, bins=50, density=True,
                              alpha=0.7, color='steelblue', edgecolor='black')
        
        mu, sigma = np.mean(flat_flux_spectra), np.std(flat_flux_spectra)
        x = np.linspace(flat_flux_spectra.min(), flat_flux_spectra.max(), 200)
        plt.plot(x, stats.norm.pdf(x, mu, sigma), color='maroon', lw=2.5, 
                 label=f'μ={mu:2e}\nσ={sigma:.2e}')
        plt.xlabel('Flux Density (Jy)')
        plt.ylabel('Probability Density')
        plt.title(f'{num_regions} Apertures of Radius {bkg_radius}"')
        plt.grid(alpha=0.3)
        plt.legend()

        plt.subplot(1,2,2)
        ax = plt.gca()
        stats.probplot(flat_flux_spectra, dist="norm", plot=ax)
        plt.title('Normal Q-Q Plot')
        plt.grid(alpha=0.3)

        plt.tight_layout()
        plt.savefig(str(dirs["fig"] / diagnostic), dpi=150, bbox_inches="tight")
        plt.close()

    return rms_per_channel

# Sigma estimation from moment 0 map
def sigma_estimate(data):
    flat = data.flatten()
    mean = np.nanmean(flat)
    sigma = np.nanstd(flat)
    return mean, sigma

# Aperture size from moment 0 map
def sigma_aperture(fitsfile, ra, dec, nsigma=nsigma):
    hdu = fits.open(fitsfile)[0]
    data = np.squeeze(hdu.data).astype(float)

    if data.ndim > 2:
        data = data[0, :, :]

    wcs = WCS(hdu.header).celestial

    coord = SkyCoord(ra, dec, unit=(u.deg, u.deg))
    x0, y0 = coord.to_pixel(wcs)
    x0, y0 = float(x0), float(y0)

    mean, sigma = sigma_estimate(data)
    thresh = mean + nsigma * sigma

    mask = data >= thresh
    labeled, _ = label(mask)

    ix = int(round(x0))
    iy = int(round(y0))

    mom_id = labeled[iy, ix]
    if mom_id == 0:
        raise ValueError("Center pixel not inside 1σ island")
    
    mom_mask = (labeled == mom_id)

    ys, xs = np.where(mom_mask)
    distances = np.sqrt((xs - x0)**2 + (ys - y0)**2)
    r_pix = distances.max()

    pixscale = np.mean(np.abs(wcs.pixel_scale_matrix.diagonal())) * 3600
    r_arcsec = r_pix * pixscale

    return data, wcs, mom_mask, (x0, y0), r_pix, r_arcsec, thresh

# Plot aperture and moment map 
def plot_aperture(data, wcs, mom_mask, center_pix, r_pix, r_arcsec):
    fig = plt.figure(figsize=(7, 6))
    ax = plt.subplot(projection=wcs)

    im = ax.imshow(data, origin='lower', cmap='inferno')
    plt.colorbar(im, ax=ax, label="Moment 0")

    ax.contour(mom_mask.astype(int),
               levels=[0.5],
               colors='cyan',
               linewidths=0.8)
    cx, cy = center_pix
    circ = Circle((cx,cy), r_pix, edgecolor='white', facecolor='none', lw=1.2)
    ax.add_patch(circ)
    plt.tight_layout()
    plt.show()

# Extract flux and frequency
def measurement(imagename, region):
    im = image()
    im.open(imagename)
    coords = im.coordsys()
    num_channels = im.shape()[2] # Number of channels
    frequencies = [coords.toworld([0, 0, freq, 0])['numeric'][2] for freq in range(num_channels)] # extract frequencies from channels - [RA, DEC, FREQ, STOKES]
    im.close()

    fluxes = []

    for channel in range(num_channels):
        stats_result = imstat(imagename=imagename, chans=str(channel), region=region)

        flux = stats_result["flux"][0] if "flux" in stats_result and stats_result["flux"].size > 0 else 0 

        fluxes.append(flux)
    
    fluxes = np.array(fluxes)
    frequencies = np.array(frequencies)

    return frequencies, fluxes

# Gaussian fit function 
def gaussian(x, amp, mean, std):
    return amp * np.exp(-(x-mean)**2 / (2*std**2))

# One gaussian fit
def one_gaussian(x, amp, mean, std):
    return gaussian(x, amp, mean, std)

# Two gaussian fit
def two_gaussian(x, amp1, mean1, std1, amp2, line_ratio, std2):
    mean2 = mean1 * line_ratio
    return gaussian(x, amp1, mean1, std1) + gaussian(x, amp2, mean2, std2)

model = lambda x, amp1, mean1, std1, amp2, std2: two_gaussian(x, amp1, mean1, std1, amp2, line_ratio, std2)

def two_gaussian_free(x, amp1, mean1, std1, amp2, mean2, std2):
    return gaussian(x, amp1, mean1, std1) + gaussian(x, amp2, mean2, std2)

# Reduced chi-squared
def reduced_chi(y_data, y_fit, dof):
    residual = y_data - y_fit
    return np.sum(residual**2) / dof

def freq_to_vel(nu):
    return 3e5 * (restfreq - nu) / restfreq

def vel_to_freq(v):
    return restfreq * (1 - v/3e5) 

# Write results to CSV file
def decimals(val, n=4):
    if isinstance(val, (int, float)):
        return f"{val:.{n}e}"
    return val

# Load coordinates
coords_data = np.loadtxt(coords_file)
if coords_data.ndim == 1:
    coords_data = coords_data.reshape(1, -1)

id = coords_data[:, 0].astype(int)

ra = coords_data[:, 1].astype(float)
dec = coords_data[:, 2].astype(float)
z_phot = coords_data[:, 3].astype(float)

region_files = []
mom_region_files = []

aperture_dict = {}
region_id_map = {}

r_arcsec_list = []

for i in range(len(id)):
    print(i)
    ra0 = ra[i]
    dec0 = dec[i]

    ra_coord = str(ra0)
    dec_coord = str(dec0)

    if "&" in line_list:
        lines_to_process = [line_list]
    else:
        lines_to_process = line_name

    for line in lines_to_process:
        line_tag = line.replace(" ", "").replace("&", "_")

        

        for spw, fits_f in enumerate(fits_files[line]):
            print(f" ID {id[i]} | {line} | spw{spw}")

            imagename = fits_f.replace(".fits", "") + ".image"
            importfits(fitsimage=fits_f, imagename=imagename, overwrite=True) 

            # Moment map region file
            mom0_img = str(dirs["inimommap"] / f"{int(id[i])}_{line_tag}_spw{spw}_mom0_{timestamp}.image")

            tmp_region = str(dirs["momregion"] / f"{int(id[i])}_{line_tag}_spw{spw}.crtf")

            with open(tmp_region, "w") as file:
                file.write("#CRTFv0 CASA Region Text Format version 0 \n")
                file.write(f"circle [[{ra_coord}deg, {dec_coord}deg], {mom_reg_size}arcsec] "
                        "coord=ICRS, corr=[I], linewidth=2, linestyle=-, symsize=1, "
                        "symthick=1, color=2ee6d6, font=Helvetica, fontsize=10, "
                        "fontstyle=bold, usetex=false")
                
            immoments(imagename=imagename, moments=[0], axis="spectral", region=tmp_region, outfile=mom0_img)
            exportfits(imagename=mom0_img, fitsimage=mom0_img.replace(".image",".fits"), overwrite=True)

            mom0_fits = mom0_img.replace(".image", ".fits")

            if aperture_selection == 'auto':
                data, wcs, mom_mask, center, r_pix, r_arcsec, _ = sigma_aperture(mom0_fits, ra0, dec0)

                r_arcsec_list.append(r_arcsec)

                print(r_arcsec_list)

                plot_aperture(data, wcs, mom_mask, center, r_pix, r_arcsec)
            else:
                r_arcsec = aperture_selection
                r_arcsec_list.append(r_arcsec)
        
best_r = max(r_arcsec_list) # common aperture size across all SPW

print("BEST_R:", best_r)

for i in range(len(id)):
    print(i)
    ra0 = ra[i]
    dec0 = dec[i]

    ra_coord = str(ra0)
    dec_coord = str(dec0)

    if "&" in line_list:
        lines_to_process = [line_list]
    else:
        lines_to_process = line_name

    for line in lines_to_process:

        line_tag = line.replace(" ", "").replace("&", "_")

        if aperture_selection == 'auto':
            radius = best_r
        elif isinstance(aperture_selection, (int, float)):
            radius = float(aperture_selection)
        else:
            raise ValueError(f"Invalud aperture_selection: {aperture_selection}. Use 'auto' or a numeric value in arcsec")

        best_region_file = str(dirs["exregion"] / f"{int(id[i])}_{line_tag}_best_{radius:.3f}arcsec.crtf")

        with open(best_region_file, "w") as file:
            file.write("#CRTFv0 CASA Region Text Format version 0 \n")
            file.write(f"circle [[{ra_coord}deg, {dec_coord}deg], {radius}arcsec] "
                        "coord=ICRS, corr=[I], linewidth=2, linestyle=-, symsize=1, "
                        "symthick=1, color=2ee6d6, font=Helvetica, fontsize=10, "
                        "fontstyle=bold, usetex=false")
            
        aperture_dict[(int(id[i]), line)] = radius
        region_id_map[(int(id[i]), line)] = best_region_file

# RMS estimation
rms_chan_dict = {}
rms_results = defaultdict(lambda: defaultdict(list))

for ln, fnames in fits_files.items():
    split_lines = [s.strip() for s in ln.split("&")] if "&" in ln else [ln]
    for f, pbfile, nopbcor_file in zip(
        fnames,
        pb_list[field][ln],
        noPBcor_list[field][ln]
    ):
        if nopbcor_file.endswith(".fits"):
            nopbcor_imagename = nopbcor_file.replace(".fits", ".image")
            importfits(
                fitsimage=nopbcor_file,
                imagename=nopbcor_imagename,
                overwrite=True
            )
        else:
            nopbcor_imagename = nopbcor_file

        if rms_diagnose == True: 
            diagnose_file = f"RMS_Diagnostic_{ln.replace(' ', '_').replace('&', 'and')}_r{radius}_{timestamp}.png"
            rms_chan = rms_noise(nopbcor_imagename, pbfile, radius, diagnostic=diagnose_file)
        else:
            rms_chan = rms_noise(nopbcor_imagename, pbfile, radius)
        rms_med = np.nanmedian(rms_chan)

        if not np.isfinite(rms_med) or rms_med <= 0:
            continue
        for ln in split_lines:
            rms_chan_dict[ln] = rms_chan

        if f.endswith(".fits"):
            flux_imagename = f.replace(".fits", ".image")
            importfits(
                fitsimage=f,
                imagename=flux_imagename,
                overwrite=True
            )
        else:
            flux_imagename = f

        ra_min, ra_max, dec_min, dec_max = image_bounds(flux_imagename)


        for tid, ra_deg, dec_deg in zip(id, ra, dec):

            if not coord_bounds(
                ra_deg, dec_deg,
                ra_min, ra_max,
                dec_min, dec_max
            ):
                continue

            for sub_line in split_lines:
                rms_results[sub_line][tid].append(rms_med)

mean_rms_dict = {
    ln: {tid: np.mean(vals) * 1e3 for tid, vals in d.items()}
    for ln, d in rms_results.items()
}

# Dictionaries
region_z_phot_map = dict(zip(id, z_phot))

region_line_z = {}
region_line_z_err = {}
gauss_flux = {}
gauss_flux_err = {}
int_flux = {}
int_flux_err = {}
gauss_lum = {}
gauss_lum_err = {}
int_lum = {}
int_lum_err = {}
line_width = {}
line_width_err = {}

for i in range(len(id)):
    region_line_z[f'{id[i]}'] = {}
    region_line_z_err[f'{id[i]}'] = {}

    gauss_flux[f'{id[i]}'] = {}
    gauss_flux_err[f'{id[i]}'] = {}

    int_flux[f'{id[i]}'] = {}
    int_flux_err[f'{id[i]}'] = {}

    line_width[f'{id[i]}'] = {}
    line_width_err[f'{id[i]}'] = {}

    gauss_lum[f'{id[i]}'] = {}
    gauss_lum_err[f'{id[i]}'] = {}

    int_lum[f'{id[i]}'] = {}
    int_lum_err[f'{id[i]}'] = {}

for (region_label, line), region_file in region_id_map.items():
    fits_file = fits_files[line] 
    
    num_spws = len(fits_file)

    plot_col = min(num_spws, 2)
    plot_row = math.ceil(num_spws / plot_col)

    plt.figure(figsize=(6 * plot_col, 4 * plot_row))

    if "&" in line:
        lines_split = [l.strip() for l in line.split("&")]
    else:
        lines_split = [line]

    best_score = np.inf
    best_fit = []
    best_err = []
    total_frequencies = []

    for i, fits_f in enumerate(fits_file, start=1):
        # Convert fits to image files
        imagename = fits_f.replace(".fits", "") + ".image"
        importfits(fitsimage=fits_f, imagename=imagename, overwrite=True)

        frequencies, fluxes = measurement(imagename, region_file)

        freq_GHz = frequencies / 1e9 # frequency in GHz

        total_frequencies.extend(frequencies.tolist())

        try:
            # Guess for Gaussian fit
            guess_amp = np.max(fluxes)
            guess_mean = frequencies[np.argmax(fluxes)]
            guess_std = 500 * (guess_mean / 3e5)

            # One Gaussian Fit
            try:
                popt1, pcov1 = curve_fit(one_gaussian, frequencies, fluxes, p0=[guess_amp, guess_mean, guess_std], maxfev=10000)
                fit1 = one_gaussian(frequencies, *popt1)
                chi1 = reduced_chi(fluxes, fit1, len(frequencies) - 3)
            except RuntimeError:
                popt1 = pcov1 = fit1 = None
                chi1 = np.inf
            
            two_fit = np.argsort(fluxes)[-2:]
            peak1_idx, peak2_idx = sorted(two_fit, key=lambda i: frequencies[i])
            amp1, mean1 = fluxes[peak1_idx], frequencies[peak1_idx]
            amp2 = fluxes[peak2_idx]

            if line_ratio is not None:
                # Two Gaussian Fit
                guess_std2 = 500 * ((mean1 * line_ratio) / 3e5)

                try:
                    popt2, pcov2 = curve_fit(model, frequencies, fluxes, p0=[amp1, mean1, guess_std, amp2, guess_std2], maxfev=10000)
                    amp1, mean1, std1, amp2, std2 = popt2
                    mean2 = mean1 * line_ratio
                    fit2 = two_gaussian(frequencies, amp1, mean1, std1, amp2, line_ratio, std2)

                    chi2 = reduced_chi(fluxes, fit2, len(frequencies) - 5)
                
                except Exception as e:
                    print("Two-Gaussian fit failed")
                    popt2 = pcov2 = fit2 = None
                    chi2 = np.inf
            else:
                guess_mean2 = frequencies[np.argmax(fluxes)] + guess_std

                try:
                    popt2, pcov2 = curve_fit(two_gaussian_free, frequencies, fluxes, p0=[amp1, mean1, guess_std, amp2, guess_mean2, guess_std], 
                                                bounds = ([0, min(frequencies), 0, 0, min(frequencies), 0], [np.inf, max(frequencies), np.inf, np.inf, max(frequencies), np.inf]),
                                                maxfev=10000)
                    amp1, mean1, std1, amp2, mean2, std2 = popt2
                    fit2 = two_gaussian_free(frequencies, amp1, mean1, std1, amp2, mean2, std2)

                    chi2 = reduced_chi(fluxes, fit2, len(frequencies) - 6)

                except Exception as e:
                    print(f"Two-Gaussian fit failed: {e}")
                    popt2 = pcov2 = fit2 = None
                    chi2 = np.inf

            # Choose best fit
            best_num_gaussian = chi2 < chi1 # least better (one gaussian preferred)

            if best_num_gaussian and popt2 is not None:
                if line_ratio is not None:
                    amp1, mean1, std1, amp2, std2 = popt2
                    mean2 = mean1 * line_ratio
                    amp1_err = np.sqrt(pcov2[0, 0])
                    mean1_err = np.sqrt(pcov2[1, 1])
                    std1_err = np.sqrt(pcov2[2, 2])
                    amp2_err = np.sqrt(pcov2[3, 3])
                    mean2_err = mean1_err * line_ratio
                    std2_err = np.sqrt(pcov2[4, 4])
                else:
                    amp1, mean1, std1, amp2, mean2, std2 = popt2
                    amp1_err = np.sqrt(pcov2[0, 0])
                    mean1_err = np.sqrt(pcov2[1, 1])
                    std1_err = np.sqrt(pcov2[2, 2])
                    amp2_err = np.sqrt(pcov2[3, 3])
                    mean2_err = np.sqrt(pcov2[4, 4])
                    std2_err = np.sqrt(pcov2[5, 5])

                chi = chi2
                fit = fit2
                line1 = 2 * np.sqrt(2 * np.log(2)) * abs(std1)
                line2 = 2 * np.sqrt(2 * np.log(2)) * abs(std2)
                label_fit = 'Double Gaussian'
            
            else:
                amp1, mean1, std1 = popt1
                amp1_err = np.sqrt(pcov1[0, 0])
                mean1_err = np.sqrt(pcov1[1, 1])
                std1_err = np.sqrt(pcov1[2, 2])

                amp2 = mean2 = std2 = amp2_err = mean2_err = std2_err = line2 = None
                chi = chi1
                fit = fit1
                line1 = 2 * np.sqrt(2 * np.log(2)) * abs(std1)
                label_fit = 'Single Gaussian'

            snr1 = amp1 / (mean_rms_dict[lines_split[0]][region_label]/1e3)
            snr2 = amp2 / (mean_rms_dict[lines_split[0]][region_label]/1e3) if amp2 else 0
            snr = max(snr1, snr2)

            fwhm_total = line1 + (line2 if line2 else 0)
            score = chi / (snr * fwhm_total**1.5)

            # Choose best out of all SPW
            if snr > snr_threshold and score < best_score:
                best_score = score
                best_fit = [amp1, mean1, std1, amp2, mean2, std2]
                best_err = [amp1_err, mean1_err, std1_err, amp2_err, mean2_err, std2_err]
                best_fluxes = fluxes
                best_freqs = frequencies
                best_label = label_fit
                best_region = region_label
                best_fits_file = fits_f

                if 'spw' in fits_f:
                    spw_label = fits_f.split('spw')[1][:1]
                else:
                    spw_label = 'SPW'
            sigma = 2
            if best_label == 'Single Gaussian':
                lower1 = mean1 - abs(sigma * std1)
                upper1 = mean1 + abs(sigma * std1)

                mask1_mom = (best_freqs >= lower1) & (best_freqs <= upper1)

                chan_inds1 = np.where(mask1_mom)[0]

                chans1 = f"{chan_inds1[0]}~{chan_inds1[-1]}" if len(chan_inds1) > 0 else None

                immoments(imagename=imagename, moments=[0,1,2], axis='spectral', region=region_file, chans=chans1, outfile=str(dirs["finalmommap"] / f"{region_label}_spw{spw_label}_{line}_G1_{mom_reg_size}arcsec_{timestamp}")) # create moment maps

            elif best_label == 'Double Gaussian':
                lower1 = mean1 - abs(sigma * std1)
                upper1 = mean1 + abs(sigma * std1)
                lower2 = mean2 - abs(sigma * std2)
                upper2 = mean2 + abs(sigma * std2)

                mask1_mom = (best_freqs >= lower1) & (best_freqs <= upper1)
                mask2_mom = (best_freqs >= lower2) & (best_freqs <= upper2)

                chan_inds1 = np.where(mask1_mom)[0]
                chan_inds2 = np.where(mask2_mom)[0]

                chans1 = f"{chan_inds1[0]}~{chan_inds1[-1]}" if len(chan_inds1) > 0 else None
                chans2 = f"{chan_inds2[0]}~{chan_inds2[-1]}" if len(chan_inds2) > 0 else None

                immoments(imagename=imagename, moments=[0,1,2], axis='spectral', region=region_file, chans=chans1, outfile=str(dirs["finalmommap"] / f"{region_label}_spw{spw_label}_{line}_G1_{mom_reg_size}arcsec_{timestamp}")) # create moment maps
                immoments(imagename=imagename, moments=[0,1,2], axis='spectral', region=region_file, chans=chans2, outfile=str(dirs["finalmommap"] / f"{region_label}_spw{spw_label}_{line}_G2_{mom_reg_size}arcsec_{timestamp}")) # create moment maps   
            
            color = ['maroon', 'darkolivegreen'] # ADD MORE COLORS IF MORE LINES

            hdul = fits.open(fits_f)
            hdr = hdul[0].header
            restfreq = hdr.get('RESTFRQ') / 1e9
            
            # Plot
            plt.subplot(plot_row, plot_col, i)
            plt.step(freq_GHz, fluxes*1e3, where='mid', color='gray', label='Data', alpha=0.8)
            plt.fill_between(freq_GHz, fluxes*1e3, step='mid', facecolor='dimgray', alpha=0.2)
            plt.plot(freq_GHz, fit*1e3, 'k', label='Gaussian Fit')
            plt.axhline(y=0, color='k', linestyle='--')
        

            for i, lines in enumerate(lines_split):
                rest_frequencies_z = rest_frequencies[lines] / (1 + region_z_phot_map[region_label])
                if np.min(freq_GHz) <= rest_frequencies_z <= np.max(freq_GHz):
                    plt.axvline(x=rest_frequencies_z, color=color[i], linestyle='--', label=f'{lines} at $z_{{\\mathrm{{SED}}}}$ = {region_z_phot_map[region_label]:.3f}')
            
            if 'spw' in fits_f:
                plt.title(f"SPW {fits_f.split('spw')[1][:1]}")
            else:
                pass

            ax = plt.gca()

            secax = ax.secondary_xaxis('top', functions=(freq_to_vel, vel_to_freq))
            secax.set_xlabel("Velocity [km/s]")
            

            plt.suptitle(f"{region_label} ({line})", fontsize=16)
            plt.xlabel("Frequency [GHz]")
            plt.ylabel("Flux Density [mJy]")
            plt.legend(loc='upper right')
        except RuntimeError:
            print(f"No fit for {fits_f}")

            # Plot
            plt.subplot(plot_row, plot_col, i)
            plt.step(freq_GHz, fluxes*1e3, where='mid', color='dimgray', label='Data', alpha=0.8)
            plt.fill_between(freq_GHz, fluxes*1e3, step='mid', facecolor='dimgray', alpha=0.2)
            plt.axhline(y=0, color='k', linestyle='--')

            for i, lines in enumerate(lines_split):
                rest_frequencies_z = rest_frequencies[lines] / (1 + region_z_phot_map[region_label])
                if np.min(freq_GHz) <= rest_frequencies_z <= np.max(freq_GHz):
                    plt.axvline(x=rest_frequencies_z, color=color[i], linestyle='--', label=f'{lines} at $z_{{\\mathrm{{SED}}}}$ = {region_z_phot_map[region_label]:.3f}')
            
            if 'spw' in fits_f:
                plt.title(f"SPW {fits_f.split('spw')[1][:1]}")
            else:
                pass
            plt.suptitle(f"{region_label} ({line})", fontsize=16)
            plt.xlabel("Frequency [GHz]")
            plt.ylabel("Flux Density [mJy]")
            plt.legend(loc='upper right')
            continue

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    #plt.savefig(str(dirs["fig"] / f"{region_label}_{timestamp}.png"), dpi='figure', bbox_inches="tight")
    plt.show()

    if best_fit:
        best_amp1, best_mean1, best_std1, best_amp2, best_mean2, best_std2 = best_fit
        best_amp1_err, best_mean1_err, best_std1_err, best_amp2_err, best_mean2_err, best_std2_err = best_err
        
        sigma_factor = 2
        niter = 1000

        if best_label == 'Single Gaussian':
            # Expected mean using photometric redshift
            expected_freq = [rest_frequencies[line] / (1 + region_z_phot_map[region_label]) for line in line_name]
            # Calculate relative distance between peak and expected means
            distance = [abs(best_mean1 - expected * 1e9) for expected in expected_freq]
            # Choose smallest distance
            min_index = [i for i, d in enumerate(distance) if d == min(distance)]
            best_index = random.choice(min_index) if len(min_index) > 1 else min_index[0]
            best_line = line_name[best_index]

            lower_bound1 = best_mean1 - abs(sigma_factor * best_std1)
            upper_bound1 = best_mean1 + abs(sigma_factor * best_std1)
            mask1 = (best_freqs >= lower_bound1) & (best_freqs <= upper_bound1)

            gaussian_fluxes1 = []
            integrated_fluxes1 = []
            fwhm1 = []
            z1 = []
            flux_density1 = []
            mean1_list = []

            gaussian_lum1 = []
            integrated_lum1 = []

            for i in range(niter):
                noise_flux = fluxes + np.random.normal(0, rms_chan_dict[best_line], size=fluxes.shape)
                #noise_flux = fluxes + np.random.normal(0, (mean_rms_dict[best_line][region_label]/1e3), size=fluxes.shape)

                try:
                    popt1_noise, pcov1_noise = curve_fit(one_gaussian, best_freqs[mask1], noise_flux[mask1], p0=[best_amp1, best_mean1, best_std1], maxfev=10000)
                    amp1_noise, mean1_noise, std1_noise = popt1_noise 
                    flux_density1.append(amp1_noise)
                    mean1_list.append(mean1_noise)

                    # Find Gaussian Line Flux
                    gaussian_flux1 = amp1_noise * std1_noise * math.sqrt(2.0 * math.pi) * erf(sigma_factor / math.sqrt(2.0))
                    gaussian_flux_vel1 = gaussian_flux1 * (3e5 / mean1_noise)

                    # Find Integrated Line Flux 
                    valid = ~np.isnan(noise_flux[mask1])
                    integrated_flux1 = np.trapz(noise_flux[mask1][valid], best_freqs[mask1][valid])
                    integrated_flux_vel1 = integrated_flux1 * (3e5 / mean1_noise)


                    # Find Line Width
                    fwhm1_noise = (2 * np.sqrt(2 * np.log(2)) * abs(std1_noise))
                    fwhm1_noise_vel = fwhm1_noise * (3e5 / (mean1_noise))

                    # Find Redshift
                    calc_z1 = (rest_frequencies[best_line] / (mean1_noise / 1e9)) - 1

                    # Luminosity
                    dl1 = cosmo.luminosity_distance(calc_z1).value
                    gauss_lum1 = 3.25 * 10**7 * (dl1**2 / ((1 + calc_z1) * (rest_frequencies[best_line])**2) ) * gaussian_flux_vel1
                    int_lum1 = 3.25 * 10**7 * (dl1**2 / ((1 + calc_z1) * (rest_frequencies[best_line])**2) ) * integrated_flux_vel1

                    gaussian_fluxes1.append(gaussian_flux_vel1)
                    gaussian_lum1.append(gauss_lum1)
                    if not np.isnan(integrated_flux_vel1):
                        integrated_fluxes1.append(integrated_flux_vel1)
                        integrated_lum1.append(int_lum1)
                    fwhm1.append(fwhm1_noise_vel)
                    z1.append(calc_z1)
                    
                except RuntimeError:
                    continue

            # Gaussian Line Flux 
            gaussian_fluxes1 = np.array(gaussian_fluxes1)
            gaussian_flux1_mean = np.mean(gaussian_fluxes1)
            gaussian_flux1_std = np.std(gaussian_fluxes1)

            gauss_flux[f"{region_label}"][f"{best_line}"] = gaussian_flux1_mean * 1e3
            gauss_flux_err[f"{region_label}"][f"{best_line}"] = gaussian_flux1_std * 1e3

            # Integrated Line Flux
            integrated_fluxes1 = np.array(integrated_fluxes1)
            integrated_flux1_mean = np.mean(integrated_fluxes1)
            integrated_flux1_std = np.std(integrated_fluxes1)

            int_flux[f"{region_label}"][f"{best_line}"] = integrated_flux1_mean * 1e3
            int_flux_err[f"{region_label}"][f"{best_line}"] = integrated_flux1_std * 1e3

            # Luminosity
            gaussian_lum1 = np.array(gaussian_lum1)
            gaussian_lum1_mean = np.mean(gaussian_lum1)
            gaussian_lum1_std = np.std(gaussian_lum1)

            gauss_lum[f"{region_label}"][f"{best_line}"] = gaussian_lum1_mean
            gauss_lum_err[f"{region_label}"][f"{best_line}"] = gaussian_lum1_std

            integrated_lum1 = np.array(integrated_lum1)
            integrated_lum1_mean = np.mean(integrated_lum1)
            integrated_lum1_std = np.std(integrated_lum1)

            int_lum[f"{region_label}"][f"{best_line}"] = integrated_lum1_mean
            int_lum_err[f"{region_label}"][f"{best_line}"] = integrated_lum1_std

            # Line Width
            fwhm1 = np.array(fwhm1)
            fwhm1_mean = np.mean(fwhm1)
            fwhm1_std = np.std(fwhm1)
            line_width[f"{region_label}"][f"{best_line}"] = fwhm1_mean
            line_width_err[f"{region_label}"][f"{best_line}"] = fwhm1_std

            # Redshift
            z1 = np.array(z1)
            z1_mean = np.mean(z1)
            z1_std = np.std(z1)

            region_line_z[f"{region_label}"][f"{best_line}"] = z1_mean
            region_line_z_err[f"{region_label}"][f"{best_line}"] = z1_std

            params = np.array([gaussian_fluxes1, integrated_fluxes1, fwhm1, z1]).T 
            labels = ["Gaussian Flux (mJy km/s)", "Integrated Flux (mJy km/s)", "FWHM (km/s)", "z"] 
            if params.ndim == 2 and params.shape[0] > params.shape[1]:
                fig = corner.corner(params, labels=labels, show_titles=True, title_fmt=".3f", quantiles=[0.25, 0.5, 0.75], title_kwargs={"fontsize": 10})
            else:
                print("Skipping corner plot")
                fig = None
            
            plt.show()

        elif best_label == 'Double Gaussian':
            if line_ratio is not None:
                expected_freq = [rest_frequencies[line] / (1 + region_z_phot_map[region_label]) for line in line_name]
                ## Peak 1
                # Calculate relative distance between peak and expected means
                distance1 = [abs(best_mean1 - expected * 1e9) for expected in expected_freq]
                # Choose smallest distance
                min_index1 = [i for i, d in enumerate(distance1) if d == min(distance1)]
                best_index1 = random.choice(min_index1) if len(min_index1) > 1 else min_index1[0]
                best_line1 = line_name[best_index1]
                ## Peak 2
                # Calculate relative distance between peak and expected means
                distance2 = [abs(best_mean2 - expected * 1e9) for expected in expected_freq]
                # Choose smallest distance
                min_index2 = [i for i, d in enumerate(distance2) if d == min(distance2)]
                best_index2 = random.choice(min_index2) if len(min_index2) > 1 else min_index2[0]
                best_line2 = line_name[best_index2]

            else:
                expected_freq = [rest_frequencies[line] / (1 + region_z_phot_map[region_label]) for line in line_name]

                candidates = []

                for i, expected in enumerate(expected_freq):
                    expected_Hz = expected * 1e9
                    distances_peak1 = abs(best_mean1 - expected_Hz)
                    distances_peak2 = abs(best_mean2 - expected_Hz)
                    
                    candidates.append((distances_peak1, "peak1", line_name[i]))
                    candidates.append((distances_peak2, "peak2", line_name[i]))

                    best_distance, best_peak, best_line = min(candidates, key=lambda x: x[0])
                    
                    if best_peak == "peak1":
                        best_line1 = best_line
                        best_line2 = None
                    else:
                        best_line1 = None
                        best_line2 = best_line
                    

            lower_bound1 = best_mean1 - abs(sigma_factor * best_std1)
            upper_bound1 = best_mean1 + abs(sigma_factor * best_std1)
            mask1 = (best_freqs >= lower_bound1) & (best_freqs <= upper_bound1)

            lower_bound2 = best_mean2 - abs(sigma_factor * best_std2)
            upper_bound2 = best_mean2 + abs(sigma_factor * best_std2)
            mask2 = (best_freqs >= lower_bound2) & (best_freqs <= upper_bound2)

            mask = mask1 | mask2

            flux_density1 = []
            mean1_list = []
            gaussian_fluxes1 = []
            integrated_fluxes1 = []
            gaussian_lum1 = []
            integrated_lum1 = []
            fwhm1 = []
            z1 = []

            flux_density2 = []
            mean2_list = []
            gaussian_fluxes2 = []
            integrated_fluxes2 = []
            gaussian_lum2 = []
            integrated_lum2 = []
            fwhm2 = []
            z2 = []

            for i in range(niter):
                if best_line1:
                    #noise_flux = fluxes + np.random.normal(0, (mean_rms_dict[best_line1][region_label]/1e3), size=fluxes.shape)
                    noise_flux = fluxes + np.random.normal(0, rms_chan_dict[best_line1], size=fluxes.shape)
                else:
                    #noise_flux = fluxes + np.random.normal(0, (mean_rms_dict[best_line2][region_label]/1e3), size=fluxes.shape)
                    noise_flux = fluxes + np.random.normal(0, rms_chan_dict[best_line2], size=fluxes.shape)
                try:
                    if line_ratio is not None:
                        popt2_noise, pcov2_noise = curve_fit(lambda x, a1, m1, s1, a2, s2: two_gaussian(x, a1, m1, s1, a2, s2, line_ratio), best_freqs[mask], noise_flux[mask], p0=[best_amp1, best_mean1, best_std1, best_amp2, best_std2], maxfev=10000)
                        amp1_noise, mean1_noise, std1_noise, amp2_noise, std2_noise = popt2_noise
                        mean2_noise = mean1_noise * line_ratio
                        flux_density1.append(amp1_noise)
                        flux_density2.append(amp2_noise)
                        mean1_list.append(mean1_noise)
                        mean2_list.append(mean2_noise)
                    else: 
                        popt2_noise, pcov2_noise = curve_fit(two_gaussian_free, best_freqs[mask], noise_flux[mask], p0=[best_amp1, best_mean1, best_std1, best_amp2, best_mean2, best_std2], maxfev=10000)
                        amp1_noise, mean1_noise, std1_noise, amp2_noise, mean2_noise, std2_noise = popt2_noise
                        flux_density1.append(amp2_noise)
                        mean1_list.append(mean2_noise)
                        

                    ## Peak 1
                    # Find Gaussian Line Flux
                    if best_line1 is not None and mean1_noise is not None:

                        gaussian_flux1 = amp1_noise * std1_noise * math.sqrt(2.0 * math.pi) * erf(sigma_factor / math.sqrt(2.0))
                        gaussian_flux_vel1 = gaussian_flux1 * (3e5 / mean1_noise)

                        # Find Integrated Line Flux 
                        valid = ~np.isnan(noise_flux[mask1])
                        integrated_flux1 = np.trapz(noise_flux[mask1][valid], best_freqs[mask1][valid])
                        integrated_flux_vel1 = integrated_flux1 * (3e5 / mean1_noise)

                        # Find Line Width
                        fwhm1_noise = (2 * np.sqrt(2 * np.log(2)) * abs(std1_noise))
                        fwhm1_noise_vel = fwhm1_noise * (3e5 / (mean1_noise))

                        # Find Redshift
                        calc_z1 = (rest_frequencies[best_line1] / (mean1_noise / 1e9)) - 1

                        # Luminosity
                        dl1 = cosmo.luminosity_distance(calc_z1).value
                        gauss_lum1 = 3.25 * 10**7 * (dl1**2 / ((1 + calc_z1) * (rest_frequencies[best_line1])**2) ) * gaussian_flux_vel1
                        int_lum1 = 3.25 * 10**7 * (dl1**2 / ((1 + calc_z1) * (rest_frequencies[best_line1])**2) ) * integrated_flux_vel1

                        gaussian_fluxes1.append(gaussian_flux_vel1)
                        gaussian_lum1.append(gauss_lum1)
                        if not np.isnan(integrated_flux_vel1):
                            integrated_fluxes1.append(integrated_flux_vel1)
                            integrated_lum1.append(int_lum1)
                        fwhm1.append(fwhm1_noise_vel)
                        z1.append(calc_z1)

                    else:
                        calc_z1 = gaussian_flux1 = gaussian_flux_vel1 = fwhm1_noise = fwhm1_noise_vel = None


                    ## Peak 2
                    if best_line2 is not None and mean2_noise is not None:

                        # Find Gaussian Line Flux
                        gaussian_flux2 = amp2_noise * std2_noise * math.sqrt(2.0 * math.pi) * erf(sigma_factor / math.sqrt(2.0))
                        gaussian_flux_vel2 = gaussian_flux2 * (3e5 / mean2_noise)

                        # Find Integrated Line Flux 
                        valid = ~np.isnan(noise_flux[mask2])
                        integrated_flux2 = np.trapz(noise_flux[mask2][valid], best_freqs[mask2][valid])
                        integrated_flux_vel2 = integrated_flux2 * (3e5 / mean2_noise)

                        # Find Line Width
                        fwhm2_noise = (2 * np.sqrt(2 * np.log(2)) * abs(std2_noise))
                        fwhm2_noise_vel = fwhm2_noise * (3e5 / (mean2_noise))

                        # Find Redshift
                        calc_z2 = (rest_frequencies[best_line2] / (mean2_noise / 1e9)) - 1

                        # Luminosity
                        dl2 = cosmo.luminosity_distance(calc_z2).value
                        gauss_lum2 = 3.25 * 10**7 * (dl2**2 / ((1 + calc_z2) * (rest_frequencies[best_line2])**2) ) * gaussian_flux_vel2
                        int_lum2 = 3.25 * 10**7 * (dl2**2 / ((1 + calc_z2) * (rest_frequencies[best_line2])**2) ) * integrated_flux_vel2

                        gaussian_fluxes2.append(gaussian_flux_vel2)
                        gaussian_lum2.append(gauss_lum2)
                        if not np.isnan(integrated_flux_vel2):
                            integrated_fluxes2.append(integrated_flux_vel2)
                            integrated_lum2.append(int_lum2)
                        fwhm2.append(fwhm2_noise_vel)
                        z2.append(calc_z2)
                    else:
                        calc_z2 = gaussian_flux2 = gaussian_flux_vel2 = fwhm2_noise = fwhm2_noise_vel = None
                    
                except RuntimeError:
                    continue
            
            ## Peak 1 
            # Gaussian Line Flux 
            gaussian_fluxes1 = np.array(gaussian_fluxes1)
            gaussian_flux1_mean = np.mean(gaussian_fluxes1)
            gaussian_flux1_std = np.std(gaussian_fluxes1)

            gauss_flux[f"{region_label}"][f"{best_line1}"] = gaussian_flux1_mean * 1e3
            gauss_flux_err[f"{region_label}"][f"{best_line1}"] = gaussian_flux1_std * 1e3

            # Integrated Line Flux
            integrated_fluxes1 = np.array(integrated_fluxes1)
            integrated_flux1_mean = np.mean(integrated_fluxes1)
            integrated_flux1_std = np.std(integrated_fluxes1)

            int_flux[f"{region_label}"][f"{best_line1}"] = integrated_flux1_mean * 1e3
            int_flux_err[f"{region_label}"][f"{best_line1}"] = integrated_flux1_std * 1e3

            # Luminosity
            gaussian_lum1 = np.array(gaussian_lum1)
            gaussian_lum1_mean = np.mean(gaussian_lum1)
            gaussian_lum1_std = np.std(gaussian_lum1)

            gauss_lum[f"{region_label}"][f"{best_line1}"] = gaussian_lum1_mean
            gauss_lum_err[f"{region_label}"][f"{best_line1}"] = gaussian_lum1_std

            integrated_lum1 = np.array(integrated_lum1)
            integrated_lum1_mean = np.mean(integrated_lum1)
            integrated_lum1_std = np.std(integrated_lum1)

            int_lum[f"{region_label}"][f"{best_line1}"] = integrated_lum1_mean
            int_lum_err[f"{region_label}"][f"{best_line1}"] = integrated_lum1_std

            # Line Width
            fwhm1 = np.array(fwhm1)
            fwhm1_mean = np.mean(fwhm1)
            fwhm1_std = np.std(fwhm1)
            line_width[f"{region_label}"][f"{best_line1}"] = fwhm1_mean
            line_width_err[f"{region_label}"][f"{best_line1}"] = fwhm1_std

            # Redshift
            z1 = np.array(z1)
            z1_mean = np.mean(z1)
            z1_std = np.std(z1)

            region_line_z[f"{region_label}"][f"{best_line1}"] = z1_mean
            region_line_z_err[f"{region_label}"][f"{best_line1}"] = z1_std

            params1 = np.array([gaussian_fluxes1, integrated_fluxes1, fwhm1, z1]).T 
            labels1 = ["Gaussian Flux (mJy km/s)", "Integrated Flux (mJy km/s)", "FWHM (km/s)", "z"] 
            if params1.ndim == 2 and params1.shape[0] > params1.shape[1]:
                fig = corner.corner(params1, labels=labels1, show_titles=True, title_fmt=".3f", quantiles=[0.25, 0.5, 0.75], title_kwargs={"fontsize": 10})
            else:
                print("Skipping corner plot")
                fig = None
            
            plt.show()

            ## Peak 2
            # Gaussian Line Flux 
            gaussian_fluxes2 = np.array(gaussian_fluxes2)
            gaussian_flux2_mean = np.mean(gaussian_fluxes2)
            gaussian_flux2_std = np.std(gaussian_fluxes2)

            gauss_flux[f"{region_label}"][f"{best_line2}"] = gaussian_flux2_mean * 1e3
            gauss_flux_err[f"{region_label}"][f"{best_line2}"] = gaussian_flux2_std * 1e3

            # Integrated Line Flux
            integrated_fluxes2 = np.array(integrated_fluxes2)
            integrated_flux2_mean = np.mean(integrated_fluxes2)
            integrated_flux2_std = np.std(integrated_fluxes2)

            int_flux[f"{region_label}"][f"{best_line2}"] = integrated_flux2_mean * 1e3
            int_flux_err[f"{region_label}"][f"{best_line2}"] = integrated_flux2_std * 1e3

            # Luminosity
            gaussian_lum2 = np.array(gaussian_lum2)
            gaussian_lum2_mean = np.mean(gaussian_lum2)
            gaussian_lum2_std = np.std(gaussian_lum2)

            gauss_lum[f"{region_label}"][f"{best_line2}"] = gaussian_lum2_mean
            gauss_lum_err[f"{region_label}"][f"{best_line2}"] = gaussian_lum2_std

            integrated_lum2 = np.array(integrated_lum2)
            integrated_lum2_mean = np.mean(integrated_lum2)
            integrated_lum2_std = np.std(integrated_lum2)

            int_lum[f"{region_label}"][f"{best_line2}"] = integrated_lum2_mean
            int_lum_err[f"{region_label}"][f"{best_line2}"] = integrated_lum2_std

            # Line Width
            fwhm2 = np.array(fwhm2)
            fwhm2_mean = np.mean(fwhm2)
            fwhm2_std = np.std(fwhm2)
            line_width[f"{region_label}"][f"{best_line2}"] = fwhm2_mean
            line_width_err[f"{region_label}"][f"{best_line2}"] = fwhm2_std

            # Redshift
            z2 = np.array(z2)
            z2_mean = np.mean(z2)
            z2_std = np.std(z2)

            region_line_z[f"{region_label}"][f"{best_line2}"] = z2_mean
            region_line_z_err[f"{region_label}"][f"{best_line2}"] = z2_std

            if best_line2 is not None and mean2_noise is not None:
                params2 = np.array([gaussian_fluxes2, integrated_fluxes2, fwhm2, z2]).T 
                labels2 = ["Gaussian Flux (mJy km/s)", "Integrated Flux (mJy km/s)", "FWHM (km/s)", "z"] 
                if params2.ndim == 2 and params2.shape[0] > params2.shape[1]:
                    fig = corner.corner(params2, labels=labels2, show_titles=True, title_fmt=".3f", quantiles=[0.25, 0.5, 0.75], title_kwargs={"fontsize": 10})
                else:
                    print("Skipping corner plot")
                    fig = None

                plt.show()

        else:
            print(f"No fit for {region_label} ({line})")
    
    else:
        print(f"No fit for {region_label} ({line})")
    
    csv_file = str(dirs["results"] / f"results_{timestamp}.csv")
    file_exists = os.path.isfile(csv_file)

    with open(csv_file, "w", newline='') as f:
        writer = csv.writer(f)

        header = ["ID", "RA (deg)", "DEC (deg)", "z (SED)"]
        for line in line_name:
            header.extend([
                f"{line} Aperture Size (arcsec)",
                f"z ({line})",
                f"z Error ({line})",
                f"{line} Gaussian Flux (mJy km/s)",
                f"{line} Gaussian Flux Error (mJy km/s)",
                f"{line} Line Flux (mJy km/s)",
                f"{line} Line Flux Error (mJy km/s)",
                f"{line} Gaussian Luminosity (K km/s pc^2)",
                f"{line} Gaussian Luminosity Error (K km/s pc^2)",
                f"{line} Line Luminosity (K km/s pc^2)",
                f"{line} Line Luminosity Error (K km/s pc^2)",
                f"{line} FWHM (km/s)",
                f"{line} FWHM Error (km/s)",
                f"{line} Mean RMS (mJy)",
            ])
        writer.writerow(header)

        for i in range(len(id)):
            row = [
                f"{float(id[i]):.0f}",
                f"{float(ra[i]):.7f}",
                f"{float(dec[i]):.7f}",
                f"{float(z_phot[i]):.8f}",
            ]

            for line in line_name:
                aperture_val = (
                    aperture_dict.get((int(id[i]), line_list), 'ND')
                    if "&" in line_list
                    else aperture_dict.get((int(id[i]), line), 'ND')
                )
                row.extend([
                    # Write Aperture size
                    decimals(aperture_val),

                    # Write z
                    decimals(region_line_z.get(str(id[i]), {}).get(line, 'ND')),
                    decimals(region_line_z_err.get(str(id[i]), {}).get(line, 'ND')),
                        
                    # Write Gaussian flux
                    decimals(gauss_flux.get(str(id[i]), {}).get(line, 'ND')),
                    decimals(gauss_flux_err.get(str(id[i]), {}).get(line, 'ND')),

                    #Write Integrated flux
                    decimals(int_flux.get(str(id[i]), {}).get(line, 'ND')),
                    decimals(int_flux_err.get(str(id[i]), {}).get(line, 'ND')),

                    # Write Gaussian flux
                    decimals(gauss_lum.get(str(id[i]), {}).get(line, 'ND')),
                    decimals(gauss_lum_err.get(str(id[i]), {}).get(line, 'ND')),

                    #Write Integrated flux
                    decimals(int_lum.get(str(id[i]), {}).get(line, 'ND')),
                    decimals(int_lum_err.get(str(id[i]), {}).get(line, 'ND')),

                    # Write Line Width
                    decimals(line_width.get(str(id[i]), {}).get(line, 'ND')),
                    decimals(line_width_err.get(str(id[i]), {}).get(line, 'ND')),

                    # Write RMS
                    decimals(mean_rms_dict.get(line, {}).get(id[i], 'ND'))
                ])

            writer.writerow(row)