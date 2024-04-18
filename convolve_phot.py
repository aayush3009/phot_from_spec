import numpy as np
import matplotlib.pyplot as plt
import seaborn as sn

from astropy.io import fits
import astropy.units as u
import glob
import pandas as pd

from multiprocessing import Pool
import os

# plt.rcParams['font.family'] = 'stix'
# plt.rcParams['mathtext.fontset'] = 'stix'

plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10

plt.rcParams['axes.labelsize'] = 12

filters_hst = ["filters/HST_ACS/HST_ACS_WFC.F606W.dat",
            "filters/HST_ACS/HST_ACS_WFC.F775W.dat",
            "filters/HST_ACS/HST_ACS_WFC.F850LP.dat"]

filters_hst_names = np.array(["F606W", "F775W", "F850LP"])

filters_hst_eff = np.array([5809.26, 7652.44, 9004.99])


### JWST NIRCam wideband filters
filters_jwst = ["filters/NIRCam/JWST_NIRCam.F070W.dat",
           "filters/NIRCam/JWST_NIRCam.F090W.dat",
           "filters/NIRCam/JWST_NIRCam.F115W.dat",
           "filters/NIRCam/JWST_NIRCam.F140M.dat",
           "filters/NIRCam/JWST_NIRCam.F150W.dat",
           "filters/NIRCam/JWST_NIRCam.F182M.dat",
           "filters/NIRCam/JWST_NIRCam.F200W.dat",
           "filters/NIRCam/JWST_NIRCam.F210M.dat",
           "filters/NIRCam/JWST_NIRCam.F277W.dat",
           "filters/NIRCam/JWST_NIRCam.F300M.dat",
           "filters/NIRCam/JWST_NIRCam.F335M.dat",
           "filters/NIRCam/JWST_NIRCam.F356W.dat",
           "filters/NIRCam/JWST_NIRCam.F360M.dat",
           "filters/NIRCam/JWST_NIRCam.F410M.dat",
           "filters/NIRCam/JWST_NIRCam.F410M.dat",
           "filters/NIRCam/JWST_NIRCam.F444W.dat",
           "filters/NIRCam/JWST_NIRCam.F460M.dat",
           "filters/NIRCam/JWST_NIRCam.F480M.dat"
          ]

filters_jwst_names = np.array(["F070W", "F090W", "F115W", "F140M", "F150W", "F182M", "F200W", "F210M", "F277W", 
        "F300M", "F335M", "F356W", "F360M", "F410M", "F430M", "F444W", "F460M", "F480M"])

filters_jwst_eff = np.array([7039.12, 9021.53, 11542.61, 14053.23, 15007.44, 18451.67, 19886.48, 20954.51, 27617.40, 
        29891.21, 33620.67, 35683.62, 36241.76, 40822.38, 42812.58, 44043.15, 46299.28, 48181.95])

# Concatenate HST and JWST filter information into a single array
filters = np.concatenate((filters_hst, filters_jwst), axis=0)
filt_eff = np.concatenate((filters_hst_eff, filters_jwst_eff), axis=0)
filt_names = np.concatenate((filters_hst_names, filters_jwst_names), axis=0)


SURVEY = "GS-HST-MEDIUM"
VERSION = "3.1"
ROOTDIR = f"/Users/aayushsaxena/Desktop/Oxford/JADES/{SURVEY}/Final_products_v{VERSION}/prism_clear/"
TABLEDIR = ROOTDIR+"predicted_phot/"
PLOTDIR = TABLEDIR+"plots/"


### Function for plotting spectra
def make_plot(filt_fluxes, filt_errors, spec_wave, spec_flux, source_ID):

    ### Choose blue to red color palette
    palette = sn.color_palette("Spectral", len(filters))
    palette.reverse()
    sn.set_palette(palette)

    fig = plt.figure(figsize=(6,4))
    plt.grid(alpha=0.4, zorder=0)

    plt.plot(spec_wave, spec_flux, c='grey', alpha=1, zorder=10)

    ### Filter phot
    # plt.scatter(filt_eff, filt_fluxes_ujy, marker='o', s=80, c='C2', ec='k', zorder=11)

    ### Plot the filter transmission curves
    for i in range(len(filters)):
        plt.scatter(filt_eff[i], filt_fluxes[i], marker='o', s=80, ec='k', zorder=11)
        plt.errorbar(filt_eff[i], filt_fluxes[i], yerr=filt_errors[i], fmt="None", capsize=4, c='k', zorder=10)

        ### Plot transmission curve
        filt_dat = np.loadtxt(filters[i])
        plt.plot(filt_dat[:,0], filt_dat[:,1], alpha=0.6, zorder=2)

        # plt.ylim(-0.7e-23, 5.5e-23)

        plt.figtext(0.2, 0.85, source_ID, fontsize=12)

    plt.xlabel(r"Wavelength ($\mu$m)")
    plt.ylabel(r"$F_\nu$ ($\mu$Jy)")

    plt.xscale('log')
    plt.yscale('log')

    plt.ylim(1e-3, 1)

    plt.xticks([7000, 10000, 15000, 20000, 30000, 40000, 50000],
              ["0.7", "1.0", "1.5", "2.0", "3.0", "4.0", "5.0"])

    plt.tight_layout()

    if not os.path.isdir(PLOTDIR):
        os.mkdir(PLOTDIR)

    plt.savefig(PLOTDIR+f"{source_ID}-predicted-phot-uJy.png", dpi=300)

    plt.close()

    return


### function for convolution
def convolve_spec(spec):
    ### Load in a R100 spectrum and define all the arrays
    hdu = fits.open(spec)
    # hdu.info()

    wave = hdu["WAVELENGTH"].data * 1e10 * u.AA ### Angstrom
    flux = hdu["DATA"].data * (u.W/u.m**2/u.m)
    flux = flux.to(u.W/u.m**2/u.AA)
    err = hdu["ERR"].data * (u.W/u.m**2/u.m)
    err = err.to(u.W/u.m**2/u.AA)

    ### Remove NaN values
    flux = np.nan_to_num(flux)
    err = np.nan_to_num(err)

    ### Split path to obtain file name
    head_tail = os.path.split(spec)
    specname = head_tail[1]

    ### Extract the source ID. Make distinction between HST and JWST selected source IDs
    source_ID = specname.split("_")[0]

    print("Analyzing source ID: ", source_ID)

    ### Convert the fluxes to uJy before convolution
    flux_ujy = []
    err_ujy = []
    for j in range(len(flux)):
        ujy_flux = flux[j].to(u.Jy, equivalencies=u.spectral_density(wave[j])) * 1e6

        ujy_err = err[j].to(u.Jy, equivalencies=u.spectral_density(wave[j])) * 1e6

        flux_ujy.append(ujy_flux.value)
        err_ujy.append(ujy_err.value)

    flux_ujy = np.array(flux_ujy)
    err_ujy = np.array(err_ujy)

    filt_fluxes_w = []
    filt_fluxes_ujy = []

    filt_errs_w = []
    filt_errs_ujy = []

    for i in range(len(filters)):
        filt = filters[i]
        eff_wavelength = filt_eff[i]

        filt_dat = np.loadtxt(filt)
        filt_w = filt_dat[:,0] * u.AA
        filt_t = filt_dat[:,1]

        
        ### Convolution

        ### Interpolate the filter transmission on to the wavelength grid
        filt_conv_t = np.interp(wave, filt_w, filt_t)

        ### Do a weighted convolution
        ### Normalize the filter profile so it sums to one
        filt_flux_w = np.sum(np.dot(flux, (filt_conv_t/np.sum(filt_conv_t))))
        filt_flux_ujy = np.sum(np.dot(flux_ujy, (filt_conv_t/np.sum(filt_conv_t))))

        ### Errors
        filt_err_w = np.sum(np.dot(err, (filt_conv_t/np.sum(filt_conv_t))))
        filt_err_ujy = np.sum(np.dot(err_ujy, (filt_conv_t/np.sum(filt_conv_t))))

        ### Append
        filt_fluxes_w.append(filt_flux_w.value)
        filt_fluxes_ujy.append(filt_flux_ujy)

        filt_errs_w.append(filt_err_w.value)
        filt_errs_ujy.append(filt_err_ujy)

    ### Convert to np arrays
    filt_fluxes_w = np.array(filt_fluxes_w)
    filt_fluxes_ujy = np.array(filt_fluxes_ujy)

    filt_errs_w = np.array(filt_errs_w)
    filt_errs_ujy = np.array(filt_errs_ujy)

    ### Store fluxes in uJy
    bb_fluxes = {filt_names[i]: filt_fluxes_ujy[i] for i in range(len(filt_names))}
    bb_errors = {filt_names[i]+"_e": filt_errs_ujy[i] for i in range(len(filt_names))}

    ### Create a dict with all the info for this source
    ID_dict = {"ID": source_ID}

    info_col = {**ID_dict, **bb_fluxes, **bb_errors}

    ### Make plot for the object
    make_plot(filt_fluxes_ujy, filt_errs_ujy, wave, flux_ujy, source_ID)

    return(info_col)

def main():
    speclist = glob.glob(ROOTDIR+"*_1D.fits")

    # Check if predicted photometry directory exists
    if not os.path.isdir(TABLEDIR):
        os.mkdir(TABLEDIR)
    

    ### Initialize a list to keep track of photometry
    # conv_phot = []

    # for spec in speclist:
    #     ### Calculate photometry
    #     conv_phot.append(convolve_spec(spec, filters))

    ### Parallelize
    with Pool() as pool:
        conv_phot = pool.map(convolve_spec, speclist)

    ### Save as pd dataframe
    phot_table = pd.DataFrame.from_dict(conv_phot)


    phot_table.to_csv(TABLEDIR+f"{SURVEY}-v{VERSION}-convolved_nircam_photometry.csv")

    print("All photometry done and plots created!")

if __name__=="__main__":
    main()
