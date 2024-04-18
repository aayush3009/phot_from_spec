# Spectral Convolution for Photometric Analysis

This Python script is designed to perform spectral convolution for photometric analysis using filters from the Hubble Space Telescope (HST) and the James Webb Space Telescope (JWST) NIRCam. It takes 1D spectra as input, convolves them with filter transmission curves, and calculates photometry.

## Requirements

- `numpy`
- `matplotlib`
- `seaborn`
- `astropy`
- `pandas`

You can install these dependencies using `pip` and the provided `requirements.txt` file:

```bash
pip install -r requirements.txt
```

## Usage

1. Clone this repository or download the Python script.
2. Install the required dependencies.
3. Update the paths to filter transmission curve files and input spectra in the script (convolution_photometry.py).
4. Run the script using the following command:

```bash
python convolution_photometry.py
```

## Description

- The script first loads the necessary libraries and defines default plotting parameters.
- Filter transmission curve files and other input parameters are specified.
- Functions for plotting spectra (make_plot) and spectral convolution (convolve_spec) are defined.
- The main function (main) processes input spectra files, performs convolution, generates photometry, and saves the results as a CSV file.
- The script utilizes multiprocessing for parallel processing to improve performance.

## File Structure

- `convolution_photometry.py`: Main Python script.
- `requirements.txt`: List of required Python libraries.
- `README.md`: This documentation file.

## Authors

- Aayush Saxena (Oxford)