# Optimized Periodic Plus Smooth Decomposition (OPSD)

**Periodic-Smooth Decomposition (PSD) routine for General 2D Images**, optimized in C with Intel MKL for efficient handling of large datasets. This repository provides Python scripts for visualizing the results and utilities for manipulating binary data.

## Overview

The **Optimized Periodic Plus Smooth Decomposition (OPSD)** algorithm is designed to separate a 2D image into a periodic and smooth component. This method is particularly effective in removing boundary artifacts that occur when applying the 2D Discrete Fourier Transform (DFT) on non-periodic images treated as periodic. By decomposing the image, OPSD mitigates these edge artifacts, resulting in a cleaner frequency domain representation.

This implementation is optimized for large images using Intel's MKL (Math Kernel Library) for high efficiency.

## Citations

- Original algorithm by Lionel Moisan: [Periodic Plus Smooth Image Decomposition (2009)](https://helios2.mi.parisdescartes.fr/~moisan/papers/2009-11r.pdf)
- Optimization explanation by Faisal Mahmood, Märt Toots, Lars-Göran Öfverstedt and Ulf Skoglund: [Algorithm and Architecture Optimization for 2D Discrete
Fourier Transforms with Simultaneous Edge Artifact Removal](https://onlinelibrary.wiley.com/doi/10.1155/2018/1403181)

## Requirements

Ensure that Intel MKL is installed on your system before running the routine. You can download the MKL library from the Intel website:  
[Click here to access the Intel MKL download page](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html).

Before running, source the MKL environment variables using:
```bash
source /opt/intel/oneapi/setvars.sh
```

## Running the OPSD Algorithm in C

To use the OPSD algorithm, first obtain a `.bin` file of the image you wish to filter. Place this binary file in `bin/` and rename it to `data.bin`. Once the file is in place, you can run the algorithm.

The `.bin` and `.png` files will be saved with the same names in the following directories:

- Binary files: `bin/{dirname}/`
- Image files: `img/{dirname}/`

To execute the OPSD algorithm in C, navigate to `Routine_OPSD` and run the following command:

```bash
make run ARGS="<rows> <columns> <routine> <precision> <save_vector> <input> <dirname> <seed>"
```
### Argument Reference

- **Routine**:
  - `ccr`: compute complete routine.
  - `css`: compute shifted spectrums.
  - `cts`: compute traditional spectrums.

- **Input**:
  - `rb`: read binary.
  - `fm`: fill matrix.

- **Precision**:
  - `single`
  - `double`

- **Save_vector**:
  - `yes`
  - `no`

- **Dirname**:
  - `example`
  - `other dirname in bin`

## Plotting Spectrums and Images

The binary files for the images and frequency spectrums can be found in `bin/{dirname}/`, and their corresponding image files in `img/{dirname}/`.

C utility functions in the `utils` folder allow reading and saving of binary files. Additionally, the `plot_all.py` script in `plot/` can be used to plot images and spectrums. To use this script, provide the directory name, the number of rows, the number of columns, the color map, and whether to transpose the matrix, as shown below:
```bash
python3 plot_all.py Lenna 312 312 grey no

