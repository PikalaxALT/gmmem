# gmmem

## Purpose
This utility will read a GSL-compliant matrix or vector of size NxP (Nx1) and fit a K-component Gaussian Mixture Model (GMM) using the Expectation Maximization (EM) algorithm.

## Prerequisites
The following packages are required:

    autoconf, autoconf-archive, automake, make, libgsl, libblas, gcc

Most of these can be installed using your distribution's package manager.

Additionally, your version of GSL must be at least 2.4.

## To build
    autoreconf --install --force
    ./configure
    make

## To run
    ./gmmem {-h/-g/-r} [-f FILENAME] [-o OUTFILENAME] [-s SEED] [-n NSAMPS] [-p NDIMS] [-k NCOMPS] [-t TOL] [-m MAXITER]

`-h` - Print this help and exit.

`-g` - Run in Generate mode.  This mode randomly samples a set of means and (co)variances from which to generate the data mixture.  It then saves the data mixture to the file specified by `-f`.

`-r` - Run in Read mode.  This mode initializes the data matrix from the file specified by `-f`.  This is useful for if you want to rerun EM on the same data to find a new local optimum.

`-f FILENAME` - File to use for reading or writing the data matrix.  By default, the file `gsl.mat` will be used.

`-o OUTFILENAME` - Text file to which the final fitted params are to be written.  By default, the file `output.txt` will be used

`-s SEED` - A number with which to seed the RNG, for reproducible results.  By default, this seed is 0 on Windows machines, and derived from `/dev/urandom` on UNIX machines.

`-n NSAMPS` - The number of samples for the data matrix.  By default, this is 1 million (1000000).

`-p NDIMS` - The number of data features.  By default, this is 1.

`-k NCOMPS` - The number of Gaussians to fit.  By default, this is 4.

`-t TOL` - This parameter governs the absolute convergence criterion: If the total squared difference between the means in successive steps falls below TOL, the EM considers itself converged.  By default, this is 1e-4.

`-m MAXITER` - If the above convergence criterion has not been satisfied within MAXITER iterations of EM, the iterator will exit and return its current state.  By default, this is 1000.

## File format description

### Data matrix
The data matrix file contains the following components:

- Null-terminated magic string of length 16 (must be "gmmem " followed by the version number)
- Simple checksum of length equal to the system's integer size.  For 64-bit CPUs, this is 8 bytes; for 32-bit systems, this is 4 bytes.
- Integer number of samples (NSAMPS)
- Integer number of dimensions (NDIMS)
- Integer number of Gaussians (NCOMPS)
- The samples as an NSAMPS-by-NDIMS double-precision float matrix.
