# gmmem

## Purpose
This utility will read a GSL-compliant matrix or vector of size NxP (Nx1) and fit a K-component Gaussian Mixture Model (GMM) using the Expectation Maximization (EM) algorithm.

## To build
    make

## To run
    ./gmmem {-g/-r} [-f FILENAME] [-s SEED] [-n NSAMPS] [-p NDIMS] [-k NCOMPS] [-t TOL] [-m MAXITER]

`-g` - Run in Generate mode.  This mode randomly samples a set of means and (co)variances from which to generate the data mixture.  It then saves the data mixture to the file specified by `-f`.

`-r` - Run in Read mode.  This mode initializes the data matrix from the file specified by `-f`.

`-f FILENAME` - File to use for reading or writing the data matrix.  By default, the file `gsl.mat` will be used.

`-s SEED` - A number with which to seed the RNG, for reproducible results.  By default, this seed is 0 on Windows machines, and derived from `/dev/urandom` on UNIX machines.

`-n NSAMPS` - The number of samples for the data matrix.  By default, this is 1 million (1000000).

`-p NDIMS` - The number of data features.  By default, this is 1.

`-k NCOMPS` - The number of Gaussians to fit.  By default, this is 4.

`-t TOL` - This parameter governs the absolute convergence criterion: If the total squared difference between the means in successive steps falls below TOL, the EM considers itself converged.  By default, this is 1e-4.

`-m MAXITER` - If the above convergence criterion has not been satisfied within MAXITER iterations of EM, the iterator will exit and return its current state.  By default, this is 1000.
