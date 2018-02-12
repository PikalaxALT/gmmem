#ifndef GUARD_GSL_EXTD_UTIL_H
#define GUARD_GSL_EXTD_UTIL_H

#include <gsl/gsl_vector.h>
#include <stdnoreturn.h>

#define FILEOP_READ  (1 << 0)
#define FILEOP_WRITE (1 << 1)

void gsl_vector_exp(gsl_vector *x);
double gsl_vector_logsumexp(const gsl_vector *x);

static inline noreturn void usage() {
    fprintf(stderr, "usage:\n"
        "\t./gmmem {-g/-r} [-f FILENAME] [-s SEED] [-n NSAMPS] [-p NDIMS] [-k NCOMPS] [-t TOL] [-m MAXITER]\n");
    exit(1);
}

static inline noreturn void help() {
    printf("usage:\n"
               "    ./gmmem {-g/-r} [-f FILENAME] [-s SEED] [-n NSAMPS] [-p NDIMS] [-k NCOMPS]\n"
               "                    [-t TOL] [-m MAXITER]\n"
               "\n"
               "-h - Print this help and exit.\n"
               "\n"
               "-g - Run in Generate mode.  This mode randomly samples a set of means and\n"
               "(co)variances from which to generate the data mixture.  It then saves the data\n"
               "mixture to the file specified by -f.\n"
               "\n"
               "-r - Run in Read mode.  This mode initializes the data matrix from the file\n"
               "specified by -f.\n"
               "\n"
               "-f FILENAME - File to use for reading or writing the data matrix.  By default,\n"
               "the file gsl.mat will be used.\n"
               "\n"
               "-s SEED - A number with which to seed the RNG, for reproducible results.  By\n"
               "default, this seed is 0 on Windows machines, and derived from /dev/urandom on\n"
               "UNIX machines.\n"
               "\n"
               "-n NSAMPS - The number of samples for the data matrix.  By default, this is 1\n"
               "million (1000000).\n"
               "\n"
               "-p NDIMS - The number of data features.  By default, this is 1.\n"
               "\n"
               "-k NCOMPS - The number of Gaussians to fit.  By default, this is 4.\n"
               "\n"
               "-t TOL - This parameter governs the absolute convergence criterion: If the\n"
               "total squared difference between the means in successive steps falls below TOL,\n"
               "the EM considers itself converged.  By default, this is 1e-4.\n"
               "\n"
               "-m MAXITER - If the above convergence criterion has not been satisfied within\n"
               "MAXITER iterations of EM, the iterator will exit and return its current state.\n"
               "By default, this is 1000.\n"
               "\n");
    exit(0);
}

#endif //GUARD_GSL_EXTD_UTIL_H
