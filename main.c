#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <memory.h>
#include <getopt.h>
#include <time.h>
#include "gsl_extd_util.h"
#include "em1dim.h"
#include "emndim.h"

int main(int argc, char * const argv[]) {
    unsigned int NSAMPS = 1000000;
    unsigned int NDIMS = 1;
    unsigned int NCOMPS = 4;
    unsigned int MAXITER = 1000;
    double TOL = 1.0e-4;

    gsl_rng *rng = gsl_rng_alloc(gsl_rng_ranlxd2);
#ifdef __WINDOWS__
    gsl_rng_set(rng, (unsigned long)time(NULL));
#else
    FILE * urandom = fopen("/dev/urandom", "r");
    unsigned long seed = fgetc(urandom) | (fgetc(urandom) << 8) | (fgetc(urandom) << 16) | (fgetc(urandom) << 24);
    fclose(urandom);
    gsl_rng_set(rng, (unsigned long)time(NULL) ^ seed);
#endif // __WINDOWS__
    int ch;
    char fname[256] = "./gsl.mat";
    unsigned int load = 0;

    if (argc < 2)
        usage();

    while ((ch = getopt(argc, argv, "hgrk:n:p:t:m:f:s:")) != -1) {
        switch (ch) {
            case 'f':
                strcpy(fname, optarg);
                break;
            case 'k':
                NCOMPS = strtoul(optarg, NULL, 0);
                break;
            case 'n':
                NSAMPS = strtoul(optarg, NULL, 0);
                break;
            case 'p':
                NDIMS = strtoul(optarg, NULL, 0);
                break;
            case 't':
                TOL = strtod(optarg, NULL);
                break;
            case 'm':
                MAXITER = strtoul(optarg, NULL, 0);
                break;
            case 'g':
                load |= FILEOP_WRITE;
                break;
            case 'r':
                load |= FILEOP_READ;
                break;
            case 's':
                gsl_rng_set(rng, strtoul(optarg, NULL, 0));
                break;
            case 'h':
                help();
            default:
                usage();
        }
    }

    if (NDIMS == 1) {
        norm_em_wrapper(load, NSAMPS, NCOMPS, MAXITER, TOL, fname, rng);
    } else {
        mvn_em_wrapper(load, NSAMPS, NDIMS, NCOMPS, MAXITER, TOL, fname, rng);
    }
    gsl_rng_free(rng);
    return 0;
}
