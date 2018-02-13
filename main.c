#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <sys/time.h>
#include <stdbool.h>
#include "gsl_extd_util.h"
#include "em1dim.h"
#include "emndim.h"
#include "main.h"

unsigned int NSAMPS = 1000000;
unsigned int NDIMS = 1;
unsigned int NCOMPS = 4;
unsigned int MAXITER = 1000;
double TOL = 1.0e-4;
gsl_rng *rng = NULL;
void (*alloc_fn)(void) = NULL;
void (*gen_fn)(void) = NULL;
void (*read_fn)(FILE *) = NULL;
void (*write_fn)(FILE *) = NULL;
void (*run_fn)(void) = NULL;
void (*print_fn)(void) = NULL;
void (*free_fn)(void) = NULL;

static inline void read_params(FILE * file) {
    fread(&NSAMPS, sizeof NSAMPS, 1, file);
    fread(&NDIMS,  sizeof NDIMS,  1, file);
    fread(&NCOMPS, sizeof NCOMPS, 1, file);
}

static inline void write_params(FILE * file) {
    fwrite(&NSAMPS, sizeof NSAMPS, 1, file);
    fwrite(&NDIMS,  sizeof NDIMS,  1, file);
    fwrite(&NCOMPS, sizeof NCOMPS, 1, file);
}

static inline void get_funcs(void) {
    NDIMS == 1 ? norm_init_funcs() : mvn_init_funcs();
}

int main(int argc, char * const argv[]) {
    rng = gsl_rng_alloc(gsl_rng_ranlxd2);
#ifdef __WINDOWS__
    gsl_rng_set(rng, (unsigned long)time(NULL));
#else
    FILE * urandom = fopen("/dev/urandom", "r");
    unsigned long seed = fgetc(urandom) | (fgetc(urandom) << 8) | (fgetc(urandom) << 16) | (fgetc(urandom) << 24);
    fclose(urandom);
    gsl_rng_set(rng, (unsigned long)time(NULL) ^ seed);
#endif // __WINDOWS__
    int ch;
    FILE * file;
    struct timeval stop, start;
    char *fname = "./gsl.mat";
    bool fname_alloced = false;
    unsigned int load = 0;

    printf("%s\nContact: %s\n\n", PACKAGE_STRING, PACKAGE_BUGREPORT);

    if (argc < 2)
        usage();

    while ((ch = getopt(argc, argv, "hgrk:n:p:t:m:f:s:")) != -1) {
        switch (ch) {
            case 'f':
                if (fname_alloced)
                    free(fname);
                fname = strdup(optarg);
                if (fname)
                    fname_alloced = true;
                else {
                    fprintf(stderr, "out of memory\n");
                    exit(EXIT_FAILURE);
                }
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

    switch (load) {
        case FILEOP_WRITE:
            file = fopen(fname, "w");
            if (!file) {
                fprintf(stderr, "file io error: could not open %s for writing\n", fname);
                return EXIT_FAILURE;
            }

            write_params(file);
            get_funcs();
            alloc_fn();
            gen_fn();
            write_fn(file);
            break;
        case FILEOP_READ:
            file = fopen(fname, "r");
            if (!file) {
                fprintf(stderr, "file io error: could not open %s for reading\n", fname);
                return EXIT_FAILURE;
            }

            read_params(file);
            get_funcs();
            alloc_fn();
            read_fn(file);
            break;
        default:
            usage();
    }
    fclose(file);
    if (fname_alloced)
        free(fname);

    printf("\nFitting %d %d-dimensional samples in a %d mixture... ", NSAMPS, NDIMS, NCOMPS);
    fflush(stdout);
    gettimeofday(&start, NULL);
    run_fn();
    gettimeofday(&stop, NULL);

    printf("Elapsed: %.6f\n\n", (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) * 1.0e-6);
    print_fn();
    free_fn();
    gsl_rng_free(rng);
    return EXIT_SUCCESS;
}
