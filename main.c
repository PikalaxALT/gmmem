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

#define MAGIC_SIZE 16

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
void (*output_fn)(FILE *) = NULL;
void (*run_fn)(void) = NULL;
void (*print_fn)(void) = NULL;
void (*free_fn)(void) = NULL;
unsigned int (*checksum_fn)(void) = NULL;

static unsigned int read_params(FILE * file) {
    char magic[MAGIC_SIZE];
    unsigned int checksum;
    fread(magic, sizeof(char), MAGIC_SIZE, file);
    if (memcmp(magic, PACKAGE_NAME, sizeof(PACKAGE_NAME) - 1) != 0) {
        fprintf(stderr, "error: file format unrecognized\n");
        exit(1);
    }
    if (strcmp(magic + sizeof(PACKAGE_NAME), "0.1") == 0) {
        fread(&checksum, sizeof checksum, 1, file);
        fread(&NSAMPS, sizeof NSAMPS, 1, file);
        fread(&NDIMS,  sizeof NDIMS,  1, file);
        fread(&NCOMPS, sizeof NCOMPS, 1, file);
    } else {
        fprintf(stderr, "error: invalid %s file version %s\n", PACKAGE_NAME, magic + sizeof(PACKAGE_NAME));
        exit(1);
    }
    return checksum;
}

static void write_params(FILE * file) {
    unsigned int checksum = checksum_fn();
    fwrite(PACKAGE_STRING, sizeof(char), MAGIC_SIZE, file);
    fwrite(&checksum, sizeof(checksum), 1, file);
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
    char *outfname = "./output.txt";
    bool fname_alloced = false;
    bool outfname_alloced = false;
    uint8_t load = 0;

    printf("%s\nContact: %s\n\n", PACKAGE_STRING, PACKAGE_BUGREPORT);

    if (argc < 2)
        usage();

    while ((ch = getopt(argc, argv, "hgrk:n:p:t:m:f:o:s:")) != -1) {
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
            case 'o':
                if (outfname_alloced)
                    free(outfname);
                outfname = strdup(optarg);
                if (outfname)
                    outfname_alloced = true;
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

            get_funcs();
            alloc_fn();
            gen_fn();
            write_params(file);
            write_fn(file);
            break;
        case FILEOP_READ:
            file = fopen(fname, "r");
            if (!file) {
                fprintf(stderr, "file io error: could not open %s for reading\n", fname);
                return EXIT_FAILURE;
            }

            unsigned int checksum = read_params(file);
            get_funcs();
            alloc_fn();
            read_fn(file);
            unsigned int checksum2 = checksum_fn();
            if (checksum != checksum2) {
                fprintf(stderr, "error: checksum mismatch. matrix might be corrupted!\n");
                exit(1);
            }
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
    file = fopen(outfname, "w");
    if (file) {
        output_fn(file);
        fclose(file);
    } else {
        fprintf(stderr, "warning: unable to output to file %s\n", outfname);
    }
    if (outfname_alloced)
        free(outfname);
    free_fn();
    gsl_rng_free(rng);
    return EXIT_SUCCESS;
}
