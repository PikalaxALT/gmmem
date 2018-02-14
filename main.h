#ifndef GUARD_MAIN_H
#define GUARD_MAIN_H

#include <stdnoreturn.h>
#include <string.h>

extern unsigned int NSAMPS;
extern unsigned int NCOMPS;
extern unsigned int NDIMS;
extern double TOL;
extern unsigned int MAXITER;
extern gsl_rng *rng;
extern void (*alloc_fn)(void);
extern void (*gen_fn)(void);
extern void (*read_fn)(FILE *);
extern void (*write_fn)(FILE *);
extern void (*output_fn)(FILE *);
extern void (*run_fn)(void);
extern void (*print_fn)(void);
extern void (*free_fn)(void);

static inline noreturn void usage() {
    fprintf(stderr, "usage:\n");

    FILE * readme = fopen("README", "r");
    char * buffer = NULL;
    size_t size = 0;

    while (getline(&buffer, &size, readme) != -1 && strcmp(buffer, "## To run\n") != 0);
    getline(&buffer, &size, readme);
    fprintf(stderr, "%s", buffer);
    free(buffer);
    exit(1);
}

static inline noreturn void help() {
    FILE * readme = fopen("README", "r");
    char * buffer = NULL;
    size_t size = 0;

    while (getline(&buffer, &size, readme) != -1 && strcmp(buffer, "## To run\n") != 0);
    while (getline(&buffer, &size, readme) != -1) {
        printf("%s", buffer);
    }
    free(buffer);
    exit(0);
}

#endif //GUARD_MAIN_H
