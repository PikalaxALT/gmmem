#ifndef GUARD_MAIN_H
#define GUARD_MAIN_H

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
extern void (*run_fn)(void);
extern void (*print_fn)(void);
extern void (*free_fn)(void);

#endif //GUARD_MAIN_H
