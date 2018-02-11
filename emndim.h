#ifndef GUARD_EMNDIM_H
#define GUARD_EMNDIM_H

void mvn_em_wrapper(unsigned int load, unsigned int NSAMPS, unsigned int NDIMS, unsigned int NCOMPS, unsigned int MAXITER, double TOL, const char *fname, gsl_rng *rng);

#endif //GUARD_EMNDIM_H
