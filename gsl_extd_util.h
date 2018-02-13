#ifndef GUARD_GSL_EXTD_UTIL_H
#define GUARD_GSL_EXTD_UTIL_H

#include <gsl/gsl_vector.h>

#define FILEOP_READ  (1 << 0)
#define FILEOP_WRITE (1 << 1)

void gsl_vector_exp(gsl_vector *x);
double gsl_vector_logsumexp(const gsl_vector *x);

#endif //GUARD_GSL_EXTD_UTIL_H
