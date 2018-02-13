#include <math.h>
#include <float.h>
#include "gsl_extd_util.h"

void gsl_vector_exp(gsl_vector *x) {
    for (unsigned int i = 0; i < x->size; i++) {
        gsl_vector_set(x, i, exp(gsl_vector_get(x, i)));
    }
}

double gsl_vector_logsumexp(const gsl_vector *x) {
    double amax = gsl_vector_max(x);
    double res = 0.0;
    for (unsigned int i = 0; i < x->size; i++) {
        res += exp(gsl_vector_get(x, i) - amax);
    }
    return log(res + DBL_MIN) + amax;
}
