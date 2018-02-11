#include <stdbool.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <assert.h>
#include "main.h"
#include "em1dim.h"
#include "gsl_extd_util.h"

static gsl_vector *data = NULL;
static gsl_matrix *indic = NULL;
static double *mu = NULL;
static double *sigma = NULL;
static double *prev_mu = NULL;
static double *prev_sigma = NULL;

static void e_step(gsl_vector *work) {
    // work: NCOMPS
    for (unsigned int i = 0; i < data->size; i++) {
        double xval = gsl_vector_get(data, i);
        for (unsigned int k = 0; k < indic->size2; k++) {
            gsl_vector_set(work, k, log(gsl_ran_gaussian_pdf(xval - mu[k], sigma[k]) + DBL_MIN));
        }
        gsl_vector_add_constant(work, -gsl_vector_logsumexp(work));
        gsl_vector_exp(work);
        gsl_matrix_set_row(indic, i, work);
    }
}

static void m_step(void) {
    for (unsigned int k = 0; k < indic->size2; k++) {
        mu[k] = 0.0;
        sigma[k] = 0.0;
        double denom = 0.0;

        // denom = 1.0 / (sum_i I_ik)
        // mu_k = denom * sum_i I_ik * x_i
        for (unsigned int i = 0; i < data->size; i++) {
            double tmp = gsl_matrix_get(indic, i, k);
            denom += tmp;
            mu[k] += tmp * gsl_vector_get(data, i);
        }
        mu[k] /= denom;

        // sigmasq_k = denom * sum_i I_ik * (x_i - mu_k)^2
        for (unsigned int i = 0; i < data->size; i++) {
            double tmp = gsl_matrix_get(indic, i, k);
            double xval = gsl_vector_get(data, i) - mu[k];
            sigma[k] += xval * xval * tmp;
        }
        sigma[k] = sqrt(sigma[k] / denom);
    }
}

static bool converged(void) {
    double diff = 0.0;
    for (unsigned int k = 0; k < NCOMPS; k++) {
        prev_mu[k] -= mu[k];
        diff += prev_mu[k] * prev_mu[k];
    }
    return diff < TOL;
}

static void fit_em(void) {
    // Initialize work buffers
    gsl_vector *work = gsl_vector_calloc(NCOMPS);
    gsl_vector *xdiff = gsl_vector_alloc(data->size);
    gsl_vector_set_all(xdiff, 1.0 / data->size);

    // Initialize buffers for the previous parameters state
    prev_mu = calloc(NCOMPS, sizeof(double));
    prev_sigma = calloc(NCOMPS, sizeof(double));

    assert(prev_mu && prev_sigma); // catch calloc failure

    // Initialize means
    for (unsigned int k = 0; k < NCOMPS; k++) {
        unsigned int i;
        unsigned int kk;
        double sum = 0.0;

        for (i = 0; i < data->size; i++) {
            double *ptr = gsl_vector_ptr(xdiff, i);
            double tmp = gsl_vector_get(data, i);
            for (kk = 0; kk < k; kk++) {
                double tmp2 = (tmp - mu[kk]) * (tmp - mu[kk]);
                if (tmp2 < *ptr)
                    *ptr = tmp2;
            }
            sum += *ptr;
        }
        gsl_vector_scale(xdiff, 1.0 / sum);
        for (i = 0; i < data->size - 1; i++) {
            *gsl_vector_ptr(xdiff, i + 1) += gsl_vector_get(xdiff, i);
        }
        double rval = gsl_rng_uniform(rng);
        for (i = 0; i < data->size; i++) {
            if (rval <= gsl_vector_get(xdiff, i)) {
                break;
            }
        }
        mu[k] = gsl_vector_get(data, i);
        sigma[k] = 1.0;
    }

    // Initialize I, the indicator matrix
    indic = gsl_matrix_calloc(data->size, NCOMPS);

    // Main EM Loop
    unsigned int iter;
    for (iter = 0; iter < MAXITER; iter++) {
#ifdef DEBUG
        printf("\nIteration %d:", iter);
        for (unsigned int k = 0; k < NCOMPS; k++) {
            printf("\t%.4f", mu[k]);
        }
        fflush(stdout);
#endif
        memcpy(prev_mu, mu, NCOMPS * sizeof(double));
        memcpy(prev_sigma, sigma, NCOMPS * sizeof(double));
        e_step(work);
        m_step();
        if (converged()) {
            break;
        }
    }
    if (iter == MAXITER)
        printf("\nEM did not converge within %d iterations\n", MAXITER);
    else
        printf("\nEM mvn_em_converged in %d iterations\n", iter + 1);

    // Free work buffers
    gsl_vector_free(xdiff);
    free(prev_sigma);
    free(prev_mu);
    gsl_matrix_free(indic);
    gsl_vector_free(work);
}

void norm_em_wrapper (void) {
    data = gsl_vector_calloc(NSAMPS);
    mu = calloc(NCOMPS, sizeof(double));
    sigma = calloc(NCOMPS, sizeof(double));

    switch (load) {
        case FILEOP_WRITE: {
            FILE * file = fopen(fname, "w");
            if (!file)
                usage();

            // Generate parameters
            for (unsigned int k = 0; k < NCOMPS; k++) {
                mu[k] = gsl_ran_flat(rng, -100.0, 100.0);
                printf("Mean %d:\t%.4f\n", k, mu[k]);
                sigma[k] = sqrt(gsl_ran_flat(rng, 1.0, 25.0));
            }

            // Generate data
            for (unsigned int i = 0; i < NSAMPS; i++) {
                unsigned long k = gsl_rng_uniform_int(rng, NCOMPS);
                gsl_vector_set(data, i, gsl_ran_gaussian(rng, sigma[k]) + mu[k]);
            }

            // Write data
            gsl_vector_fwrite(file, data);
            fclose(file);
            break;
        }
        case FILEOP_READ: {
            FILE * file = fopen(fname, "r");
            if (!file)
                usage();
            gsl_vector_fread(file, data);
            fclose(file);
            break;
        }
        default:
            usage();
    }

    printf("\nFitting...");
    struct timeval stop, start;
    gettimeofday(&start, NULL);
    fit_em();
    gettimeofday(&stop, NULL);

    printf("Elapsed: %.6f\n\n", (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) * 1.0e-6);
    for (unsigned int k = 0; k < NCOMPS; k++) {
        printf("New mean %d:\t%.4f\n", k, mu[k]);
    }
    free(sigma);
    free(mu);
    gsl_vector_free(data);
}
