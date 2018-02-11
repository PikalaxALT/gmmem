#include <stdbool.h>
#include <sys/time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics.h>
#include "emndim.h"
#include "gsl_extd_util.h"
#include "main.h"

static gsl_matrix *data = NULL;
static gsl_matrix *indic = NULL;
static gsl_vector **mu = NULL;
static gsl_matrix **sigma = NULL;
static gsl_vector **prev_mu = NULL;
static gsl_matrix **prev_sigma = NULL;

static void e_step(gsl_vector *work1, gsl_vector *work2, gsl_vector *work3) {
    // work1: p
    // work2: p
    // work3: NCOMPS
    for (unsigned int i = 0; i < data->size1; i++) {
        gsl_matrix_get_row(work1, data, i);
        for (unsigned int k = 0; k < indic->size2; k++) {
            gsl_ran_multivariate_gaussian_log_pdf(work1, mu[k], sigma[k], gsl_vector_ptr(work3, k), work2);
        }
        gsl_vector_add_constant(work3, -gsl_vector_logsumexp(work3));
        gsl_vector_exp(work3);
        gsl_matrix_set_row(indic, i, work3);
    }
}

static void m_step(gsl_vector *work1) {
    for (unsigned int k = 0; k < indic->size2; k++) {
        gsl_vector_set_zero(mu[k]);
        gsl_matrix_set_zero(sigma[k]);
        double denom = 0.0;

        // denom = 1.0 / (sum_i I_ik)
        // mu_k = denom * sum_i I_ik * x_i
        for (unsigned int i = 0; i < data->size1; i++) {
            double tmp = gsl_matrix_get(indic, i, k);
            denom += tmp;
            gsl_matrix_get_row(work1, data, i);
            gsl_vector_scale(work1, tmp);
            gsl_vector_add(mu[k], work1);
        }
        denom = 1.0 / denom;
        gsl_vector_scale(mu[k], denom);

        // sigmasq_k = denom * sum_i I_ik * (x_i - mu_k) * (x_i - mu_k)^T
        for (unsigned int i = 0; i < data->size1; i++) {
            double tmp = gsl_matrix_get(indic, i, k);
            gsl_matrix_get_row(work1, data, i);
            gsl_vector_sub(work1, mu[k]);
            gsl_blas_dger(tmp, work1, work1, sigma[k]);
        }
        gsl_matrix_scale(sigma[k], denom);
        gsl_linalg_cholesky_decomp(sigma[k]); // because the multivariate gaussian requires a Cholesky decomposition
    }
}

static bool converged(void) {
    double diff = 0.0;
    for (unsigned int k = 0; k < NCOMPS; k++) {
        gsl_vector_sub(prev_mu[k], mu[k]);
        diff += gsl_stats_tss_m(prev_mu[k]->data, prev_mu[k]->stride, prev_mu[k]->size, 0.0);
    }
    return diff < TOL;
}

static void fit_em(void) {
    // Initialize work buffers
    gsl_vector *work1 = gsl_vector_calloc(data->size2);
    gsl_vector *work2 = gsl_vector_calloc(data->size2);
    gsl_vector *work3 = gsl_vector_calloc(NCOMPS);
    gsl_vector *xdiff = gsl_vector_calloc(data->size1);
    gsl_vector_set_all(xdiff, 1.0 / data->size1);

    // Initialize buffers for the previous parameters state
    prev_mu = calloc(NCOMPS, sizeof(gsl_vector *));
    prev_sigma = calloc(NCOMPS, sizeof(gsl_matrix *));

    // Initialize means
    for (unsigned int k = 0; k < NCOMPS; k++) {
        prev_mu[k] = gsl_vector_calloc(data->size2);
        prev_sigma[k] = gsl_matrix_calloc(data->size2, data->size2);

        unsigned int i;
        unsigned int j;
        unsigned int kk;
        double sum = 0.0;

        for (i = 0; i < data->size1; i++) {
            double *ptr = gsl_vector_ptr(xdiff, i);
            gsl_matrix_get_row(work1, data, i);
            for (kk = 0; kk < k; kk++) {
                gsl_vector_memcpy(work2, work1);
                gsl_vector_sub(work2, mu[kk]);
                double tmp2 = gsl_stats_tss_m(work2->data, work2->stride, work2->size, 0);
                if (tmp2 < *ptr)
                    *ptr = tmp2;
            }
            sum += *ptr;
        }
        gsl_vector_scale(xdiff, 1.0 / sum);
        for (i = 0; i < data->size1 - 1; i++) {
            *gsl_vector_ptr(xdiff, i + 1) += gsl_vector_get(xdiff, i);
        }
        double rval = gsl_rng_uniform(rng);
        for (i = 0; i < data->size1; i++) {
            double tmp = gsl_vector_get(xdiff, i);
            if (rval <= gsl_vector_get(xdiff, i)) {
                break;
            }
        }
        gsl_matrix_get_row(mu[k], data, i);
        gsl_matrix_set_identity(sigma[k]);
    }

    // Initialize I, the indicator matrix
    indic = gsl_matrix_calloc(data->size1, NCOMPS);
    for (unsigned int i = 0; i < data->size1; i++) {
        gsl_matrix_set(indic, i, gsl_rng_uniform_int(rng, NCOMPS), 1.0);
    }

    // Main EM Loop
    unsigned int iter;
    for (iter = 0; iter < MAXITER; iter++) {
#ifdef DEBUG
        printf("\nIteration %d:", iter);
        for (unsigned int k = 0; k < NCOMPS; k++) {
            printf("\n\t[");
            for (unsigned int j = 0; j < data->size2; j++) {
                printf("%s%.4f", j == 0 ? "" : "\t", gsl_vector_get(mu[k], j));
            }
            printf("]");
        }
        fflush(stdout);
#endif
        for (unsigned int k = 0; k < NCOMPS; k++) {
            gsl_vector_memcpy(prev_mu[k], mu[k]);
            gsl_matrix_memcpy(prev_sigma[k], sigma[k]);
        }
        e_step(work1, work2, work3);
        m_step(work1);
        if (converged())
            break;
    }
    if (iter == MAXITER)
        printf("\nEM did not converge within %d iterations\n", MAXITER);
    else
        printf("\nEM converged in %d iterations\n", iter + 1);

    // Free work buffers
    for (unsigned int k = 0; k < NCOMPS; k++) {
        gsl_matrix_free(prev_sigma[k]);
        gsl_vector_free(prev_mu[k]);
    }
    free(prev_sigma);
    free(prev_mu);
    gsl_matrix_free(indic);
    gsl_vector_free(work3);
    gsl_vector_free(work2);
    gsl_vector_free(work1);
}

void mvn_em_wrapper(void) {
    data = gsl_matrix_calloc(NSAMPS, NDIMS);
    mu = calloc(NCOMPS, sizeof(gsl_vector *));
    sigma = calloc(NCOMPS, sizeof(gsl_matrix *));
    gsl_vector *tmp = gsl_vector_calloc(NDIMS);

    switch (load) {
        case FILEOP_WRITE: {
            FILE * file = fopen(fname, "w");
            if (!file)
                usage();
            gsl_error_handler_t *handler = gsl_set_error_handler_off();

            // Generate parameters
            for (unsigned int k = 0; k < NCOMPS; k++) {
                do {
                    mu[k] = gsl_vector_calloc(NDIMS);
                    sigma[k] = gsl_matrix_calloc(NDIMS, NDIMS);
                    for (unsigned int j = 0; j < NDIMS; j++) {
                        gsl_vector_set(mu[k], j, gsl_ran_flat(rng, -100.0, 100.0));
                        gsl_matrix_set(sigma[k], j, j, gsl_ran_flat(rng, 1.0, 25.0));
                        for (unsigned int jj = j + 1; jj < NDIMS; jj++) {
                            double val = gsl_ran_gaussian(rng, 1.0);
                            gsl_matrix_set(sigma[k], j, jj, val);
                            gsl_matrix_set(sigma[k], jj, j, val);
                        }
                    }
                } while (gsl_linalg_cholesky_decomp(sigma[k]));
                printf("Mean %d:", k);
                for (unsigned int j = 0; j < NDIMS; j++) {
                    printf("\t%.4f", gsl_vector_get(mu[k], j));
                };
                printf("\n");
            }
            gsl_set_error_handler(handler);

            // Generate data
            for (unsigned int i = 0; i < NSAMPS; i++) {
                unsigned long k = gsl_rng_uniform_int(rng, NCOMPS);
                gsl_ran_multivariate_gaussian(rng, mu[k], sigma[k], tmp);
                gsl_matrix_set_row(data, i, tmp);
            }

            // Write data
            gsl_matrix_fwrite(file, data);
            fclose(file);
            break;
        }
        case FILEOP_READ: {
            FILE * file = fopen(fname, "r");
            if (!file)
                usage();
            gsl_matrix_fread(file, data);
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
    gsl_vector_free(tmp);
    for (unsigned int k = 0; k < NCOMPS; k++) {
        printf("New mean %d:", k);
        for (unsigned int j = 0; j < NDIMS; j++) {
            printf("\t%.4f", gsl_vector_get(mu[k], j));
        }
        printf("\n");
        gsl_matrix_free(sigma[k]);
        gsl_vector_free(mu[k]);
    }
    free(sigma);
    free(mu);
    gsl_matrix_free(data);
}
