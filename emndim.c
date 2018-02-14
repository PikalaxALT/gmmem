#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <sys/time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>
#include "emndim.h"
#include "gsl_extd_util.h"
#include "main.h"

static gsl_matrix *data = NULL;
static gsl_matrix *indic = NULL;
static gsl_vector **mu = NULL;
static gsl_matrix **sigma = NULL;
static double *theta = NULL;
static gsl_vector **prev_mu = NULL;

static void e_step(gsl_vector *work1, gsl_vector *work2, gsl_vector *work3) {
    // work1: p
    // work2: p
    // work3: NCOMPS
    for (unsigned int i = 0; i < NSAMPS; i++) {
        gsl_matrix_get_row(work1, data, i);
        for (unsigned int k = 0; k < NCOMPS; k++) {
            double *ptr = gsl_vector_ptr(work3, k);
            gsl_ran_multivariate_gaussian_log_pdf(work1, mu[k], sigma[k], ptr, work2);
            *ptr += log(theta[k] + DBL_MIN);
        }
        gsl_vector_add_constant(work3, -gsl_vector_logsumexp(work3));
        gsl_vector_exp(work3);
        gsl_matrix_set_row(indic, i, work3);
    }
}

static void m_step(gsl_vector *work1) {
    for (unsigned int k = 0; k < NCOMPS; k++) {
        gsl_vector_set_zero(mu[k]);
        gsl_matrix_set_zero(sigma[k]);
        double denom = 0.0;

        // denom = 1.0 / (sum_i I_ik)
        // mu_k = denom * sum_i I_ik * x_i
        for (unsigned int i = 0; i < NSAMPS; i++) {
            double tmp = gsl_matrix_get(indic, i, k);
            denom += tmp;
            gsl_matrix_get_row(work1, data, i);
            gsl_vector_scale(work1, tmp);
            gsl_vector_add(mu[k], work1);
        }
        theta[k] = denom / NSAMPS;
        denom = 1.0 / denom;
        gsl_vector_scale(mu[k], denom);

        // sigmasq_k = denom * sum_i I_ik * (x_i - mu_k) * (x_i - mu_k)^T
        for (unsigned int i = 0; i < NSAMPS; i++) {
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
    gsl_vector *work1 = gsl_vector_calloc(NDIMS);
    gsl_vector *work2 = gsl_vector_calloc(NDIMS);
    gsl_vector *work3 = gsl_vector_calloc(NCOMPS);
    gsl_vector *xdiff = gsl_vector_calloc(NSAMPS);

    // Initialize buffers for the previous parameters state
    prev_mu = calloc(NCOMPS, sizeof(gsl_vector *));
    assert(prev_mu);

    // Initialize means
    prev_mu[0] = gsl_vector_calloc(NDIMS);
    gsl_matrix_get_row(mu[0], data, gsl_rng_uniform_int(rng, NSAMPS));
    gsl_matrix_set_identity(sigma[0]);
    theta[0] = 1.0 / NCOMPS;
    for (unsigned int k = 1; k < NCOMPS; k++) {
        prev_mu[k] = gsl_vector_calloc(NDIMS);

        unsigned int i;
        unsigned int kk;
        double sum = 0.0;

        for (i = 0; i < NSAMPS; i++) {
            double *ptr = gsl_vector_ptr(xdiff, i);
            *ptr = INFINITY;
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
        assert(sum && !isnan(sum));
        gsl_vector_scale(xdiff, 1.0 / sum);
        double rval = gsl_rng_uniform(rng);
        sum = 0.0;
        for (i = 0; i < NSAMPS - 1; i++) {
            sum += gsl_vector_get(xdiff, i);
            if (rval <= sum) {
                break;
            }
        }
        gsl_matrix_get_row(mu[k], data, i);
        gsl_matrix_set_identity(sigma[k]);
        theta[k] = 1.0 / NCOMPS;
    }

    // Initialize I, the indicator matrix
    indic = gsl_matrix_calloc(NSAMPS, NCOMPS);
    for (unsigned int i = 0; i < NSAMPS; i++) {
        gsl_matrix_set(indic, i, gsl_rng_uniform_int(rng, NCOMPS), 1.0);
    }

    // Main EM Loop
    unsigned int iter;
    for (iter = 0; iter < MAXITER; iter++) {
#ifdef DEBUG
        printf("\nIteration %d:", iter);
        for (unsigned int k = 0; k < NCOMPS; k++) {
            printf("\n\t[");
            for (unsigned int j = 0; j < NDIMS; j++) {
                printf("%s%.4f", j == 0 ? "" : "\t", gsl_vector_get(mu[k], j));
            }
            printf("]");
        }
        printf("\n");
        for (unsigned int k = 0; k < NCOMPS; k++) {
            printf("\t%.4f", theta[k]);
        }
        fflush(stdout);
#else
        printf("%c\b", "/-\\|"[iter % 4]);
        fflush(stdout);
#endif
        for (unsigned int k = 0; k < NCOMPS; k++) {
            gsl_vector_memcpy(prev_mu[k], mu[k]);
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
        gsl_vector_free(prev_mu[k]);
    }
    free(prev_mu);
    gsl_matrix_free(indic);
    gsl_vector_free(work3);
    gsl_vector_free(work2);
    gsl_vector_free(work1);
}

static void alloc_data(void) {
    data = gsl_matrix_calloc(NSAMPS, NDIMS);
    mu = calloc(NCOMPS, sizeof(gsl_vector *));
    sigma = calloc(NCOMPS, sizeof(gsl_matrix *));
    theta = calloc(NCOMPS, sizeof(double));
    assert(mu && sigma && theta);
    for (unsigned int k = 0; k < NCOMPS; k++) {
        mu[k] = gsl_vector_calloc(NDIMS);
        sigma[k] = gsl_matrix_calloc(NDIMS, NDIMS);
    }
}

static void generate_data(void) {
    gsl_vector *tmp = gsl_vector_calloc(NDIMS);
    gsl_error_handler_t *handler = gsl_set_error_handler_off();

    // Generate parameters
    double thetasum = 0.0;
    for (unsigned int k = 0; k < NCOMPS; k++) {
        do {
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
        thetasum += theta[k] = gsl_rng_uniform_pos(rng);
    }
    gsl_set_error_handler(handler);
    for (unsigned int k = 0; k < NCOMPS; k++) {
        theta[k] /= thetasum;
        printf("Mean %d:", k);
        for (unsigned int j = 0; j < NDIMS; j++) {
            printf("\t%.4f", gsl_vector_get(mu[k], j));
        };
        printf("\t(p = %.4f)\n", theta[k]);
    }

    // Generate data
    for (unsigned int i = 0; i < NSAMPS; i++) {
        double tmp2 = gsl_rng_uniform(rng);
        double cumsum = 0.0;
        unsigned int k;
        for (k = 0; k < NCOMPS; k++) {
            cumsum += theta[k];
            if (tmp2 <= cumsum) {
                break;
            }
        }
        gsl_ran_multivariate_gaussian(rng, mu[k], sigma[k], tmp);
        gsl_matrix_set_row(data, i, tmp);
    }
    gsl_vector_free(tmp);
}

static void data_write(FILE * file) {
    gsl_matrix_fwrite(file, data);
}

static void data_read(FILE * file) {
    gsl_matrix_fread(file, data);
}

static void print_fit(void) {
    for (unsigned int k = 0; k < NCOMPS; k++) {
        printf("New mean %d:", k);
        for (unsigned int j = 0; j < NDIMS; j++) {
            printf("\t%.4f", gsl_vector_get(mu[k], j));
        }
        printf("\t(p = %.4f)\n", theta[k]);
    }
}

static void free_data(void) {
    for (unsigned int k = 0; k < NCOMPS; k++) {
        gsl_matrix_free(sigma[k]);
        gsl_vector_free(mu[k]);
    }
    free(theta);
    free(sigma);
    free(mu);
    gsl_matrix_free(data);
}

static void params_output(FILE * file) {
    gsl_matrix *square = gsl_matrix_alloc(NDIMS, NDIMS);
    for (unsigned int k = 0; k < NCOMPS; k++) {
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, sigma[k], sigma[k], 0.0, square);
        fprintf(file, "COMPONENT %d:\n", k + 1);
        fprintf(file, "mu:\t[");
        for (unsigned int j = 0; j < NDIMS; j++) {
            fprintf(file, "%.4f%c", gsl_vector_get(mu[k], j), j == NDIMS - 1 ? ']' : '\t');
        }
        fputc('\n', file);
        fprintf(file, "sigma:\n\t[");
        for (unsigned int j = 0; j < NDIMS; j++) {
            fputc('[', file);
            for (unsigned int jj = 0; jj < NDIMS; jj++) {
                fprintf(file, "%.4f%c", gsl_matrix_get(square, j, jj), jj == NDIMS - 1 ? ']' : '\t');
            }
            if (j != NDIMS - 1) {
                fprintf(file, "\n\t ");
            }
        }
        fprintf(file, "]\n");
        fprintf(file, "theta:\t%.4f\n", theta[k]);
        fputc('\n', file);
    }
}

void mvn_init_funcs(void) {
    alloc_fn = alloc_data;
    gen_fn = generate_data;
    read_fn = data_read;
    write_fn = data_write;
    output_fn = params_output;
    run_fn = fit_em;
    print_fn = print_fit;
    free_fn = free_data;
}
