#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <sys/time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include "main.h"
#include "em1dim.h"
#include "gsl_extd_util.h"

static gsl_vector *data = NULL;
static gsl_matrix *indic = NULL;
static double *mu = NULL;
static double *sigma = NULL;
static double *theta = NULL;
static double *prev_mu = NULL;

static void e_step(gsl_vector *work) {
    // work: NCOMPS
    for (unsigned int i = 0; i < NSAMPS; i++) {
        double xval = gsl_vector_get(data, i);
        for (unsigned int k = 0; k < NCOMPS; k++) {
            gsl_vector_set(work, k, log(gsl_ran_gaussian_pdf(xval - mu[k], sigma[k]) + DBL_MIN) + log(theta[k] + DBL_MIN));
        }
        gsl_vector_add_constant(work, -gsl_vector_logsumexp(work));
        gsl_vector_exp(work);
        gsl_matrix_set_row(indic, i, work);
    }
}

static void m_step(void) {
    for (unsigned int k = 0; k < NCOMPS; k++) {
        mu[k] = 0.0;
        sigma[k] = 0.0;
        double denom = 0.0;

        // denom = 1.0 / (sum_i I_ik)
        // mu_k = denom * sum_i I_ik * x_i
        for (unsigned int i = 0; i < NSAMPS; i++) {
            double tmp = gsl_matrix_get(indic, i, k);
            denom += tmp;
            mu[k] += tmp * gsl_vector_get(data, i);
        }
        theta[k] = denom / NSAMPS;
        mu[k] /= denom;

        // sigmasq_k = denom * sum_i I_ik * (x_i - mu_k)^2
        for (unsigned int i = 0; i < NSAMPS; i++) {
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
    gsl_vector *xdiff = gsl_vector_alloc(NSAMPS);

    // Initialize buffers for the previous parameters state
    prev_mu = calloc(NCOMPS, sizeof(double));

    assert(prev_mu); // catch calloc failure

    // Initialize means
    mu[0] = gsl_vector_get(data, gsl_rng_uniform_int(rng, NSAMPS));
    sigma[0] = 1.0;
    theta[0] = 1.0 / NCOMPS;
    for (unsigned int k = 1; k < NCOMPS; k++) {
        unsigned int i;
        unsigned int kk;
        double sum = 0.0;

        for (i = 0; i < NSAMPS; i++) {
            double *ptr = gsl_vector_ptr(xdiff, i);
            *ptr = INFINITY;
            double tmp = gsl_vector_get(data, i);
            for (kk = 0; kk < k; kk++) {
                double tmp2 = (tmp - mu[kk]) * (tmp - mu[kk]);
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
        mu[k] = gsl_vector_get(data, i);
        sigma[k] = 1.0;
        theta[k] = 1.0 / NCOMPS;
    }

    // Initialize I, the indicator matrix
    indic = gsl_matrix_calloc(NSAMPS, NCOMPS);

    // Main EM Loop
    unsigned int iter;
    for (iter = 0; iter < MAXITER; iter++) {
#ifdef DEBUG
        printf("\nIteration %d:", iter);
        for (unsigned int k = 0; k < NCOMPS; k++) {
            printf("\t%.4f", mu[k]);
        }
        printf("\t|");
        for (unsigned int k = 0; k < NCOMPS; k++) {
            printf("\t%.4f", theta[k]);
        }
        fflush(stdout);
#else
        printf("%c\b", "/-\\|"[iter % 4]);
        fflush(stdout);
#endif
        memcpy(prev_mu, mu, NCOMPS * sizeof(double));
        e_step(work);
        m_step();
        if (converged()) {
            break;
        }
    }
    if (iter == MAXITER)
        printf("\nEM did not converge within %d iterations\n", MAXITER);
    else
        printf("\nEM converged in %d iterations\n", iter + 1);

    // Free work buffers
    gsl_vector_free(xdiff);
    free(prev_mu);
    gsl_matrix_free(indic);
    gsl_vector_free(work);
}

static unsigned int checksum(void) {
    double *ptr;
    unsigned int sum = 0;
    for (unsigned int i = 0; i < NSAMPS; i++) {
        ptr = gsl_vector_ptr(data, i);
        sum += *(unsigned int *)ptr ^ XORHASH;
    }
    return sum;
}

static void alloc_data(void) {
    data = gsl_vector_calloc(NSAMPS);
    mu = calloc(NCOMPS, sizeof(double));
    sigma = calloc(NCOMPS, sizeof(double));
    theta = calloc(NCOMPS, sizeof(double));
    assert(mu && sigma && theta);
}

static void generate_data(void) {
    // Generate parameters
    double thetasum = 0.0;
    for (unsigned int k = 0; k < NCOMPS; k++) {
        mu[k] = gsl_ran_flat(rng, -100.0, 100.0);
        sigma[k] = sqrt(gsl_ran_flat(rng, 1.0, 25.0));
        thetasum += theta[k] = gsl_rng_uniform_pos(rng);
    }
    for (unsigned int k = 0; k < NCOMPS; k++) {
        theta[k] /= thetasum;
        printf("Mean %d:\t%.4f\t(p = %.4f)\n", k, mu[k], theta[k]);
    }

    // Generate data
    for (unsigned int i = 0; i < NSAMPS; i++) {
        double tmp = gsl_rng_uniform(rng);
        double cumsum = 0.0;
        unsigned int k;
        for (k = 0; k < NCOMPS; k++) {
            cumsum += theta[k];
            if (tmp <= cumsum) {
                break;
            }
        }
        gsl_vector_set(data, i, gsl_ran_gaussian(rng, sigma[k]) + mu[k]);
    }
}

static void data_write(FILE * file) {
    gsl_vector_fwrite(file, data);
}

static void data_read(FILE * file) {
    gsl_vector_fread(file, data);
}

static void print_fit(void) {
    for (unsigned int k = 0; k < NCOMPS; k++) {
        printf("New mean %d:\t%.4f\t(p = %.4f)\n", k, mu[k], theta[k]);
    }
}

static void free_data(void) {
    free(theta);
    free(sigma);
    free(mu);
    gsl_vector_free(data);
}

static void params_output(FILE * file) {
    for (unsigned int k = 0; k < NCOMPS; k++) {
        fprintf(file, "COMPONENT %d:\n", k + 1);
        fprintf(file, "mu:\t%.4f\n", mu[k]);
        fprintf(file, "sigma:\t%.4f\n", sigma[k] * sigma[k]);
        fprintf(file, "theta:\t%.4f\n", theta[k]);
        fputc('\n', file);
    }
}

void norm_init_funcs(void) {
    alloc_fn = alloc_data;
    gen_fn = generate_data;
    read_fn = data_read;
    write_fn = data_write;
    output_fn = params_output;
    run_fn = fit_em;
    print_fn = print_fit;
    free_fn = free_data;
    checksum_fn = checksum;
}
