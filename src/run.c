#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <string.h>
#include <time.h>
#include "tools.h"
#include <stdbool.h>

#define ALPHA 0.15

typedef struct {
    double* walkers; 
    int n;
    double mean_E_t;
} result_dmc;

typedef struct {
    double* walkers;
    int n;
} result_dmc_one_step;

void task1(void);
void task2(void);
void task3(void);
void task3b(void);
void task4(void);
gsl_rng* get_rand(void);
double displace_x(double x, double delta_tau, gsl_rng* r);
double weight_1d(double x, double E_t, double dt);
double update_E_t(double E_t, double gamma, int n, int n0);
void diffusion_monte_carlo_1d(double* walkers, int n0, double E_t, double gamma, double dt, int n_iter, int n_eq);
result_dmc diffusion_monte_carlo_6d(double** walkers, int n0, double E_t, double gamma, double dt, int n_iter, int n_eq, int decomposition, int task, bool print);
double** init_walkers_6d(int n);
void diffusion_monte_carlo_task2(double** walkers, int n0, double E_t, double gamma, double dt, int n_iter, int n_eq);
double** polar_to_cart(double** polar, int n);
void first_order_drift_6d(double* walker, double alpha, double dt);
void second_order_drift_6d(double* walker, double alpha, double dt);
double radius(double* walker);
double hamiltonian_potential(double* walker);
double get_E_l(double* walker, double alpha);
double** reactive_part(double** walkers, int* n_ptr, double dt, double E_t, gsl_rng* r);
void diffusive_part(double** walkers, int n, double dt, gsl_rng* r);

int
run(
    int argc,
    char *argv[]
   )
{
    //task1();
    //task2();
    //task3();
    task3b();
    //task4();
    return 0;
}

void task1(void)
{
    int n = 200;
    double* walkers = linspace(-5, 5, n, false);
    double E_t = 0.5;
    double delta_tau = 0.02;
    double gamma = 0.5;
    int n_iter = 50000;
    int n_eq = 1500;
    diffusion_monte_carlo_1d(walkers, n, E_t, gamma, delta_tau, n_iter, n_eq);
}

void task2(void)
{
    int n = 1000;
    double** walkers = init_walkers_6d(n);
    double** walkers_cartesian = polar_to_cart(walkers, n);
    destroy_2D_array(walkers);
    double E_t = - 3;
    double delta_tau = 0.01;
    double gamma = 0.5;
    int n_iter = 1000000;
    int n_eq = 1500;
    diffusion_monte_carlo_task2(walkers_cartesian, n, E_t, gamma, delta_tau, n_iter, n_eq);
}

void task3(void)
{
    int n = 1000;
    double** walkers = init_walkers_6d(n);
    double** walkers_cartesian = polar_to_cart(walkers, n);
    destroy_2D_array(walkers);
    double E_t = - 3;
    double delta_tau = 0.1;
    double gamma = 0.5;
    int n_iter = 25000;
    int n_eq = 1500;
    diffusion_monte_carlo_6d(walkers_cartesian, n, E_t, gamma, delta_tau, n_iter, n_eq, 1, 3, true);
}

void task3b(void)
{
    int n = 1000;
    double** walkers = init_walkers_6d(n);
    double** walkers_cartesian = polar_to_cart(walkers, n);
    destroy_2D_array(walkers);
    double E_t = - 3;
    double delta_tau = 0.1;
    double gamma = 0.5;
    int n_iter = 20000;
    int n_eq = 1500;
    diffusion_monte_carlo_6d(walkers_cartesian, n, E_t, gamma, delta_tau, n_iter, n_eq, 2, 4, true);   
}

void task4(void)
{
    int n = 1000;
    double gamma = 0.5;
    int t = 2000;
    int n_eq = 1500;
    double E_t = - 3;
    int n_runs = 20;
    double* dts = linspace(0.01, 0.4, n_runs, true);
    FILE* file = fopen("data/task5.csv", "w+");

    for (int i = 0; i < n_runs; i++)
    {   
        printf("dt: %.3lf\n", dts[i]);
        int n_iter = t / dts[i];
        double** walkers = init_walkers_6d(n);
        result_dmc result1 = diffusion_monte_carlo_6d(polar_to_cart(walkers, n), n, E_t, gamma, dts[i], n_iter, n_eq, 1, 5, false);  
        result_dmc result2 = diffusion_monte_carlo_6d(polar_to_cart(walkers, n), n, E_t, gamma, dts[i], n_iter, n_eq, 2, 5, false);  
        destroy_2D_array(walkers);
        fprintf(file, "%lf, %lf\n", result1.mean_E_t, result2.mean_E_t);
    }
    fclose(file);
}

void diffusion_monte_carlo_1d(double* walkers, int n0, double E_t, double gamma, double dt, int n_iter, int n_eq)
{
    gsl_rng* r = get_rand();
    int n = n0;
    double E_t_sum = 0;
    int denominator = 0;
    FILE* file = fopen("data/task1.csv", "w+");
    FILE* positions = fopen("data/positions_task1.csv", "w+");
    int* walker_multiplier;
    double* new_walkers;

    for (int i = 0; i < n_iter; i++)
    { // For every DMC iteration
        for (int j = 0; j < n; j++)
        {
            fprintf(positions, "%lf\n", walkers[j]);
        }
        fprintf(file, "%d, %lf\n", n, E_t);

        walker_multiplier = (int*) malloc(sizeof(int) * n);
        int multiplier_sum = 0; 

        for (int j = 0; j < n; j++)
        {
            walkers[j] = displace_x(walkers[j], dt, r);
            double w = weight_1d(walkers[j], E_t, dt);
            int m = (int) (w + gsl_rng_uniform(r));
            walker_multiplier[j] = m;
            multiplier_sum += m;
        }

        new_walkers = (double*) malloc(sizeof(double) * multiplier_sum);
        int l = 0;    
        
        for (int j = 0; j < n; j++)
        { // For every previous walker
            for (int k = 0; k < walker_multiplier[j]; k++)
            { // For every walker child    
                new_walkers[l] = walkers[j];
                l++;
            }
        }
        free(walker_multiplier);
        free(walkers);
        walkers = (double*) malloc(sizeof(double) * multiplier_sum);
        
        for(int j = 0; j < multiplier_sum; ++j)
        {
            walkers[j] = new_walkers[j];
        }

        free(new_walkers);
        n = multiplier_sum;
        if (i == n_eq)
        {
            E_t_sum = 0;
            denominator = 0;
        }
        E_t_sum += E_t;
        denominator++;
        E_t = update_E_t(E_t_sum / denominator, gamma, n, n0);
    }
    free(walkers);
    fclose(file);
    fclose(positions);
}

void diffusion_monte_carlo_task2(double** walkers, int n0, double E_t, double gamma, double dt, int n_iter, int n_eq)
{
    gsl_rng* r = get_rand();
    int n = n0;
    double E_t_sum = 0;
    FILE* file = fopen("data/task2.csv", "w+");
    int* walker_multiplier;
    double** new_walkers;
    int dim = 6;
    int denominator = 0;

    for (int i = 0; i < n_iter; i++)
    {

        fprintf(file, "%d, %lf\n", n, E_t);

        walker_multiplier = (int*) malloc(sizeof(int) * n);
        int multiplier_sum = 0; 

        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < dim; k++)
            {
                walkers[j][k] = displace_x(walkers[j][k], dt, r);
            }
            double v = hamiltonian_potential(walkers[j]);
            double w = exp(- (v - E_t) * dt);
            int m = (int) (w  + gsl_rng_uniform(r));
            m = m > 3 ? 3 : m;
            walker_multiplier[j] = m;
            multiplier_sum += m;
        }

        new_walkers = create_2D_array(multiplier_sum, dim);
        int l = 0;    

        for (int j = 0; j < n; j++)
        { //For every previous walker
            for (int k = 0; k < walker_multiplier[j]; k++)
            { //For every walker child    
                for (int m = 0; m < dim; m++)
                {
                    new_walkers[l][m] = walkers[j][m];
                }
                l++;
            }
        }
        free(walker_multiplier);
        destroy_2D_array(walkers);
        walkers = create_2D_array(multiplier_sum, dim);
        
        for(int j = 0; j < multiplier_sum; j++)
        {
            for (int k = 0; k < dim; k++)
            {
                walkers[j][k] = new_walkers[j][k];
            }
        }

        destroy_2D_array(new_walkers);
        n = multiplier_sum;
        if (i == n_eq)
        {
            E_t_sum = 0;
            denominator = 0;
        }
        E_t_sum += E_t;
        denominator++;
        E_t = update_E_t(E_t_sum / denominator, gamma, n, n0);
    }
    destroy_2D_array(walkers);
    fclose(file);
}

result_dmc diffusion_monte_carlo_6d(double** walkers, int n0, double E_t, double gamma, double dt, int n_iter, int n_eq, int decomposition, int task, bool print)
{
    gsl_rng* r = get_rand();
    int n = n0;
    double E_t_sum = 0;
    char buffer1[50];
    char buffer2[50];
    sprintf(buffer1, "data/task%d.csv", task);
    sprintf(buffer2, "data/positions_task%d.csv", task);
    FILE* file; 
    FILE* positions;
    if (print)
    {
        file = fopen(buffer1, "w+");
        positions = fopen(buffer2, "w+");
    }
    int denominator = 0;

    for (int i = 0; i < n_iter; i++)
    { // For every DMC iteration
        if (print)
        {
            for (int j = 0; j < n; j++)
            {
                fprintf(positions, "%lf, %lf\n", radius(walkers[j]), radius(walkers[j] + 3));
            }
            fprintf(file, "%d, %lf\n", n, E_t);
        }

        if (decomposition == 1) // Basic decomposition.
        { // Reactive part -> Diffusive part -> Drift
            walkers = reactive_part(walkers, &n, dt, E_t, r);
            diffusive_part(walkers, n, dt, r);
            for (int j = 0; j < n; j++)
            {
                first_order_drift_6d(walkers[j], ALPHA, dt);
            }
        }
        else if (decomposition == 2) // Advanced decomposition.
        {   // Half drift -> Half diffusive -> Reactive -> Half diffusive -> Half drift
            for (int j = 0; j < n; j++)
            {
                second_order_drift_6d(walkers[j], ALPHA, dt / 2);
            }
            diffusive_part(walkers, n, dt / 2, r);
            walkers = reactive_part(walkers, &n, dt, E_t, r);
            diffusive_part(walkers, n, dt / 2, r);
            for (int j = 0; j < n; j++)
            {
                second_order_drift_6d(walkers[j], ALPHA, dt / 2);
            }
        }
        
        // Updating E_t
        if (i == n_eq)
        {
            E_t_sum = 0;
            denominator = 0;
        }
        E_t_sum += E_t;
        denominator++;
        E_t = update_E_t(E_t_sum / denominator, gamma, n, n0);
    }
    destroy_2D_array(walkers);
    if (print)
    {
        fclose(file);
        fclose(positions);
    }
    result_dmc result;
    result.mean_E_t = E_t_sum / denominator;
    return result;
}

double** polar_to_cart(double** polar, int n)
{
    int dim = 6;
    double** rect = create_2D_array(n, dim);

    for(int i = 0; i < n; i++)
    {
        rect[i][0] = polar[i][0] * sin(polar[i][1]) * cos(polar[i][2]);
        rect[i][1] = polar[i][0] * sin(polar[i][1]) * sin(polar[i][2]);
        rect[i][2] = polar[i][0] * cos(polar[i][1]);
        rect[i][3] = polar[i][3] * sin(polar[i][4]) * cos(polar[i][5]);
        rect[i][4] = polar[i][3] * sin(polar[i][4]) * sin(polar[i][5]);
        rect[i][5] = polar[i][3] * cos(polar[i][4]);
    }
    return rect;
}

double displace_x(double x, double delta_tau, gsl_rng* r)
{
    double G = gsl_ran_gaussian(r, 1);
    return x + sqrt(delta_tau) * G;
}

double potential_1d(double x)
{
    double y = 1 - exp(- x);
    return 0.5 * y * y;
}

double weight_1d(double x, double E_t, double dt)
{
    return exp(- (potential_1d(x) - E_t) * dt);
} 

double hamiltonian_potential(double* walker)
{
    double r1[] = {walker[0], walker[1], walker[2]};
    double r2[] = {walker[3], walker[4], walker[5]};
    double r12_len = distance_between_vectors(r1, r2, 3);
    double r1_len = vector_norm(r1, 3);
    double r2_len = vector_norm(r2, 3);
    return - 2 / (r1_len + 1e-4) - 2 / (r2_len + 1e-4)+ 1 / (r12_len + 1e-4);
}

void v_F(double* v, double* walker, double alpha)
{
    // r_1, r_2
    double r1[] = {walker[0], walker[1], walker[2]};
    double r2[] = {walker[3], walker[4], walker[5]};

    // r_1/|r_1|, r_2/|r_2|
    double r1_norm[] = {walker[0], walker[1], walker[2]}; 
    double r2_norm[] = {walker[3], walker[4], walker[5]};
    normalize_vector(r1_norm, 3);
    normalize_vector(r2_norm, 3);

    // r_12/|r_12|
    double r12_norm[3];
    elementwise_subtraction(r12_norm, r1, r2, 3);
    normalize_vector(r12_norm, 3);

    // |r_12|
    double r12_len = distance_between_vectors(r1, r2, 3);

    // 2(1 + alpha |r_12|)^2
    double shared_denominator = 2. * (1 + alpha * r12_len) * (1 + alpha * r12_len);

    for(int i = 0; i < 3; i++)
    {
        v[i]   = -2 * r1_norm[i] - r12_norm[i] / shared_denominator;
        v[i+3] = -2 * r2_norm[i] + r12_norm[i] / shared_denominator;
    }
}

void first_order_drift_6d(double* walker, double alpha, double dt)
{
    
    int dim = 6;
    double v[] = {0., 0., 0., 0., 0., 0.};
    v_F(v, walker, alpha);
    
    for(int i = 0; i < dim; i++)
    {
        walker[i] += v[i]*dt;
    }
}

void second_order_drift_6d(double* walker, double alpha, double dt)
{
    int dim = 6;

    double v[] = {0., 0., 0., 0., 0., 0.};
    
    double tmp_walker[dim];
    memcpy(tmp_walker, walker, sizeof(tmp_walker)); 
    
    // R_{1/2}
    v_F(v, walker, alpha);
    for(int i = 0; i < dim; i++)
    {
        tmp_walker[i] += v[i]*(dt/2);
    }

    // R(dt)
    v_F(v, tmp_walker, alpha);
    for(int i = 0; i < dim; i++)
    {
        walker[i] += v[i]*dt;
    }    
}

double get_E_l(double* walker, double alpha)
{
    // r_1, r_2
    double r1[] = {walker[0], walker[1], walker[2]};
    double r2[] = {walker[3], walker[4], walker[5]};
    
    // |r_12|
    double r12_len = distance_between_vectors(r1, r2, 3);

    // (1 + alpha*r_12)
    double denominator = 1 + alpha * r12_len;

    // r^_1, r^_2 = r_1/|r_1|, r_2/|r_2|
    double r1_norm[] = {walker[0], walker[1], walker[2]};
    double r2_norm[] = {walker[3], walker[4], walker[5]};
    normalize_vector(r1_norm, 3);
    normalize_vector(r2_norm, 3);

    // (r^_1 - r^_2)
    double r_norm_diff[] = {0., 0., 0.};
    elementwise_subtraction(r_norm_diff, r1_norm, r2_norm, 3);

    // (r_1 - r_2)
    double r_diff[] = {0., 0., 0.};
    elementwise_subtraction(r_diff, r1, r2, 3);
    
    // (r^_1 - r^_2)•(r_1 - r_2)
    double dot = (dot_product(r_norm_diff, r_diff, 3));

    return - 4 + dot / (r12_len * pow(denominator, 2)) -
        1 / (r12_len * pow(denominator, 3)) - 1 / (4 * pow(denominator, 4)) + 1 / r12_len;
}

double update_E_t(double E_t, double gamma, int n, int n0)
{
    return E_t - gamma * log((double) n / n0);
}

double** init_walkers_6d(int n)
{
    gsl_rng* r = get_rand();
    int dim = 6;
    double** walkers = create_2D_array(n, dim);
    for (int i = 0; i < n; i++)
    {
        walkers[i][0] = 0.7 + gsl_rng_uniform(r);
        walkers[i][1] = acos(2 * gsl_rng_uniform(r) - 1);
        walkers[i][2] = 2 * M_PI * gsl_rng_uniform(r);
        walkers[i][3] = 0.7 + gsl_rng_uniform(r);
        walkers[i][4] = acos(2 * gsl_rng_uniform(r) - 1);
        walkers[i][5] = 2 * M_PI * gsl_rng_uniform(r);
    }
    return walkers;
}

double radius(double* walker)
{
    return sqrt(walker[0] * walker[0] + walker[1] * walker[1] +
        walker[2] * walker[2]);
}

gsl_rng* get_rand(void){
    const gsl_rng_type* T;
    gsl_rng* r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    time_t seed = time(NULL);
    gsl_rng_set(r, seed);
    return r;
}

double** reactive_part(double** walkers, int* n_ptr, double dt, double E_t, gsl_rng* r)
{
    int n = *n_ptr;
    int multiplier_sum = 0;
    int* walker_multiplier = (int*) malloc(sizeof(int) * n);
    int dim = 6;

    for (int j = 0; j < n; j++)
    {
        double E_l = get_E_l(walkers[j], ALPHA);
        double w = exp(- (E_l - E_t) * dt);
        int m = (int) (w  + gsl_rng_uniform(r));
        walker_multiplier[j] = m;
        multiplier_sum += m;
    }

    double** new_walkers = create_2D_array(multiplier_sum, dim);
    int l = 0;    

    for (int j = 0; j < n; j++)
    { //For every previous walker
        for (int k = 0; k < walker_multiplier[j]; k++)
        { //For every walker child    
            for (int m = 0; m < dim; m++)
            {
                new_walkers[l][m] = walkers[j][m];
            }
            l++;
        }
    }

    destroy_2D_array(walkers);
    walkers = create_2D_array(multiplier_sum, dim);
    
    for(int j = 0; j < multiplier_sum; j++)
    {
        for (int k = 0; k < dim; k++)
        {
            walkers[j][k] = new_walkers[j][k];
        }
    }

    free(walker_multiplier);
    destroy_2D_array(new_walkers);
    *n_ptr = multiplier_sum;
    //printf("%d\n", multiplier_sum);
    return walkers;
}

void diffusive_part(double** walkers, int n, double dt, gsl_rng* r)
{
    for (int j = 0; j < n; j++)
    {
        for (int k = 0; k < 6; k++)
        {
            walkers[j][k] = displace_x(walkers[j][k], dt, r);
        }
    }
}
