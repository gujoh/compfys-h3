#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <string.h>
#include <time.h>
#include "tools.h"
#include <stdbool.h>

typedef struct {
    double* walkers; 
    int n;
} result_dmc;

typedef struct {
    double* walkers;
    int n;
} result_dmc_one_step;

void task1(void);
void task2(void);
void task3(void);
gsl_rng* get_rand(void);
double displace_x(double x, double delta_tau, gsl_rng* r);
double weight_1d(double x, double E_t, double dt);
double weight_6d(double* walker, double E_t, double dt);
double update_E_t(double E_t, double gamma, int n, int n0);
void diffusion_monte_carlo_1d(double* walkers, int n0, double E_t, double gamma, double dt, int n_iter, int n_eq);
void diffusion_monte_carlo_6d(double** walkers, int n0, double E_t, double gamma, double dt, int n_iter, int n_eq, int add_drift, int task);
double** init_walkers_6d(int n);
double** polar_to_cart(double** polar, int n);
void drift1(double* walker, double alpha, double dt);
double radius(double* walker);

int
run(
    int argc,
    char *argv[]
   )
{
    //task1();
    //task2();
    task3();
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
    double E_t = - 3;
    double delta_tau = 0.01;
    double gamma = 0.5;
    int n_iter = 100000;
    int n_eq = 1500;
    diffusion_monte_carlo_6d(walkers_cartesian, n, E_t, gamma, delta_tau, n_iter, n_eq, 0, 2);
}

void task3(void)
{
    int n = 1000;
    double** walkers = init_walkers_6d(n);
    double** walkers_cartesian = polar_to_cart(walkers, n);
    double E_t = - 3;
    double delta_tau = 0.01;
    double gamma = 0.5;
    int n_iter = 20000;
    int n_eq = 1500;
    diffusion_monte_carlo_6d(walkers_cartesian, n, E_t, gamma, delta_tau, n_iter, n_eq, 1, 3);
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

void diffusion_monte_carlo_6d(double** walkers, int n0, double E_t, double gamma, double dt, int n_iter, int n_eq, int add_drift, int task)
{
    gsl_rng* r = get_rand();
    int n = n0;
    double E_t_sum = 0;
    char buffer1[50];
    char buffer2[50];
    sprintf(buffer1, "data/task%d.csv", task);
    sprintf(buffer2, "data/positions_task%d.csv", task);
    FILE* file = fopen(buffer1, "w+");
    FILE* positions = fopen(buffer2, "w+");
    int* walker_multiplier;
    double** new_walkers;
    int dim = 6;
    int denominator = 0;

    for (int i = 0; i < n_iter; i++)
    { // For every DMC iteration
        // for (int j = 0; j < n; j++)
        // {
        //     fprintf(positions, "%lf, %lf\n", radius(walkers[j]), radius(walkers[j] + 3));
        // }
        fprintf(file, "%d, %lf\n", n, E_t);

        walker_multiplier = (int*) malloc(sizeof(int) * n);
        int multiplier_sum = 0; 

        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < dim; k++)
            {
                walkers[j][k] = displace_x(walkers[j][k], dt, r);
            }
            if (add_drift == 1)
            {
                drift1(walkers[j], 0.15, dt);
            }
            double w = weight_6d(walkers[j], E_t, dt);
            int m = (int) (w  + gsl_rng_uniform(r));
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
    fclose(positions);
    // result_dmc result;
    // result.n = n;
    // return result;
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

double weight_6d(double* walker, double E_t, double dt)
{
    double r1[] = {walker[0], walker[1], walker[2]};
    double r2[] = {walker[3], walker[4], walker[5]};
    double r12_len = distance_between_vectors(r1, r2, 3);
    double r1_len = vector_norm(r1, 3);
    double r2_len = vector_norm(r2, 3);
    double v = - 2 / r1_len - 2 / r2_len + 1 / r12_len;
    return exp(- (v - E_t) * dt);
} 

void drift1(double* walker, double alpha, double dt)
{
    double r1[] = {walker[0], walker[1], walker[2]};
    double r2[] = {walker[3], walker[4], walker[5]};
    double r1_norm[3]; 
    double r2_norm[3];
    memcpy(r1_norm, r1, sizeof(r1_norm)); 
    memcpy(r2_norm, r2, sizeof(r2_norm)); 
    normalize_vector(r1_norm, 3);
    normalize_vector(r2_norm, 3);
    double r12_norm[3];
    elementwise_subtraction(r12_norm, r1, r2, 3);
    normalize_vector(r12_norm, 3);
    double r12_len = distance_between_vectors(r1, r2, 3);
    for (int i = 0; i < 3; i++)
    {
        walker[i] += (- 2. * r1_norm[i] - r12_norm[i] / (2. * pow(1 + alpha * r12_len, 2))) * dt;
        walker[i + 3] += (- 2. * r2_norm[i] + r12_norm[i] / (2. * pow(1 + alpha * r12_len, 2))) * dt;
    }
}

// double E_l(double* walker, double E_t, double alpha, double dt)
// {
//        double r1[] = {walker[0], walker[1], walker[2]};
//     double r2[] = {walker[3], walker[4], walker[5]};
//     double r12_len = distance_between_vectors(r1, r2, 3);
//     double r_norm_diff[] = {0, 0, 0};
//     double r_diff[] = {0, 0, 0};
//     double r1_norm[3]; 
//     double r2_norm[3];
//     memcpy(r1_norm, r1, sizeof(r1_norm)); 
//     memcpy(r2_norm, r2, sizeof(r2_norm)); 
//     double denominator = 1 + alpha * r12_len;
//     normalize_vector(r1_norm, 3);
//     normalize_vector(r2_norm, 3);
//     elementwise_subtraction(r_norm_diff, r1_norm, r2_norm, 3);
//     elementwise_subtraction(r_diff, r1, r2, 3);
//     double E_l = - 4 + (dot_product(r_norm_diff, r_diff, 3)) / (r12_len * pow(denominator, 2)) -
//         1 / (r12_len * pow(denominator, 3)) - 1 / (4 * pow(denominator, 4)) + 1 / r12_len;
// }

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

