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
gsl_rng* get_rand(void);
double displace_x(double x, double delta_tau, gsl_rng* r);
double weight_1d(double x, double E_t, double dt);
double update_E_t(double E_t, double gamma, int n, int n0);
result_dmc diffusion_monte_carlo_1d(double* walkers, int n0, double E_t, double gamma, double dt, int n_iter, int n_eq);
result_dmc diffusion_monte_carlo_6d(double** walkers, int n0, double E_t, double gamma, double dt, int n_iter, int n_eq);
double** init_walkers_6d(int n);

int
run(
    int argc,
    char *argv[]
   )
{
    //task1();
    task2();
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
    int n_eq = 10;
    result_dmc result = diffusion_monte_carlo_1d(walkers, n, E_t, gamma, delta_tau, n_iter, n_eq);
}

void task2(void)
{
    int n = 1000;
    double** walkers = init_walkers_6d(n);
    double E_t = 0.5;
    double delta_tau = 0.01;
    double gamma = 0.5;
    int n_iter = 10000;
    int n_eq = 10;
    result_dmc result = diffusion_monte_carlo_6d(walkers, n, E_t, gamma, delta_tau, n_iter, n_eq);
}

result_dmc diffusion_monte_carlo_1d(double* walkers, int n0, double E_t, double gamma, double dt, int n_iter, int n_eq)
{
    gsl_rng* r = get_rand();
    int n = n0;
    double E_t_sum = 0;
    FILE* file = fopen("data/task1.csv", "w+");
    FILE* positions = fopen("data/positions_task1.csv", "w+");
    int* walker_multiplier;
    double* new_walkers;

    for (int i = 0; i < n_iter; i++)
    { //For every DMC iteration
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
        { //For every previous walker
            for (int k = 0; k < walker_multiplier[j]; k++)
            { //For every walker child    
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
        if (i >= n_eq)
        {
            E_t_sum += E_t;
            E_t = update_E_t(E_t_sum / (i - n_eq + 1), gamma, n, n0);
        }
    }
    free(walkers);
    fclose(file);
    fclose(positions);
    result_dmc result;
    result.n = n;
    return result;
}

result_dmc diffusion_monte_carlo_6d(double** walkers, int n0, double E_t, double gamma, double dt, int n_iter, int n_eq)
{
    gsl_rng* r = get_rand();
    int n = n0;
    double E_t_sum = 0;
    FILE* file = fopen("data/task2.csv", "w+");
    FILE* positions = fopen("data/positions_task2.csv", "w+");
    int* walker_multiplier;
    double** new_walkers;


    free(walkers);
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

double update_E_t(double E_t, double gamma, int n, int n0)
{
    return E_t - gamma * log((double) n / n0);
}

double** init_walkers_6d(int n)
{
    gsl_rng* r = get_rand();
    int dim = 6;
    double** walkers = create_2D_array(dim, n);
    for (int i = 0; i < n; i++)
    {
        walkers[0][i] = 0.7 + gsl_rng_uniform(r);
        walkers[1][i] = acos(2 * gsl_rng_uniform(r) - 1);
        walkers[2][i] = 2 * M_PI * gsl_rng_uniform(r);
        walkers[3][i] = 0.7 + gsl_rng_uniform(r);
        walkers[4][i] = acos(2 * gsl_rng_uniform(r) - 1);
        walkers[5][i] = 2 * M_PI * gsl_rng_uniform(r);
    }
    return walkers;
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

