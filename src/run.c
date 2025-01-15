// #include <stdio.h>
// #include <math.h>
// #include <stdlib.h>
// #include <gsl/gsl_rng.h>
// #include <gsl/gsl_randist.h>
// #include <string.h>
// #include <time.h>
// #include "tools.h"
// #include <stdbool.h>

// typedef struct {
//     double* walkers; 
//     int n;
// } result_dmc;

// typedef struct {
//     double* walkers;
//     int n;
// } result_dmc_one_step;

// void task1(void);
// gsl_rng* get_rand(void);
// result_dmc diffusion_monte_carlo(double* walkers, int n0, double E_t, double gamma, double dt, int t, int t_eq);
// result_dmc_one_step diffusion_monte_carlo_one_step(double* walkers, int n, double E_t, double dt, gsl_rng* r);
// double displace_x(double x, double delta_tau, gsl_rng* r);
// double weight(double x, double E_t, double dt);
// double update_E_t(double E_t, double gamma, int n, int n0);

// int
// run(
//     int argc,
//     char *argv[]
//    )
// {
//     task1();
//     return 0;
// }

// void task1(void){

//     int n = 200;
//     double* walkers = linspace(-5, 5, n, false);
//     double E_t = 0.5;
//     double delta_tau = 0.02;
//     double gamma = 0.5;
//     int t = 10000;
//     int t_eq = 0;
//     result_dmc result = diffusion_monte_carlo(walkers, n, E_t, gamma, delta_tau, t, t_eq);

//     free(result.walkers);
// }

// result_dmc diffusion_monte_carlo(double* walkers, int n0, double E_t, double gamma, double dt, int t, int t_eq)
// {
//     gsl_rng* r = get_rand();
//     int n = n0;
//     FILE* file = fopen("data/task1.csv", "w+");
//     FILE* positions = fopen("data/positions_task1.csv", "w+");
//     double E_t_sum = 0;

//     for (int i = 0; i < t; i++)
//     {
//         result_dmc_one_step result_one_step = diffusion_monte_carlo_one_step(walkers, n, E_t, dt, r);
//         n = result_one_step.n;
//         walkers = (double*) realloc(walkers, sizeof(double) * n);
//         memcpy(walkers, result_one_step.walkers, sizeof(double) * n);
//         free(result_one_step.walkers);
//         if (i >= t_eq)
//         { // Updating E_t
//             E_t_sum += E_t;
//             E_t = update_E_t(E_t_sum / (i - t_eq + 1), gamma, n, n0);
//         }

//         //Writing to file 
//         for (int j = 0; j < n; j++)
//         {
//             fprintf(positions, "%.10lf\n", walkers[j]);
//         }
//         //fprintf(positions, "\n");
//         fprintf(file, "%d, %lf\n", n, E_t);
//     }
//     fclose(file);
//     fclose(positions);
//     result_dmc result;
//     result.walkers = walkers;
//     result.n = n;
//     return result;
// }

// result_dmc_one_step diffusion_monte_carlo_one_step(double* walkers, int n, double E_t, double dt, gsl_rng* r)
// {
//     int* walker_multiplier = (int*) malloc(sizeof(int) * n);
//     int n_multiplier = 0; // Total number of new walkers.
//     for (int i = 0; i < n; i++)
//     { // Displaces walkers, and counts how many we should create. 
//         walkers[i] = displace_x(walkers[i], dt, r);
//         double w = weight(walkers[i], E_t, dt);
//         int m = (int) (w + gsl_rng_uniform(r));
//         walker_multiplier[i] = m;
//         n_multiplier += m;
//     }
//     double* new_walkers = (double*) malloc(sizeof(double) * n_multiplier);
//     int k = 0;
//     for (int i = 0; i < n; i++)
//     { // Creates new walkers.
//         bool inner_executed = false;
//         for (int j = 0; j < walker_multiplier[i]; j++)
//         { // The k variable is used so we do not skip places in the array.
//             new_walkers[k + j] = walkers[i];
//             inner_executed = true;
//         }
//         if (inner_executed)
//         {
//             k++;
//         }
//     }
//     result_dmc_one_step result;
//     result.walkers = (double*) malloc(sizeof(double) * n_multiplier);
//     memcpy(result.walkers, new_walkers, sizeof(double) * n_multiplier);
//     free(new_walkers);
//     result.n = n_multiplier;
//     return result;
// }

// double displace_x(double x, double delta_tau, gsl_rng* r)
// {
//     double G = gsl_ran_gaussian(r, 1);
//     return x + sqrt(delta_tau) * G;
// }

// double potential_1d(double x)
// {
//     double y = 1 - exp(- x);
//     return 0.5 * y * y;
// }

// double weight(double x, double E_t, double dt)
// {
//     return exp(- (potential_1d(x) - E_t) * dt);
// } 

// double update_E_t(double E_t, double gamma, int n, int n0)
// {
//     return E_t - gamma * log((double) n / n0);
// }

// gsl_rng* get_rand(void){
//     const gsl_rng_type* T;
//     gsl_rng* r;
//     gsl_rng_env_setup();
//     T = gsl_rng_default;
//     r = gsl_rng_alloc(T);
//     time_t seed = time(NULL);
//     gsl_rng_set(r, seed);
//     return r;
// }
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
gsl_rng* get_rand(void);
result_dmc diffusion_monte_carlo(double* walkers, int n0, double E_t, double gamma, double dt);
result_dmc_one_step diffusion_monte_carlo_one_step(double* walkers, int n, double E_t, double dt, gsl_rng* r);
double displace_x(double x, double delta_tau, gsl_rng* r);
double weight(double x, double E_t, double dt);
double update_E_t(double E_t, double gamma, int n, int n0);
result_dmc diffusion_monte_carlo_1d(double* walkers, int n0, double E_t, double gamma, double dt, int n_iter);

int
run(
    int argc,
    char *argv[]
   )
{
    task1();
    return 0;
}

void task1(void){

    int n = 200;
    double* walkers = linspace(-5, 5, n, false);
    double E_t = 0.5;
    double delta_tau = 0.02;
    double gamma = 0.5;
    double n_iter = 10000;
    diffusion_monte_carlo_1d(walkers, n, E_t, gamma, delta_tau, n_iter);
}

result_dmc diffusion_monte_carlo_1d(double* walkers, int n0, double E_t, double gamma, double dt, int n_iter)
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
            fprintf(positions, "%.10lf\n", walkers[j]);
        }
        fprintf(file, "%d, %lf\n", n, E_t);

        walker_multiplier = (int*) malloc(sizeof(int) * n);
        int multiplier_sum = 0; 

        for (int j = 0; j < n; j++)
        {

            walkers[j] = displace_x(walkers[j], dt, r);
            double w = weight(walkers[j], E_t, dt);
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
        E_t_sum += E_t;
        E_t = update_E_t(E_t_sum / (i + 1), gamma, multiplier_sum, n0);
        n = multiplier_sum;
    }
    fclose(file);
    fclose(positions);
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

double weight(double x, double E_t, double dt)
{
    return exp(- (potential_1d(x) - E_t) * dt);
} 

double update_E_t(double E_t, double gamma, int n, int n0)
{
    return E_t - gamma * log((double) n / n0);
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

