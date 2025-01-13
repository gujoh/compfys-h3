#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <string.h>
#include <time.h>
#include "tools.h"
#include <stdbool.h>

void task1();
gsl_rng* get_rand(void);

int
run(
    int argc,
    char *argv[]
   )
{
    task1();
    return 0;
}

void task1(){

    double* walkers = linspace(-5, 5, 200, false);
    double E_t = 0.5;
    double delta_tau = 0.02;
    gsl_rng* r = get_rand();



    free(walkers);

}

float displace_x(float x, float delta_tau, gsl_rng* r)
{
    float G = gsl_ran_gaussian(r, 1);
    return x + sqrt(delta_tau)*G;
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