#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <string.h>
#include <time.h>
#include "tools.h"
#include <stdbool.h>

void task1();

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



    free(walkers);

}


