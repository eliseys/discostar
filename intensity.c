#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "def.h"

double B_lambda(double T, double lambda)
{
  return 2.0 * H_PLANCK * C * C * pow(lambda,-5.0) / (exp(H_PLANCK * C/(lambda * K * T)) - 1.0);
}
