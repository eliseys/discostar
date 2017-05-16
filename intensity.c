#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "def.h"

double F_lambda(double T, double lambda)
{
  /* C1 and C2 constants defined in def.h */
  
  return C1 * pow(lambda,-5.0) / (exp(C2/(lambda * T)) - 1.0);
}
