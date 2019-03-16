#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "def.h"

double F_lambda(double T, double lambda)
{

  /* intencity at specified lambda */


  /* C1 and C2 constants defined in def.h */
  return C1 * pow(lambda,-5.0) / (exp(C2/(lambda * T)) - 1.0);
}


double F_filter(double T, char filter)
{

  /* flux through the filter */

  
  double summa = 0.0;
  int i;

  double A_cm = 1E-8;
  
  /* normalized transmission*/
  /* ftp://obsftp.unige.ch/pub/mermio/filters/ */
  /* ApJ 141, 923 (1965) */
  
  double lambda_B[] = {3600, 3650, 3700, 3750, 3800, 3850, 3900, 3950, 4000, 4050, 4100, 4150, 4200, 4250, 4300, 4350, 4400, 4450, 4500, 4550, 4600, 4650, 4700, 4750, 4800, 4850, 4900, 4950, 5000, 5050, 5100, 5150, 5200, 5250, 5300, 5350, 5400, 5450, 5500, 5550};

  double transmisson_B[] = {0.000, 0.000, 0.020, 0.050, 0.110, 0.180, 0.350, 0.550, 0.920, 0.950, 0.980, 0.990, 1.000, 0.990, 0.980, 0.960, 0.940, 0.910, 0.870, 0.830, 0.790, 0.740, 0.690, 0.630, 0.580, 0.520, 0.460, 0.410, 0.360, 0.300, 0.250, 0.200, 0.150, 0.120, 0.090, 0.060, 0.040, 0.020, 0.010, 0.000};

  double lambda_V[] = {4600, 4650, 4700, 4750, 4800, 4850, 4900, 4950, 5000, 5050, 5100, 5150, 5200, 5250, 5300, 5350, 5400, 5450, 5500, 5550, 5600, 5650, 5700, 5750, 5800, 5850, 5900, 5950, 6000, 6050, 6100, 6150, 6200, 6250, 6300, 6350, 6400, 6450, 6500, 6550, 6600, 6650, 6700, 6750, 6800, 6850, 6900, 6950, 7000, 7050, 7100, 7150, 7200, 7250, 7300, 7350};
  
  double transmisson_V[] = {0.000, 0.000, 0.010, 0.010, 0.020, 0.051, 0.112, 0.204, 0.388, 0.684, 0.796, 0.867, 0.929, 0.959, 0.980, 1.000, 1.000, 0.969, 0.888, 0.806, 0.735, 0.724, 0.704, 0.663, 0.633, 0.592, 0.531, 0.469, 0.408, 0.347, 0.296, 0.245, 0.204, 0.173, 0.143, 0.112, 0.082, 0.061, 0.051, 0.031, 0.020, 0.020, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.000};
 
  if (filter == 'B')
    {
      double I_lambda[40] = {0.0};

      for (i = 0; i < 40; i++)
	{
	  I_lambda[i] = C1 * pow(lambda_B[i]*A_cm,-5.0) / (exp(C2/(lambda_B[i]*A_cm * T)) - 1.0);
	}
      for (i = 1; i < 40; i++)
	{
	  summa = summa + ( (lambda_B[i]*A_cm - lambda_B[i-1]*A_cm) * (I_lambda[i-1] * transmisson_B[i-1] + I_lambda[i] * transmisson_B[i]) * 0.5 );
	}
    }
  else if (filter == 'V')
    {

      double I_lambda[56];

      for (i = 0; i < 56; i++)
	{
	  I_lambda[i] = C1 * pow(lambda_V[i],-5.0) / (exp(C2/(lambda_V[i] * T)) - 1.0);
	}
      for (i = 0; i < 56; i++)
	{
	  summa = summa + I_lambda[i] * transmisson_V[i];
	}
       
    }

  
  return summa;
  
}
