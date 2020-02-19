#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "def.h"
#include <omp.h>
#include <errno.h>


double F_filter(double T, char * filter)
{

  /* flux through the filter */

  double lambda;

  char * filter_name;
  
  lambda = strtod(filter, &filter_name); 


  /* printf("++++++++ %lf\n", lambda); */


  
  if ( lambda == 0.0 && strcmp(filter_name, "") == 0 )
    {
      printf("filter is empty\n");
      return EINVAL;
    }
  else if ( lambda != 0.0 && strcmp(filter_name, "") == 0 )
    {
      return 1.0 * A_cm * C1 * pow(lambda * A_cm,-5.0) / (exp(C2/(lambda * A_cm * T)) - 1.0);
    }
  else if ( lambda == 0.0 && filter_name != "" )
    {
      if ( strcmp(filter_name, "B") != 0 && strcmp(filter_name, "V") != 0 && strcmp(filter_name, "WASP") != 0 )
	{
	  printf("444 %s filter does not exist\n", filter_name );
	  return EINVAL;  
	}
      
      double summa = 0.0;
      int i;
  
      double I_lambda_0;
      double I_lambda_1;

  

      double lambda_B[] = {3765.0,
			   3915.0,
			   4065.0,
			   4215.0,
			   4365.0,
			   4515.0,
			   4665.0,
			   4815.0,
			   4965.0,
			   5115.0,
			   5265.0,
			   5415.0,
			   5565.0,
			   5715.0,
			   5865.0};


      double transmisson_B[] = {0.002609,
				0.683457,
				0.779870,
				0.796901,
				0.779850,
				0.734233,
				0.632581,
				0.415888,
				0.216261,
				0.088407,
				0.023910,
				0.010236,
				0.018927,
				0.011319,
				0.001020};


		      

      double lambda_V[] = {4845.0,
			   4995.0,
			   5145.0,
			   5295.0,
			   5445.0,
			   5595.0,
			   5745.0,
			   5895.0,
			   6045.0,
			   6195.0,
			   6345.0,
			   6495.0,
			   6645.0,
			   6795.0,
			   6945.0};


      double transmisson_V[] = {0.049660,
				0.710850,
				0.931175,
				0.936013,
				0.899921,
				0.829204,
				0.719113,
				0.574333,
				0.413209,
				0.262001,
				0.143645,
				0.067348,
				0.026601,
				0.008895,
				0.002514};


      double lambda_WASP[] = {4100.0,
			      4395.0,
			      4690.0,
			      4985.0,
			      5280.0,
			      5575.0,
			      5870.0,
			      6165.0,
			      6460.0,
			      6755.0,
			      7050.0};


      double transmisson_WASP[] = {0.85,
				   0.85,
				   0.85,
				   0.85,
				   0.85,
				   0.85,
				   0.85,
				   0.85,
				   0.85,
				   0.85,
				   0.85};
  
      if (strcmp(filter_name, "B") == 0)
	{

	  I_lambda_0 = C1 * pow(lambda_B[0]*A_cm,-5.0) / (exp(C2/(lambda_B[0] * A_cm * T)) - 1.0);
	  i = 1;
      
	  while (i < 15)
	    {
	      I_lambda_1 = C1 * pow(lambda_B[i]*A_cm,-5.0) / (exp(C2/(lambda_B[i] * A_cm * T)) - 1.0);
	  
	      summa = summa + ( ((lambda_B[i]-lambda_B[i-1])*A_cm) * (I_lambda_1 * transmisson_B[i] + I_lambda_0 * transmisson_B[i-1]) * 0.5 );
	  
	      I_lambda_0 = I_lambda_1;
	  
	      i++;
	    }

	}
  
      else if (strcmp(filter_name, "V") == 0)
	{
      
	  I_lambda_0 = C1 * pow(lambda_V[0]*A_cm,-5.0) / (exp(C2/(lambda_V[0] * A_cm * T)) - 1.0);
	  i = 1;
      
	  while (i < 15)
	    {
	      I_lambda_1 = C1 * pow(lambda_V[i]*A_cm,-5.0) / (exp(C2/(lambda_V[i] * A_cm * T)) - 1.0);
	  
	      summa = summa + ( ((lambda_V[i]-lambda_V[i-1])*A_cm) * (I_lambda_1 * transmisson_V[i] + I_lambda_0 * transmisson_V[i-1]) * 0.5 );
	  
	      I_lambda_0 = I_lambda_1;
	  
	      i++;
	    }

       
	}

      else if (strcmp(filter_name, "WASP") == 0)
	{
      
	  I_lambda_0 = C1 * pow(lambda_WASP[0]*A_cm,-5.0) / (exp(C2/(lambda_WASP[0] * A_cm * T)) - 1.0);
	  i = 1;
      
	  while (i < 11)
	    {
	      I_lambda_1 = C1 * pow(lambda_WASP[i]*A_cm,-5.0) / (exp(C2/(lambda_WASP[i] * A_cm * T)) - 1.0);
	  
	      summa = summa + ( ((lambda_WASP[i]-lambda_WASP[i-1])*A_cm) * (I_lambda_1 * transmisson_WASP[i] + I_lambda_0 * transmisson_WASP[i-1]) * 0.5 );
	  
	      I_lambda_0 = I_lambda_1;
	  
	      i++;
	    }
	}
  
      return summa;
    }


}
