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


  //printf("filter %s\n", filter);
  
  /* printf("++++++++ %lf\n", lambda); */


  /* printf("strcmp(filter_name, \"V\") %d\n", strcmp(filter_name, "V")); */
  /* printf("filter_name %s\n", filter_name); */
  /* printf("%s\n", ""); */
  
  
  
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
      if ( strcmp(filter_name, "B") != 0 && strcmp(filter_name, "V") != 0 && strcmp(filter_name, "WASP") != 0 && strcmp(filter_name, "vis") != 0)
	{
	  //printf("%s filter does not exist\n", filter_name );
	  return EINVAL;  
	}
      
      double summa = 0.0;
      int i;
  
      double I_lambda_0;
      double I_lambda_1;
      double I_energy_0;
      double I_energy_1;


  

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


      double lambda_vis[] = {3500.0,
			     3600.0,
			     3700.0,
			     3800.0,
			     3900.0,
			     4000.0,
			     4100.0,
			     4200.0,
			     4300.0,
			     4400.0,
			     4500.0,
			     4600.0,
			     4700.0,
			     4800.0,
			     4900.0,
			     5000.0,
			     5100.0,
			     5200.0,
			     5300.0,
			     5400.0,
			     5500.0,
			     5600.0,
			     5700.0,
			     5800.0,
			     5900.0,
			     6000.0,
			     6100.0,
			     6200.0,
			     6300.0,
			     6400.0,
			     6500.0,
			     6600.0,
			     6700.0,
			     6800.0,
			     6900.0,
			     7000.0};


      double transmisson_vis[] = {1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0,
				  1.0};






      
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

      else if (strcmp(filter_name, "vis") == 0)
	{
      
	  I_lambda_0 = C1 * pow(lambda_vis[0]*A_cm,-5.0) / (exp(C2/(lambda_vis[0] * A_cm * T)) - 1.0);
	  i = 1;
      
	  while (i < 11)
	    {
	      I_lambda_1 = C1 * pow(lambda_vis[i]*A_cm,-5.0) / (exp(C2/(lambda_vis[i] * A_cm * T)) - 1.0);
	  
	      summa = summa + ( ((lambda_vis[i]-lambda_vis[i-1])*A_cm) * (I_lambda_1 * transmisson_vis[i] + I_lambda_0 * transmisson_vis[i-1]) * 0.5 );
	  
	      I_lambda_0 = I_lambda_1;
	  
	      i++;
	    }
	}

   
      return summa;
      
    }


}



double kappa_mm(double energy)
{
  /* Robert Morrison and Dan McCammon, 1983 */
  /* energy in keV */

  
  double kappa;
  
  double c0[] = {
		 17.3,
		 34.6,
		 78.1,
		 71.4,
		 95.5,
		 308.9,
		 120.6,
		 141.3,
		 202.7,
		 342.7,
		 352.2,
		 433.9,
		 629.0,
		 701.2};

  double c1[] = {
		 608.1,
		 267.9,
		 18.8,
		 66.8,
		 145.8,
		 -380.6,
		 169.3,
		 146.8,
		 104.7,
		 18.7,
		 18.7,
		 -2.4,
		 30.9,
		 25.2};

  double c2[] = {
		 -2150.0,
		 -476.1,
		 4.3,
		 -51.4,
		 -61.1,
		 294.0,
		 -47.7,
		 -31.5,
		 -17.0,
		 0.0,
		 0.0,
		 0.75,
		 0.0,
		 0.0};
    


  double energy_ranges[] = {
			  0.030,
			  0.100,
			  0.284,
			  0.400,
			  0.532,
			  0.707,
			  0.867,
			  1.303,
			  1.840,
			  2.471,
			  3.210,
			  4.038,
			  7.111,
			  8.331,
			  10.000};


  int i;
  
  for (i = 0; i < 14; i++)
    {
      if ((energy - energy_ranges[i]) >= 0.0 && (energy_ranges[i+1] - energy) > 0.0)
	{
	  kappa = (c0[i] + c1[i]*energy + c2[i]*energy*energy)*pow(energy, -3.0);
	}
      else if ( fabs(energy - 10.0) < eps )
	{
	  kappa = (c0[13] + c1[13]*energy + c2[13]*energy*energy)*pow(energy, -3.0);
	}
      else
	{
	  continue;
	
	}

    }

  
  
  return kappa * 1.0e-24;

}



double x_ray_reflected_fraction(double zeta, double energy)
{
  /* Sobolev approximation */

  //double lambda = SIGMA_THOMSON/(SIGMA_THOMSON + kappa_mm(energy));

  double kappa = 1.0 * SIGMA_THOMSON;
  double lambda = SIGMA_THOMSON/(SIGMA_THOMSON + kappa);

  //printf("lambda %f\n", );
  
  double k = sqrt(3.0*(1.0-lambda));

  double D = (3.0*lambda*zeta*zeta)/(1.0 - k*k*zeta*zeta);

  double D2 = D * (1.0 + 2.0/(3.0*zeta))/(1.0 + (2.0*k)/3.0);


  //printf("kappa %e\n", kappa_mm(energy));
  return D/(3.0*zeta*zeta) - k*D2/(3.0*zeta);

  
}



double normalized_x_ray_spectrum(double energy)
{

  /* Her X-1 spectral parameters from */
  /* Vrtilek and Halpern, ApJ 1985, table 2 */

  /* energy in keV */

  double A_BB = 417.0;
  double A_PL = 0.126;
  double alpha_1 = 0.93;
  double kT = 0.12; /* keV */


  /* enegy range 0.3 - 2.0 keV */
  //double normalization = 1.0/0.6167598821359219;

  /* enegy range 0.3 - 1.0 keV */
  //double normalization = 1.0/0.4697978454640511;

  /*energy range 2.0 - 6.0 keV */
  //double normalization = 1.0/0.5537607286596092;

  /* enegy range 1.0 - 6.0 keV */
  double normalization = 1.0/0.70072276533148;


  
  /* F(energy) */ 
  return normalization*energy*(A_BB*energy*energy/(exp(energy/kT) - 1.0) + A_PL*pow(energy, -alpha_1));
  
}





double x_ray_flux_integrated(double zeta)
{

  int i = 0;
  double F0, F1;


  double n_steps = 400.0;

  //double energy_step = (2.0 - 0.3)/n_steps;
  //double energy_step = (1.0 - 0.3)/n_steps;
  //double energy_step = (6.0 - 2.0)/n_steps;
  double energy_step = (6.0 - 1.0)/n_steps;

  //double energy_0 = 0.3;
  //double energy_0 = 2.0;
  double energy_0 = 1.0;


  
  double energy_1;

  double summa = 0;

  
  while (i < 401)
    {

      energy_1 = energy_0 + energy_step;


 
      F0 = normalized_x_ray_spectrum(energy_0) * x_ray_reflected_fraction(zeta, energy_0);
     
      F1 = normalized_x_ray_spectrum(energy_1) * x_ray_reflected_fraction(zeta, energy_1);
      
 
      /* F0 = Lx * normalized_x_ray_spectrum(energy_0); */
     
      /* F1 = Lx * normalized_x_ray_spectrum(energy_1); */

      summa = summa + energy_step * (F0+F1)/2.0;


      //printf("%f\t%e\n", (energy_0 + energy_1)/2.0, (F0 + F1)/2.0);

      energy_0 = energy_1;
      i++; 

    }

  return summa;
  

}

/* double x_ray_flux(double zeta, ) */
/* { */
/*   /\* integrated by energy range *\/ */


/* 	  I_energy_0 = Lx_spektrum(energy_03_2keV[0]); */
	  

/* 	  i = 1; */
      
/* 	  while (i < 11) */
/* 	    { */
/* 	      I_energy_1 = Lx_spektrum(energy_03_2keV[1]); */
	  
/* 	      summa = summa + ( (energy_03_2keV[i]-energy_03_2keV[i-1]) * (x_rays_reflected_fraction(energy_03_2keV[1]) + x_rays_reflected_fraction(energy_03_2keV[0])) * 0.5 ); */
	  
/* 	      I_energy_0 = I_energy_1; */
	  
/* 	      i++; */
/* 	    } */

/* } */
