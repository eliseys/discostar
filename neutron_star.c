#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "def.h"

sp arcgen(double theta, double delta_theta, double N_theta, double N_phi, int i)
{
  /* i is the index of point */
  sp result;
  
  if (theta != 0.0)
    {
      /* arc */
      
      double epsilon_theta = delta_theta/N_theta;
      double epsilon_phi = (2.0 * M_PI)/N_phi;
            
      result.phi = epsilon_phi * (i % (int) (N_phi + 1.0));
      result.theta = theta - 0.5 * delta_theta + 0.5 * epsilon_theta + epsilon_theta * (i % (int) (N_theta - 1.0));
      result.r = 1.0;
    }
  else
    {
      /* magnetic pole */

      double epsilon_theta = 0.5 * delta_theta/N_theta;
      double epsilon_phi = (2.0 * M_PI)/N_phi;

      if (i == 0)
	{}
      else
	{
	  result.phi = epsilon_phi * (i % (int) (N_phi + 1.0));
	  result.theta = 0.5 * epsilon_theta + epsilon_theta * (i % (int) (N_theta - 1.0));
	  result.r = 1.0;
	}
    }
      
  return result;
}


double * x_ray_direction_diagram(double PSI_pr)
{

  //double PSI_pr = 330.0 * M_PI/180.0 ; /* NS precession angle */
  double angle_JI = 50.0 * M_PI/180.0; /* angle between J and I_3 vectors */

  double N_theta = 50.0; /* arc pixelization */
  double N_phi = 1500.0; /* arc pixelization */

  double N = N_theta * N_phi;

  int arc_n; /* arc index */
  int i, f; /* point index, directional diagram index */
  vec3 decrt;
  sp sphr;

  double F[180]; /* flux in f-th solid angle */
  double I[180]; /* mean intensity in f-th solid angle */

  for (f = 0; f < 180; f++)
    {
      F[f] = 0.0;
    }
  
  double b[18], l[18], R[18], BEG[18], END[18], delta_theta[18], intensity[18];

  /* !!!!!! I`d changed b[2] from 65.0 to 75.0 */
  /* !!!!!! I`d changed END[8] from 350.0 to 380.0 */
  
  b[0] = 70.0,    l[0] = 160.0,  R[0] = 50.0,   BEG[0] = 205.0,  END[0] = 250.0,  delta_theta[0] = 14.0,  intensity[0] = 54.3;  
  b[1] = 75.0,    l[1] = 180.0,  R[1] = 43.0,   BEG[1] = 250.0,  END[1] = 300.0,  delta_theta[1] = 13.0,  intensity[1] = 53.1;  
  b[2] = 75.0,    l[2] = 180.0,  R[2] = 43.0,   BEG[2] = 302.0,  END[2] = 350.0,  delta_theta[2] = 12.0,  intensity[2] = 51.5;  
  b[3] = 65.0,    l[3] = 170.0,  R[3] = 53.0,   BEG[3] = 355.0,  END[3] = 400.0,  delta_theta[3] = 11.0,  intensity[3] = 19.3;  
  b[4] = 60.0,    l[4] = 160.0,  R[4] = 55.0,   BEG[4] = 45.0,   END[4] = 60.0,   delta_theta[4] = 11.0,  intensity[4] = 4.5;  
  b[5] = 60.0,    l[5] = 150.0,  R[5] = 50.0,   BEG[5] = 65.0,   END[5] = 95.0,   delta_theta[5] = 11.0,  intensity[5] = 20.2;  
  b[6] = 58.0,    l[6] = 160.0,  R[6] = 55.0,   BEG[6] = 95.0,   END[6] = 120.0,  delta_theta[6] = 11.0,  intensity[6] = 25.8;  
  b[7] = 70.0,    l[7] = 140.0,  R[7] = 60.0,   BEG[7] = 120.0,  END[7] = 160.0,  delta_theta[7] = 13.0,  intensity[7] = 29.0;  
  b[8] = -97.0,   l[8] = 180.0,  R[8] = 44.0,   BEG[8] = 330.0,  END[8] = 380.0,  delta_theta[8] = 11.0,  intensity[8] = 4.1;  
  b[9] = -97.0,   l[9] = 180.0,  R[9] = 44.0,   BEG[9] = 20.0,   END[9] = 50.0,   delta_theta[9] = 11.0,  intensity[9] = 7.6;  
  b[10] = -97.0,  l[10] = 180.0, R[10] = 44.0,  BEG[10] = 50.0,  END[10] = 100.0, delta_theta[10] = 11.0, intensity[10] = 10.2;  
  b[11] = -95.0,  l[11] = 200.0, R[11] = 42.0,  BEG[11] = 110.0, END[11] = 170.0, delta_theta[11] = 8.0,  intensity[11] = 17.0;  
  b[12] = -100.0, l[12] = 200.0, R[12] = 65.0,  BEG[12] = 205.0, END[12] = 260.0, delta_theta[12] = 11.0, intensity[12] = 3.6;  
  b[13] = 60.0,   l[13] = 180.0, R[13] = 80.0,  BEG[13] = 250.0, END[13] = 280.0, delta_theta[13] = 11.0, intensity[13] = 11;  
  b[14] = 75.0,   l[14] = 180.0, R[14] = 105.0, BEG[14] = 90.0,  END[14] = 120.0, delta_theta[14] = 11.0, intensity[14] = 2.7;  

  /* magnetic poles: N, S, P1 */
  b[15] = 60.0,   l[15] = 180.0, R[15] = 0.0, BEG[15] = 0.0,  END[15] = 360.0, delta_theta[15] = 16.0, intensity[15] = 100.0;  
  b[16] = -85.0,  l[16] = 0.0,   R[16] = 0.0, BEG[16] = 0.0,  END[16] = 360.0, delta_theta[16] = 16.0, intensity[16] = 32.6;  
  b[17] = -17.0,  l[17] = 80.0,  R[17] = 0.0, BEG[17] = 0.0,  END[17] = 360.0, delta_theta[17] = 14.0, intensity[17] = 0.94;  
  
  double theta_n, delta_theta_n, epsilon_theta_n;
  double N_phi_n, epsilon_phi_n;
  double ds; /* surface element */

  for (arc_n = 0; arc_n < 18; arc_n ++)
    {
  
      theta_n = R[arc_n] * M_PI/180.0;
      delta_theta_n = delta_theta[arc_n] * M_PI/180.0;
      
      epsilon_theta_n = delta_theta_n/N_theta;
      epsilon_phi_n = (2.0 * M_PI)/N_phi;


      for(i = 0; i <= N; i++)
	{
	  sphr = arcgen(theta_n, delta_theta_n, N_theta, N_phi, i);
	  
	  ds = epsilon_theta_n * epsilon_phi_n * sin(sphr.theta); /* solid angle element */
	  
	  decrt = sp2dec(sphr);

	  decrt = rotate(decrt, (-(90.0 - b[arc_n]) * M_PI/180.0), (-l[arc_n] * M_PI/180.0));
	  
	  sphr = dec2sp(decrt);

	  if ( END[arc_n] <= 360.0 && sphr.phi >= (BEG[arc_n] * M_PI/180.0) && sphr.phi <= (END[arc_n] * M_PI/180.0) )
	    {

	      decrt = sp2dec(sphr);
	      decrt = rotate(decrt, 0.0, -PSI_pr);
	      decrt = rotate(decrt, angle_JI, 0.0);

	      //printf("%.20f\t %.20f\t %.20f\t %.20f\n", decrt.x, decrt.y, decrt.z, intensity[arc_n]);
	      
	      sphr = dec2sp(decrt);
	      
	      
	      f = (int) floor( sphr.theta * 180.0/M_PI );
	      
	      F[f] = F[f] + ds * intensity[arc_n];
	      
	    }
	  else if ( END[arc_n] > 360.0 && sphr.phi >= (BEG[arc_n] * M_PI/180.0) || sphr.phi <= (END[arc_n] * M_PI/180.0 - 2.0 * M_PI) )
	    {
	      
	      decrt = sp2dec(sphr);
	      decrt = rotate(decrt, 0.0, -PSI_pr);
	      decrt = rotate(decrt, angle_JI, 0.0);

	      //printf("%.20f\t %.20f\t %.20f\t %.20f\n", decrt.x, decrt.y, decrt.z, intensity[arc_n]);
	      
	      sphr = dec2sp(decrt);

	      //ds = epsilon_theta_n * epsilon_phi_n * sin(sphr.theta); /* surface element, r = 1.0 */

	      f = (int) floor( sphr.theta * 180.0/M_PI );
	      
	      F[f] = F[f] + ds * intensity[arc_n];

	    }

	}

      /* printf("\n\n"); */
      /* printf("0.00000\t 0.00000\t 1.00000\n0.00000\t 0.00000\t -1.00000\n"); */
      /* printf("0.00000\t 1.00000\t 0.00000\n0.00000\t -1.00000\t 0.00000\n"); */
      /* printf("1.00000\t 0.00000\t 0.00000\n-1.00000\t 0.00000\t 0.00000\n"); */
    }

  double * result = (double *) malloc(sizeof(double) * 180);
  
  double I_sum = 0.0;
  double f_angle;

  for (f = 0; f < 180; f++)
    {
      //printf("%i\t %.20f\t %.20f\n", j, F[j], F_sum);

      f_angle = (double) (f * M_PI/180.0 + 0.5 * M_PI/180.0); 
	
      I[f] = F[f]/(2.0 * M_PI * sin(f_angle) * M_PI/180.0); /* mean intensity in f-th zone */

      //printf("%.20f\n", result[f]);
    }
  
  for (f = 0; f < 180; f++)
    {
      I_sum = I_sum + I[f];
      //printf("%.20f\n", F_sum);
    }
  
  for (f = 0; f < 180; f++)
    {
      result[f] = I[f]/I_sum;
      
      //printf("%d\t %f\t %f\n", f, result[f], I_sum);
    }

  
  return result;
}
