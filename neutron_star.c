#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "def.h"

sp arcgen(double theta, double N_phi, int i)
{
  /* i is the index of point */

  sp result;
      
  double epsilon_phi = (2.0 * M_PI)/N_phi;

  if (theta != 0.0)
    {
      result.phi = epsilon_phi * (i % (int) (N_phi + 1.0));
      result.theta = theta;
      result.r = 1.0;
    }
  else if (theta == 0.0 && i == 0)
    {
      result.phi = 0.0;
      result.theta = 0.0;
      result.r = 1.0;
    }
  else
    {}
  
  return result;
}




double * x_ray_direction_diagram(double PSI_pr)
{

  //double angle_JI = 50.0 * M_PI/180.0; /* angle between J and I_3 vectors */
  double angle_JI = 50.0 * M_PI/180.0; /* angle between J and I_3 vectors */
    
  /* sky pixelization */
  double n_theta;
  double n_phi;
  
  double N_phi_n = 100.0; /* arc pixelization */

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

  /* !!!!!! I`ve changed b[2] from 65.0 to 75.0 */
  /* !!!!!! I`ve changed END[8] from 350.0 to 380.0 */
  
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

  /* b[0] = 70.0,    l[0] = 160.0,  R[0] = 50.0,   BEG[0] = 205.0,  END[0] = 250.0,  delta_theta[0] = 14.0,  intensity[0] = 54.3; */
  /* b[1] = 75.0,    l[1] = 180.0,  R[1] = 43.0,   BEG[1] = 250.0,  END[1] = 300.0,  delta_theta[1] = 13.0,  intensity[1] = 53.1; */
  /* b[2] = 75.0,    l[2] = 180.0,  R[2] = 43.0,   BEG[2] = 302.0,  END[2] = 350.0,  delta_theta[2] = 12.0,  intensity[2] = 51.5; */
  /* b[3] = 65.0,    l[3] = 170.0,  R[3] = 53.0,   BEG[3] = 355.0,  END[3] = 400.0,  delta_theta[3] = 11.0,  intensity[3] = 19.3; */
  /* b[4] = 60.0,    l[4] = 160.0,  R[4] = 55.0,   BEG[4] = 45.0,   END[4] = 60.0,   delta_theta[4] = 11.0,  intensity[4] = 4.5; */
  /* b[5] = 60.0,    l[5] = 150.0,  R[5] = 50.0,   BEG[5] = 65.0,   END[5] = 95.0,   delta_theta[5] = 11.0,  intensity[5] = 20.2; */
  /* b[6] = 58.0,    l[6] = 160.0,  R[6] = 55.0,   BEG[6] = 95.0,   END[6] = 120.0,  delta_theta[6] = 11.0,  intensity[6] = 25.8; */
  /* b[7] = 70.0,    l[7] = 140.0,  R[7] = 60.0,   BEG[7] = 120.0,  END[7] = 160.0,  delta_theta[7] = 13.0,  intensity[7] = 29.0; */

  /* b[8] = -97.0,   l[8] = 180.0,  R[8] = 44.0,   BEG[8] = 330.0,  END[8] = 380.0,  delta_theta[8] = 11.0,  intensity[8] = 4.1; */
  /* b[9] = -97.0,   l[9] = 180.0,  R[9] = 44.0,   BEG[9] = 20.0,   END[9] = 50.0,   delta_theta[9] = 11.0,  intensity[9] = 7.6; */
  /* b[10] = -97.0,  l[10] = 180.0, R[10] = 44.0,  BEG[10] = 50.0,  END[10] = 100.0, delta_theta[10] = 11.0, intensity[10] = 10.2; */
  /* b[11] = -95.0,  l[11] = 200.0, R[11] = 42.0,  BEG[11] = 110.0, END[11] = 170.0, delta_theta[11] = 8.0,  intensity[11] = 17.0; */

  /* b[12] = -100.0, l[12] = 200.0, R[12] = 65.0,  BEG[12] = 205.0, END[12] = 260.0, delta_theta[12] = 11.0, intensity[12] = 3.6; */
  /* b[13] = 60.0,   l[13] = 180.0, R[13] = 80.0,  BEG[13] = 250.0, END[13] = 280.0, delta_theta[13] = 11.0, intensity[13] = 11; */
  /* b[14] = 75.0,   l[14] = 180.0, R[14] = 105.0, BEG[14] = 90.0,  END[14] = 120.0, delta_theta[14] = 11.0, intensity[14] = 2.7; */

  /* /\* magnetic poles: N, S, P1 *\/ */
  /* b[15] = 60.0,   l[15] = 180.0, R[15] = 0.0, BEG[15] = 0.0,  END[15] = 360.0, delta_theta[15] = 16.0, intensity[15] = 100.0; */
  /* b[16] = -85.0,  l[16] = 0.0,   R[16] = 0.0, BEG[16] = 0.0,  END[16] = 360.0, delta_theta[16] = 16.0, intensity[16] = 32.6; */
  /* b[17] = -17.0,  l[17] = 80.0,  R[17] = 0.0, BEG[17] = 0.0,  END[17] = 360.0, delta_theta[17] = 14.0, intensity[17] = 0.94; */

  
  double theta_n, delta_theta_n, epsilon_theta_n;
  double epsilon_phi_n;
  
  double dl; /* line element */


  int j=0;
  int k=0;
    
  int Nj = 360;
  int Nk = 180;
  double alpha[Nj][Nk]; /* sky coordinate */ 
  double delta[Nj][Nk]; /* sky coordinate */
  double flux[Nj][Nk];
  double ds_sky[Nj][Nk]; /* surface element on the sky */
  double flux_sum[Nk]; /* pre - directional diagram */

  double flux_through_ds = 0.0;
  double cos_p_decrt;
  double p_decrt;

  sp sp_sky; /* point on the sky spherical coordinates */
  vec3 p_sky; /* point on the sky cartesian coordinates */

  int arc_len_count = 0;
  
  
  for (j = 0; j < Nj; j++)
    {      
      for (k = 1; k < Nk; k++)
	{
	  alpha[j][k] = ((double) j)/((double) Nj) * 2.0 * M_PI; 
          delta[j][k] = ((double) k)/((double) Nk) * M_PI; 
	  ds_sky[j][k] = ((2.0 * M_PI)/((double) Nj)) * (M_PI/((double) Nk)) * sin(delta[j][k]); /* surface element on the sky */
	  //ds_sky[j][k] = ((2.0 * M_PI)/((double) Nj)) * (M_PI/((double) Nk));
	  flux[j][k] = 0.0;
	  flux_sum[k] = 0.0; /* pre - directional diagram */

	  //flux_through_ds = 0.0;

	  //sp_sky.phi = delta;
	  //sp_sky.theta = alpha;
	  //sp_sky.r = 1.0;

	  //p_sky = sp2dec(sp_sky);
	  
	  //printf("%.20f\t%.20f\t1.000000\n", alpha[j][k], delta[j][k]);
	}

    }





	  
  for (arc_n = 0; arc_n < 18; arc_n ++)
    {

      arc_len_count = 0;
      
      theta_n = R[arc_n] * M_PI/180.0;
      delta_theta_n = delta_theta[arc_n] * M_PI/180.0;
      
      if (arc_n == 15 || arc_n == 16 || arc_n == 17)
	{
	  N_phi_n = 1.0;
	}
      
      epsilon_phi_n = (2.0 * M_PI)/N_phi_n;
      
      for(i = 0; i < N_phi_n; i++)
	{
	  sphr = arcgen(theta_n, N_phi_n, i);

	  if (sphr.theta != 0.0)
	    {
	      dl = epsilon_phi_n * sin(sphr.theta)/(2.0 * M_PI); /* line element */
	    }
	  else if (sphr.theta == 0.0)
	    {
	      dl = 1.0; /* line element */
	    }
	  
	  decrt = sp2dec(sphr);	  
	  decrt = rotate(decrt, (-(90.0 - b[arc_n]) * M_PI/180.0), (-l[arc_n] * M_PI/180.0));		  
	  sphr = dec2sp(decrt);

	  if ( END[arc_n] <= 360.0 && sphr.phi >= (BEG[arc_n] * M_PI/180.0) && sphr.phi <= (END[arc_n] * M_PI/180.0) )
	    {
	      arc_len_count++;
	    }
	  else if ( END[arc_n] > 360.0 && sphr.phi >= (BEG[arc_n] * M_PI/180.0) || sphr.phi <= (END[arc_n] * M_PI/180.0 - 2.0 * M_PI) )
	    {
	      arc_len_count++;
	    }
	}

      //printf("%d\t%d\n", arc_n, arc_len_count);

      for(i = 0; i < N_phi_n; i++)
	{
	  sphr = arcgen(theta_n, N_phi_n, i);

	  if (sphr.theta != 0.0)
	    {
	      dl = epsilon_phi_n * sin(sphr.theta)/(2.0 * M_PI); /* line element */
	    }
	  else if (sphr.theta == 0.0)
	    {
	      dl = 1.0; /* line element */
	    }
	  
	  decrt = sp2dec(sphr);	  
	  decrt = rotate(decrt, (-(90.0 - b[arc_n]) * M_PI/180.0), (-l[arc_n] * M_PI/180.0));		  
	  sphr = dec2sp(decrt);
	  
	  if ( END[arc_n] <= 360.0 && sphr.phi >= (BEG[arc_n] * M_PI/180.0) && sphr.phi <= (END[arc_n] * M_PI/180.0) )
	    {
	      	      
	      decrt = sp2dec(sphr);
	      decrt = rotate(decrt, 0.0, -PSI_pr);
	      /* printf("%f\t%f\t%f\n", decrt.x, decrt.y, decrt.z); */
	      
	      decrt = rotate(decrt, angle_JI, 0.0);

	      for (j = 0; j < Nj; j++)
		{      
		  for (k = 1; k < Nk; k++)
		    {
		      sp_sky.theta = delta[j][k];
		      sp_sky.phi = alpha[j][k];
		      sp_sky.r = 1.0;

		      p_sky = sp2dec(sp_sky);
		      
		      cos_p_decrt = dot(p_sky, decrt);
		      p_decrt = acos(cos_p_decrt);
	      
		      if (cos_p_decrt > 0.0)
			{
			  flux[j][k] = flux[j][k] + (1.0/((double) arc_len_count)) * intensity[arc_n] * exp((-1.0)*(p_decrt/(delta_theta[arc_n]*M_PI/180.0))*(p_decrt/(delta_theta[arc_n]*M_PI/180.0))) * cos_p_decrt * ds_sky[j][k];
			}
		      else if (cos_p_decrt < 0.0)
			{
			  flux[j][k] = flux[j][k] + 0.0;
			}
		    }
		}
	      
	      //printf("%.20f\t %.20f\t %.20f\t %.20f\n", decrt.x, decrt.y, decrt.z, intensity[arc_n]);
	      
		      
	      //sphr = dec2sp(decrt);
	      
	      //f = (int) floor( sphr.theta * 180.0/M_PI );
	      
	      //F[f] = F[f] + ds * intensity[arc_n];
	      
	    }
	  else if ( END[arc_n] > 360.0 && sphr.phi >= (BEG[arc_n] * M_PI/180.0) || sphr.phi <= (END[arc_n] * M_PI/180.0 - 2.0 * M_PI) )
	    {
	      
	      decrt = sp2dec(sphr);
	      decrt = rotate(decrt, 0.0, -PSI_pr);
	      decrt = rotate(decrt, angle_JI, 0.0);
	      
	      for (j = 0; j < Nj; j++)
		{      
		  for (k = 1; k < Nk; k++)
		    {
		      sp_sky.theta = delta[j][k];
		      sp_sky.phi = alpha[j][k];
		      sp_sky.r = 1.0;

		      p_sky = sp2dec(sp_sky);
		      
		      cos_p_decrt = dot(p_sky, decrt);
		      p_decrt = acos(cos_p_decrt);
	      
		      if (cos_p_decrt > 0.0)
			{
			  //flux[j][k] = flux[j][k] + (1.0/((double) arc_len_count)) * intensity[arc_n] * exp((-1.0)*(p_decrt/(delta_theta[arc_n]*M_PI/180.0))*(p_decrt/(delta_theta[arc_n]*M_PI/180.0))) * cos_p_decrt * ds_sky[j][k];

			  flux[j][k] = flux[j][k] + (1.0/((double) arc_len_count)) * intensity[arc_n] * exp((-1.0)*(p_decrt/(delta_theta[arc_n]*M_PI/180.0))*(p_decrt/(delta_theta[arc_n]*M_PI/180.0))) * cos_p_decrt * ds_sky[j][k];
			    
			}
		      else if (cos_p_decrt < 0.0)
			{
			  flux[j][k] = flux[j][k] + 0.0;
			}
		    }
		}
	      
	      //printf("%.20f\t %.20f\t %.20f\t %.20f\n", decrt.x, decrt.y, decrt.z, intensity[arc_n]);
	      
	      //sphr = dec2sp(decrt);
	      
	      //f = (int) floor( sphr.theta * 180.0/M_PI );
	      
	      //F[f] = F[f] + ds * intensity[arc_n];
	      
	    }
	  

	  //printf("%.20f\t %.20f\t %.20f\t %.20f\n", p_sky.x, p_sky.y, p_sky.z, flux_through_ds);
		  
	  
	  //printf("%.20f\t %.20f\t %.20f\t %.20f\n", decrt.x, decrt.y, decrt.z, intensity[arc_n]);
	  

		  
	}

      /* printf("\n\n"); */
      /* printf("0.00000\t 0.00000\t 1.00000\n0.00000\t 0.00000\t -1.00000\n"); */
      /* printf("0.00000\t 1.00000\t 0.00000\n0.00000\t -1.00000\t 0.00000\n"); */
      /* printf("1.00000\t 0.00000\t 0.00000\n-1.00000\t 0.00000\t 0.00000\n"); */
    }



  double full_sum = 0.0;
  
  for (j = 0; j < Nj; j++)
    {      
      for (k = 1; k < Nk; k++)
	{
	  //sp_sky.theta = delta[j][k];
	  //sp_sky.phi = alpha[j][k];
	  //sp_sky.r = 1.0;
	  
	  //flux_sum[k] = flux_sum[k] + flux[j][k];
	  
	  if (isnan(flux[j][k]))
	    {
	      flux[j][k] = 0.0;
	    }

	  flux_sum[k] = flux_sum[k] + flux[j][k];
	  
	  full_sum = full_sum + flux[j][k];
	  
	  /* printf("PHI THETA FLUX(j,k) full_sum \t%.20f\t %.20f\t %.20f\t %.20f\n", sp_sky.phi, sp_sky.theta, flux[j][k], full_sum); */
	  /* if (flux[j][k] == NAN) */
	  /*   { */
	  /*     printf("%.20f\t %d\t %d\n", flux[j][k], j, k); */
	  /*   } */
	}

    }

  
  /* for (j = 0; j < Nj; j++) */
  /*   {       */

  /*     printf("%d\t %.20f\n", j, flux[j][86]); */

  /*   } */









  
  /* for (k = 1; k < Nk; k++) */
  /*   { */
  /*     printf("%.20f\t %.20f\t %.20f\n", (PSI_pr * 180.0)/M_PI, ((double) k * M_PI)/((double) Nk), flux_mean[k]); */
  /*   } */


  
	  
  double * result = (double *) malloc(sizeof(double) * 180);
  
  double I_sum = 0.0;
  double f_angle;

  /* for (f = 0; f < 180; f++) */
  /*   { */

  /*     f_angle = (double) (f * M_PI/180.0 + 0.5 * M_PI/180.0); */
	
  /*     I[f] = flux_sum[f]/(2.0 * M_PI * sin(f_angle) * M_PI/180.0); /\* mean intensity in f-th zone *\/ */

  /*   } */
  
  /* for (f = 0; f < 180; f++) */
  /*   { */
  /*     I_sum = I_sum + I[f]; */
  /*     //printf("%.20f\n", F_sum); */
  /*   } */


  /* double test = 0.0; */
  
  /* for (f = 0; f < 180; f++) */
  /*   { */
  /*     result[f] = Lx * 4.0 * M_PI * I[f]/I_sum; */

  /*     test = test + result[f]; */
      
  /*     //printf("%d\t %f\t %f\n", f, result[f], I_sum); */
  /*   } */

  /* printf("%.10f\n", test); */

  

  //FILE *diagram;

  //diagram = fopen("DIAGRAM", "w");

  for (k = 1; k < 180; k++)
    {

      //result[k] = Lx * 4.0 * M_PI * flux_sum[k]/full_sum;
      
      result[k] = 4.0 * M_PI * flux_sum[k]/full_sum;
      
      /* printf("%.20f\t", result[k]); */
      
            
      //test for directions
      /* if (k <= 90) */
      /* 	{ */
      /* 	  result[k] = Lx * 4.0 * M_PI * 1.0/90.0;   */
      /* 	} */
      /* else */
      /* 	{ */
      /* 	  result[k] = 0.0; */
      /* 	} */
      
      //printf("%f\t %d\t %f\n", (PSI_pr * 180.0/M_PI), k, result[k]);
      
    }

  //fclose(diagram);

  /* printf("full_sum %.20f\t", full_sum); */

  
  result[0] = 0.0;
  

  //result[0] = 1.0/90.0;
  
  return result;

}
