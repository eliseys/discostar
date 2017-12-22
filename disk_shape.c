#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "def.h"


double radius_disk(disk disk, double phi, double theta)
{
  double r_out = disk.R;
  double h = 2.0 * len(disk.h);
  
  double r;
  
  if ((theta <= atan(2.0*r_out/h)) || (theta >= (M_PI - atan(2.0*r_out/h))))
    r = 0.5 * h/(sqrt(cos(theta) * cos(theta)));
  if ((theta > atan(2.0*r_out/h)) && (theta < (M_PI - atan(2.0*r_out/h))))
    {
      r = r_out / sin(theta);
    }
  
  return r;
}






double * phi_func_disk(int steps_phi_disk)
{
  int i;

  double * result = (double*) malloc(sizeof(double) * steps_phi_disk);
   
  for (i = 0; i < steps_phi_disk; i++)
    {
      result[i] = (double) i * 2.0 * M_PI/steps_phi_disk + 0.5 * 2.0 * M_PI/steps_phi_disk;
    }

  return result;
  
}


double * theta_func_disk(int disk_tiles)
{

  //int * result = (int *) malloc(sizeof(int));
  
  //result = 0;
  
  return 0;
  
  /* double h = len(disk.h); */
  /* double R = disk.R; */

  /* int steps = sqrt(disk_tiles/2.0); */
  
  /* int N = steps * R / (2.0 * (R + h)); */
  /* int M = steps * (R + 2.0 * h) / (2.0 * (R + h)); */
  /* int steps_theta = N + M - 1; */


  /* double delta_N = R/N; */
  /* double delta_M = 2.0 * h/(M - N); */

  /* double phi; */

  /* int j; */

  
  /* for (j = 0; j <= steps_theta; j++) */
  /*   { */

	  
  /* 	  if (j <= N - 1) */
  /* 	    { */
  /* 	      theta = atan( (j + 0.5) * (delta_N/h) ); */
  /* 	      theta_0 = atan(j*(delta_N/h)); */
  /* 	      theta_1 = atan((j + 1)*(delta_N/h)); */
  /* 	      delta_theta = theta_1 - theta_0; */
  /* 	      /\**\/ */
  /* 	      //cos_on = dot(o,n1); */
		  
  /* 	    } */
  /* 	  else if (j >= N && j <= M - 1) */
  /* 	    { */
  /* 	      delta = delta_M * (j - N) + 0.5 * delta_M; */
  /* 	      delta_0 = delta_M * (j - N); */
  /* 	      delta_1 = delta_M * (j + 1 - N); */
  /* 	      theta = atan(R/h) + acos((h * h + R * R - h * delta)/sqrt((h * h + R * R)*(R * R + (h - delta)*(h - delta)))); */
  /* 	      theta_0 = atan(R/h) + acos((h * h + R * R - h * delta_0)/sqrt((h * h + R * R)*(R * R + (h - delta_0)*(h - delta_0)))); */
  /* 	      theta_1 = atan(R/h) + acos((h * h + R * R - h * delta_1)/sqrt((h * h + R * R)*(R * R + (h - delta_1)*(h - delta_1)))); */
  /* 	      delta_theta = theta_1 - theta_0; */
	      
  /* 	      /\* disc`s normal vector for side *\/ */
  /* 	      n3.x = sin(phi) * sin(z_tilt - phi_orb + M_PI) + cos(phi) * cos(y_tilt) * cos(z_tilt - phi_orb + M_PI); */
  /* 	      n3.y = sin(phi) * cos(z_tilt - phi_orb + M_PI) - cos(phi) * cos(y_tilt) * sin(z_tilt - phi_orb + M_PI); */
  /* 	      n3.z = cos(phi) * sin(y_tilt); */
  /* 	      /\**\/ */
  /* 	      //cos_on = dot(o,n3); */

  /* 	    } */
  /* 	  else if (j >= M && j <= steps_theta) */
  /* 	    { */
  /* 	      delta = delta_N * (j - M) + 0.5 * delta_N; */
  /* 	      delta_0 = delta_N * (j - M); */
  /* 	      delta_1 = delta_N * (j + 1 - M); */
  /* 	      theta = 0.5 * M_PI + atan(h/R) + acos((h * h + R * R - R * delta)/sqrt((h * h + R * R)*(h * h + (R - delta)*(R - delta)))); */
  /* 	      theta_0 = 0.5 * M_PI + atan(h/R) + acos((h * h + R * R - R * delta_0)/sqrt((h * h + R * R)*(h * h + (R - delta_0)*(R - delta_0)))); */
  /* 	      theta_1 = 0.5 * M_PI + atan(h/R) + acos((h * h + R * R - R * delta_1)/sqrt((h * h + R * R)*(h * h + (R - delta_1)*(R - delta_1)))); */
  /* 	      delta_theta = theta_1 - theta_0; */
  /* 	      /\**\/ */
  /* 	      //cos_on = dot(o,n2); */
  /* 	    } */


}
