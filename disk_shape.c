#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "def.h"


double * disk_geometry(disk disk, parameters parameters)
{
  int N_r = parameters.N_r;
  
  int N_phi[N_r];
    
  double delta_r = disk.R/N_r;

  int N_z = (int) round(disk.h/delta_r);

  double delta_z = (double) disk.h/N_z;

  double r;

  double delta_phi;
  double phi;
  double phi_random_shift;

  vec3 p0, p1;
  vec3 n0, n1;
  double nr;
  double s;

  int N = 0;
  for (int i = 0; i < N_r; i++)
    {
      N_phi[i] = (int) round(2.0 * M_PI * (i + 0.5));
      N = N + N_phi[i];
    }

  N = 2 * N + 2 * N_phi[N_r-1] * N_z ; 

  /* printf(">>>> %d\n", N); */
    
  double * output = (double * ) malloc(sizeof(double) * N * 7 + 1);

  output[0] = N;
  
  int k = 0;
  for (int i = 0; i < N_r; i++)
    {
      r = ((double) i + 0.5) * delta_r;
      delta_phi = (double) (2.0 * M_PI)/N_phi[i];
      phi_random_shift = ((double) rand()/RAND_MAX) * delta_phi;
      
      for (int j = 0; j < N_phi[i]; j++)
  	{

  	  phi = (double) j * delta_phi + phi_random_shift;

  	  p0.x = r * cos(phi);
  	  p0.y = r * sin(phi);
  	  p0.z = disk.h * pow(((r - disk.r_in)/(disk.R-disk.r_in)), disk.gamma);

	  p1.x = p0.x;
  	  p1.y = p0.y;
  	  p1.z = - p0.z;
	  
	  nr = - sin(atan(disk.h * disk.gamma * pow(((r - disk.r_in)/(disk.R - disk.r_in)), disk.gamma)/(r - disk.r_in)));
	  
	  n0.x = nr * cos(phi);
	  n0.y = nr * sin(phi);
	  n0.z = sqrt(1.0 - pow(nr, 2));
	  
	  n1.x = n0.x;
	  n1.y = n0.y;
	  n1.z = - n0.z;

	  s = delta_phi * (2.0 * r * delta_r/n0.z + pow(delta_r/n0.z,2))/2.0;

	  p0 = R_y(p0, disk.theta_out);
	  p0 = R_z(p0, disk.phi_out);
	  n0 = R_y(n0, disk.theta_out);
	  n0 = R_z(n0, disk.phi_out);
	  p1 = R_y(p1, disk.theta_out);
	  p1 = R_z(p1, disk.phi_out);
	  n1 = R_y(n1, disk.theta_out);
	  n1 = R_z(n1, disk.phi_out);

	  output[14*k + 1] = p0.x;
	  output[14*k + 2] = p0.y;
	  output[14*k + 3] = p0.z;
	  output[14*k + 4] = n0.x;
	  output[14*k + 5] = n0.y;
	  output[14*k + 6] = n0.z;
	  output[14*k + 7] = s;
	  output[14*k + 8] = p1.x;
	  output[14*k + 9] = p1.y;
	  output[14*k + 10] = p1.z;
	  output[14*k + 11] = n1.x;
	  output[14*k + 12] = n1.y;
	  output[14*k + 13] = n1.z;
	  output[14*k + 14] = s;
	  k++;
	  /* printf(">>> %d\n", k); */
  	  /* printf("%f\t%f\t%f\n", p0.x, p0.y, p0.z); */
	  /* printf("%f\t%f\t%f\n", p1.x, p1.y, p1.z); */

  	}
    }


  
  for (int i = 0; i < N_z; i++)
    {
      delta_phi = (double) (2.0 * M_PI)/N_phi[N_r-1];
      
      phi_random_shift = ((double) rand()/RAND_MAX) * delta_phi;
      
      for (int j = 0; j < N_phi[N_r-1]; j++)
	{

  	  phi = (double) j * delta_phi + phi_random_shift;

	  p0.x = disk.R * cos(phi);
	  p0.y = disk.R * sin(phi);
	  p0.z = ((double) i + 0.5) * delta_z;

	  p1.x = p0.x;
	  p1.y = p0.y;
	  p1.z = - p0.z;
	  
	  n0.x = cos(phi);
	  n0.y = sin(phi);
	  n0.z = 0.0;

	  n1.x = n0.x;
	  n1.y = n0.y;
	  n1.z = - n0.z;

	  s = disk.R * delta_phi * delta_z;
	    
	  p0 = R_y(p0, disk.theta_out);
	  p0 = R_z(p0, disk.phi_out);
	  n0 = R_y(n0, disk.theta_out);
	  n0 = R_z(n0, disk.phi_out);
	  p1 = R_y(p1, disk.theta_out);
	  p1 = R_z(p1, disk.phi_out);
	  n1 = R_y(n1, disk.theta_out);
	  n1 = R_z(n1, disk.phi_out);

	  output[14*k + 1] = p0.x;
	  output[14*k + 2] = p0.y;
	  output[14*k + 3] = p0.z;
	  output[14*k + 4] = n0.x;
	  output[14*k + 5] = n0.y;
	  output[14*k + 6] = n0.z;
	  output[14*k + 7] = s;
	  output[14*k + 8] = p1.x;
	  output[14*k + 9] = p1.y;
	  output[14*k + 10] = p1.z;
	  output[14*k + 11] = n1.x;
	  output[14*k + 12] = n1.y;
	  output[14*k + 13] = n1.z;
	  output[14*k + 14] = s;
	  k++;
	  /* printf(">>> %d\n", k); */

	  /* printf("%f\t%f\t%f\n", p0.x, p0.y, p0.z); */
	  /* printf("%f\t%f\t%f\n", p1.x, p1.y, p1.z); */


	}

    }

  return output;

}












double radius_disk(disk disk, double phi, double theta)
{
  double r_out = disk.R;
  double h = 2.0 * disk.h;
  
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
