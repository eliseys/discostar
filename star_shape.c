#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "def.h"

double fr(double r, double phi, double theta, double q, double omega)
{
  double l = sin(theta) * cos(phi);
  double n = cos(theta);
  return omega - 1.0/r - q * (1.0/sqrt(1.0 + r * r - 2 * r * l) - r * l) - 0.5 * (1.0 + q) * r * r * (1 - n * n);
}


double dfr(double r, double phi, double theta, double q, double omega)
{
  double l = sin(theta) * cos(phi);
  double n = cos(theta);
  return q * ((r - l)/((r * r - 2.0 * l * r + 1.0) * sqrt(r * r - 2.0 * l * r + 1.0)) + l) - (1.0 - n * n) * (q + 1.0) * r + 1.0/(r * r);
}


double radius_star(double phi, double theta, double q, double omega)
{
  double r0 = 0.001; /* zeroth approximation */
  double r1 = r0 - fr(r0, phi, theta, q, omega)/dfr(r0, phi, theta, q, omega);

  while (fabs(r1 - r0) > eps)
    {
      r0 = r1;
      r1 = r1 - fr(r1, phi, theta, q, omega)/dfr(r1, phi, theta, q, omega);
    }
  return r1;  
}

double * polar(double q, double omega)
{
  /* Computes polar_g_abs and polar_r of the star */

  double polar_r = radius_star(0.0, 0.0, q, omega);
  double max_r = radius_star(0.0, 0.5*M_PI, q, omega);
  
  double z = polar_r;
  
  double d_omega_polar_x = q/(sqrt(z*z + 1.0)*(z*z + 1.0)) - q;

  double d_omega_polar_y = 0.0;
  double d_omega_polar_z = - q*z/(sqrt(z*z + 1.0)*(z*z + 1.0)) - 1.0/(z*z);

  double polar_g_abs =
    sqrt(
	 d_omega_polar_x * d_omega_polar_x
	 + d_omega_polar_y * d_omega_polar_y
	 + d_omega_polar_z * d_omega_polar_z
	 );

  double *result = (double*) malloc(sizeof(double) * 3);

  result[0] = polar_g_abs;
  result[1] = polar_r;
  result[2] = max_r;
    
  return result;

}


double * gradient(double phi, double theta, double q, double omega)
{

  double d_omega_x, d_omega_y, d_omega_z;
  double g_x, g_y, g_z;
  double g_abs; 
  double n_x, n_y, n_z;
  
  double r = radius_star(phi, theta, q, omega);

  double x = r * sin(theta) * cos(phi);
  double y = r * sin(theta) * sin(phi);
  double z = r * cos(theta);

  
  d_omega_x =
    q * ((1.0 - x)/(sqrt(x*x + y*y + z*z - 2.0*x + 1.0)*(x*x + y*y + z*z - 2.0*x + 1.0)) - 1.0)
    - x/(sqrt(x*x + y*y + z*z)*(x*x + y*y + z*z))
    + (q + 1.0)*x;
  
  d_omega_y =
    - q*y/(sqrt(x*x + y*y + z*z - 2.0*x + 1.0)*(x*x + y*y + z*z - 2.0*x + 1.0))
    - y/(sqrt(x*x + y*y + z*z)*(x*x + y*y + z*z))
    + (q + 1.0)*y;
  
  d_omega_z =
    - q*z/(sqrt(x*x + y*y + z*z - 2.0*x + 1.0)*(x*x + y*y + z*z - 2.0*x + 1.0))
    - z/(sqrt(x*x + y*y + z*z)*(x*x + y*y + z*z));
    
  g_x = d_omega_x;
  g_y = d_omega_y;
  g_z = d_omega_z;
  
  g_abs = sqrt(g_x * g_x + g_y * g_y + g_z * g_z);

  n_x = - d_omega_x/g_abs;
  n_y = - d_omega_y/g_abs;
  n_z = - d_omega_z/g_abs;
  
  double *result = (double*) malloc(sizeof(double) * 4);

  result[0] = g_abs;
  result[1] = n_x;
  result[2] = n_y;
  result[3] = n_z;
  
  return result;
}


double fx(double x, double q)
{

  return
    pow(x,5) * (1 + q)
    + pow(x,4) * (-3*q - 2)
    + pow(x,3) * (3*q + 1)
    - x * x
    + x * 2.0
    - 1;

}

double dfx(double x, double q)
{
  
  return
    5.0 * pow(x,4) * (1 + q)
    + 4.0 * pow(x,3) * (-3*q - 2)
    + 3.0 * x * x * (3*q + 1)
    - 2.0 * x
    + 2.0;
    
}

double fomega(double r, double q, double omega_crit)
{
  return 1.0/r + q * (1.0/sqrt(1.0 + r * r)) - omega_crit;
  
}

double dfomega(double r, double q)
{
  return (-1.0) * q * r /(sqrt(r * r + 1.0) * (r * r + 1.0)) - 1.0 /(r * r);
}

double omg(double q, double mu)
{
  /* Returns dimensionless potential omega on the surface of the star */
  
  double x_crit;
  double r_polar_crit;
  double x0 = 0.001; /* zeroth approximation */
  double x1 = x0 - fx(x0, q)/dfx(x0, q);

  while (fabs(x1 - x0) > eps)
    {
      x0 = x1;
      x1 = x1 - fx(x1, q)/dfx(x1, q);
    }
  x_crit = x1;
    
  double omega_crit = 1.0/x_crit + q * (1.0/(1.0 - x_crit) - x_crit) + x_crit * x_crit * (1.0 + q) * 0.5;
  
  double r0 = 0.001; /* zeroth approximation */
  double r1 = r0 - fomega(r0, q, omega_crit)/dfomega(r0, q);
  
  while (fabs(r1 - r0) > eps)
    {
      r0 = r1;
      r1 = r1 - fomega(r1, q, omega_crit)/dfomega(r1, q);
    }
  r_polar_crit = r1;
  
  double r_polar = r_polar_crit * mu;

  double value = 1.0 / r_polar + q / sqrt(1.0 + r_polar * r_polar);

  return value;

}



double * phi_func(int steps_phi, int threads)
{
  int i;

  double * result = (double*) malloc(sizeof(double) * steps_phi);


  omp_set_dynamic(0);
  omp_set_num_threads(threads);

#pragma omp parallel for private(i)
  for (i = 0; i < steps_phi; i++)
    {
      result[i] = (double) i * 2.0 * M_PI/steps_phi + 0.5 * 2.0 * M_PI/steps_phi;
    }

  return result;
}


double * theta_func(int steps_theta, int threads)
{
  int j;

  double * result = (double*) malloc(sizeof(double) * steps_theta);

  omp_set_dynamic(0);
  omp_set_num_threads(threads);

#pragma omp parallel for private(j)
  for (j = 0; j < steps_theta; j++)
    {
      result[j] = (double) j * M_PI/steps_theta + 0.5 * M_PI/steps_theta;
    }
  
  return result;
}



double * shape_g_abs(int steps_phi, int steps_theta, double * phi_array, double * theta_array, double q, double omega, int threads)
{
  
  int N = steps_theta * steps_phi;
  double * result = (double*) malloc(sizeof(double) * N * 4);

  //double * grd;
  //double * phi_array;
  //double * theta_array;
  double * polar_array;

  //phi_array = phi_func(steps_phi, threads);
  //theta_array = theta_func(steps_theta, threads);

  polar_array = polar(q, omega);

  double polar_g_abs = polar_array[0];
  double g;
  double g_abs;


  int i, j;

  omp_set_dynamic(0);
  omp_set_num_threads(threads);

#pragma omp parallel for private(j, g_abs, g)
  for (i = 0; i < steps_phi; i++)
    {

      for (j = 0; j < steps_theta; j++)
	{ 
	  double * grd = gradient(phi_array[i], theta_array[j], q, omega);

	  g_abs = grd[0];

	  g = g_abs/polar_g_abs;

	  result[(steps_phi * j + i)*4 + 0] = g;
	  result[(steps_phi * j + i)*4 + 1] = grd[1];
	  result[(steps_phi * j + i)*4 + 2] = grd[2];
	  result[(steps_phi * j + i)*4 + 3] = grd[3];
	  free(grd);
	}

      
    }

  //free(grd);
  //free(phi_array);
  //free(theta_array);
  free(polar_array);
  
  return result;

}



double * shape_r(int steps_phi, int steps_theta, double * phi_array, double * theta_array, double q, double omega, int threads)
{
  /* star */

  int N = steps_theta * steps_phi;

  double * result = (double*) malloc(sizeof(double) * N);

  /* double * phi_array; */
  /* double * theta_array; */
  
  /* phi_array = phi_func(steps_phi, threads); */
  /* theta_array = theta_func(steps_theta, threads); */

  int i, j;
 
  omp_set_dynamic(0);
  omp_set_num_threads(threads);


#pragma omp parallel for private(j)  
  for (i = 0; i < steps_phi; i++)
    {
      for (j = 0; j < steps_theta; j++)
	{ 
	  result[steps_phi * j + i] = radius_star(phi_array[i], theta_array[j], q, omega);
	}
    }


  
  //free(phi_array);
  //free(theta_array);

  return result;
    
}




