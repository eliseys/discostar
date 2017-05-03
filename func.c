#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "def.h"

double len(vec3 p)
{
  return sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
}

double dot(vec3 a, vec3 b)
{
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

double *coordinate_transformation(double x, double y, double z)
{
  /* transforms cartesian coordinates to spherical */

  double phi;
  double theta;
  double r;

  /* r [0,inf) */
  r = sqrt(x * x + y * y + z * z);

  /* theta [0, M_PI] */
  if (x != 0.0 || y != 0.0 || z != 0.0)
    theta = acos(z/sqrt(x * x + y * y + z * z));

  else if (x == 0.0 && y == 0.0 && z == 0.0)
    theta = 0.5 * M_PI;
  
  /* phi [0, 2*M_PI) */  
  if (x > 0.0 && y > 0.0)
    phi = asin((y/x)/sqrt(1 + (y/x)*(y/x)));
  else if (x > 0.0 && y < 0.0)
    phi = 2.0 * M_PI - asin(((-y)/x)/sqrt(1 + ((-y)/x)*((-y)/x)));
  else if (x < 0.0 && y > 0.0)
    phi = M_PI - asin((y/(-x))/sqrt(1 + (y/(-x))*(y/(-x))));
  else if (x < 0.0 && y < 0.0)
    phi = M_PI + asin(((-y)/(-x))/sqrt(1 + ((-y)/(-x))*((-y)/(-x))));
  else if (x == 0.0 && y > 0.0)
    phi = 0.5 * M_PI;
  else if (x == 0.0 && y < 0.0)
    phi = 1.5 * M_PI;
  else if (x >= 0.0 && y == 0.0)
    phi = 0.0;
  else if (x < 0.0 && y == 0.0)
    phi = M_PI;

  double *result = malloc(sizeof(double) * 3);

  result[0] = r;
  result[1] = phi;
  result[2] = theta;
  
  return result;
}

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

  double *result = malloc(sizeof(double) * 3);

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
  
  double *result = malloc(sizeof(double) * 4);

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

double distance_to_star(vec3 p, double omega, double q)
{

  double *coord;
  coord = coordinate_transformation(p.x, p.y, p.z);

  double r = coord[0];
  double phi = coord[1];
  double theta = coord[2];

  free(coord);
    
  double d = r - radius_star(phi, theta, q, omega);
  /* double d = r - 0.5; */

  return d;
}


double eclipse_by_star(double omega, double q, vec3 o, vec3 p)
{
  
  /* delta for calculating derivation (now equal to eps) */
  vec3 d;
  d.x = eps * o.x;
  d.y = eps * o.y;
  d.z = eps * o.z;

  /* zeroth approximation */
  double shift0 = distance_to_star(p, omega, q);
  vec3 position;
  position.x = p.x + o.x * shift0;
  position.y = p.y + o.y * shift0;
  position.z = p.z + o.z * shift0;   
  vec3 position_p;
  position_p.x = position.x + d.x * 0.5;
  position_p.y = position.y + d.y * 0.5;
  position_p.z = position.z + d.z * 0.5;
  vec3 position_m;
  position_m.x = position.x - d.x * 0.5;
  position_m.y = position.y - d.y * 0.5;
  position_m.z = position.z - d.z * 0.5;
  double plus_delta = distance_to_star(position_p, omega, q);
  double mins_delta = distance_to_star(position_m, omega, q);
  double derivative = (plus_delta - mins_delta)/len(d);
      
  double shift1 = shift0 - distance_to_star(position, omega, q)/derivative;

  /* double shift1 = shift0 + distance_to_star(position, omega, q); */
  
  int i = 0;
  while (fabs(shift1 - shift0) > eps)
    {
      shift0 = shift1;

      position.x = p.x + o.x * shift0;
      position.y = p.y + o.y * shift0;
      position.z = p.z + o.z * shift0;
   
      position_p.x = position.x + d.x * 0.5;
      position_p.y = position.y + d.y * 0.5;
      position_p.z = position.z + d.z * 0.5;
      position_m.x = position.x - d.x * 0.5;
      position_m.y = position.y - d.y * 0.5;
      position_m.z = position.z - d.z * 0.5;
      plus_delta = distance_to_star(position_p, omega, q);
      mins_delta = distance_to_star(position_m, omega, q);
      derivative = (plus_delta - mins_delta)/len(d);
      
      shift1 = shift0 - distance_to_star(position, omega, q)/derivative;
      
      /* shift1 = shift0 + distance_to_star(position, omega, q); */
      i++;
      
      if (i > 10)
	{
	  break;
	}
    }
    
  double result = distance_to_star(position, omega, q);

  vec3 shift; 
  shift.x = position.x - p.x;
  shift.y = position.y - p.y;
  shift.z = position.z - p.z;

  double direction_test = dot(shift, o);

  double ray;
  
  if (result < eps && direction_test < 0.0)
    ray = 1.0;
  else if (result < eps && direction_test > 0.0 )
    ray = 0.0;
  else if (result > eps)
    ray = 1.0;

  return ray;

}

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


double distance_to_disk(vec3 p, disk disk)
{
  double R = disk.R;
  double h = len(disk.h);
  double r = len(p);
  double hr = dot(disk.h,p);
  double r1 = sqrt(fabs(r * r - (hr * hr)/(h * h)));
  double r2 = fabs(hr/h); 

  double d;

  if (r2 >= h && r1 < R || r2 < h && r1 < R - r2)
    {
      d = r2 - h;
    }
  else if (r2 < h && r1 >= R - r2)
    {
      d = r1 - R;
    }
  else if (r2 >= h && r1 >= R)
    {
      d = sqrt((r1 - R) * (r1 - R) + (r2 - h) * (r2 - h));
    }

  return d;
}


double eclipse_by_disk(disk disk, vec3 o, vec3 p)
{
    
  /* delta for calculating derivation (now equal to eps) */
  vec3 d;
  d.x = eps * o.x;
  d.y = eps * o.y;
  d.z = eps * o.z;

  /* zeroth approximation */
  double shift0 = distance_to_disk(p, disk);
  vec3 position;
  position.x = p.x + o.x * shift0;
  position.y = p.y + o.y * shift0;
  position.z = p.z + o.z * shift0;
   
  vec3 position_p;
  position_p.x = position.x + d.x * 0.5;
  position_p.y = position.y + d.y * 0.5;
  position_p.z = position.z + d.z * 0.5;  
  vec3 position_m;	           
  position_m.x = position.x - d.x * 0.5;
  position_m.y = position.y - d.y * 0.5;
  position_m.z = position.z - d.z * 0.5;
  double plus_delta = distance_to_disk(position_p, disk);
  double mins_delta = distance_to_disk(position_m, disk);
  double derivative = (plus_delta - mins_delta)/len(d);
      
  double shift1 = shift0 - distance_to_disk(position, disk)/derivative;

  
  int i = 0;
  while (fabs(shift1 - shift0) > eps)
    {
      shift0 = shift1;

      position.x = p.x + o.x * shift0;
      position.y = p.y + o.y * shift0;
      position.z = p.z + o.z * shift0;
   
      position_p.x = position.x + d.x * 0.5;
      position_p.y = position.y + d.y * 0.5;
      position_p.z = position.z + d.z * 0.5;  
      position_m.x = position.x - d.x * 0.5;
      position_m.y = position.y - d.y * 0.5;
      position_m.z = position.z - d.z * 0.5;
      plus_delta = distance_to_disk(position_p, disk);
      mins_delta = distance_to_disk(position_m, disk);
      derivative = (plus_delta - mins_delta)/len(d);
      
      shift1 = shift0 - distance_to_disk(position, disk)/derivative;

      i++;

      if (i > 10)
	{
	  break;
	} 
    }
  
  double result = distance_to_disk(position, disk);

  /* test if the direction of the ray opposite to the direction of trace */

  vec3 shift;
  shift.x = position.x - p.x;
  shift.y = position.y - p.y;
  shift.z = position.z - p.z;

  double direction_test = dot(shift, o);

  double ray;
  
  if (result < eps && direction_test < 0.0)
    ray = 1.0;
  else if (result < eps && direction_test > 0.0 )
    ray = 0.0;
  else if (result > eps)
    ray = 1.0;

  return ray;
}



double flux_star(vec3 o, double q, double omega, double beta, double u, disk disk, double Lx, double albedo, int tiles, double T, double lambda, double a)
{
  /* */
  int steps = sqrt(tiles/2.0);
  int steps_phi = 2 * steps;
  int steps_theta = steps;

  double delta_phi = 2.0 * M_PI / steps_phi;
  double delta_theta = M_PI / steps_theta;
    
  /* polar values */
  double *plr;
  plr = polar(q, omega);
  double polar_g_abs = plr[0];
  double polar_r = plr[1];
  double max_r = plr[2]; 
  free(plr);
  
  /* */
  double r;
  double pl;
  vec3 p;
  vec3 n;
  vec3 pn;
  vec3 ps;
  vec3 psn;
  double lps;
  
  double phi, theta;

  /* gradient omega */
  double g_abs, g;
    
  /* dot products */
  double cos_on, cos_rn;
  
  /* overlapping */
  double ray;

  /* irradiation */
  double cos_irr;
  double cos_in;
  double h = len(disk.h);
  double R = disk.R;
  double disk_shadow_semi_angle = atan(h/R);
  double cos_disk_shadow_semi_angle = cos(0.5 * M_PI - disk_shadow_semi_angle);
  double S;
  double Fx;
  double T_star_4 = pow(T,4);
  double T_irr_4;
  double T_sum;
  double F_0;

  
  double color;
  
  /* unity vectior along disk.h */
  vec3 hn;
  hn.x = disk.h.x/h;
  hn.y = disk.h.y/h;
  hn.z = disk.h.z/h;
  
  /* flux from star */
  double result = 0.0;

  /* */
  /* double F_0 = F_lambda(T,lambda); */
  
  int i, j;
  
  for (i = 0; i < steps_phi; i++)
    {

      phi = (double) i * 2.0 * M_PI/steps_phi + 0.5 * 2.0 * M_PI/steps_phi;

      for (j = 0; j < steps_theta; j++)
	{

	  theta = (double) j * M_PI/steps_theta + 0.5 * M_PI/steps_theta;
	    
	  /* gradient omega */
	  double *grd;
	  grd = gradient(phi, theta, q, omega);
	  g_abs = grd[0];
	  g = g_abs/polar_g_abs;
	  /* surface normal vector */	  
	  n.x = grd[1];
	  n.y = grd[2];
	  n.z = grd[3];
	  free(grd);

	  /* star`s dot products */
	  cos_on = dot(o,n);

	  if (cos_on < - eps)
	    {
	      continue;
	    }

	  /* star */
	  r = radius_star(phi, theta, q, omega);
	  p.x = r * sin(theta) * cos(phi);
	  p.y = r * sin(theta) * sin(phi);
	  p.z = r * cos(theta);

	  /* shifted points */
	  ps.x = p.x - 1.0;
	  ps.y = p.y;
	  ps.z = p.z;

	  /* ray = eclipse_by_disk(disk, o, ps); */
	  if (o.x < cos(R + max_r))
	    {
	      ray = 1.0;
	    }
	  else if (o.x > cos(R + max_r))
	    {
	      ray = eclipse_by_disk(disk, o, ps);
	    }
	  
	  /* normalized vector */
	  pl = len(p);
	  pn.x = p.x/pl;
	  pn.y = p.y/pl;
	  pn.z = p.z/pl;
	  
	  cos_rn = dot(pn,n);
	  /* if (cos_rn < 0.0) */
	  /*   { */
	  /*     printf("%f\n", cos_rn); */
	  /*     printf("%f %f %f\t %f %f %f\n", pn.x, pn.y, pn.z, n.x, n.y, n.z); */
		      
	  /*   } */


	  lps = len(ps); /* distance from the secondary to the point p */
	  psn.x = ps.x/lps;
	  psn.y = ps.y/lps;
	  psn.z = ps.z/lps;

	  /* Irradiation */	  
	  /* */
	  cos_irr = dot(psn, hn);
	  cos_in  = dot(psn, n);
	  
	  /* Surface element */
	  S = a * a * r * r * sin(theta) * delta_phi * delta_theta / cos_rn;

	  /* printf("%f\n", cos_rn); */

	  /* X-ray flux incident of the surface element */
	  if ( cos_in < - eps && fabs(cos_irr) > cos_disk_shadow_semi_angle)
	    {
	      /* Fx = Lx * S * fabs(cos_in) / (4.0 * M_PI * lps * lps); */
	      
	      Fx = (1.0 - albedo) * Lx * fabs(cos_in) / (4.0 * M_PI * lps * lps * a * a);

	      /* printf("%f\n", Fx); */
	    }
	  else
	    {
	      Fx = 0.0;
	    }
	  
	  if (ray == 1.0 && cos_on > 0.0 + eps)
	    {

	      T_irr_4 = Fx / SIGMA;
	      
	      T_sum = pow((T_star_4 + T_irr_4),0.25); 

	      F_0 = F_lambda(T_sum, lambda);
	      
	      result = result + F_0 * (1 - u + u * cos_on) * pow(g,beta) * cos_on * S;

	      /* result = result + F_0 * (Fx + (1 - u + u * cos_on) * pow(g,beta)) * cos_on * S; */
	    
	      /* color = (Fx * (1.0 - albedo) + (1 - u + u * cos_on) * pow(g,beta)) * cos_on * S; */
	      /* printf("%f\t %f\t %f\t %f\n", p.x, p.y, p.z, color); */

	    }

	}
      
    }      

  return result;
  
}


double flux_disk(vec3 o, disk disk, double y_tilt, double z_tilt, double omega, double q, double b, int disk_tiles, double phi_orb, double T, double lambda, double a)
{

  double R = disk.R;
  double h = len(disk.h); /* semithickness of the disk */

  /* */
  double theta, phi;

  int steps = sqrt(disk_tiles/2.0);
  
  int steps_phi = 2 * steps;

  int N = steps * R / (2.0 * (R + h));  

  int M = steps * (R + 2.0 * h) / (2.0 * (R + h));

  int steps_theta = N + M - 1;
  
  double delta_phi = 2.0 * M_PI / steps_phi;
  double delta_theta; /* calculated on each j step */
  
  double delta_N = R/N;
  
  double delta_M = 2.0 * h/(M - N);

  double delta;
  double delta_0, delta_1;
  double theta_0, theta_1;
  
  double r;
  vec3 p;
  vec3 pt;
  vec3 pn;
  vec3 ps;

  /* unity normal vectors to top and bottom surface of the disk */
  vec3 n1;
  n1.x = disk.h.x/h;
  n1.y = disk.h.y/h;
  n1.z = disk.h.z/h;

  vec3 n2;
  n2.x = - disk.h.x/h;
  n2.y = - disk.h.y/h;
  n2.z = - disk.h.z/h; 

  vec3 n3;

  /* printf("%f %f %f\n", n1.x, n1.y, n1.z); */

  vec3 nu, nd, ns;
  
  double lp;
  
  /* dot products */
  double cos_rn_u;
  double cos_rn_d;
  double cos_rn_s;
    
  /* overlapping */
  double ray;

  /* */
  double S;
  
  /* brightness coefficient of the disk */
  double B;
  double rho;
    
  double color;

  /* */
  double F_0 = F_lambda(T, lambda);
  
  /* from top side of the disk */
  double result_1 = 0.0;
  /* from ridge of the disk */
  double result_2 = 0.0;
  /* from bottom side of the disk */
  double result_3 = 0.0;

  double cos_on;

  double *plr;
  plr = polar(q, omega);
  double max_r = plr[2]; 
  free(plr);

  int i, j;

  for (i = 0; i < steps_phi; i++)
    {
      phi = (double) i * 2.0 * M_PI/steps_phi + 0.5 * 2.0 * M_PI/steps_phi;

      for (j = 0; j <= steps_theta; j++)
	{	  

	  /* theta = (double) j * M_PI / steps_theta + 0.5 * M_PI/steps_theta; */
	  
	  if (j <= N - 1)
	    {
	      theta = atan( (j + 0.5) * (delta_N/h) );
	      theta_0 = atan(j*(delta_N/h));
	      theta_1 = atan((j + 1)*(delta_N/h));
	      delta_theta = theta_1 - theta_0;
	      /**/
	      cos_on = dot(o,n1);
		  
	    }
	  else if (j >= N && j <= M - 1)
	    {
	      delta = delta_M * (j - N) + 0.5 * delta_M;
	      delta_0 = delta_M * (j - N);
	      delta_1 = delta_M * (j + 1 - N);
	      theta = atan(R/h) + acos((h * h + R * R - h * delta)/sqrt((h * h + R * R)*(R * R + (h - delta)*(h - delta))));
	      theta_0 = atan(R/h) + acos((h * h + R * R - h * delta_0)/sqrt((h * h + R * R)*(R * R + (h - delta_0)*(h - delta_0))));
	      theta_1 = atan(R/h) + acos((h * h + R * R - h * delta_1)/sqrt((h * h + R * R)*(R * R + (h - delta_1)*(h - delta_1))));
	      delta_theta = theta_1 - theta_0;
	      
	      /* disc`s normal vector for side */
	      n3.x = sin(phi) * sin(z_tilt - phi_orb + M_PI) + cos(phi) * cos(y_tilt) * cos(z_tilt - phi_orb + M_PI);
	      n3.y = sin(phi) * cos(z_tilt - phi_orb + M_PI) - cos(phi) * cos(y_tilt) * sin(z_tilt - phi_orb + M_PI);
	      n3.z = cos(phi) * sin(y_tilt);
	      /**/
	      cos_on = dot(o,n3);

	    }
	  else if (j >= M && j <= steps_theta)
	    {
	      delta = delta_N * (j - M) + 0.5 * delta_N;
	      delta_0 = delta_N * (j - M);
	      delta_1 = delta_N * (j + 1 - M);
	      theta = 0.5 * M_PI + atan(h/R) + acos((h * h + R * R - R * delta)/sqrt((h * h + R * R)*(h * h + (R - delta)*(R - delta))));
	      theta_0 = 0.5 * M_PI + atan(h/R) + acos((h * h + R * R - R * delta_0)/sqrt((h * h + R * R)*(h * h + (R - delta_0)*(R - delta_0))));
	      theta_1 = 0.5 * M_PI + atan(h/R) + acos((h * h + R * R - R * delta_1)/sqrt((h * h + R * R)*(h * h + (R - delta_1)*(R - delta_1))));
	      delta_theta = theta_1 - theta_0;
	      /**/
	      cos_on = dot(o,n2);
	      
	    }

	  if (cos_on < -eps)
	    {
	      continue;
	    }

	  
	  r = radius_disk(disk, phi, theta);

	  /* cartesian coordinates of the point */
	  p.x = r * sin(theta) * cos(phi);
	  p.y = r * sin(theta) * sin(phi);
	  p.z = r * cos(theta);


	  /* brightness profile of the disk */
	  /* rho = r * sin(theta); */
	  /* B = b * exp(-(rho*rho)/(1.0 * R * 1.0 * R)); */
	  B = b;
	  
	  /* shifted coordinates of the disk */
	  ps.x = p.x + 1.0;
	  ps.y = p.y;
	  ps.z = p.z;
	  
	  /* normalized radius vector */
	  lp = len(p);
	  pn.x = p.x/lp;
	  pn.y = p.y/lp;
	  pn.z = p.z/lp;

	  /* surface normal vectors */
	  /* top of the disc */
	  nu.x = 0.0;
	  nu.y = 0.0;
	  nu.z = 1.0;
	  
	  /* bottom of the disc */
	  nd.x = 0.0;
	  nd.y = 0.0;
	  nd.z = -1.0;

	  /* side of the disc */
	  ns.x = cos(phi);
	  ns.y = sin(phi);
	  ns.z = 0.0;

	  cos_rn_u = dot(pn,nu);
	  cos_rn_d = dot(pn,nd);
	  cos_rn_s = dot(pn,ns);

	  /* printf("%f %f %f\n", len(n1), len(n2), len(n3)); */
	  
	  /* cos_on_u = dot(o,n1); */
	  /* cos_on_d = dot(o,n2);  */
	  /* cos_on_s = dot(o,n3); */

	  /* printf("%f %f %f\n", o.x, o.y, o.z); */

	  /* sphere surface element */
	  S = a * a * r * r * sin(theta) * delta_phi * delta_theta;
	  	  
	  /* tilt and shift disc */
	  pt.x =   p.x * cos(y_tilt) * cos(z_tilt - phi_orb + M_PI) + p.y * sin(z_tilt - phi_orb + M_PI) - p.z * sin(y_tilt) * cos(z_tilt - phi_orb + M_PI) + 1.0;
	  pt.y = - p.x * cos(y_tilt) * sin(z_tilt - phi_orb + M_PI) + p.y * cos(z_tilt - phi_orb + M_PI) + p.z * sin(y_tilt) * sin(z_tilt - phi_orb + M_PI);
	  pt.z =   p.x * sin(y_tilt) + p.z * cos(y_tilt);
	  
	  /* ray = eclipse_by_star(omega, q, o, pt); */
	  if (o.x < - cos(R + max_r))
	    {
	      ray = eclipse_by_star(omega, q, o, pt);
	    }
	  else if (o.x > - cos(R + max_r))
	    {
	      ray = 1.0;
	    }
	  
	  /* ray = 1.0; */
	  /* printf("%d\n", j); */
	  
	  if ( j <= N - 1 && ray == 1.0 )
	    {
	      /* top surface */
	      /* printf("%f\t %f\t %f\n", pt.x, pt.y, pt.z); */

	      result_1 = result_1 + (F_0 * cos_on * S)/cos_rn_u;
	      /* result_1 = result_1 + cos_rn_u; */

	      /* color = B * cos_on_u * S / cos_rn_u ; */
	      /* printf("%f\t %f\t %f\t %f\n", pt.x, pt.y, pt.z, cos_on_u); */
	      /* printf("%f\n", color); */

	    }
	  else if ( j >= M && j <= steps_theta && ray == 1.0 )
	    {
	      /* bottom surface */
	      /* printf("%f\t %f\t %f\n", pt.x, pt.y, pt.z); */
	      
	      result_2 = result_2 + (F_0 * cos_on * S)/cos_rn_d;
	      
	      /* result_2 = result_2 + cos_rn_d; */

	      /* color = B * cos_on_d * S / cos_rn_d; */
	      /* printf("%f\t %f\t %f\t %f\n", pt.x, pt.y, pt.z, cos_on_d); */
	      /* printf("%f\n", color); */

	    }
	  else if ( j >= N && j <= M - 1 && ray == 1.0 )
	    {
	      /* side surface */
	      /* printf("%f\t %f\t %f\n", pt.x, pt.y, pt.z); */
	    	  
	      result_3 = result_3 + (F_0 * cos_on * S)/cos_rn_s;
	      /* result_3 = result_3 + cos_rn_s; */

	      /* color = 1.0 * cos_on_s * S / cos_rn_s; */
	      /* printf("%f\t %f\t %f\t %f\n", pt.x, pt.y, pt.z, cos_on_s); */
	      /* printf("%f\n", color); */
	    }



	  /* printf("%d\n", j); */
	}
      
    }

  
  /* double *coord; */
  /* coord = coordinate_transformation(o.x, o.y, o.z); */

  /* double r_o = coord[0]; */
  /* double phi_o = coord[1]; */
  /* double theta_o = coord[2]; */

  /* free(coord); */

  /* printf("%f %f %f %f %f\n", phi_o/(2.0*M_PI), result_1, result_2, result_3, result_1 + result_2 + result_3); */

  return result_1 + result_2 + result_3;
}



