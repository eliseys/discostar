#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "def.h"

double distance_to_star(vec3 p, double omega, double q)
{
  sp d = dec2sp(p);
    
  return d.r - radius_star(d.phi, d.theta, q, omega);
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


double distance_to_disk_inside(vec3 p, disk disk)
{
  // distance function for rays outgoing from disk 

  double R = disk.R;
  double h = len(disk.h);
  double r = len(p);
  double hr = dot(disk.h,p);
  double r1 = sqrt(fabs(r * r - (hr * hr)/(h * h)));
  double r2 = fabs(hr/h); 

  double d;
  
  if (r2 < h)
    {
      d = - (r1 - R);
    }
  else if (r2 >= h)
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

double eclipse_by_disk_inside(disk disk, vec3 o, vec3 p)
{
    
  /* delta for calculating derivation (now equal to eps) */
  vec3 d;
  d.x = eps * o.x;
  d.y = eps * o.y;
  d.z = eps * o.z;

  /* zeroth approximation */
  double shift0 = distance_to_disk_inside(p, disk);
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
  double plus_delta = distance_to_disk_inside(position_p, disk);
  double mins_delta = distance_to_disk_inside(position_m, disk);
  double derivative = (plus_delta - mins_delta)/len(d);
      
  double shift1 = shift0 - distance_to_disk_inside(position, disk)/derivative;

  
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
      plus_delta = distance_to_disk_inside(position_p, disk);
      mins_delta = distance_to_disk_inside(position_m, disk);
      derivative = (plus_delta - mins_delta)/len(d);
      
      shift1 = shift0 - distance_to_disk_inside(position, disk)/derivative;

      i++;

      if (i > 10)
	{
	  break;
	} 
    }
  
  double result = distance_to_disk_inside(position, disk);

  /* test if the direction of the ray opposite to the direction of trace */

  vec3 shift;
  shift.x = position.x - p.x;
  shift.y = position.y - p.y;
  shift.z = position.z - p.z;

  double direction_test = dot(shift, o);

  double ray;
  
  if (result < eps && direction_test < 0.0 )
    ray = 1.0;
  else if (result < eps && direction_test > 0.0 )
    ray = 0.0;
  else if (result > eps )
    ray = 1.0;


  return ray;
}



double flux_star(vec3 o, double q, double omega, double beta, double u, disk disk, vec3 d2, double Lx, double Lx_disk, double Lx_iso, double Lx_disk_2, double albedo, int tiles, double T, double lambda, double a, vec3 neutron_star, double PSI_pr, int picture, int isotrope, sp disk_reflection_diagr, double * r_array, double * g_array, double * phi_array, double * theta_array)
{
  /* */
  int steps = sqrt(tiles/2.0);
  int steps_phi = 2 * steps;
  int steps_theta = steps;

  double delta_phi = 2.0 * M_PI / steps_phi;
  double delta_theta = M_PI / steps_theta;

  sp observer_sph = dec2sp(o);
  double phase_orb = observer_sph.phi;
  
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
  double cos_irr2;
  double cos_irr_min;
  double cos_irr2_min;
  double cos_in;
  double h = len(disk.h);
  double R = disk.R;
  double disk_shadow_semi_angle = atan(h/R);
  double cos_disk_shadow_semi_angle = cos(0.5 * M_PI - disk_shadow_semi_angle);
  double cos_drd;
  
  double S;
  double Fx;
  double T_star_4 = pow(T,4);
  double T_irr_4;
  double T_sum;
  double F_0;
  double T_Lagrange_point = 0.0;
  
  double color;
  
  /* unity vectior along disk.h */
  vec3 hn;
  hn.x = disk.h.x/h;
  hn.y = disk.h.y/h;
  hn.z = disk.h.z/h;

  vec3 hn_min;
  hn_min.x = -hn.x;
  hn_min.y = -hn.y;
  hn_min.z = -hn.z;

  vec3 hn2;
  hn2.x = d2.x/h;
  hn2.y = d2.y/h;
  hn2.z = d2.z/h;

  vec3 hn2_min;
  hn2_min.x = -hn2.x;
  hn2_min.y = -hn2.y;
  hn2_min.z = -hn2.z;

  static double *Ix_dd = NULL;

  if (Ix_dd == NULL)
    {
      Ix_dd = x_ray_direction_diagram(PSI_pr, Lx);
    }
  else
    {}
  
  /* double *Fx_dd; */
  /* Fx_dd = x_ray_direction_diagram(PSI_pr); */

  vec3 drd_vec3 = sp2dec(disk_reflection_diagr);
    
  double Fx_sum = 0.0;
  sp star_surface_spherical_coordinates;
  
  int diagr_index;

  
  /* flux from star */
  double result = 0.0;


  /* */
  /* double F_0 = F_lambda(T,lambda); */


  
  int true_n_tiles = steps_theta * steps_phi;

  
  /* static double *phi_array = NULL; */

  /* if (phi_array == NULL) */
  /*   { */
  /*     phi_array = phi_func(steps_phi); */
  /*   } */
  /* else */
  /*   {} */

  
  /* static double *theta_array = NULL; */

  /* if (theta_array == NULL) */
  /*   { */
  /*     theta_array = theta_func(steps_theta); */
  /*   } */
  /* else */
  /*   {} */


  /* static double *r_array = NULL; */

  /* if (r_array == NULL) */
  /*   { */
  /*     r_array = shape_r(steps_phi, steps_theta, q, omega); */
  /*   } */
  /* else */
  /*   {} */

  
  /* static double *g_array = NULL; */

  /* if (g_array == NULL) */
  /*   { */
  /*     g_array = shape_g_abs(steps_phi, steps_theta, q, omega); */
  /*   } */
  /* else */
  /*   {} */

  
  int i, j;
  
  for (i = 0; i < steps_phi; i++)
    {

      //phi = (double) i * 2.0 * M_PI/steps_phi + 0.5 * 2.0 * M_PI/steps_phi;
      phi = phi_array[i];
      
      for (j = 0; j < steps_theta; j++)
	{

	  //theta = (double) j * M_PI/steps_theta + 0.5 * M_PI/steps_theta;
	  theta = theta_array[j];

	  /* gradient omega */
	  //double *grd;
	  //grd = gradient(phi, theta, q, omega);
	  //g_abs = grd[0];
	  //g = g_abs/polar_g_abs;

	  /* surface normal vector */	  
	  g   = g_array[(steps_phi * j + i)*4 + 0];
	  n.x = g_array[(steps_phi * j + i)*4 + 1];
	  n.y = g_array[(steps_phi * j + i)*4 + 2];
	  n.z = g_array[(steps_phi * j + i)*4 + 3];
	  //free(grd);

	  /* star`s dot products */
	  cos_on = dot(o,n);

	  if (cos_on < - eps)
	    {
	      continue;
	    }

	  /* star */
	  //r = radius_star(phi, theta, q, omega);

	  r = r_array[steps_phi * j + i];
	    
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


	  cos_drd = dot(psn, drd_vec3);

	  //printf("PSN %f %f %f\t DRD %f %f %f\n", psn.x, psn.y, psn.z, drd_vec3.x, drd_vec3.y, drd_vec3.z);

	  //printf("COS_DRD %f\n", cos_drd);
	  
	  /* Irradiation */	  
	  /* */
	  cos_irr = dot(psn, hn);
	  cos_irr2 = dot(psn, hn2);

	  cos_in  = dot(psn, n);
	  
	  cos_irr_min = dot(psn, hn_min);
	  cos_irr2_min = dot(psn, hn2_min);

	  diagr_index = (int) floor( acos(dot(neutron_star, psn)) * 180.0/M_PI );

	  //printf("%f\t %f\t %f\n", neutron_star.x, neutron_star.y, neutron_star.z);
	  
	  /* Surface element */
	  
	  S = a * a * r * r * sin(theta) * delta_phi * delta_theta / cos_rn;

	  /* printf("%f\n", cos_rn); */

	  /* X-ray flux incident of the surface element */
	  /* if ( cos_in < - eps && fabs(cos_irr) > cos_disk_shadow_semi_angle) */

	  if ( cos_in < - eps && cos_irr > cos_disk_shadow_semi_angle && cos_irr2 > cos_disk_shadow_semi_angle || cos_in < - eps && cos_irr_min > cos_disk_shadow_semi_angle && cos_irr2_min > cos_disk_shadow_semi_angle)
	    {

	      if (isotrope == 0)
		{
		  Fx = Ix_dd[diagr_index] * (1.0 - albedo) * fabs(cos_in) / (lps * lps * a * a) + (1.0 - albedo) * fabs(cos_drd) * Lx_disk_2 * fabs(cos_in) / (2.0 * M_PI * lps * lps * a * a) + (1.0 - albedo) * fabs(cos_irr2) * Lx_disk * fabs(cos_in) / (2.0 * M_PI * lps * lps * a * a) + (1.0 - albedo) * Lx_iso * fabs(cos_in) / (4.0 * M_PI * lps * lps * a * a);
		}
	      else if (isotrope == 1)
		{
		  Fx = (1.0 - albedo) * Lx * fabs(cos_in) / (4.0 * M_PI * lps * lps * a * a);
		}

	    }
	  else
	    {
	      Fx = 0.0;
	    }

	  if (j == steps_theta/2 && i == 0) /* Lagrange point */
		{
		  T_irr_4 = Fx / SIGMA;
		  T_sum = pow((T_star_4 + T_irr_4),0.25); 

		  T_Lagrange_point = T_sum;
		}	  
	  
	  if (ray == 1.0 && cos_on > 0.0 + eps)
	    {

	      T_irr_4 = Fx / SIGMA;
	      
	      T_sum = pow((T_star_4 + T_irr_4),0.25); 

	      F_0 = F_lambda(T_sum, lambda);
	      
	      result = result + F_0 * (1 - u + u * cos_on) * pow(g,beta) * cos_on * S;

	      
	      //color = F_0 * (1 - u + u * cos_on) * pow(g,beta) * cos_on * S;
	      if (picture == 1)
		{
		  star_surface_spherical_coordinates = dec2sp(p);
		  printf("%f\t %f\t %f\t %f\t %f\t %f\t %f\n", phase_orb*180.0/M_PI , p.x, p.y, p.z, T_sum, star_surface_spherical_coordinates.phi, star_surface_spherical_coordinates.theta);
		}
	      else if (picture == 0)
		{}
	    }
	  	  
	}
      
    }

  
  //free(Ix_dd);
  //double *output = (double*) malloc(sizeof(double) * 2);
 
  //output[0] = result; 
  //output[1] = T_Lagrange_point; 

  return result;
  
}



double B(disk disk, double A, double rho_in, double T)
{
  /* auxuliary function for determining T(rho)  */

  double R = disk.R;
  double h = 2.0 * len(disk.h);
  
  return (pow(T,4) * (R*R + R*h) - 2.0 * A * (1.0/rho_in - 1.0/R + h/(2.0*R*R)))/(2.0*(log(R) - log(rho_in) + h/(2.0*R)));

}



  
double flux_disk(vec3 o, disk disk, double rho_in, double A, double uniform_disk, double y_tilt, double z_tilt, double omega, double q, int disk_tiles, double phi_orb, double T, double lambda, double a, int picture, int spot_disk, double T_spot, double spot_beg, double spot_end, double spot_rho_in, double spot_rho_out)
{

  sp coord = dec2sp(o);
  double phase_orb = coord.phi;

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

  vec3 pt_non_shifted;
  
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

  vec3 z_u, z_d;
  z_u.x = 0.0;
  z_u.y = 0.0;
  z_u.z = 1.0;
  z_d.x = 0.0;
  z_d.y = 0.0;
  z_d.z = -1.0;

  /* printf("%f %f %f\n", n1.x, n1.y, n1.z); */

  vec3 nu, nd, ns;
  vec3 nut, ndt, nst;
  
  double lp;
  double len_nu, len_nd;
  
  /* dot products */
  double cos_rn_u;
  double cos_rn_d;
  double cos_rn_s;
    
  /* overlapping */
  double ray;
  double ray_2;

  /* */
  double S_sp, S_cy;
  
  /* Temperature profile of the disk */
  double rho;
  /* double rho_in = 0.000236; inner radius where T becomes 0 for rho < rho_in (approx. 100 R_NS, R_NS = 15 km) */
  //double rho_scale = 0.5 * R; /* rho scale for T */

  //double T_scale = (sqrt(R)*T)/pow(2.0*(-1.0/R + 1.0/rho_in), 1.0/4.0);

  // 2.0 * h because h is semithickness
  //double T_scale = T * pow((R * R + R * 2.0 * h)/(1.0/(R*R) - 1.0/R + 1.0/rho_in), 1.0/4.0);
  //double T_scale = T * pow((R * R + R * 2.0 * h)/(2.0*h/R + log(R) - log(rho_in)), 1.0/4.0);

  double T_scale_A = (A * pow(T,4.0) * (R * R + R * 2.0*h)) / (2.0 * (1.0/rho_in - 1.0/R + 2.0*h/(2.0*R*R)));
  double T_scale_B = ((1.0 - A) * pow(T,4.0) * (R * R + R * 2.0*h)) / (2.0 * (log(R) - log(rho_in) + 2.0*h/(2.0*R)));

  //printf("T_scale_A %f\nT_scale_B %f\n", T_scale_A, T_scale_B);

  
  //printf("T_scale_B %f\n", T_scale_B);

  //printf("T_scale\t%f\tlog(R)\t%f\tlog(rho_in)\t%f\n", T_scale, log(R), log(rho_in));
  
  double T_rho;
  double T_color; /* picture color */
  
  double color;

  sp disk_spherical_coord;

  
  /* */
  /* double F_0 = F_lambda(T, lambda); temperature may depend on rho, in that case set it in the cycle below */
  double F_0; /* temperature may depend on rho, in that case set it in the cycle below */

  double F_spot = F_lambda(T_spot, lambda); /* side of the disk */

  /* from top side of the disk */
  double result_1 = 0.0;
  /* from ridge of the disk */
  double result_2 = 0.0;
  /* from bottom side of the disk */
  double result_3 = 0.0;

  double cos_on;
  double cos_on1;
  double cos_on2;

  double cos_onut;
  double cos_ondt;
  double cos_onst;

  double *plr;
  plr = polar(q, omega);
  double max_r = plr[2]; 
  free(plr);


  static double *phi_array_disk = NULL;

  if (phi_array_disk == NULL)
    {
      phi_array_disk = phi_func_disk(steps_phi);
    }
  else
    {}




  
  int i, j;

  for (i = 0; i < steps_phi; i++)
    {
      //phi = (double) i * 2.0 * M_PI/steps_phi + 0.5 * 2.0 * M_PI/steps_phi;

      phi = phi_array_disk[i];
      
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
	      //cos_on = dot(o,n1);
		  
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
	      //cos_on = dot(o,n3);

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
	      //cos_on = dot(o,n2);  
	    }

	  /* if (cos_on < -eps) */
	  /*   { */
	  /*     continue; */
	  /*   } */

	  
	  r = radius_disk(disk, phi, theta);

	  /* cartesian coordinates of the point */
	  p.x = r * sin(theta) * cos(phi);
	  p.y = r * sin(theta) * sin(phi);
	  p.z = r * cos(theta);


	  /* disk with complex profile h = h(rho) */

	  if (j <= N - 1)
	    {
	      p.z = p.z + (h * sqrt(p.x*p.x + p.y*p.y))/R - h;
	    }
	  if (j >= M && j <= steps_theta)
	    {
	      p.z = p.z - (h * sqrt(p.x*p.x + p.y*p.y))/R + h;
	    }

	  /* Temperature profile of the disk */
	  rho = r * sin(theta);


	  
	  if (rho > rho_in)
	    {
	      /* T_rho = T_scale * pow(rho, -3./4.); */
	      if (uniform_disk == 0.0)
		{
		  T_rho = pow(T_scale_A * pow(rho, -3.0) + T_scale_B * pow(rho, -2.0), 1.0/4.0);
		}
	      else if (uniform_disk == 1.0)
		{
		  T_rho = T;
		}
	    }
	  else
	    {
	      T_rho = 0.0;
	    }

	  //printf("%f\t%f\n", uniform_disk, T_rho);
	  
	  //F_0 = F_lambda(T_rho, lambda);

	  //printf("T_rho %f \t F_0 %f \n", T_rho, F_0);
	  
	  /* flux on on wavelenght lambda */
	  F_0 = F_lambda(T_rho, lambda);
	  

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
	  /* simple constant profile disk */
	  nu.x = 0.0;
	  nu.y = 0.0;
	  nu.z = 1.0;

	  /* normal vector for complex profile disk */
	  nu.x = -(h/R) * p.x/sqrt(p.x*p.x + p.y*p.y);
	  nu.y = -(h/R) * p.y/sqrt(p.x*p.x + p.y*p.y);
	  nu.z = 1.0;

	  len_nu = len(nu);
	  nu.x = nu.x/len_nu;
	  nu.y = nu.y/len_nu;
	  nu.z = nu.z/len_nu;
	  /**/
	  
	  /* bottom of the disc */
	  /* simple constant profile disk */	  
	  nd.x = -(h/R) * p.x/sqrt(p.x*p.x + p.y*p.y);
	  nd.y = -(h/R) * p.y/sqrt(p.x*p.x + p.y*p.y);
	  nd.z = -1.0;


	  len_nd = len(nd);
	  nd.x = nd.x/len_nd;
	  nd.y = nd.y/len_nd;
	  nd.z = nd.z/len_nd;

	  
	  /* side of the disc */
	  ns.x = cos(phi);
	  ns.y = sin(phi);
	  ns.z = 0.0;
	  
	  cos_rn_u = dot(z_u,nu);
	  cos_rn_d = dot(z_d,nd);
	  cos_rn_s = dot(pn,ns);

	  /* printf("%f %f %f\n", len(n1), len(n2), len(n3)); */
	  
	  /* cos_on_u = dot(o,n1); */
	  /* cos_on_d = dot(o,n2);  */
	  /* cos_on_s = dot(o,n3); */

	  /* printf("%f %f %f\n", o.x, o.y, o.z); */

	  /* sphere surface element */
	  S_sp = a * a * r * r * sin(theta) * delta_phi * delta_theta;

	  /* cylindrical surface element */
	  S_cy = a * a * r * delta_phi * delta_N;
 
	  
	  /* tilt and shift disc */
	  pt.x =   p.x * cos(y_tilt) * cos(z_tilt - phi_orb + M_PI) + p.y * sin(z_tilt - phi_orb + M_PI) - p.z * sin(y_tilt) * cos(z_tilt - phi_orb + M_PI) + 1.0;
	  pt.y = - p.x * cos(y_tilt) * sin(z_tilt - phi_orb + M_PI) + p.y * cos(z_tilt - phi_orb + M_PI) + p.z * sin(y_tilt) * sin(z_tilt - phi_orb + M_PI);
	  pt.z =   p.x * sin(y_tilt) + p.z * cos(y_tilt);

	  nut.x =   nu.x * cos(y_tilt) * cos(z_tilt - phi_orb + M_PI) + nu.y * sin(z_tilt - phi_orb + M_PI) - nu.z * sin(y_tilt) * cos(z_tilt - phi_orb + M_PI);
	  nut.y = - nu.x * cos(y_tilt) * sin(z_tilt - phi_orb + M_PI) + nu.y * cos(z_tilt - phi_orb + M_PI) + nu.z * sin(y_tilt) * sin(z_tilt - phi_orb + M_PI);
	  nut.z =   nu.x * sin(y_tilt) + nu.z * cos(y_tilt);

	  ndt.x =   nd.x * cos(y_tilt) * cos(z_tilt - phi_orb + M_PI) + nd.y * sin(z_tilt - phi_orb + M_PI) - nd.z * sin(y_tilt) * cos(z_tilt - phi_orb + M_PI);
	  ndt.y = - nd.x * cos(y_tilt) * sin(z_tilt - phi_orb + M_PI) + nd.y * cos(z_tilt - phi_orb + M_PI) + nd.z * sin(y_tilt) * sin(z_tilt - phi_orb + M_PI);
	  ndt.z =   nd.x * sin(y_tilt) + nd.z * cos(y_tilt);

	  nst.x =   ns.x * cos(y_tilt) * cos(z_tilt - phi_orb + M_PI) + ns.y * sin(z_tilt - phi_orb + M_PI) - ns.z * sin(y_tilt) * cos(z_tilt - phi_orb + M_PI);
	  nst.y = - ns.x * cos(y_tilt) * sin(z_tilt - phi_orb + M_PI) + ns.y * cos(z_tilt - phi_orb + M_PI) + ns.z * sin(y_tilt) * sin(z_tilt - phi_orb + M_PI);
	  nst.z =   ns.x * sin(y_tilt) + ns.z * cos(y_tilt);
	  
	  cos_onut = dot(o,nut);
	  cos_ondt = dot(o,ndt);
	  cos_onst = dot(o,nst);
	  
	  /* ray = eclipse_by_star(omega, q, o, pt); */
	  if (o.x < - cos(R + max_r))
	    {
	      ray = eclipse_by_star(omega, q, o, pt);
	      //ray = 1.0;
	    }
	  else if (o.x > - cos(R + max_r))
	    {
	      ray = 1.0;
	    }
	  
	  /* ray = 1.0; */
	  /* printf("%d\n", j); */

	  
	  pt_non_shifted.x = pt.x - 1.0;
	  pt_non_shifted.y = pt.y;
	  pt_non_shifted.z = pt.z;
	  
	  disk_spherical_coord = dec2sp(pt_non_shifted);

	  //cos_on1 = dot(n1,o);
	  //cos_on2 = dot(n2,o);

	  
	  ray_2 = eclipse_by_disk_inside(disk, o, pt_non_shifted);
	  //ray_2 = 1.0;
	  
	  if ( j <= N - 1 && ray == 1.0  &&  ray_2 == 1.0  && cos_onut > eps )
	    {
	      /* top surface */
	      /* printf("%f\t %f\t %f\n", pt.x, pt.y, pt.z); */

	      if (spot_disk == 1)
		{
		  if ( disk_spherical_coord.phi >= spot_beg && disk_spherical_coord.phi <= spot_end &&
		       disk_spherical_coord.r >= spot_rho_in && disk_spherical_coord.r <= spot_rho_out)
		    {
		      result_1 = result_1 + (F_spot * cos_onut * S_cy)/cos_rn_u;
		      T_color = T_spot;
		      /* T_color = (F_spot * cos_onut * S_cy)/cos_rn_u; */
		    }
		  else
		    {
		      result_1 = result_1 + (F_0 * cos_onut * S_cy)/cos_rn_u;
		      T_color = T_rho;
		      /* T_color = (F_0 * cos_onut * S_cy)/cos_rn_u; */
		    }
		}
	      else if (spot_disk == 0 || spot_disk == 2 || spot_disk == 3)
		{
		  result_1 = result_1 + (F_0 * cos_onut * S_cy)/cos_rn_u;
		  T_color = T_rho;
		  /* T_color = (F_0 * cos_onut * S_cy)/cos_rn_u; */
		}
	      
	      //result_1 = result_1 + (F_0 * cos_on * S)/cos_rn_u;
	      /* result_1 = result_1 + cos_rn_u; */

	      //color = (F_0 * cos_on * S)/cos_rn_u;
	      if (picture == 1)
		{
		  printf("%f\t %f\t %f\t %f\t %f\n", phase_orb*180.0/M_PI, pt.x, pt.y, pt.z, T_color);
		}
	      else if (picture == 0)
		{}

	    }
	  else if ( j >= M && j <= steps_theta && ray == 1.0 && ray_2 == 1.0 && cos_ondt > eps)
	    {
	      /* bottom surface */
	      /* printf("%f\t %f\t %f\n", pt.x, pt.y, pt.z); */

	      if (spot_disk == 3)
		{
		  if ( disk_spherical_coord.phi >= spot_beg && disk_spherical_coord.phi <= spot_end &&
		       disk_spherical_coord.r >= spot_rho_in && disk_spherical_coord.r <= spot_rho_out)
		    {
		      result_2 = result_2 + (F_spot * cos_ondt * S_cy)/cos_rn_d;
		      T_color = T_spot; 
		    }
		  else
		    {
		      result_2 = result_2 + (F_0 * cos_ondt * S_cy)/cos_rn_d;
		      T_color = T_rho;
		    }
		}
	      else if (spot_disk == 0 || spot_disk == 1 || spot_disk == 2)
		{
		  result_2 = result_2 + (F_0 * cos_ondt * S_cy)/cos_rn_d;
		  T_color = T_rho;
		}
	      
	      //result_2 = result_2 + (F_0 * cos_on * S)/cos_rn_d;
	      
	      /* result_2 = result_2 + cos_rn_d; */
	      
	      if (picture == 1)
		{
		  printf("%f\t %f\t %f\t %f\t %f\n", phase_orb*180.0/M_PI, pt.x, pt.y, pt.z, T_color);
		}
	      else if (picture == 0)
		{}
	      /* printf("%f\n", color); */

	    }
	  else if ( j >= N && j <= M - 1 && ray == 1.0 && cos_onst > eps)
	    {
	      /* side surface */
	      /* printf("%f\t %f\t %f\n", pt.x, pt.y, pt.z); */
	      
	      //printf("%d\t%f\t\t%f\n", i, (phi_orb * 180.0/M_PI + 180.0), disk_spherical_coord.phi);

	      if (spot_disk == 2)
		{
		  if ( disk_spherical_coord.phi >= spot_beg && disk_spherical_coord.phi <= spot_end )
		    {
		      result_3 = result_3 + (F_spot * cos_onst * S_sp)/cos_rn_s;
		      T_color = T_spot;
		      /* T_color = (F_spot * cos_onst * S_sp)/cos_rn_s; */
		      
		    }
		  else
		    {
		      result_3 = result_3 + (F_0 * cos_onst * S_sp)/cos_rn_s;
		      T_color = T_rho;
		      /* T_color = (F_0 * cos_onst * S_sp)/cos_rn_s; */
		    }
		}
	      else if (spot_disk == 0 || spot_disk == 1 || spot_disk == 3)
		{
		  result_3 = result_3 + (F_0 * cos_onst * S_sp)/cos_rn_s;
		  T_color = T_rho;
		  /* T_color = (F_0 * cos_onst * S_sp)/cos_rn_s; */
		}
	      //result_3 = result_3 + (F_side * cos_on * S)/cos_rn_s;

	      //color = (F_0 * cos_on * S)/cos_rn_u;
	      if (picture == 1)
		{
		  printf("%f\t %f\t %f\t %f\t %f\n", phase_orb*180.0/M_PI, pt.x, pt.y, pt.z, T_color);
		}
	      else if (picture == 0)
		{}
	    }



	  /* printf("%d\n", j); */
	}
      
    }

  
  return result_1 + result_2 + result_3;
}



