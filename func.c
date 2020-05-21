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
      
      //printf("Ray tracing iteration %d\n", i);
      
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
  double h = disk.h;
  double r = len(p);
  double hr = dot(scale(disk.n, h), p);
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



double disk_h_profile(disk disk, double r, double gamma)
{

  double r_in = disk.r_in;
  
  return disk.h * pow(((r - r_in)/(disk.R-r_in)), gamma);

}


double disk_h_diff_profile(disk disk, double r, double gamma)
{
  /* differential by r */
  double r_in = disk.r_in;
  
  return disk.h * gamma * pow(((r - r_in)/(disk.R-r_in)), gamma)/(r - r_in);

}





double star_F(double * star_elements, parameters parameters, disk disk, vec3 observer)
{
  
  int N = star_elements[0];

  double T;

  vec3 p;
  vec3 n;

  double s, g;

  double zeta;
  double fx;


  vec3 ns;
  ns.x = - 1.0;
  ns.y = 0.0;
  ns.z = 0.0;
  
  
  for (int i = 0; i < N; i++)
    {
      p.x = star_elements[8*i + 1];
      p.y = star_elements[8*i + 2];
      p.z = star_elements[8*i + 3];

      n.x = star_elements[8*i + 4];
      n.y = star_elements[8*i + 5];
      n.z = star_elements[8*i + 6];

      s = star_elements[8*i + 7];
      g = star_elements[8*i + 8];

      
      zeta = dot(n, scale(sum(p, ns), len(sum(p,ns))));


      
      if ( ray_disk(disk, observer, sum(p, ns)) && dot(n, observer) > 0 )
	{
	  fx = 0.0;

	  T = parameters.T_star_polar * pow(g, parameters.beta);


	  
	  if (zeta < 0)
	    {
	      fx = (parameters.Lx/(4.0 * M_PI * pow(len(p)*parameters.a,2))) * fabs(zeta) * (double) disk_shadow(sum(p, ns), disk);
	      T = pow((pow(T,4) + fx/SIGMA), 1.0/4.0);

	      //printf(">>>> %f\t%f\n", (double) disk_shadow(p, disk), fx);
	    }
	  else
	    {}

	  //F_0 = F_filter(T, filter);
	  //F_0 = photometric_filter(T, filter);

	  //printf("%f\t%f\t%f\t%f\t%f\n", w.phi * (180.0/M_PI), p.x, p.y, p.z, T);
	  /* printf("%f\t%f\t%f\t%f\n", p.x, p.y, p.z, T); */
	  
	  //printf("%e\n", s);
	}
    }

  return 0;

}




double flux_star(vec3 o, double q, double omega, double beta, double u, disk disk, double Lx_noniso, double Lx_disk, double Lx_iso, double Lx_disk_2, double albedo, int tiles, double T, double a, vec3 neutron_star, double PSI_pr, int picture, sp disk_reflection_diagr, double * r_array, double * g_array, double * phi_array, double * theta_array, double * Ix_dd, double y_tilt, double y_tilt2, double z_tilt, double z_tilt2, double phi_orb, char * filter)
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
  

  /* irradiation */

  double cos_in;
  double h = disk.h;
  double R = disk.R;
  double r_in = 0.0002;
  //double disk_shadow_semi_angle = atan(h/R);
  //double cos_disk_shadow_semi_angle = cos(0.5 * M_PI - disk_shadow_semi_angle);
  //double cos_drd;
  
  double S;
  double Fx;
  double T_star_4 = pow(T,4);
  double T_irr_4;
  double T_sum;
  double F_0;
  double T_Lagrange_point = 0.0;
  
  double color;
  


  /* vec3 drd_vec3 = sp2dec(disk_reflection_diagr); */
    
  double Fx_sum = 0.0;
  sp star_surface_spherical_coordinates;
  
  int diagr_index;



  double zeta; // cos incident X-ray and observer
  
  double energy = 1.0; // energy in keV


  
  /* flux from star */
  double result = 0.0;
  
  int true_n_tiles = steps_theta * steps_phi;


  bool X_shadow;
  bool ray;

  int number_of_rings_in_disk = 20;

  vec3 ring_k;
  vec3 ring_kplus1;

  double y_tilt_delta, y_tilt2_delta;
  double z_tilt2_delta;

  double y_tilt_delta_plus1, y_tilt2_delta_plus1;
  double z_tilt2_delta_plus1;

  double Z;
  
  int i, j, k;
  
  for (i = 0; i < steps_phi; i++)
    {

      //phi = (double) i * 2.0 * M_PI/steps_phi + 0.5 * 2.0 * M_PI/steps_phi;
      phi = phi_array[i];
      
      for (j = 0; j < steps_theta; j++)
	{

	  //theta = (double) j * M_PI/steps_theta + 0.5 * M_PI/steps_theta;
	  theta = theta_array[j];


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

	  //ray = eclipse_by_disk(disk, o, ps);

	  /* ray = 1.0; */
	  
	  /* if (o.x < cos(R + max_r)) */
	  /*   { */
	  /*     ray = 1.0; */
	  /*   } */
	  /* else if (o.x > cos(R + max_r)) */
	  /*   { */
	  /*     ray = eclipse_by_disk(disk, o, ps); */
	  /*     //ray = 1.0; */
	  /*   } */

	  //ray = 1.0;
	  
	  /* normalized vector */
	  pl = len(p);
	  pn.x = p.x/pl;
	  pn.y = p.y/pl;
	  pn.z = p.z/pl;
	  
	  cos_rn = dot(pn,n);


	  lps = len(ps); /* distance from the secondary to the point p */
	  psn.x = ps.x/lps;
	  psn.y = ps.y/lps;
	  psn.z = ps.z/lps;


	  /* Irradiation */	  

	  zeta = dot(psn, n); 
	  
	  
	  diagr_index = (int) floor( acos(dot(neutron_star, psn)) * 180.0/M_PI );

	  
	  /* Surface element of the star */	  
	  S = a * a * r * r * sin(theta) * delta_phi * delta_theta / cos_rn;
	  
	  	  

	  X_shadow = disk_shadow(ps, disk);
	  ray = ray_disk(disk, o, ps);



	  
	  //if ( (cos_in < - eps && cos_irr > cos_disk_shadow_semi_angle && cos_irr2 > cos_disk_shadow_semi_angle || cos_in < - eps && cos_irr_min > cos_disk_shadow_semi_angle && cos_irr2_min > cos_disk_shadow_semi_angle) && shadow_condition )
	  if (zeta < - eps && X_shadow )
	    {
	      
	      /* Fx = */
	      /* 	Lx_noniso * Ix_dd[diagr_index] * (1.0 - albedo) * fabs(zeta) / (lps * lps * a * a) + */
	      /* 	Lx_iso * (1.0 - albedo) * fabs(zeta) / (4.0 * M_PI * lps * lps * a * a); */

	      Fx = Lx_iso / (4.0 * M_PI * lps * lps * a * a);
	      //Fx = Lx_noniso * Ix_dd[diagr_index] / (lps * lps * a * a);
	      
	    }
	  else
	    {
	      Fx = 0.0;
	    }

	  
	  if (ray && cos_on > eps)
	    {

	      T_irr_4 = Fx / SIGMA;
	      
	      T_sum = pow((T_star_4 + T_irr_4),0.25); 
	      F_0 = F_filter(T_sum, filter);


	      /* optical flux */
	      //result = result + F_0 * (1 - u + u * cos_on) * pow(g,beta) * cos_on * S;

	      
	      /* reflected X-ray flux */
	      //result = result + Fx * S * cos_on * x_rays_reflected_fraction(zeta, energy);
	      /* result = result + Fx * S * cos_on * x_rays_reflected_fraction(fabs(zeta), energy); */
	      result = result + Fx * S * cos_on * x_ray_flux_integrated(-zeta);

	      	      
	      //color = F_0 * (1 - u + u * cos_on) * pow(g,beta) * cos_on * S;
	      if (picture == 1)
	      	{
	      	  printf("%f\t %f\t %f\t %f\t %f\t \n", phase_orb*180.0/M_PI , p.x, p.y, p.z, T_sum);
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

  //printf("result %f\n", result);
  return result;
  
}






double disk_F(double * disk_elements, parameters parameters, disk disk, vec3 observer)
{

  int N = (int) disk_elements[0];
  //printf(">> %d\n", N);
  
  vec3 p;
  vec3 n;
  double s;

  double fx;
  
  double T = 10000.0;
  
  double zeta;

  sp w = dec2sp(observer);
  
  for (int i = 0; i < N; i++)
    {
      p.x = disk_elements[7*i + 1];
      p.y = disk_elements[7*i + 2];
      p.z = disk_elements[7*i + 3];

      n.x = disk_elements[7*i + 4];
      n.y = disk_elements[7*i + 5];
      n.z = disk_elements[7*i + 6];

      s = disk_elements[7*i + 7];

      zeta = dot(n, scale(p, len(p)));

      
      
      if ( ray_disk(disk, observer, p) && dot(n, observer) > 0 )
	{
	  fx = 0.0;
	  T = 10000.0;

	  if (zeta < 0)
	    {
	      fx = (parameters.Lx/(4.0 * M_PI * pow(len(p)*parameters.a,2))) * fabs(zeta) * (double) disk_shadow(p, disk);
	      T = pow((pow(T,4) + fx/SIGMA), 1.0/4.0);

	      //printf(">>>> %f\t%f\n", (double) disk_shadow(p, disk), fx);
	    }
	      

	  //F_0 = F_filter(T, filter);
	  //F_0 = photometric_filter(T, filter);

	  /* printf("%f\t%f\t%f\t%f\t%f\n", w.phi * (180.0/M_PI), p.x + 1.0, p.y, p.z, T); */
	  /* printf("%f\t%f\t%f\t%f\n", p.x + 1.0, p.y, p.z, T); */
	  
	  //printf("%e\n", s);

	  
	}
    }


  return 0.0;


}






  
double * flux_disk(vec3 o, disk disk, double rho_in, double A, double uniform_disk, double y_tilt, double z_tilt, double omega, double q, int disk_tiles, double phi_orb, double T, double a, int picture, int spot_disk, double T_spot, double spot_beg, double spot_end, double spot_rho_in, double spot_rho_out, double h_warp, double * Ix_dd, double Lx_iso, double Lx_noniso, vec3 neutron_star, char * filter)
{

  sp coord = dec2sp(o);
  double phase_orb = coord.phi;

  double R = disk.R;
  double h = disk.h; /* semithickness of the disk */
  /* double h = h_warp; /\* semithickness of the disk *\/ */

  /* */
  double theta, phi;
  
  int i, j;
  
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
  bool ray_2;

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

  double T_rho = 0.0;

  double T_rho_up;
  double T_rho_down;

  double T_color; /* picture color */
  
  double color;

  sp disk_spherical_coord;



  /* temperature of the disk due to irradiation */


  sp s_disk_irradiation;
  s_disk_irradiation.theta = 0.0;
  s_disk_irradiation.phi = 0.0;
  s_disk_irradiation.r = 1.0;

  double theta_disk_irradiation_num_of_steps = 50;
  double phi_disk_irradiation_num_of_steps = 2.0 * theta_disk_irradiation_num_of_steps;

  vec3 p_disk_irradiation;


  //vec3 d = n1;

  /* double d2_len = len(d2); */

  /* d2.x = d2.x/d2_len; */
  /* d2.y = d2.y/d2_len; */
  /* d2.z = d2.z/d2_len; */

  //vec3 d_in_NS;
  //vec3 d2_in_NS;

  sp s_neutron_star = dec2sp(neutron_star);

  //d = rotate(d, 0.0, -phi_orb);
  //d2 = rotate(d2, 0.0, -phi_orb);

  //double kappa = 15.0 * M_PI/180.0;
  //double ns_theta = -3.0 * M_PI/180.0;

  // check sign of the phi_orb!!! DK 8Feb2019
  //d_in_NS = rotate(d, -s_neutron_star.theta, -s_neutron_star.phi + phi_orb);
  //d2_in_NS = rotate(d2, -s_neutron_star.theta, -s_neutron_star.phi + phi_orb);

  double L_sum_up = 0.0;
  double L_sum_down = 0.0;

  int diagr_index;

  //double alpha = acos(dot(d,d2));

  



	 
  /* double T_up = pow(((L_sum_up/(R*R*a*a) + (Lx*h)/(14*M_PI*R*R*R*a*a))/SIGMA), 0.25); */

  /* double T_down = pow(((L_sum_down/(R*R*a*a) + (Lx*h)/(14*M_PI*R*R*R*a*a))/SIGMA), 0.25); */

  double T_up = T;
  double T_down = T;

  
  /* */
  double F_0 = 0.0; /* temperature may depend on rho, in that case set it in the cycle below */
  double F_0_up; /* temperature may depend on rho, in that case set it in the cycle below */
  double F_0_down; /* temperature may depend on rho, in that case set it in the cycle below */

  double Fx;

  double energy = 1.0; // keV
  
  double F_spot = F_filter(T_spot, filter); /* side of the disk */

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


  double pl;
  double cos_pnu, cos_pnd; 

  /* disk's h profile */
  double gamma = 1.0/8.0;

  
  static double *phi_array_disk = NULL;

  if (phi_array_disk == NULL)
    {
      phi_array_disk = phi_func_disk(steps_phi);
    }
  else
    {}


  for (i = 0; i < steps_phi; i++)
    {
      phi = (double) i * 2.0 * M_PI/steps_phi + 0.5 * 2.0 * M_PI/steps_phi;

      //phi = phi_array_disk[i];
      
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
	      /* n3.x = sin(phi) * sin(z_tilt + phi_orb) + cos(phi) * cos(y_tilt) * cos(z_tilt + phi_orb); */
	      /* n3.y = sin(phi) * cos(z_tilt + phi_orb) - cos(phi) * cos(y_tilt) * sin(z_tilt + phi_orb); */
	      /* n3.z = cos(phi) * sin(y_tilt); */
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

	  
	  //p.z = disk_h_profile(disk, r, gamma, rho_in);


	  /* normalized p vector */
	  pl = len(p);
	  pn.x = p.x/pl;
	  pn.y = p.y/pl;
	  pn.z = p.z/pl;


	  
	  /* disk with complex profile h = h(rho) */
	  

	  
	  
	  if (j <= N - 1)
	    {
	      //p.z = p.z + (h * sqrt(p.x*p.x + p.y*p.y))/R - h;
	      p.z = disk_h_profile(disk, r, gamma);
	      
	    }
	  if (j >= M && j <= steps_theta)
	    {
	      //p.z = p.z - (h * sqrt(p.x*p.x + p.y*p.y))/R + h;
	      p.z = - disk_h_profile(disk, r, gamma);

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
		  T_rho_up = T_up;
		  T_rho_down = T_down;
		}
	    }
	  else
	    {
	      T_rho = 0.0;
	    }

	  
	  
	  /* F_0_up   = F_filter(T_rho_up, filter); */
	  /* F_0_down = F_filter(T_rho_down, filter); */

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

	  /* normal vector for complex profile disk */


	  nu.x = - sin(atan(disk_h_diff_profile(disk, r, gamma))) * cos(phi);
	  nu.y = - sin(atan(disk_h_diff_profile(disk, r, gamma))) * sin(phi);
	  nu.z = cos(atan(disk_h_diff_profile(disk, r, gamma)));

	  
	  /* bottom of the disc */
	  /* simple constant profile disk */	  

	  nd.x = - sin(atan(disk_h_diff_profile(disk, r, gamma))) * cos(phi);
	  nd.y = - sin(atan(disk_h_diff_profile(disk, r, gamma))) * sin(phi);
	  nd.z = - cos(atan(disk_h_diff_profile(disk, r, gamma)));
	  
	  
	  /* side of the disc */
	  ns.x = cos(phi);
	  ns.y = sin(phi);
	  ns.z = 0.0;
	  
	  cos_rn_u = dot(z_u,nu);
	  cos_rn_d = dot(z_d,nd);
	  cos_rn_s = dot(pn,ns);


	  /* sphere surface element */
	  S_sp = a * a * r * r * sin(theta) * delta_phi * delta_theta;

	  /* cylindrical surface element */
	  S_cy = a * a * r * delta_phi * delta_N;


	  diagr_index = (int) floor( acos(dot(neutron_star, pn)) * 180.0/M_PI );


	  cos_pnu = dot(pn, nu);
	  cos_pnd = dot(pn, nd);

     	    
	  Fx = Lx_iso / (4.0 * M_PI * r * r * a * a);


	  //Lx_noniso * Ix_dd[diagr_index] * fabs(cos_pnu) / (r * r * a * a)

	  //printf("F_0_up\t %f\n", F_0_up);
	  
	  //Lx_noniso * Ix_dd[diagr_index] * fabs(cos_pnd) / (r * r * a * a)
	    
	  //printf("F_0_down\t %f\n", F_0_down);
	  

	  pt = rotate(p, y_tilt, z_tilt + phi_orb);
	  pt.x = pt.x + 1;
	  

	  nut = rotate(nu, y_tilt, z_tilt + phi_orb);
	  

	  ndt = rotate(nd, y_tilt, z_tilt + phi_orb);
	  

	  nst = rotate(ns, y_tilt, z_tilt + phi_orb);


	  	  
	  cos_onut = dot(o, nut);
	  cos_ondt = dot(o, ndt);
	  cos_onst = dot(o, nst);
	  
	  ray = eclipse_by_star(omega, q, o, pt);
	  	  
	  pt_non_shifted.x = pt.x - 1.0;
	  pt_non_shifted.y = pt.y;
	  pt_non_shifted.z = pt.z;
	  
	  disk_spherical_coord = dec2sp(pt_non_shifted);
	  
	  ray_2 = ray_disk(disk, o, p);
	  
	  //if ( j <= N - 1 && ray == 1.0  &&  ray_2 == 1.0  && cos_onut > eps && n1to > eps)
	  if ( j <= N - 1 && ray == 1.0  && ray_2 && cos_onut > eps )
	    {
	      /* top surface */

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
		      result_1 = result_1 + (F_0_up * cos_onut * S_cy)/cos_rn_u;
		      T_color = T_rho_up;
		      /* T_color = (F_0 * cos_onut * S_cy)/cos_rn_u; */
		    }
		}
	      else if (spot_disk == 0 || spot_disk == 2 || spot_disk == 3)
		{
		  result_1 = result_1 + (F_0_up * cos_onut * S_cy)/cos_rn_u;
		  T_color = T_rho_up;
		  /* T_color = (F_0 * cos_onut * S_cy)/cos_rn_u; */
		}
	      
	      //result_1 = result_1 + (F_0_up * cos_onut * S_cy)/cos_rn_u;
	      // result_1 = result_1 + Fx * S_cy * cos_onut * x_rays_reflected_fraction(fabs(cos_pnu), energy);
	      result_1 = result_1 + Fx * S_cy * cos_onut * x_ray_flux_integrated(fabs(cos_pnu));
	      //result_1 = result_1 + (F_0_up * S_cy)/cos_rn_u;

	      
	      /* result_1 = result_1 + cos_rn_u; */
	      
	      //color = (F_0 * cos_on * S)/cos_rn_u;
	      if (picture == 1)
		{
		  printf("%f\t %f\t %f\t %f\t %f\n", phase_orb*180.0/M_PI, pt.x, pt.y, pt.z, 0.0);
		}
	      else if (picture == 0)
		{}

	    }
	  else if ( j >= M && j <= steps_theta && ray == 1.0 && ray_2 == 1.0 && cos_ondt > eps )
	    {
	      /* bottom surface */

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
		      result_2 = result_2 + (F_0_down * cos_ondt * S_cy)/cos_rn_d;
		      T_color = T_rho_down;
		    }
		}
	      else if (spot_disk == 0 || spot_disk == 1 || spot_disk == 2)
		{
		  result_2 = result_2 + (F_0_down * cos_ondt * S_cy)/cos_rn_d;
		  T_color = T_rho_down;
		}
	      
	      //result_2 = result_2 + (F_0_down * cos_ondt * S_cy)/cos_rn_d;
	      //result_2 = result_2 + Fx * S_cy * cos_ondt * x_rays_reflected_fraction(fabs(cos_pnd), energy);
	      result_2 = result_2 + Fx * S_cy * cos_ondt * x_ray_flux_integrated(fabs(cos_pnd));

	      /* result_2 = result_2 + cos_rn_d; */

	      if (picture == 1)
		{
		  printf("%f\t %f\t %f\t %f\t %f\n", phase_orb*180.0/M_PI, pt.x, pt.y, pt.z, 0.0);
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
		  printf("%f\t %f\t %f\t %f\t %f\n", phase_orb*180.0/M_PI, pt.x, pt.y, pt.z, 0.0);
		}
	      else if (picture == 0)
		{}
	    }



	  /* printf("%d\n", j); */
	}
      
    }

  result_3 = 0.0;

  double * result = (double *) malloc(sizeof(double) * 3);
  
  /* return result_1 + result_2 + result_3; */


  result[0] = result_1 + result_2 + result_3;
  result[1] = T_up;
  result[2] = T_down;

  return result;

}
