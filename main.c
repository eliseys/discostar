#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "def.h"

int main(int argc, char **argv)
{

  /* star parameters */
  double q; /* mass ratio q = m_x/m_opt */
  double mu; /* roche lobe filling */
  double beta;    
  double u;   
  double albedo;  /* 1 - (X-ray photons reprocessing efficiency) */	 

  /* neutron star parameter */
  double Lx_noniso;     /* in units of L of the optical star */ 

  double h;      /* semi-thickness */ 
  double R;      /* radius */
  double y_tilt; 
  double z_tilt;  

  double y_tilt2; 
  double z_tilt2;  

  int picture;      /* print 3D picture */
  
  /* observer */
  double inclination;

  int lc_num; /* number of light curve points */

  /* star and disk surface partition */
  int star_tiles;
  int disk_tiles;

  int threads; /* OpenMP threads */

  /* temperature of star and disk and spectral band */
  double T_disk;
  double T_star;

  char filter[5];

  double a_cm;

  double PSI_pr;
  double kappa;

  double Lx_disk;
  double Lx_disk_2;
  double Lx_iso;

  int spot_disk;
  double T_spot;
  double spot_beg;
  double spot_end;

  double ns_theta;

  double spot_rho_in;
  double spot_rho_out;

  double drd_phi;
  double drd_theta;

  double rho_in;

  double A;

  double uniform_disk;

  double disk_flux;
  double h_warp;

  double lc_start;
  double lc_end;

  
  /**/
  
  sscanf(argv[1], "%lf", &q);
  sscanf(argv[2], "%lf", &mu);
  sscanf(argv[3], "%lf", &beta);
  sscanf(argv[4], "%lf", &u);
  sscanf(argv[5], "%lf", &albedo);  
  sscanf(argv[6], "%lf", &Lx_noniso);  
  sscanf(argv[7], "%lf", &h); 
  sscanf(argv[8], "%lf", &R); 
  sscanf(argv[9], "%lf", &y_tilt); 
  sscanf(argv[10], "%lf", &z_tilt);
  sscanf(argv[11], "%d", &picture);
  sscanf(argv[12], "%lf", &inclination);
  sscanf(argv[13], "%d", &lc_num);
  sscanf(argv[14], "%d", &star_tiles);
  sscanf(argv[15], "%d", &disk_tiles);
  sscanf(argv[16], "%d", &threads);
  sscanf(argv[17], "%lf", &T_disk);
  sscanf(argv[18], "%lf", &T_star);
  sscanf(argv[19], "%lf", &a_cm);
  sscanf(argv[20], "%lf", &y_tilt2);
  sscanf(argv[21], "%lf", &z_tilt2);
  sscanf(argv[22], "%lf", &PSI_pr);
  sscanf(argv[23], "%lf", &kappa);
  sscanf(argv[24], "%lf", &Lx_disk);
  sscanf(argv[25], "%d", &spot_disk);
  sscanf(argv[26], "%lf", &T_spot);
  sscanf(argv[27], "%lf", &spot_beg);
  sscanf(argv[28], "%lf", &spot_end);
  sscanf(argv[29], "%lf", &ns_theta);
  sscanf(argv[30], "%lf", &spot_rho_in);
  sscanf(argv[31], "%lf", &spot_rho_out);
  sscanf(argv[32], "%lf", &drd_phi);
  sscanf(argv[33], "%lf", &drd_theta);
  sscanf(argv[34], "%lf", &Lx_disk_2);
  sscanf(argv[35], "%lf", &Lx_iso);
  sscanf(argv[36], "%lf", &rho_in);
  sscanf(argv[37], "%lf", &A);
  sscanf(argv[38], "%lf", &uniform_disk);
  sscanf(argv[39], "%lf", &disk_flux);
  sscanf(argv[40], "%lf", &h_warp);
  sscanf(argv[41], "%s", &filter);
  sscanf(argv[42], "%lf", &lc_start);
  sscanf(argv[43], "%lf", &lc_end);

  
  /* h in the parameters list is the full width of the disk */
  h = h * 0.5; /* here h is the semiwidth of the disk */

  /* convert angles to radians */
  y_tilt = y_tilt * (M_PI/180.0);
  z_tilt = z_tilt * (M_PI/180.0);
  inclination = inclination * (M_PI/180.0);
  
  y_tilt2 = y_tilt2 * (M_PI/180.0);
  z_tilt2 = z_tilt2 * (M_PI/180.0);

  PSI_pr = 2.0 * M_PI - PSI_pr * (M_PI/180.0);
  kappa = kappa * (M_PI/180.0);

  spot_beg = spot_beg * (M_PI/180.0);
  spot_end = spot_end * (M_PI/180.0);

  drd_phi = drd_phi * (M_PI/180.0);
  drd_theta = drd_theta * (M_PI/180.0);
  
  
  double omega = omg(q, mu);  /* Dimentionless potential */

  disk d;

  vec3 d2;
  
  /* observer */
  /* double phi = 90.0 * (M_PI/180.0); */
  double phi;
  double flux_from_the_star;
  
  vec3 o;
  
  double phase[lc_num];
  double phase_array[lc_num];
  double flx[lc_num];
  double T_Lagrange_point[lc_num];

  //double disk_flx = flux_disk(o, d, y_tilt, z_tilt, omega, q, b, disk_tiles, phi, T_disk, lambda_cm, a_cm);
  //printf ("disk_flx %f", disk_flx);
  /* double star_flx = flux_star(o, q, omega, beta, u, d, d2, Lx, albedo, star_tiles, T_star, lambda_cm, a_cm); */

  sp neutron_star_sp;
  
  vec3 neutron_star;

  sp neutron_star_sp_greenwich;
  vec3 neutron_star_greenwich;

  int i; /* light curve step */

  vec3 w;

  double star;

  sp disk_reflection_diagr;

  vec3 drd_vec3;

  int steps = sqrt(star_tiles/2.0);
  int steps_phi = 2 * steps;
  int steps_theta = steps;

  double *plr;
  plr = polar(q, omega);
  double max_r = plr[2]; 
  free(plr);

  phi = 0.0 * M_PI/180.0;
  o.x = sin(inclination) * cos(phi);
  o.y = sin(inclination) * sin(phi);
  o.z = cos(inclination);

  vec3 DN;
  
  d.h.x = h * sin(y_tilt) * cos(z_tilt - phi);
  d.h.y = h * sin(y_tilt) * sin(z_tilt - phi);
  d.h.z = h * cos(y_tilt);
  d.R = R;
  
  DN.x = d.h.x;
  DN.y = d.h.y;
  DN.z = d.h.z;

  double DNL = len(DN);

  DN.x = d.h.x/DNL;
  DN.y = d.h.y/DNL;
  DN.z = d.h.z/DNL;

  d2.x = h * sin(y_tilt2) * cos(z_tilt2 - phi);
  d2.y = h * sin(y_tilt2) * sin(z_tilt2 - phi);
  d2.z = h * cos(y_tilt2);

  disk_reflection_diagr.phi = drd_phi;
  disk_reflection_diagr.theta = drd_theta;
  disk_reflection_diagr.r = 1.0;

  int diagram_size = 180*360;

  double diagram[diagram_size];


  FILE *file;
  file = fopen("DIAGRAM", "r");
  fread(&diagram, sizeof(double), diagram_size, file);
  fclose(file);

  int PSI_pr_i = (int) (360.0 * 180.0 * PSI_pr)/(360 * M_PI);

  double *Ix_dd;
  Ix_dd = &diagram[180 * PSI_pr_i + 0];

  //double F_disk_const = flux_disk(o, d, rho_in, A, uniform_disk, y_tilt, z_tilt, omega, q, disk_tiles, phi, T_disk, lambda_cm, a_cm, picture, spot_disk, T_spot, spot_beg, spot_end, spot_rho_in, spot_rho_out, h_warp, Ix_dd, Lx);

  //double F_disk_const = 0.0;


  double * phi_array;
  double * theta_array;

  double * r_array;
  double * g_array;

  double F_disk;

  phi_array = phi_func(steps_phi, threads);
  theta_array = theta_func(steps_theta, threads);
  
  r_array = shape_r(steps_phi, steps_theta, phi_array, theta_array, q, omega, threads);
  g_array = shape_g_abs(steps_phi, steps_theta, phi_array, theta_array, q, omega, threads);
    

  double disk1_normal_vector_x[lc_num]; 
  double disk1_normal_vector_y[lc_num];
  double disk1_normal_vector_z[lc_num];
	
  double disk2_normal_vector_x[lc_num];
  double disk2_normal_vector_y[lc_num];
  double disk2_normal_vector_z[lc_num];

  double observer_vector_x[lc_num];
  double observer_vector_y[lc_num];
  double observer_vector_z[lc_num];

  double epsilon1[lc_num];
  double epsilon2[lc_num];

  double disk_precession_over_phase = 0.0 * M_PI/180.0;
  
  double disk_delta_phase;
  double * disk_flux_evaluation;

  omp_set_dynamic(0);
  omp_set_num_threads(threads);

#pragma omp parallel for private(i, phi, o, d, d2, drd_vec3, disk_reflection_diagr, neutron_star_sp, neutron_star, star, flux_from_the_star, disk_delta_phase, disk_flux_evaluation)
  for(i = 0; i < lc_num; i++)
    {
      //phi = (double) i * 2.0 * M_PI/(lc_num - 1) - M_PI;

      phi = (double) i * 2.0 * M_PI * (lc_end - lc_start)/(lc_num - 1) - M_PI + 2.0 * M_PI * lc_start;
      
      //phi = phase_array[i];
      
      o.x = sin(inclination) * cos(phi);
      o.y = sin(inclination) * sin(phi);
      o.z = cos(inclination);
      
      disk_delta_phase = ((double) i / (lc_num - 1) - 0.5) * disk_precession_over_phase * (-1.0);

      /* "z_tilt - phi" because disk doesn`t rotatates with respect to observer */
      d.h.x = h * sin(y_tilt) * cos(z_tilt + phi + disk_delta_phase);
      d.h.y = h * sin(y_tilt) * sin(z_tilt + phi + disk_delta_phase);
      d.h.z = h * cos(y_tilt);
      d.R = R;

      d2.x = h * sin(y_tilt2) * cos(z_tilt2 + phi + disk_delta_phase);
      d2.y = h * sin(y_tilt2) * sin(z_tilt2 + phi + disk_delta_phase);
      d2.z = h * cos(y_tilt2);
      
      disk1_normal_vector_x[i] = sin(y_tilt) * cos(z_tilt + phi + disk_delta_phase);
      disk1_normal_vector_y[i] = sin(y_tilt) * sin(z_tilt + phi + disk_delta_phase);
      disk1_normal_vector_z[i] = cos(y_tilt);
	
      disk2_normal_vector_x[i] = sin(y_tilt2) * cos(z_tilt2 + phi + disk_delta_phase);
      disk2_normal_vector_y[i] = sin(y_tilt2) * sin(z_tilt2 + phi + disk_delta_phase);
      disk2_normal_vector_z[i] = cos(y_tilt2);

      observer_vector_x[i] = o.x;
      observer_vector_y[i] = o.y;
      observer_vector_z[i] = o.z;


      epsilon1[i] = 90.0 - (180.0/M_PI)*acos(disk1_normal_vector_x[i] * observer_vector_x[i] + disk1_normal_vector_y[i] * observer_vector_y[i] + disk1_normal_vector_z[i] * observer_vector_z[i]);
      epsilon2[i] = 90.0 - (180.0/M_PI)*acos(disk2_normal_vector_x[i] * observer_vector_x[i] + disk2_normal_vector_y[i] * observer_vector_y[i] + disk2_normal_vector_z[i] * observer_vector_z[i]);
      
      disk_reflection_diagr.phi = drd_phi;
      disk_reflection_diagr.theta = drd_theta;
      disk_reflection_diagr.r = 1.0;

      
      //printf("%f %f %f\n", drd_vec3.x, drd_vec3.y, drd_vec3.z);
      drd_vec3 = sp2dec(disk_reflection_diagr);

      drd_vec3 = rotate(drd_vec3, y_tilt, z_tilt + phi);

      //printf("%f %f %f\n", drd_vec3.x, drd_vec3.y, drd_vec3.z);
      
      disk_reflection_diagr = dec2sp(drd_vec3);

      //printf("%f %f %f\n", disk_reflection_diagr.phi, disk_reflection_diagr.theta, disk_reflection_diagr.r);

      neutron_star_sp_greenwich.phi = 0.0 * M_PI/180.0;
      neutron_star_sp_greenwich.theta = ns_theta * M_PI/180.0;

      
      neutron_star_sp.phi = 0.0 * M_PI/180.0;
      
      //neutron_star_sp.theta = 3.0 * M_PI/180.0;

      /* new angle !!! */
      neutron_star_sp.theta = ns_theta * M_PI/180.0;
      
      neutron_star_sp.r = 1.0;

      neutron_star = sp2dec(neutron_star_sp);

      neutron_star = rotate(neutron_star, 0.0, -phi);
      neutron_star = axrot(neutron_star, o, kappa);

      //printf("%f\t%f\t%f\n", neutron_star.x, neutron_star.y, neutron_star.z);

      
      star = flux_star(o, q, omega, beta, u, d, d2, Lx_noniso, Lx_disk, Lx_iso, Lx_disk_2, albedo, star_tiles, T_star, a_cm, neutron_star, PSI_pr, picture, disk_reflection_diagr, r_array, g_array, phi_array, theta_array, Ix_dd, y_tilt, y_tilt2, z_tilt+disk_delta_phase, z_tilt2+disk_delta_phase, phi, filter);

      flux_from_the_star = star;
      
      if (o.x < - cos(R + max_r) && picture == 0)
	{
	  F_disk = 0;
	}
      else if (o.x > - cos(R + max_r) && picture == 0)
	{ 
	  //F_disk = disk_flux;
	  F_disk = 0;

	}
      else if (picture == 1)
	{
	  //F_disk = disk_flux;
	  disk_flux_evaluation = flux_disk(o, d, rho_in, A, uniform_disk, y_tilt, z_tilt, omega, q, disk_tiles, phi, T_disk, a_cm, picture, spot_disk, T_spot, spot_beg, spot_end, spot_rho_in, spot_rho_out, h_warp, Ix_dd, Lx_noniso, neutron_star, d2, filter);
	}

      
      /* flx[i] = F_disk + flux_from_the_star; */
      flx[i] = flux_from_the_star;

      /* T_disk_from_Lx_up[i] = F_disk[1]; */
      /* T_disk_from_Lx_down[i] = F_disk[2]; */

      phase[i] = phi;
      //free(Ix_dd);



    }
  
  free(r_array);
  free(g_array);
  free(phi_array);
  free(theta_array);

  double min = flx[0]; /* searching minimum of the light curve */  
  for(i = 1; i < lc_num; i++)
    {
      if (flx[i] < min)
  	{
  	  min = flx[i];
  	}
    }

  
  /* double disk_flux; */
  /* /\* double lc_slope; *\/ */
  
  /* if (UBV_filter == 'B') */
  /*   { */
  /*     disk_flux = disk_flux_B; */
  /*     /\* lc_slope = lc_slope_B; *\/ */
  /*   } */
  /* else if (UBV_filter == 'V') */
  /*   { */
  /*     disk_flux = disk_flux_V; */
  /*     /\* lc_slope = lc_slope_V; *\/ */
  /*   } */
  /* else if (UBV_filter == 'W') */
  /*   { */
  /*     disk_flux = disk_flux_WASP; */
  /*     /\* lc_slope = lc_slope_V; *\/ */
  /*   } */
  /* else if (UBV_filter == '') */
  /*   { */
  /*     disk_flux = disk_flux_0 */
  /*     /\* lc_slope = lc_slope_V; *\/ */
  /*   } */

  
  /* printf ("%c %f\n", UBV_filter, disk_flux); */
  
  /* for(i = 0; i < lc_num; i++) */
  /*   { */
  /*     if (picture == 0) */
  /* 	{ */
  /* 	  if (((phase[i]/(2.0 * M_PI) + 0.5) >= 0.13 - eps) && ((phase[i]/(2.0 * M_PI) + 0.5) <= 0.87 + eps)) */
  /* 	    { */
  /* 	      printf("%.10f\t %.10f\n", phase[i]/(2.0 * M_PI), flx[i]/min + disk_flux + (phase[i]/(2.0 * M_PI))/tan(lc_slope)); */
  /* 	    } */
  /* 	  else if (((phase[i]/(2.0 * M_PI) + 0.5) < 0.13) || ((phase[i]/(2.0 * M_PI) + 0.5) > 0.87)) */
  /* 	    { */
  /* 	      printf("%.10f\t %.10f\n", phase[i]/(2.0 * M_PI), 1.0); */
  /* 	    } */
  /* 	} */
  /*     else if (picture == 1) */
  /* 	{} */
  /*   } */



  /* normalized flux */
  for(i = 0; i < lc_num; i++)
    {
      if (picture == 0)
  	{
  	  if (((phase[i]/(2.0 * M_PI) + 0.5) >= lc_start - eps) && ((phase[i]/(2.0 * M_PI) + 0.5) <= lc_end + eps))
  	    {
  	      printf("%.10f\t %.10f\n", phase[i]/(2.0 * M_PI), flx[i]/min + disk_flux);
  	    }
  	  else if (((phase[i]/(2.0 * M_PI) + 0.5) < lc_start) || ((phase[i]/(2.0 * M_PI) + 0.5) > lc_end))
  	    {
  	      printf("%.10f\t %.10f\n", phase[i]/(2.0 * M_PI), 1.0);
  	    }
  	}
      else if (picture == 1)
  	{}
    }



  /* /\* absolute flux *\/ */
  /* /\* there is no disk *\/ */
  /* for(i = 0; i < lc_num; i++) */
  /*   { */
  /*     if (picture == 0) */
  /* 	{ */
  /* 	  printf("%.10f\t %.10f\n", phase[i]/(2.0 * M_PI), flx[i]); */
  /* 	} */
  /*     else if (picture == 1) */
  /* 	{} */
  /*   } */


  
  return 0;
}
