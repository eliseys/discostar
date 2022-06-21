#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>
#include "def.h"

int main(int argc, char **argv)
{

  parameters parameters;
  
  /**/
  
  sscanf(argv[1], "%lf", &parameters.q);
  sscanf(argv[2], "%lf", &parameters.mu);
  sscanf(argv[3], "%lf", &parameters.beta);
  sscanf(argv[4], "%lf", &parameters.u);
  sscanf(argv[5], "%lf", &parameters.X_albedo);  
  sscanf(argv[6], "%lf", &parameters.T_star_polar);  
  sscanf(argv[7], "%lf", &parameters.a); 
  sscanf(argv[8], "%lf", &parameters.inclination); 
  sscanf(argv[9], "%lf", &parameters.Lx); 
  sscanf(argv[10], "%lf", &parameters.NS_phi);
  sscanf(argv[11], "%lf", &parameters.NS_kappa);
  sscanf(argv[12], "%lf", &parameters.NS_theta);
  sscanf(argv[13], "%lf", &parameters.h_out);
  sscanf(argv[14], "%lf", &parameters.r_out);
  sscanf(argv[15], "%lf", &parameters.r_in);
  sscanf(argv[16], "%lf", &parameters.gamma);
  sscanf(argv[17], "%lf", &parameters.theta_out);
  sscanf(argv[18], "%lf", &parameters.phi_out);
  sscanf(argv[19], "%lf", &parameters.theta_in);
  sscanf(argv[20], "%lf", &parameters.phi_in);
  sscanf(argv[21], "%d",  &parameters.do_lc);
  sscanf(argv[22], "%d", &parameters.N_lc);
  sscanf(argv[23], "%d", &parameters.N_theta);
  sscanf(argv[24], "%d", &parameters.N_r);
  sscanf(argv[25], "%d", &parameters.OMP_threads);
  sscanf(argv[26], "%d",  &parameters.do_corona);
  sscanf(argv[27], "%d",  &parameters.N_corona);
  sscanf(argv[28], "%lf", &parameters.h_corona);
  sscanf(argv[29], "%d",  &parameters.rs_corona);
  sscanf(argv[30], "%c", &parameters.filter);

  
  parameters.omega = omg(parameters.q, parameters.mu);  /* Dimentionless potential */

  double phi;
  
  vec3 observer;

  disk disk;

  
  double orbital_phase[parameters.N_lc];
  double F[parameters.N_lc];
  double corona_Fx[parameters.N_lc];

  /******/
  /**/
  sp neutron_star_sp;
  vec3 neutron_star;
  /**/
  /******/

  int diagram_size = 180*360;
  double diagram[diagram_size];
  
  double * star_elements = star_geometry(parameters);

  //double * disk_elements;

  double * corona_elements;

  
  if (parameters.do_corona)
    {
      corona_elements = rand_p(parameters);
    }
  
  FILE * file;
  file = fopen("DIAGRAM", "r");
  

  fread(&diagram, sizeof(double), diagram_size, file);
  fclose(file);


  
  int NS_phi_i = (int) (360.0 * 180.0 * parameters.NS_phi)/(360 * M_PI);

  double * Ix_dd;
  Ix_dd = &diagram[180 * NS_phi_i + 0];

  
  int i;
  
  omp_set_dynamic(0);
  omp_set_num_threads(parameters.OMP_threads);


  
  if (parameters.do_lc)
    {
#pragma omp parallel for private(i, phi, disk, observer, neutron_star)
      for(i = 0; i < parameters.N_lc; i++)
	{ 
	  phi = ((double) i / (parameters.N_lc - 1) - 0.5) * 2.0 * M_PI;
	  orbital_phase[i] = phi;


	  disk.R = parameters.r_out;
	  disk.r_out = disk.R;
	  disk.r_in = parameters.r_in;
	  disk.h = parameters.h_out;
	  disk.gamma = parameters.gamma;


	  disk.theta_out = parameters.theta_out;
	  disk.theta_in = parameters.theta_in;

	  disk.phi_out = parameters.phi_out + phi;
	  disk.phi_in = parameters.phi_in + phi;
  
	  disk.n.x = 0.0;
	  disk.n.y = 0.0;
	  disk.n.z = 1.0;

	  disk.n = R_y(disk.n, disk.theta_out);
	  disk.n = R_z(disk.n, disk.phi_out);

	  double * disk_elements;
	  disk_elements = disk_geometry(disk, parameters);

	  vec3 observer;
	  observer.x = 0.0;
	  observer.y = 0.0;
	  observer.z = 1.0;

	  observer = R_y(observer, parameters.inclination);
	  observer = R_z(observer, phi);

	  neutron_star.x = 0.0;
	  neutron_star.y = 0.0;
	  neutron_star.z = 1.0;	  
	  neutron_star = R_y(neutron_star, - (M_PI/2 - parameters.inclination) + parameters.NS_theta);
	  neutron_star = axrot(neutron_star, observer, parameters.NS_kappa);
	  neutron_star = R_z(neutron_star, phi);
	  
	  F[i] = 
	    star_F(star_elements, parameters, disk, observer, neutron_star, Ix_dd) + 
	     disk_F(disk_elements, parameters, disk, observer); 

	  //F[i] = star_F(star_elements, parameters, disk, observer, neutron_star, Ix_dd);
	  
	  //if (parameters.do_corona && phi == 0.0)
	  //  {
	  //   corona_Fx[i] = x_ray_corona(parameters, disk, observer, corona_elements);
	  //  }

	  free(disk_elements);

	  //printf("%d\t%f\n", i, F[i]);

	}
    }

  free(star_elements);
  //free(disk_elements);
  
  if (parameters.do_corona)
    {
      free(corona_elements);
    }

  //free(Ix_dd);


  //exit(0);

  
  phi = 0.0 * (M_PI/180.0);

  /* Disk setup */
  disk.R = parameters.r_out;
  disk.r_out = disk.R;
  disk.r_in = parameters.r_in;
  disk.h = parameters.h_out;
  disk.gamma = parameters.gamma;

  disk.theta_out = parameters.theta_out;
  disk.theta_in = parameters.theta_in;

  disk.phi_out = parameters.phi_out + phi;
  disk.phi_in = parameters.phi_in + phi;
  
  disk.n.x = 0.0;
  disk.n.y = 0.0;
  disk.n.z = 1.0;

  disk.n = R_y(disk.n, disk.theta_out);
  disk.n = R_z(disk.n, disk.phi_out);

  /* Observer setup  */
  observer.x = 0.0;
  observer.y = 0.0;
  observer.z = 1.0;

  observer = R_y(observer, parameters.inclination);

  observer = R_z(observer, phi);


  /* double qq, aa, E; */
  
  /* for (i = 0; i<990; i++) */
  /*   { */
  /*     E = 0.1 + (double) i*9.9/(990.0-1.0); // energy in keV */
  /*     qq = star_X(star_elements, parameters, disk, observer, neutron_star, Ix_dd, E); */
  /*     printf("%f %f\n", E, qq); */
  /*   } */


  /* for (i = 0; i<91; i++) */
  /*   { */
  /*     E = 1.0 + (double) i*9.0/(91.0-1.0); // energy in keV */
  /*     aa = star_X(star_elements, parameters, disk, observer, neutron_star, Ix_dd, E); */
  /*     printf("%f %f\n", E, aa); */
  /*   } */

  
  /* E = 10.0; // energy in keV */
  /* qq = star_X(star_elements, parameters, disk, observer, neutron_star, Ix_dd, E); */

  /* printf("qq %f\n", qq); */

  double min = F[0]; /* searching minimum of the light curve */
  for(i = 1; i < parameters.N_lc; i++)
    {
      if (F[i] < min)
  	{
  	  min = F[i];
  	}
    }

  double disk_add = 0.0;


  double max = F[0]; /* searching maximum of the light curve */
  for(i = 1; i < parameters.N_lc; i++)
    {
      if (F[i] > max)
  	{
  	  max = F[i];
  	}
    }

  
  if (parameters.do_lc)
    { 
      for(i = 0; i < parameters.N_lc; i++)
  	{
  	  printf("%.10f\t %.10f\n", orbital_phase[i]/(2.0 * M_PI) + 0.5, F[i]/max);
  	}
    }

  
  /* if (parameters.do_lc && parameters.do_corona) */
  /*   { */
  /*     printf("\n\n"); */
  /*     for(i = 0; i < parameters.N_lc; i++) */
  /* 	{ */
  /* 	  printf("%.10f\t %.10f\n", orbital_phase[i]/(2.0 * M_PI), corona_Fx[i]/parameters.N_corona); */
  /* 	} */
  /*   } */


  /* if (!parameters.do_lc) */
  /*   { */
  /*     phi = 0.0 * (360.0/50) * (M_PI/180.0); */
  
  /*     disk.R = parameters.r_out; */
  /*     disk.r_out = disk.R; */
  /*     disk.r_in = parameters.r_in; */
  /*     disk.h = parameters.h_out; */
  /*     disk.gamma = parameters.gamma; */


  /*     disk.theta_out = parameters.theta_out; */
  /*     disk.theta_in = parameters.theta_in; */

  /*     disk.phi_out = parameters.phi_out + phi; */
  /*     disk.phi_in = parameters.phi_in + phi; */
  
  /*     disk.n.x = 0.0; */
  /*     disk.n.y = 0.0; */
  /*     disk.n.z = 1.0; */

  /*     disk.n = R_y(disk.n, disk.theta_out); */
  /*     disk.n = R_z(disk.n, disk.phi_out); */
	  
  /*     double * disk_elements; */
  /*     disk_elements = disk_geometry(disk, parameters); */

  /*     observer.x = 0.0; */
  /*     observer.y = 0.0; */
  /*     observer.z = 1.0; */

  /*     observer = R_y(observer, parameters.inclination); */
  /*     observer = R_z(observer, phi); */
  
  /*     star_F(star_elements, parameters, disk, observer, neutron_star, Ix_dd); */
  /*     disk_F(disk_elements, parameters, disk, observer); */
  /*   } */
  //printf(">>>> TEST\n");


  

  /* int M = 201; */

  /* double theta, r, x, y, z; */
  
  /* for (i = 0; i<M; i++) */
  /*   { */
  /*     phi = 3.0*M_PI/2.0; */
  /*     theta = (double) i/(M-1) * M_PI; */
      
  /*     r = radius_star(phi, theta, parameters.q, parameters.omega); */

  /*     z = r * cos(theta); */
  /*     x = - r * sin(theta); */


  /*     printf("%f\t%f\n", x, z); */


  /*   } */

  /* printf("\n\n"); */

  /* for (i = M-1; i>=0; i--) */
  /*   {       */
  /*     phi = 0.0; */
  /*     theta = (double) i/(M-1) * M_PI; */

  /*     r = radius_star(phi, theta, parameters.q, parameters.omega); */

  /*     z = r * cos(theta); */
  /*     x = r * sin(theta); */

  /*     printf("%f\t%f\n", x, z); */


  /*   } */

  /* double x, y, z; */

  /* double * corona_points = rand_p(parameters); */


  //double summa = x_ray_corona(parameters, disk, observer, corona_elements);

  //printf("summa %f\n", summa);
  
  /* for (i = 0; i<parameters.N_corona; i++) */
  /*   { */
      
  /*     x = corona_points[3*i+1]; */
  /*     y = corona_points[3*i+2]; */
  /*     z = corona_points[3*i+3]; */

  /*     printf("%f\t%f\t%f\n", x, y, z); */

      
  /*   } */


  // printf("\n\n");

  /* int N_d = 50; */

  /* double theta_out = 10.0 * (M_PI/180.0); */

  /* double r_f, z_f; */
  
  /* for (i = 1; i<N_d; i++) */
  /*   { */

  /*     r = (double) i/(N_d - 1) * parameters.r_out; */

  /*     /\* printf("parameters.r_out %f\n", parameters.r_out); *\/ */
  /*     /\* printf("r %f\n", r); *\/ */
      
  /*     z = parameters.h_out * pow(((r - parameters.r_in)/(parameters.r_out-parameters.r_in)), parameters.gamma); */

      
  /*     r_f = r * cos(theta_out) - z * sin(theta_out); */
  /*     z_f = r * sin(theta_out) + z * cos(theta_out); */

  /*     printf("%f\t%f\n", r_f + 1.0, z_f); */
  /*     //printf("%f\t%f\n", r_f + 1.0, -z_f); */
  /*   } */
  
  /* for (i = N_d-1; i>=0; i--) */
  /*   { */

  /*     r = (double) i/(N_d - 1) * parameters.r_out; */

  /*     z = - parameters.h_out * pow(((r - parameters.r_in)/(parameters.r_out-parameters.r_in)), parameters.gamma); */

  /*     r_f = r * cos(theta_out) - z * sin(theta_out); */
  /*     z_f = r * sin(theta_out) + z * cos(theta_out); */

  /*     printf("%f\t%f\n", r_f + 1.0, z_f); */

  /*   } */


  /* for (i = 1; i<N_d; i++) */
  /*   { */

  /*     r = (double) i/(N_d - 1) * parameters.r_out; */

  /*     /\* printf("parameters.r_out %f\n", parameters.r_out); *\/ */
  /*     /\* printf("r %f\n", r); *\/ */
      
  /*     z = parameters.h_out * pow(((r - parameters.r_in)/(parameters.r_out-parameters.r_in)), parameters.gamma); */
 
  /*     r_f = -r * cos(theta_out) - z * sin(theta_out); */
  /*     z_f = -r * sin(theta_out) + z * cos(theta_out); */

  /*     printf("%f\t%f\n", r_f + 1.0, z_f); */
  /*   } */

  /* for (i = N_d-1; i>=0; i--) */
  /*   { */

  /*     r = (double) i/(N_d - 1) * parameters.r_out; */

  /*     z = - parameters.h_out * pow(((r - parameters.r_in)/(parameters.r_out-parameters.r_in)), parameters.gamma); */

  /*     r_f = -r * cos(theta_out) - z * sin(theta_out); */
  /*     z_f = -r * sin(theta_out) + z * cos(theta_out); */

  /*     printf("%f\t%f\n", r_f + 1.0, z_f); */


  /*   } */




  

  /* int N_d2 = (parameters.h_out/parameters.r_out) * N_d; */


  
  /* for (i = 0; i<N_d2; i++) */
  /*   { */
  /*     r = parameters.r_out; */

  /*     z = (double) i/(N_d2 - 1) * parameters.h_out; */

  /*     r_f = r * cos(theta_out) - z * sin(theta_out); */
  /*     z_f = r * sin(theta_out) + z * cos(theta_out); */

  /*     printf("%f\t%f\n", r_f + 1.0, z_f); */
  /*     printf("%f\t%f\n", -r_f + 1.0, -z_f); */

  /*     z = - (double) i/(N_d2 - 1) * parameters.h_out; */

  /*     r_f = r * cos(theta_out) - z * sin(theta_out); */
  /*     z_f = r * sin(theta_out) + z * cos(theta_out); */

  /*     printf("%f\t%f\n", r_f + 1.0, z_f); */
  /*     printf("%f\t%f\n", -r_f + 1.0, -z_f); */

     
           
  /*   } */

  //ion_init();
  
  return 0;

  
}
