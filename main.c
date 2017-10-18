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
  double Lx;     /* in units of L of the optical star */ 

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
  double lambda_A;
  double a_cm;

  double PSI_pr;
  double kappa;

  int isotrope;

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
  
  /**/
  
  sscanf(argv[1], "%lf", &q);
  sscanf(argv[2], "%lf", &mu);
  sscanf(argv[3], "%lf", &beta);
  sscanf(argv[4], "%lf", &u);
  sscanf(argv[5], "%lf", &albedo);  
  sscanf(argv[6], "%lf", &Lx);  
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
  sscanf(argv[19], "%lf", &lambda_A);
  sscanf(argv[20], "%lf", &a_cm);
  sscanf(argv[21], "%lf", &y_tilt2);
  sscanf(argv[22], "%lf", &z_tilt2);
  sscanf(argv[23], "%lf", &PSI_pr);
  sscanf(argv[24], "%lf", &kappa);
  sscanf(argv[25], "%d", &isotrope);
  sscanf(argv[26], "%lf", &Lx_disk);
  sscanf(argv[27], "%d", &spot_disk);
  sscanf(argv[28], "%lf", &T_spot);
  sscanf(argv[29], "%lf", &spot_beg);
  sscanf(argv[30], "%lf", &spot_end);
  sscanf(argv[31], "%lf", &ns_theta);
  sscanf(argv[32], "%lf", &spot_rho_in);
  sscanf(argv[33], "%lf", &spot_rho_out);
  sscanf(argv[34], "%lf", &drd_phi);
  sscanf(argv[35], "%lf", &drd_theta);
  sscanf(argv[36], "%lf", &Lx_disk_2);
  sscanf(argv[37], "%lf", &Lx_iso);



  /**/

  double A_cm = 1E-8;
  double lambda_cm = lambda_A * A_cm;


  /* h in the parameters list is the full width of the disk */
  h = h * 0.5; /* here h is the semiwidth of the disk */
  
  /* convert angles to radians */
  y_tilt = y_tilt * (M_PI/180.0);
  z_tilt = z_tilt * (M_PI/180.0);
  inclination = inclination * (M_PI/180.0);
  
  y_tilt2 = y_tilt2 * (M_PI/180.0);
  z_tilt2 = z_tilt2 * (M_PI/180.0);

  PSI_pr = PSI_pr * (M_PI/180.0);
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
  /* o.x = sin(inclination) * cos(phi); */
  /* o.y = sin(inclination) * sin(phi); */
  /* o.z = cos(inclination); */


  /* "z_tilt - phi" because disk doesn`t rotatate with respect to observer */
  /* d.h.x = - h * sin(y_tilt) * cos(z_tilt - phi + M_PI); */
  /* d.h.y = h * sin(y_tilt) * sin(z_tilt - phi + M_PI); */
  /* d.h.z = h * cos(y_tilt); */
  /* d.R = R; */

  /* d2.x = - h * sin(y_tilt2) * cos(z_tilt + z_tilt2 - phi + M_PI); */
  /* d2.y = h * sin(y_tilt2) * sin(z_tilt + z_tilt2 - phi + M_PI); */
  /* d2.z = h * cos(y_tilt2); */
  
  double phase[lc_num];
  double flx[lc_num];
  double T_Lagrange_point[lc_num];

  /* double disk_flx = flux_disk(o, d, y_tilt, z_tilt, omega, q, b, disk_tiles, phi, T_disk, lambda_cm, a_cm); */
  /* double star_flx = flux_star(o, q, omega, beta, u, d, d2, Lx, albedo, star_tiles, T_star, lambda_cm, a_cm); */

  sp neutron_star_sp;
  
  vec3 neutron_star;
  
  int i; /* light curve step */

  vec3 w;

  double * star;

  sp disk_reflection_diagr;

  vec3 drd_vec3;
  
  omp_set_dynamic(0);
  omp_set_num_threads(threads);


#pragma omp parallel for private(i, phi, o, d, d2, drd_vec3, disk_reflection_diagr, neutron_star_sp, neutron_star, star, flux_from_the_star)
  for(i = 0; i < lc_num; i++)
    {
      phi = (double) i * 2.0 * M_PI/(lc_num - 1) - M_PI;

      o.x = sin(inclination) * cos(phi);
      o.y = sin(inclination) * sin(phi);
      o.z = cos(inclination);

      /* "z_tilt - phi" because disk doesn`t rotatate with respect to observer */
      d.h.x = - h * sin(y_tilt) * cos(z_tilt - phi + M_PI);
      d.h.y = h * sin(y_tilt) * sin(z_tilt - phi + M_PI);
      d.h.z = h * cos(y_tilt);
      d.R = R;

      d2.x = - h * sin(y_tilt2) * cos(z_tilt + z_tilt2 - phi + M_PI);
      d2.y = h * sin(y_tilt2) * sin(z_tilt + z_tilt2 - phi + M_PI);
      d2.z = h * cos(y_tilt2);

      disk_reflection_diagr.phi = drd_phi;
      disk_reflection_diagr.theta = drd_theta;
      disk_reflection_diagr.r = 1.0;
      
      //printf("%f %f %f\n", drd_vec3.x, drd_vec3.y, drd_vec3.z);
      drd_vec3 = sp2dec(disk_reflection_diagr);

      drd_vec3 = rotate(drd_vec3, y_tilt, z_tilt - phi + M_PI);

      //printf("%f %f %f\n", drd_vec3.x, drd_vec3.y, drd_vec3.z);
      
      disk_reflection_diagr = dec2sp(drd_vec3);

      //printf("%f %f %f\n", disk_reflection_diagr.phi, disk_reflection_diagr.theta, disk_reflection_diagr.r);

      neutron_star_sp.phi = 0.0 * M_PI/180.0;
      //neutron_star_sp.theta = 3.0 * M_PI/180.0;

      /* new angle !!! */
      neutron_star_sp.theta = ns_theta * M_PI/180.0;
      
      neutron_star_sp.r = 1.0;

      neutron_star = sp2dec(neutron_star_sp);      
      neutron_star = rotate(neutron_star, 0.0, -phi);

      neutron_star = axrot(neutron_star, o, kappa);

      phase[i] = phi;

      star = flux_star(o, q, omega, beta, u, d, d2, Lx, Lx_disk, Lx_disk_2, Lx_iso, albedo, star_tiles, T_star, lambda_cm, a_cm, neutron_star, PSI_pr, picture, isotrope, disk_reflection_diagr);
      flux_from_the_star = star[0];
      T_Lagrange_point[i] = star[1];


      flx[i] = flux_disk(o, d, y_tilt, z_tilt, omega, q, disk_tiles, phi, T_disk, lambda_cm, a_cm, picture, spot_disk, T_spot, spot_beg, spot_end, spot_rho_in, spot_rho_out) + flux_from_the_star;

    }
  
  free(star);
  
  double min = flx[0]; /* searching minimum of the light curve */
  //double max = T_Lagrange_point[0];
  //printf("FIRST ASSIGMENT %f\n", max);

  
  for(i = 1; i < lc_num; i++)
    {
      if (flx[i] < min)
  	{
  	  min = flx[i];
  	}

    }
  
  for(i = 0; i < lc_num; i++)
    {
      if (picture == 0)
  	{
  	  printf("%.20f\t %.20f\t %.20f\n", phase[i]/(2.0 * M_PI), flx[i]/min, T_Lagrange_point[i]/8400.0 - 1.0);
  	}
      else if (picture == 1)
  	{}
    }

  /* for(i = 0; i < lc_num; i++) */
  /*   { */
  /*     /\* 2nd phase just copy of the first one *\/ */
  /*     printf("%.20f\t %.20f\n", phase[i]/(2.0 * M_PI) + 1.0, flx[i]/min); */
  /*   } */

  return 0;
}
