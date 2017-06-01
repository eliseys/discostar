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

  double b;      /* reserve parameter */ 

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
  sscanf(argv[11], "%lf", &b);
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


  /**/

  double A_cm = 1E-8;
  double lambda_cm = lambda_A * A_cm;
  
  /* convert angles to radians */
  y_tilt = y_tilt * (M_PI/180.0);
  z_tilt = z_tilt * (M_PI/180.0);
  inclination = inclination * (M_PI/180.0);
  
  y_tilt2 = y_tilt2 * (M_PI/180.0);
  z_tilt2 = z_tilt2 * (M_PI/180.0);

  double omega = omg(q, mu);  /* Dimentionless potential */

  disk d;

  vec3 d2;
  
  /* observer */
  /* double phi = 90.0 * (M_PI/180.0); */
  double phi;

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
  
  /* double disk_flx = flux_disk(o, d, y_tilt, z_tilt, omega, q, b, disk_tiles, phi, T_disk, lambda_cm, a_cm); */
  /* double star_flx = flux_star(o, q, omega, beta, u, d, d2, Lx, albedo, star_tiles, T_star, lambda_cm, a_cm); */


  int i; /* light curve step */

  omp_set_dynamic(0);
  omp_set_num_threads(threads);


#pragma omp parallel for private(i, phi, o, d, d2)
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
      
      phase[i] = phi;

      flx[i] = flux_disk(o, d, y_tilt, z_tilt, omega, q, disk_tiles, phi, T_disk, lambda_cm, a_cm) + flux_star(o, q, omega, beta, u, d, d2, Lx, albedo, star_tiles, T_star, lambda_cm, a_cm);

    }

  double min = flx[0]; /* searching minimum of the light curve */
    
  for(i = 1; i < lc_num; i++)
    {
      if (flx[i] < min)
  	{
  	  min = flx[i];
  	}
    }

  for(i = 0; i < lc_num; i++)
    {
      printf("%.20f\t %.20f\n", phase[i]/(2.0 * M_PI), flx[i]/min);
    }
  
  /* for(i = 0; i < lc_num; i++) */
  /*   { */
  /*     /\* 2nd phase just copy of the first one *\/ */
  /*     printf("%.20f\t %.20f\n", phase[i]/(2.0 * M_PI) + 1.0, flx[i]/min); */
  /*   } */

  return 0;
}
