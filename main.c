#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "def.h"

main()
{
  /* disk parameters */
  double h = 0.023;
  double R = 0.266;
  double y_tilt = 28.0 * (M_PI/180.0);
  double z_tilt = 60.0 * (M_PI/180.0);

  disk d;
  d.h.x = - h * sin(y_tilt) * cos(z_tilt);
  d.h.y = h * sin(y_tilt) * sin(z_tilt);
  d.h.z = h * cos(y_tilt);
  d.R = R;
  
  /* observer */
  double phi = 200.0 * (M_PI/180.0);
  /* double phi; */
  double inclination = 81.0 * (M_PI/180.0);

  vec3 o;
  o.x = sin(inclination) * cos(phi);
  o.y = sin(inclination) * sin(phi);
  o.z = cos(inclination);
  
  int i;
  int lc_num = 100;

  double phase[lc_num];
  double flx[lc_num];

  /* mass ratio q = m_x/m_opt */
  double q = 0.58;

  /* Roche lobe filling */
  double mu = 1.0;

  /* Dimentionless potential */
  double omega = omg(q, mu);
  
  /* double omega = 3.027; */

  double beta = 0.08;
  double u = 0.4;

  double f;
  
  double result_1 = flux_disk(o, d, y_tilt, z_tilt, omega, q);
  double result_2 = flux_star(o, q, omega, beta, u, d);

/*   omp_set_dynamic(0); */
/*   omp_set_num_threads(4); */

/* #pragma omp parallel for private(i, phi, o) */
/*   for(i = 0; i < lc_num; i++) */
/*     { */
/*       phi = (double) i * 2.0 * M_PI/(lc_num - 1) + M_PI; */
/*       o.x = sin(inclination) * cos(phi); */
/*       o.y = sin(inclination) * sin(phi); */
/*       o.z = cos(inclination); */

/*       phase[i] = phi; */

/*       flx[i] = flux_disk(o, d, y_tilt, z_tilt, omega, q) + flux_star(o, q, omega, beta, u, d); */
/*     } */

/*   FILE *file; */
/*   file = fopen( "LC", "w" ); */

/*   for(i = 0; i < lc_num; i++) */
/*     fprintf(file, "%.20f\t %.20f\n", phase[i]/(2.0 * M_PI), flx[i]); */
  
/*   fclose(file); */

/*   printf("*******STAR******\n"); */
/*   printf("beta.............%.10f\n", beta); */
/*   printf("u................%.10f\n", u); */
/*   printf("q................%.10f\n", q); */
/*   printf("Omega............%.10f\n", omega); */
/*   /\* printf("mu...............%.10f\n", mu); *\/ */
/*   printf("inclination......%.10f\n", 360.0 * inclination/(2.0 * M_PI)); */
/*   printf("*******DISC******\n"); */
/*   printf("R................%.10f\n", R); */
/*   printf("h................%.10f\n", h); */
/*   printf("y_tilt...........%.10f\n", y_tilt * (180.0/M_PI)); */
/*   printf("z_tilt...........%.10f\n", z_tilt * (180.0/M_PI)); */
/*   printf("***************\n"); */

}
