#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct decart {
  double x;
  double y;
  double z;
};

typedef struct decart vec3;


struct spherical {
  double phi;
  double theta;
  double r;
};

typedef struct spherical sp;


struct A {
  double b;
  double l;
  double r;
  double begin;
  double end;
  double delta_theta;
};

typedef struct A arc;


sp dec2sp(vec3 a)
{
  /* transforms cartesian coordinates to spherical */

  double x = a.x;
  double y = a.y;
  double z = a.z;

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

  sp result;

  result.phi = phi;
  result.theta = theta;
  result.r = r;

  return result;
}


vec3 sp2dec(sp a)
{
  vec3 result;

  result.x = a.r * sin(a.theta) * cos(a.phi);
  result.y = a.r * sin(a.theta) * sin(a.phi);
  result.z = a.r * cos(a.theta);

  return result;  
}


vec3 rotate(vec3 a, double Y, double Z)
{
  /* rotates counterclockwise around Y and Z axes */
  
  vec3 result;

  result.x = a.x * cos(Y) * cos(Z) + a.y * sin(Z) - a.z * sin(Y) * cos(Z);
  result.y = - a.x * cos(Y) * sin(Z) + a.y * cos(Z) + a.z * sin(Y) * sin(Z);
  result.z = a.x * sin(Y) + a.z * cos(Y); 

  return result;
}


sp arcgen(double theta, double delta_theta, double N_theta, double N_phi, int i)
{
  /* i is the index of points */

  double epsilon_theta = delta_theta/N_theta;
  double epsilon_phi = (2.0 * M_PI)/N_phi;

  sp result;
  
  result.phi = epsilon_phi * (i % (int) (N_phi + 1.0));
  result.theta = theta - 0.5 * delta_theta + 0.5 * epsilon_theta + epsilon_theta * (i % (int) (N_theta - 1.0));
  result.r = 1.0;
    
  return result;
}



int main()
{

  double PSI_pr = 90.0 * M_PI/180.0 ; /* NS precession angle */
  double angle_JI = 50.0 * M_PI/180.0; /* angle between J and I_3 vectors */

  double N_theta = 20.0; /* arc pixelization */
  double N_phi = 300.0; /* arc pixelization */

  double N = N_theta * N_phi;

  int arc_n; /* arc index */
  int i, f; /* point index, directional diagram index */
  vec3 decrt;
  sp sphr;

  double F[180]; /* directional diagram */

  double b[15], l[15], R[15], BEG[15], END[15], delta_theta[15], intensity[15];
  
  b[0] = 70.0,    l[0] = 160.0,  R[0] = 50.0,   BEG[0] = 205.0,  END[0] = 250.0,  delta_theta[0] = 14.0,  intensity[0] = 54.3;  
  b[1] = 75.0,    l[1] = 180.0,  R[1] = 43.0,   BEG[1] = 250.0,  END[1] = 300.0,  delta_theta[1] = 13.0,  intensity[1] = 53.1;  
  b[2] = 65.0,    l[2] = 180.0,  R[2] = 43.0,   BEG[2] = 302.0,  END[2] = 350.0,  delta_theta[2] = 12.0,  intensity[2] = 51.5;  
  b[3] = 65.0,    l[3] = 170.0,  R[3] = 53.0,   BEG[3] = 355.0,  END[3] = 400.0,  delta_theta[3] = 11.0,  intensity[3] = 19.3;  
  b[4] = 60.0,    l[4] = 160.0,  R[4] = 55.0,   BEG[4] = 45.0,   END[4] = 60.0,   delta_theta[4] = 11.0,  intensity[4] = 4.5;  
  b[5] = 60.0,    l[5] = 150.0,  R[5] = 50.0,   BEG[5] = 65.0,   END[5] = 95.0,   delta_theta[5] = 11.0,  intensity[5] = 20.2;  
  b[6] = 58.0,    l[6] = 160.0,  R[6] = 55.0,   BEG[6] = 95.0,   END[6] = 120.0,  delta_theta[6] = 11.0,  intensity[6] = 25.8;  
  b[7] = 70.0,    l[7] = 140.0,  R[7] = 60.0,   BEG[7] = 120.0,  END[7] = 160.0,  delta_theta[7] = 13.0,  intensity[7] = 29.0;  
  b[8] = -97.0,   l[8] = 180.0,  R[8] = 44.0,   BEG[8] = 330.0,  END[8] = 350.0,  delta_theta[8] = 11.0,  intensity[8] = 4.1;  
  b[9] = -97.0,   l[9] = 180.0,  R[9] = 44.0,   BEG[9] = 20.0,   END[9] = 50.0,   delta_theta[9] = 11.0,  intensity[9] = 7.6;  
  b[10] = -97.0,  l[10] = 180.0, R[10] = 44.0,  BEG[10] = 50.0,  END[10] = 100.0, delta_theta[10] = 11.0, intensity[10] = 10.2;  
  b[11] = -95.0,  l[11] = 200.0, R[11] = 42.0,  BEG[11] = 110.0, END[11] = 170.0, delta_theta[11] = 8.0,  intensity[11] = 17.0;  
  b[12] = -100.0, l[12] = 200.0, R[12] = 65.0,  BEG[12] = 205.0, END[12] = 260.0, delta_theta[12] = 11.0, intensity[12] = 3.6;  
  b[13] = 60.0,   l[13] = 180.0, R[13] = 80.0,  BEG[13] = 250.0, END[13] = 280.0, delta_theta[13] = 11.0, intensity[13] = 11;  
  b[14] = 75.0,   l[14] = 180.0, R[14] = 105.0, BEG[14] = 90.0,  END[14] = 120.0, delta_theta[14] = 11.0, intensity[14] = 2.7;  

  double theta_n, delta_theta_n, epsilon_theta_n;
  double N_phi_n, epsilon_phi_n;
  double ds; /* surface element */

  for (arc_n = 0; arc_n < 15; arc_n ++)
    {
  
      theta_n = R[arc_n] * M_PI/180.0;
      delta_theta_n = delta_theta[arc_n] * M_PI/180.0;
      
      epsilon_theta_n = delta_theta_n/N_theta;
      epsilon_phi_n = (2.0 * M_PI)/N_phi;

      for(i = 0; i <= N; i++)
	{
	  sphr = arcgen(theta_n, delta_theta_n, N_theta, N_phi, i);
	  
	  //printf("%.20f\t %.20f\t 1.00000\n", sphr.phi, sphr.theta);

	  ds = epsilon_theta_n * epsilon_phi_n * sin(sphr.theta); /* surface element, r = 1.0 */
	  
	  decrt = sp2dec(sphr);

	  //printf("%.20f\t %.20f\t %.20f\n", decrt.x, decrt.y, decrt.z);

	  decrt = rotate(decrt, (-(90.0 - b[arc_n]) * M_PI/180.0), (-l[arc_n] * M_PI/180.0));
	  
	  sphr = dec2sp(decrt);

	  //printf("%.20f\t %.20f\t %.20f\n", decrt.x, decrt.y, decrt.z);

	  if ( END[arc_n] <= 360.0 && sphr.phi >= (BEG[arc_n] * M_PI/180.0) && sphr.phi <= (END[arc_n] * M_PI/180.0) )
	    {
	      decrt = sp2dec(sphr);
	      
	      //printf("%.20f\t %.20f\t %.20f\n", decrt.x, decrt.y, decrt.z);

	      decrt = rotate(decrt, 0.0, PSI_pr);
	      decrt = rotate(decrt, angle_JI, 0.0);

	      printf("%.20f\t %.20f\t %.20f\n", decrt.x, decrt.y, decrt.z);

	      
	      sphr = dec2sp(decrt);
	      
	      f = (int) floor( sphr.theta * 180.0/M_PI );
	      
	      F[f] = F[f] + ds;
	      //printf("%.20f\t %.20f\t %.20f\n", decrt.x, decrt.y, decrt.z);
	    }
	  else if ( END[arc_n] > 360.0 && sphr.phi >= (BEG[arc_n] * M_PI/180.0) || sphr.phi <= (END[arc_n] * M_PI/180.0 - 2.0 * M_PI) )
	    {
	      decrt = sp2dec(sphr);
	      
	      //printf("%.20f\t %.20f\t %.20f\n", decrt.x, decrt.y, decrt.z);

	      decrt = rotate(decrt, 0.0, PSI_pr);
	      decrt = rotate(decrt, angle_JI, 0.0);

	      printf("%.20f\t %.20f\t %.20f\n", decrt.x, decrt.y, decrt.z);
	      
	      sphr = dec2sp(decrt);
	      
	      f = (int) floor( sphr.theta * 180.0/M_PI );
	      
	      F[f] = F[f] + ds;

	    }

	}
      printf("\n\n");
      printf("0.00000\t 0.00000\t 1.00000\n0.00000\t 0.00000\t -1.00000\n");
      printf("0.00000\t 1.00000\t 0.00000\n0.00000\t -1.00000\t 0.00000\n");
      printf("1.00000\t 0.00000\t 0.00000\n-1.00000\t 0.00000\t 0.00000\n");
    }
  
  int j;
  for (j = 0; j < 180; j++)
    {
      //printf("%i %.20f\n", j, F[j]);
    }
  
  return 0;
}
