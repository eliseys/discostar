#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "def.h"



double * X_map_2()
{

  /* sp I_sp; */
  /* I_sp.phi = 180.0 * (M_PI/180.0); */
  /* I_sp.theta = 50.0 * (M_PI/180.0); */
  /* I_sp.r = 1.0;  */

  /* vec3 I = sp2dec(I_sp); // NS's inertia axis */

  
  double b[18], l[18], R[18], BEG[18], END[18], delta_theta[18], intensity[18], N[18];

  /* ORIGINAL */
  /* !!!!!! b[2] changed: 65.0 to 75.0 */
  /* !!!!!! BEG[8] changed: 330.0 to 350.0 */
  /* !!!!!! END[8] changed: 350.0 to 380.0 */
  
  /* b[ 0] =  70.0, l[ 0] = 160.0, R[ 0] =  50.0, BEG[ 0] = 205.0, END[ 0] = 250.0, delta_theta[ 0] = 14.0, intensity[ 0] = 54.3; */
  /* b[ 1] =  75.0, l[ 1] = 180.0, R[ 1] =  43.0, BEG[ 1] = 250.0, END[ 1] = 300.0, delta_theta[ 1] = 13.0, intensity[ 1] = 53.1; */
  /* b[ 2] =  75.0, l[ 2] = 180.0, R[ 2] =  43.0, BEG[ 2] = 302.0, END[ 2] = 350.0, delta_theta[ 2] = 12.0, intensity[ 2] = 51.5; */
  /* b[ 3] =  65.0, l[ 3] = 170.0, R[ 3] =  53.0, BEG[ 3] = 355.0, END[ 3] = 400.0, delta_theta[ 3] = 11.0, intensity[ 3] = 19.3; */
  /* b[ 4] =  60.0, l[ 4] = 160.0, R[ 4] =  55.0, BEG[ 4] =  45.0, END[ 4] =  60.0, delta_theta[ 4] = 11.0, intensity[ 4] =  4.5; */
  /* b[ 5] =  60.0, l[ 5] = 150.0, R[ 5] =  50.0, BEG[ 5] =  65.0, END[ 5] =  95.0, delta_theta[ 5] = 11.0, intensity[ 5] = 20.2; */
  /* b[ 6] =  58.0, l[ 6] = 160.0, R[ 6] =  55.0, BEG[ 6] =  95.0, END[ 6] = 120.0, delta_theta[ 6] = 11.0, intensity[ 6] = 25.8; */
  /* b[ 7] =  70.0, l[ 7] = 140.0, R[ 7] =  60.0, BEG[ 7] = 120.0, END[ 7] = 160.0, delta_theta[ 7] = 13.0, intensity[ 7] = 29.0; */

  /* b[ 8] = -97.0, l[ 8] = 180.0, R[ 8] =  44.0, BEG[ 8] = 350.0, END[ 8] = 380.0, delta_theta[ 8] = 11.0, intensity[ 8] =  4.1; */
  /* b[ 9] = -97.0, l[ 9] = 180.0, R[ 9] =  44.0, BEG[ 9] =  20.0, END[ 9] =  50.0, delta_theta[ 9] = 11.0, intensity[ 9] =  7.6; */
  /* b[10] = -97.0, l[10] = 180.0, R[10] =  44.0, BEG[10] =  50.0, END[10] = 100.0, delta_theta[10] = 11.0, intensity[10] = 10.2; */
  /* b[11] = -95.0, l[11] = 200.0, R[11] =  42.0, BEG[11] = 110.0, END[11] = 170.0, delta_theta[11] =  8.0, intensity[11] = 17.0; */

  /* b[12] =-100.0, l[12] = 200.0, R[12] =  65.0, BEG[12] = 205.0, END[12] = 260.0, delta_theta[12] = 11.0, intensity[12] =  3.6; */
  /* b[13] =  60.0, l[13] = 180.0, R[13] =  80.0, BEG[13] = 250.0, END[13] = 280.0, delta_theta[13] = 11.0, intensity[13] = 11.0; */
  /* b[14] =  75.0, l[14] = 180.0, R[14] = 105.0, BEG[14] =  90.0, END[14] = 120.0, delta_theta[14] = 11.0, intensity[14] =  2.7; */
  /* /\* magnetic poles: N, S, P1 *\/ */
  /* b[15] =  60.0, l[15] = 180.0, R[15] =   0.0, BEG[15] =   0.0, END[15] = 360.0, delta_theta[15] = 16.0, intensity[15] =100.0; */
  /* b[16] = -85.0, l[16] =   0.0, R[16] =   0.0, BEG[16] =   0.0, END[16] = 360.0, delta_theta[16] = 16.0, intensity[16] = 32.6; */
  /* b[17] = -17.0, l[17] =  80.0, R[17] =   0.0, BEG[17] =   0.0, END[17] = 360.0, delta_theta[17] = 14.0, intensity[17] = 0.94; */


  /* cycle 308 */
  
  b[ 0] =  70.0, l[ 0] = 160.0, R[ 0] =  50.0, BEG[ 0] = 200.0, END[ 0] = 250.0, delta_theta[ 0] = 14.0, intensity[ 0] = 54.0;
  b[ 1] =  75.0, l[ 1] = 180.0, R[ 1] =  43.0, BEG[ 1] = 250.0, END[ 1] = 300.0, delta_theta[ 1] = 13.0, intensity[ 1] = 32.0;
  b[ 2] =  75.0, l[ 2] = 180.0, R[ 2] =  43.0, BEG[ 2] = 302.0, END[ 2] = 350.0, delta_theta[ 2] = 12.0, intensity[ 2] = 54.0;
  b[ 3] =  65.0, l[ 3] = 170.0, R[ 3] =  53.0, BEG[ 3] = 355.0, END[ 3] = 400.0, delta_theta[ 3] = 11.0, intensity[ 3] = 21.0;
  b[ 4] =  60.0, l[ 4] = 160.0, R[ 4] =  55.0, BEG[ 4] =  45.0, END[ 4] =  60.0, delta_theta[ 4] = 11.0, intensity[ 4] =  6.0;
  b[ 5] =  60.0, l[ 5] = 150.0, R[ 5] =  50.0, BEG[ 5] =  65.0, END[ 5] =  95.0, delta_theta[ 5] = 11.0, intensity[ 5] = 22.0;
  b[ 6] =  58.0, l[ 6] = 160.0, R[ 6] =  55.0, BEG[ 6] =  95.0, END[ 6] = 120.0, delta_theta[ 6] = 11.0, intensity[ 6] = 22.0;
  b[ 7] =  70.0, l[ 7] = 140.0, R[ 7] =  60.0, BEG[ 7] = 120.0, END[ 7] = 160.0, delta_theta[ 7] = 13.0, intensity[ 7] = 27.0;

  b[ 8] = -97.0, l[ 8] = 180.0, R[ 8] =  44.0, BEG[ 8] = 350.0, END[ 8] = 380.0, delta_theta[ 8] = 11.0, intensity[ 8] =  4.3;
  b[ 9] = -97.0, l[ 9] = 180.0, R[ 9] =  44.0, BEG[ 9] =  20.0, END[ 9] =  50.0, delta_theta[ 9] = 11.0, intensity[ 9] =  6.5;
  b[10] = -97.0, l[10] = 180.0, R[10] =  44.0, BEG[10] =  50.0, END[10] = 100.0, delta_theta[10] = 11.0, intensity[10] =  6.5;
  b[11] = -95.0, l[11] = 200.0, R[11] =  42.0, BEG[11] = 110.0, END[11] = 170.0, delta_theta[11] =  8.0, intensity[11] = 19.0;

  b[12] =-100.0, l[12] = 200.0, R[12] =  65.0, BEG[12] = 205.0, END[12] = 260.0, delta_theta[12] = 11.0, intensity[12] =  3.2;
  b[13] =  60.0, l[13] = 180.0, R[13] =  80.0, BEG[13] = 250.0, END[13] = 280.0, delta_theta[13] = 11.0, intensity[13] =  6.5;
  b[14] =  75.0, l[14] = 180.0, R[14] = 105.0, BEG[14] =  90.0, END[14] = 120.0, delta_theta[14] = 11.0, intensity[14] =  1.0;
  
  /* magnetic poles: N, S, P1 */
  b[15] =  60.0, l[15] = 180.0, R[15] =   0.0, BEG[15] =   0.0, END[15] = 360.0, delta_theta[15] = 16.0, intensity[15] =100.0;
  b[16] = -85.0, l[16] =   0.0, R[16] =   0.0, BEG[16] =   0.0, END[16] = 360.0, delta_theta[16] = 16.0, intensity[16] = 43.0;
  b[17] = -17.0, l[17] =  80.0, R[17] =   0.0, BEG[17] =   0.0, END[17] = 360.0, delta_theta[17] = 14.0, intensity[17] =  0.0;

  N[15] = 1.0;
  N[16] = 1.0;
  N[17] = 1.0;

  
  for (int i = 0; i < 18; i++)
    {
      intensity[i] = intensity[i]/100.0; /* normalize to unity */
    }

  

  sp v;
  v.r = 1.0;

  sp u;
  u.r = 1.0;
  
  vec3 c, a, a_psi;
  vec3 pole, pole_psi;
  
  double * output = (double * ) malloc(sizeof(double) * (18 * 360 * 7 + 1));

  int k = 0;
  for (int i = 0; i < 15; i++)
    {
      /* arcs */
      N[i] = 0.0;
      
      v.phi = l[i] * (M_PI/180.0);
      v.theta = (90.0 - b[i]) * (M_PI/180.0);

      u.phi = v.phi;
      u.theta = v.theta - R[i] * (M_PI/180.0);

      c = sp2dec(v);
      a = sp2dec(u);
      
      for (int j = 0; j<360; j++)
	{
	  a = axrot(a, c, M_PI/180.0); /* rotate by 1 degree */
	  u = dec2sp(a);
	  
	  if ((u.phi > fmod(BEG[i]*(M_PI/180.0), 2.0*M_PI) && u.phi < fmod(END[i]*(M_PI/180.0), 2.0*M_PI)) || ((END[i] > 360.0 && u.phi > BEG[i]*(M_PI/180.0)) || (END[i] > 360.0 && u.phi < fmod(END[i]*(M_PI/180.0), 2.0*M_PI))))      
	    {
	      N[i] = N[i] + 1.0;
	    }

	}

      for (int j = 0; j<360; j++)
	{
	  a = axrot(a, c, M_PI/180.0); /* rotate by 1 degree */
	  u = dec2sp(a);
	  
	  if ((u.phi > fmod(BEG[i]*(M_PI/180.0), 2.0*M_PI) && u.phi < fmod(END[i]*(M_PI/180.0), 2.0*M_PI)) || ((END[i] > 360.0 && u.phi > BEG[i]*(M_PI/180.0)) || (END[i] > 360.0 && u.phi < fmod(END[i]*(M_PI/180.0), 2.0*M_PI))))  
	    {
	      /* a_psi = R_y(a, -50.0 * (M_PI/180.0)); */
	      /* a_psi = axrot(a_psi, I, psi); */

	      output[k + 1] = a.x;
	      output[k + 2] = a.y;
	      output[k + 3] = a.z;
	      output[k + 4] = delta_theta[i];
	      output[k + 5] = intensity[i];
	      output[k + 6] = N[i];
	      output[k + 7] = (double) i;

	      k = k + 7;
	   
	    }

	  

	}

      /* printf("%d\t%f\n", i, N[i]); */
    }


  for (int i = 15; i < 18; i++)
    {
      /* poles */
      v.phi = l[i] * (M_PI/180.0);
      v.theta = (90.0 - b[i]) * (M_PI/180.0);
      
      pole = sp2dec(v);

      /* pole_psi = R_y(pole, -50.0 * (M_PI/180.0)); */
      /* pole_psi = axrot(pole_psi, I, psi); */
      
      output[k + 1] = pole.x;
      output[k + 2] = pole.y;
      output[k + 3] = pole.z;
      output[k + 4] = delta_theta[i];
      output[k + 5] = intensity[i];
      output[k + 6] = N[i];
      output[k + 7] = (double) i;

      k = k + 7;

      
    }
  
  output[0] = (double) k/7;

  //printf("output 0: %f\n", output[0]);
  
  output = (double *) realloc(output, sizeof(double) * (k + 1));
  //output = (double *) realloc(output, 1000);

  return output;
  
}


double * psi_rotator(double * map, double psi)
{

  int N = (int) map[0]; 
  
  double * output = (double * ) malloc(sizeof(double) * (7*N + 1));

  vec3 a, a_psi;

  sp I_sp;
  I_sp.phi = 180.0 * (M_PI/180.0);
  I_sp.theta = 50.0 * (M_PI/180.0);
  I_sp.r = 1.0; 
  vec3 I = sp2dec(I_sp); /* NS's inertia axis */
  
  for (int k = 0; k < N; k++)
    {
      a.x = map[7*k + 1];
      a.y = map[7*k + 2];
      a.z = map[7*k + 3];

      a_psi = R_y(a, -50.0 * (M_PI/180.0));
      a_psi = axrot(a_psi, I, psi);

      output[7*k + 1] = a_psi.x;
      output[7*k + 2] = a_psi.y;
      output[7*k + 3] = a_psi.z;
      
      output[7*k + 4] = map[7*k + 4];
      output[7*k + 5] = map[7*k + 5];
      output[7*k + 6] = map[7*k + 6];
      output[7*k + 7] = map[7*k + 7];
    }
  
  output[0] = map[0];
  
  return output;
    
}





void map_show(double * map)
{

  vec3 a;
  double n;

  double N = map[0];

  for (int k=0; k < N; k++)
    {
  
      a.x = map[1 + 7*k];
      a.y = map[2 + 7*k];
      a.z = map[3 + 7*k];
      n = map[7 + 7*k];

      printf("%f\t%f\t%f\t%d\n", a.x, a.y, a.z, (int) n);
      
    }

  vec3 o;
  sp o_sp;
  
  o_sp.r = 1.0;
  o_sp.phi = 0.0;
  o_sp.theta = (90.0 - 3.0)  * (M_PI/180.0);
  int lc_num = 360;
  
  for (int i=0; i < lc_num; i++)
    {
      o_sp.phi = 2.0 * M_PI * (1.0 - (double) i/lc_num);
      o = sp2dec(o_sp);
      //o = R_y(o, 50.0 * (M_PI/180.0));


      printf("%f\t%f\t%f\t%f\n", o.x, o.y, o.z, 18.0);

    }

}



double * pulse(double * map, double index, int lc_num)
{
  
  /* double * raw_map = X_map_2(); */
  /* double * map = psi_rotator(raw_map, psi); */
  
  int N_sum = map[0];

  double * output = (double * ) malloc(sizeof(double) * lc_num);
  
  
  sp o_sp;
  o_sp.r = 1.0;
  o_sp.phi = 0.0;
  o_sp.theta = (90.0 - 3.0) * (M_PI/180.0); /* 3 degree is NS orientation */

  vec3 o;

  vec3 a;

  double phi;

  double delta_theta;
  double intensity;
  double n;
  double theta_ao;
  
  for (int i=0; i < lc_num; i++)
    {

      output[i] = 0.0;
      
      phi = 2.0 * M_PI * ((double) i)/((double) lc_num);
      
      o_sp.phi = 2.0 * M_PI - phi;
      
      o = sp2dec(o_sp);
      
      for (int k=0; k < N_sum; k++)
	{
	      
	  if (index >= 0.0 && map[7 + 7*k] != index)
	    {
	      continue;
	    }

	      
	  a.x = map[1 + 7*k];
	  a.y = map[2 + 7*k];
	  a.z = map[3 + 7*k];
	  	      
	  theta_ao = acos(dot(a,o));

	  if (theta_ao <= 0.0)
	    {
	      continue;
	    }

	  
	  delta_theta = map[4 + 7*k] * (M_PI/180.0);
	  intensity = map[5 + 7*k];
	  n = map[6 + 7*k];
	      
	  output[i] = output[i] + (intensity/n) * exp(-pow(theta_ao/delta_theta, 2)/2.0) * cos(theta_ao);
	  //printf("%f\t%f\t %f\t %f\t %f\n", index, intensity, n, theta_ao, delta_theta);
	      
	}

      //output[i] = output[i] + B_psi(psi, parameters) * sin(phi + phi_0);
      //printf("TEST %d\t%f\n",  i, output[i]);

    }
    
  return output;
  
}





void pulsar_X_ray_flux(double index, parameters parameters)
{

  double psi = 0.0 * (M_PI/180.0);

  int N_sum;

  int lc_num = 100;
  int super_orb_num = 400;
  
  sp o_sp;

  o_sp.r = 1.0;
  o_sp.phi = 0.0;
  o_sp.theta = (M_PI/2.0) - parameters.NS_kappa;

  vec3 o;

  vec3 a;
  double delta_theta;
  double intensity;
  double n;
  double theta_ao;

  double arc_id;
  
  double phase[lc_num];
  
  
  /* double tau_eff = 0.0; */
  /* double B_psi = 0.0; */
  /* double phi_0 = 0.6 * 2.0 * M_PI; */

  double phi;

  double integrated_flux[super_orb_num];

  double * raw_map = X_map_2();
  double * map;
  double * summa;
  
  for (int j = 0; j < super_orb_num; j++)
    {
      
      psi = (2.0*M_PI) * ((double) j)/((double) super_orb_num);
      
      map = psi_rotator(raw_map, psi);
      summa = pulse(map, index, lc_num);

      integrated_flux[j] = 0.0;

      for (int i = 1; i < lc_num; i++)
	{

	  integrated_flux[j] = integrated_flux[j] + (double) (1.0/lc_num) * (summa[i-1] + summa[i])/2.0;
	  
	}

      //integrated_flux[j] = (integrated_flux[j] + B_psi(psi)) * (1.0 - exp(-tau(psi)));
      //integrated_flux[j] = integrated_flux[j] * (1.0 - exp(-tau(psi, parameters)));

      /* integrated_flux[j] = integrated_flux[j] + B_psi(psi); */

    }

  free(raw_map);
  free(map);
  free(summa);


  double s_koeff = 0.01;
  double delta_s = 0.1;

    
  for (int j = 0; j < super_orb_num; j++)
    {
      printf("%f\t%f\t%d\n", (double) j/super_orb_num, integrated_flux[j], (int) index);
    }
  printf("%f\t%f\t%d\n", 1.0, integrated_flux[0], (int) index);


  

}





int main(int argc, char **argv)
{

  parameters parameters;

  parameters.inclination = 88.0 * (M_PI/180.0);

  parameters.NS_kappa = 10.0 * (M_PI/180.0);
  parameters.NS_theta = 3.0 * (M_PI/180.0);

  double theta_i = parameters.NS_theta + (M_PI/2.0 - parameters.inclination); /* "+" if north magnetic pole towads our semispace */
  
  double kappa = acos(cos(parameters.NS_kappa) * cos(theta_i)); /* angle between normal vector orbital plane and rotation axis */

  double NS_J_psi = 3.0*M_PI/2.0 + acos((cos(theta_i) - cos(parameters.NS_kappa)*cos(kappa))/(sin(parameters.NS_kappa)*sin(kappa)));
  //double NS_J_psi = M_PI/2.0 - acos((cos(theta_i) - cos(parameters.NS_kappa)*cos(kappa))/(sin(parameters.NS_kappa)*sin(kappa)));

  
  double h = 8.0 * (M_PI/180.0); /* set h_out using degrees */
  
  parameters.r_out = 0.24; /* outer radius of the disk */
  parameters.h_out = tan(h) * parameters.r_out; /* semi-thickness of the outer edge of the disk */ 


  
  parameters.theta_out = 14.0 * (M_PI/180.0); 
  parameters.theta_in = 10.0 * (M_PI/180.0); 

    
  // tau
  parameters.epsilon_0_in =  8.0 * (M_PI/180.0);
  parameters.epsilon_0_out = 5.0 * (M_PI/180.0);


  
  int v35;
  double psi;

  sscanf(argv[1], "%d", &v35);
  sscanf(argv[2], "%lf", &psi);

  int lc_num = 400;

  double * raw_map = X_map_2();
  double * map = psi_rotator(raw_map, psi);


  if (v35)
    {
      for (int i = -1; i < 0; i++)
  	{
  	  pulsar_X_ray_flux(i, parameters);
  	}
    }
  else
    {

      for (int i = -1; i < 0; i++)
  	{
	  double * p;
  	  p = pulse(map, i, lc_num);
  	  for (int j = 0; j < lc_num; j++)
  	    {
  	      printf("%f\t%f\t%d\n", ((double) j)/((double) lc_num), p[j], i);

  	    }
	  free(p);
  	}

      //free(p);
    }



  
  /* double * raw_map = X_map_2(); */
  /* double * map = psi_rotator(raw_map, psi); */

  /* int i = -1; */
  /* double * p = pulse(map, i, lc_num); */

  /* for(int i = 0; i < lc_num; i++) */
  /*   { */
  /*     printf("%f\t%f\n", ((double) i)/((double) lc_num), p[i]); */
  /*   } */
  
  printf("\n");
  
  map_show(map);

  free(raw_map);
  free(map);


  
  return 0;
}
