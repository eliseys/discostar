#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "def.h"

double integrate(double * f, int N)
{

  double delta_x = 1.0/(N - 1); 
  
  double integral = 0.0;
  
  for (int i = 0; i < N-1; i++)
    {
      integral = integral + ((f[i+1] + f[i])/2.0)*delta_x;
    }

  return integral;

}


double * Legendre_2(int N)
{
  // 2nd Legendre polynomial

  double * output = (double *) malloc(sizeof(double) * N);

  double x;
  
  for (int i = 0; i<N; i++)
    {
      x = (double) i/(N-1);
      output[i] = 0.5*(3.0*pow(x,2) - 1.0);
    }
  
  return output;

}


double * h_0(int N)
{
  // initial guess for h
  
  double * output = (double *) malloc(sizeof(double) * N);
  
  for (int i=0; i<N; i++)
    {
      output[i] = 1.0 + (double) i/(N-1);
    }

  return output;
}



double * phi_scatter(double x_0, double x_2, double a, int N)
{
  // scattering diagram

  double * output = (double *) malloc(sizeof(double) * N);
  double * P_2 = Legendre_2(N);

  double mu;
  
  for (int i=0; i<N; i++)
    {
      mu = (double) i/(N-1);
      output[i] = (a/2.0)*(x_0 + (x_2/2.0)*(3.0*(1.0 - a)*pow(mu,2) - 1.0)*P_2[i]);
      //printf("phi %f\n", output[i]);
    }

  return output;

}



double I_h(double * h, double * phi, double mu, int N)
{
  // intermediate step for h calculation
  
  double integral = 0;

  double x, xp1;

  int i;
  
  if (mu == 0)
    {
      for (i = 0; i < N-1; i++)
	{

	  x = (double) i/(N - 1); 
	  xp1 = (double) (i + 1)/(N - 1);
	  
	  integral = integral + 0.5*(phi[i] * h[i] + phi[i+1] * h[i+1])*(1.0/(N-1));
	}
    }
  else
    {
      for (i = 0; i < N-1; i++)
	{
	  x = (double) i/(N - 1); 
	  xp1 = (double) (i + 1)/(N - 1);
	  
	  integral = integral + 0.5*(x * phi[i] * h[i]/(mu + x) + xp1 * phi[i+1] * h[i+1]/(mu + xp1))*(1.0/(N-1));
	}
    }
  
  //printf("%f ", integral);
  return integral;

}



double * H(double x_0, double x_2, double a, int N)
{
   /* N is the number of division points mu_i */
  int t = 0;
  int i;

  double * h_next = (double *) malloc(sizeof(double) * N);

  double * h_prev = h_0(N);

  double h_hat[N];

  double * phi = phi_scatter(x_0, x_2, a, N);

  
  double phi_0 = integrate(phi, N); // integral of scatteting diagram

  if (fabs(phi_0 - 0.5) < 1e-6)
    {
      phi_0 = 0.5;
    }


  
  /* printf("phi_0 %1.30f\n", phi_0); */

  
  double lambda = 0.5 * (1 + sqrt(1 - 2 * phi_0));
  

  
  /* printf("lambda %f\n", lambda); */
  /* printf("sqrt(0) %f\n", sqrt(0)); */

  
  
  double mu[N];

  for (int j=0; j < N; j++)
    {
      mu[j] = ((double) j)/((double) N - 1);
      //printf("mu[j] %f\n", mu[j]);
      
    }



  
  do {
    
    for (i=0; i<N; i++)
      {

    	double I = I_h(h_prev, phi, mu[i], N);
	
	h_hat[i] = 1.0/(sqrt(1 - 2 * phi_0) + I);
	
      }
    
    for (i=0; i<N; i++)
      {

	h_next[i] = lambda * h_hat[i] + (1.0 - lambda) * h_prev[i];

      }
    
    for (i=0; i<N; i++)
      {
	h_prev[i] = h_next[i];
      }
    
    t++;
    
  } while (t < 20);


  /* for (i = 0; i < N; i++) */
  /*   { */
  /*     printf("H: h[i] %f\n", h_next[i]); */
  /*   } */
  
  return h_next;
  
}



double * rho(double x_0, double x_2, double a, int N)
{

  int i, j;
  
  double * rho_mu_mu = (double *) malloc(sizeof(double) * N * N);

  double * h = H(x_0, x_2, a, N);
  double * P_2 = Legendre_2(N);

  double * F = (double *) malloc(sizeof(double) * N);
  double * mu = (double *) malloc(sizeof(double) * N);

  double * q_0 = (double *) malloc(sizeof(double) * N);
  double * q_2 = (double *) malloc(sizeof(double) * N);

  double * phi_0 = (double *) malloc(sizeof(double) * N);
  double * phi_2 = (double *) malloc(sizeof(double) * N);

  
  for (i = 0; i<N; i++)
    {
      mu[i] = (double) i/(N-1);
    }
 
  for (i = 0; i<N; i++)
    {
      F[i] = h[i] * (3.0*x_2*P_2[i]/2.0) * mu[i];
      /* printf("h[i] %f\n", h[i]); */
    }


  
  /* printf("integrate M_0 %f\n", integrate(F, N)); */
  
  double M_0 = - (a/2.0) * (1.0 - a) * integrate(F, N);

  
  for (i = 0; i<N; i++)
    {
      F[i] = h[i] * (1.0 - x_2*P_2[i]/2.0);
    }
 
  double M_1 = 1.0 - (a/2.0)*integrate(F, N);


  for (i = 0; i<N; i++)
    {
      F[i] = h[i] * (1.0 - x_2*P_2[i]/2.0) * mu[i];
    }
  
  double M_2 = - (a/2.0)*integrate(F, N);


  for (i = 0; i<N; i++)
    {
      F[i] = h[i] * P_2[i];
    }

  double N_0 = -(a/4.0) * x_2 * 3.0 * (1.0 - a) * integrate(F, N);


  for (i = 0; i<N; i++)
    {
      F[i] = h[i] * P_2[i] * mu[i];
    }

  double N_1 = (a/4.0) * x_2 * 3.0 * (1.0 - a) * integrate(F, N);


  double N_2 = M_1;

  
  double Delta = M_1 * N_2 - M_2 * N_1;

  /* printf("M_0 %f\n", M_0); */
  /* printf("M_1 %f\n", M_1); */
  /* printf("M_2 %f\n", M_2); */
  /* printf("N_0 %f\n", N_0); */
  /* printf("N_1 %f\n", N_1); */
  /* printf("N_2 %f\n", N_2); */
  
  /* printf("Delta %f\n", Delta); */
  
  for (i = 0; i<N; i++)
    {
      q_0[i] = 1.0 + mu[i]*(N_2*M_0-M_2*N_0)/Delta + pow(mu[i],2)*(M_1*N_0-N_1*M_0)/Delta;
    }

  for (i = 0; i<N; i++)
    {
      q_2[i] = - (1.0/2.0)*q_0[i] - (3.0*(1.0-a)/(2.0*Delta))*(M_2*mu[i] - M_1*pow(mu[i],2));
    }

  for (i = 0; i<N; i++)
    {
      phi_0[i] = h[i]*q_0[i];
      phi_2[i] = h[i]*q_2[i];

      /* printf("phi_0 %f\t phi_2 %f\n", phi_0[i], phi_0[i]); */
    }

  for (j = 0; j<N; j++)
    {
      for (i = 0; i<N; i++)
	{
	  rho_mu_mu[j*N + i] = (a/4.0)*(x_0 * (phi_0[j]*phi_0[i]/(mu[j] + mu[i])) + x_2*(phi_2[j]*phi_2[i]/(mu[j] + mu[i])));
	}
    }

  return rho_mu_mu;
  
}


double * I(double flux, double x_0, double x_2, double a, int N)
{

  double * output = (double *) malloc(sizeof(double) * N * N);
 
  int i, j;

  double * rho_mu_mu = rho(x_0, x_2, a, N); 

  double * mu = (double *) malloc(sizeof(double) * N);
  for (i = 0; i<N; i++)
    {
      mu[i] = (double) i/(N-1);
    }


  for (j = 0; j<N; j++)
    {
      for (i = 0; i<N; i++)
	{
	  output[j*N + i] = flux * rho_mu_mu[j*N + i] * mu[i];
	}
    }

 
  return output;

}





double * albedo_mu(double x_0, double x_2, double a, int N)
{

  int i, j;
  
  double * albedo = (double *) malloc(sizeof(double) * N * N);

  double * h = H(x_0, x_2, a, N);
  double * P_2 = Legendre_2(N);

  double * F = (double *) malloc(sizeof(double) * N);
  double * mu = (double *) malloc(sizeof(double) * N);

  double * q_0 = (double *) malloc(sizeof(double) * N);
  double * q_2 = (double *) malloc(sizeof(double) * N);

  double * phi_0 = (double *) malloc(sizeof(double) * N);
  double * phi_2 = (double *) malloc(sizeof(double) * N);

  
  for (i = 0; i<N; i++)
    {
      mu[i] = (double) i/(N-1);
    }
 
  for (i = 0; i<N; i++)
    {
      F[i] = h[i] * (3.0*x_2*P_2[i]/2.0) * mu[i];
      /* printf("h[i] %f\n", h[i]); */
    }


  
  /* printf("integrate M_0 %f\n", integrate(F, N)); */
  
  double M_0 = - (a/2.0) * (1.0 - a) * integrate(F, N);

  
  for (i = 0; i<N; i++)
    {
      F[i] = h[i] * (1.0 - x_2*P_2[i]/2.0);
    }
 
  double M_1 = 1.0 - (a/2.0)*integrate(F, N);


  for (i = 0; i<N; i++)
    {
      F[i] = h[i] * (1.0 - x_2*P_2[i]/2.0) * mu[i];
    }
  
  double M_2 = - (a/2.0)*integrate(F, N);


  for (i = 0; i<N; i++)
    {
      F[i] = h[i] * P_2[i];
    }

  double N_0 = -(a/4.0) * x_2 * 3.0 * (1.0 - a) * integrate(F, N);


  for (i = 0; i<N; i++)
    {
      F[i] = h[i] * P_2[i] * mu[i];
    }

  double N_1 = (a/4.0) * x_2 * 3.0 * (1.0 - a) * integrate(F, N);


  double N_2 = M_1;

  
  double Delta = M_1 * N_2 - M_2 * N_1;

  
  for (i = 0; i<N; i++)
    {
      albedo[i] = 1.0 - ((1.0 - a)/Delta)*h[i]*(N_2 - N_1*mu[i]);
    }

  return albedo;

}
