#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "def.h"



bool ray_star(double omega, double q, vec3 o, vec3 p)
{

  double lambda = distance_to_star(p, omega, q);

  do {
    p = sum(p, scale(o, lambda));
      
    lambda = distance_to_star(p, omega, q);

  } while (lambda < 1.0 && lambda > eps);

  if (lambda >= 1.0)
    {
      return true;
    }
  else if (lambda <= eps)
    {
      return false;
    }

}



bool ray_disk(disk disk, vec3 o, vec3 p)
{

  
  double lambda;
  
  
  if ( (fabs(dot(disk.n, p)) < disk.h) &&
       (sqrt(pow(len(p),2) - pow(fabs(dot(disk.n, p)),2)) < disk.R - eps)
       )
    {

      double a = pow(dot(o, disk.n),2) - dot(o,o);
      double b = 2.0 * (dot(p,disk.n)*dot(o,disk.n) - dot(p,o));
      double c = pow(disk.r_out, 2) + pow(dot(p,disk.n),2) - dot(p,p);

      lambda = (- b - sqrt(pow(b,2) - 4.0 * a * c))/(2.0*a);

      vec3 q = sum(p, scale(o, lambda)); 
        
      if ( fabs(dot(disk.n, q)) > disk.h && dot(disk.n, p)*dot(disk.n, q) > 0 )
	{
	  return true;
	}
      else
	{
	  return false;
	}
    }
  else
    {

      lambda = distance_to_disk(p, disk);
      
      if (lambda < eps)
	{
	  /* point belongs to the surface of the disk */
	  return true;
	}
      else
	{      
	  do {
	    p = sum(p, scale(o, lambda));
      
	    lambda = distance_to_disk(p, disk);

	  } while (lambda < 5.0 && lambda > eps);

	  if (lambda >= 5.0)
	    {
	      return true;
	    }
	  else if (lambda <= eps)
	    {
	      return false;
	    }
	}
    }
}



bool disk_shadow(vec3 p, disk disk)
{

  double theta_out = disk.theta_out ;
  double theta_in = disk.theta_in;
  double phi_out = disk.phi_out;
  double phi_in = disk.phi_in;
  double r_out = disk.r_out;
  double r_in = disk.r_in;

  sp v = dec2sp(p);

  
  double t = v.theta;


  sp w = dec2sp(disk.n);
  
  phi_out = M_PI/2.0 + phi_out;
  phi_in = M_PI/2.0 + phi_in;

  double A = (theta_out - theta_in)/(r_out - r_in);  
  double B = (phi_out - phi_in)/(r_out - r_in);
 

  double in = (A * (r_in - r_out) + theta_out) * sin(v.phi - (B*(r_in - r_out) + phi_out));
  double out = theta_out * sin(v.phi - phi_out);


  double r_1, r_0;

  if (len(p) > sqrt(pow(disk.h, 2) + pow(disk.R, 2)))
    {

      //printf("TEST");
      
      if ( (fabs(M_PI/2.0 - t - out) > atan(disk.h/disk.R)) && (max(theta_in, theta_out) - fabs(M_PI/2.0 - t)) < 0 )
	{
	  return true;
	}
      else if ( (fabs(M_PI/2.0 - t - out) < atan(disk.h/disk.R)) || (M_PI/2.0 - t - in)*(M_PI/2.0 - t - out) < 0 )
	{
	  return false;
	}
      else if ( (fabs(M_PI/2.0 - t - out) > atan(disk.h/disk.R)) && (M_PI/2.0 - t - in)*(M_PI/2.0 - t - out) > 0 && (max(theta_in, theta_out) - fabs(M_PI/2.0 - t)) > 0 )
	{

	  if (B == 0)
	    {
	      return true;
	    }
	  else if (B != 0)
	    {
      
	      r_0 = (r_in + r_out)*0.5 ;

	      do {
	
		r_1 = r_0 -
		  (( - B * cos(B * (r_0 - r_out) + phi_out - v.phi) * (theta_out + A * (r_0 - r_out))) - A * sin(B * (r_0 - r_out) + phi_out - v.phi))/( B * B * sin(B * (r_0 - r_out) + phi_out - v.phi) * (theta_out + A * (r_0 - r_out)) - 2 * A * B * cos(B * (r_0 - r_out) + phi_out - v.phi));
	
		r_0 = r_1;

	    
	      } while (fabs(r_0 - r_1) > eps);

	      if ((M_PI/2.0 - t - in)*(M_PI/2.0 - t - (A * ((r_1 - r_out) + theta_out) * sin(v.phi - (B*(r_1 - r_out) + phi_out)))) > 0)
		{
		  return true;
		}
	      else
		{
		  return false;
		}
	    }
	}

    }
  else
    {
  
  
      if ( (M_PI/2.0 - t - in)*(M_PI/2.0 - t - out) < 0 )
	{
	  return false;
	}
      else
	{
	  if (B == 0)
	    {
	      return true;
	    }
	  else if (B != 0)
	    {
      
	      r_0 = (r_in + r_out)*0.5 ;

	      do {
	
		r_1 = r_0 -
		  (( - B * cos(B * (r_0 - r_out) + phi_out - v.phi) * (theta_out + A * (r_0 - r_out))) - A * sin(B * (r_0 - r_out) + phi_out - v.phi))/( B * B * sin(B * (r_0 - r_out) + phi_out - v.phi) * (theta_out + A * (r_0 - r_out)) - 2 * A * B * cos(B * (r_0 - r_out) + phi_out - v.phi));
	
		r_0 = r_1;

	    
	      } while (fabs(r_0 - r_1) > eps);

	      if ((M_PI/2.0 - t - in)*(M_PI/2.0 - t - (A * ((r_1 - r_out) + theta_out) * sin(v.phi - (B*(r_1 - r_out) + phi_out)))) > 0)
		{
		  return true;
		}
	      else
		{
		  return false;
		}
	    }	
	}

    }



  
}









vec3 rand_p(int seed, double q, double omega, double delta_r)
{

  //srand(seed);
  
  vec3 p;
  p.x = 0.0;
  p.y = 0.0;
  p.z = 0.0;


  vec3 r;
  
  sp v;

  vec3 x_ray;
  vec3 n_x_ray;

  
  vec3 ns_shift;
  ns_shift.x = -1.0;
  ns_shift.y = 0.0;
  ns_shift.z = 0.0;

  //double * normal;
  vec3 n;

  
  do {

      double * normal;

      v.phi = ((double) rand()/RAND_MAX) * M_PI + 1.5 * M_PI;
      v.theta = ((double) rand()/RAND_MAX) * M_PI;
      v.r = radius_star(v.phi, v.theta, q, omega);

      r = sp2dec(v);

      //printf("r >>>>>> %f\t%f\t%f\t%f\n", r.x, r.y, r.z, len(r));
	     
      normal = gradient(v.phi, v.theta, q, omega);
      n.x = normal[1];
      n.y = normal[2];
      n.z = normal[3];

      //printf("%f\n", n.x);

      
      //delta_r = ((double) rand()/RAND_MAX) * delta_r * len(r);

      delta_r = 2.0 * eps;
      //printf("%f\n", v.r);
      
      
      p = sum(r, scale(n, delta_r)); 
      
      //printf("n %f\t%f\t%f\t%f\n", n.x, n.y, n.z, len(n));

      x_ray = sum(r, ns_shift);
      //printf("x %f\t%f\t%f\t%f\n", x_ray.x, x_ray.y, x_ray.z, len(x_ray));

      n_x_ray = scale(x_ray, 1.0/len(x_ray));
      //printf("nx %f\t%f\t%f\t%f\n", n_x_ray.x, n_x_ray.y, n_x_ray.z, len(n_x_ray));


      //printf("dot n nx >>>>>> %f\n", dot(n, n_x_ray));
      free(normal);
   } while (dot(n_x_ray, n) > eps);


  //free(normal); 
  return p;

}
  













