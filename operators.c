#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "def.h"


double len(vec3 p)
{
  return sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
}

double dot(vec3 a, vec3 b)
{
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

vec3 sp2dec(sp a)
{
  vec3 result;

  result.x = a.r * sin(a.theta) * cos(a.phi);
  result.y = a.r * sin(a.theta) * sin(a.phi);
  result.z = a.r * cos(a.theta);

  return result;  
}

sp dec2sp(vec3 a)
{
  /* transforms cartesian coordinates to spherical */
  
  sp result;
  
  double x = a.x;
  double y = a.y;
  double z = a.z;

  double phi, theta;

  /* r [0,inf) */
  result.r = sqrt(x * x + y * y + z * z);

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

  result.phi = phi;
  result.theta = theta;
  
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


vec3 axrot(vec3 a, vec3 u, double theta)
{
  /* rotates by an angle of theta about an axis in the direction of u */

  vec3 result;  
  result.x = a.x * (cos(theta) + u.x * u.x * (1.0 - cos(theta))) + a.y * (u.x * u.y * (1.0 - cos(theta)) - u.z * sin(theta)) + a.z * (u.x * u.z * (1.0 - cos(theta)) + u.y * sin(theta));
  result.y = a.x * (u.y * u.x * (1.0 - cos(theta)) + u.z * sin(theta)) + a.y * (cos(theta) + u.y * u.y * (1.0 - cos(theta))) + a.z * (u.y * u.z * (1.0 - cos(theta)) - u.x * sin(theta));
  result.z = a.x * (u.z * u.x * (1.0 - cos(theta)) - u.y * sin(theta)) + a.y * (u.z * u.y * (1.0 - cos(theta)) + u.x * sin(theta)) + a.z * (cos(theta) + u.z * u.z * (1.0 - cos(theta)));

  return result;
}

