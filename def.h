#ifndef DEF_H
#define DEF_H

#include <math.h>

#define eps 0.000001

/* Planck constant, erg * s */
#define H_PLANCK 6.626070040e-27

/* Speed on light, cm/s */
#define C 29979245800.0

/* Boltzmann constant */
#define K 1.38064852e-16

/* Stefan–Boltzmann constant */
#define SIGMA 5.6704e-5

/* constants in Planck function, see F_lambda in intensity.c */
#define C1 3.74185e-5
#define C2 1.43883


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


struct disk {
  vec3 h;
  double R;
};

typedef struct disk disk;


double len(vec3 p);
double dot(vec3 a, vec3 b);

vec3 sp2dec(sp a);
sp dec2sp(vec3 a);

vec3 rotate(vec3 a, double Y, double Z);
vec3 axrot(vec3 a, vec3 u, double theta);

double F_lambda(double T, double lambda);

sp arcgen(double theta, double N_phi, int i);
double * x_ray_direction_diagram(double PSI_pr, double Lx);

double fr(double r, double phi, double theta, double q, double omega);
double dfr(double r, double phi, double theta, double q, double omega);
double radius_star(double phi, double theta, double q, double omega);

double * polar(double q, double omega);
double * gradient(double phi, double theta, double q, double omega);

double fx(double x, double q);
double dfx(double x, double q);
double fomega(double r, double q, double omega_crit);
double dfomega(double r, double q);
double omg(double q, double mu);
double distance_to_star(vec3 p, double omega, double q);

double eclipse_by_star(double omega, double q, vec3 o, vec3 p);
double radius_disk(disk disk, double phi, double theta);

double distance_to_disk(vec3 p, disk disk);
double distance_to_disk_inside(vec3 p, disk disk);

double eclipse_by_disk(disk disk, vec3 o, vec3 p);
double eclipse_by_disk_inside(disk disk, vec3 o, vec3 p);

double * flux_star(vec3 o, double q, double omega, double beta, double u, disk disk, vec3 d2, double Lx, double Lx_disk, double Lx_iso, double Lx_disk_2, double albedo, int star_tiles, double T_star, double lambda, double a, vec3 neutron_star, double PSI_pr, int picture, int isotrope, sp disk_reflection_diagr);

double B(disk disk, double A, double rho_in, double T);

double flux_disk(vec3 o, disk disk, double rho_in, double A, double uniform_disk, double y_tilt, double z_tilt, double omega, double q, int disk_tiles, double phi_orb, double T_disk, double lambda, double a, int picture, int spot_disk, double T_spot, double spot_beg, double spot_end, double spot_rho_in, double spot_rho_out);



#endif //DEF_H
