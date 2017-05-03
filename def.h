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


struct vec3p {
  double x;
  double y;
  double z;
};

typedef struct vec3p vec3;

struct disk {
  vec3 h;
  double R;
};

typedef struct disk disk;

double len(vec3 p);
double dot(vec3 a, vec3 b);
double * coordinate_transformation(double x, double y, double z);

double F_lambda(double T, double lambda);

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
double eclipse_by_disk(disk disk, vec3 o, vec3 p);

double flux_star(vec3 o, double q, double omega, double beta, double u, disk disk, double Lx, double albedo, int star_tiles, double T_star, double lambda, double a);
double flux_disk(vec3 o, disk disk, double y_tilt, double z_tilt, double omega, double q, double b, int disk_tiles, double phi_orb, double T_disk, double lambda, double a);



