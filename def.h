#ifndef DEF_H
#define DEF_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <stdbool.h>

#define eps 0.000001

/* Planck constant, erg * s */
#define H_PLANCK 6.626070040e-27

/* Speed on light, cm/s */
#define C 29979245800.0

/* Boltzmann constant */
#define K 1.38064852e-16

/* Stefan–Boltzmann constant */
#define SIGMA 5.6704e-5

/* Thomson scattering cross-section */
#define SIGMA_THOMSON 6.65245873e-25

/* constants in Planck function, see F_lambda in intensity.c */
#define C1 1.19104295e-5
#define C2 1.43883

#define min(x, y) (((x) < (y)) ? (x) : (y))
#define max(x, y) (((x) > (y)) ? (x) : (y))


#define A_cm 1e-8

/* star parameters */
//double q;                /* mass ratio q = m_x/m_opt */
//double mu;               /* roche lobe filling */
//double beta;             /* gravity darkening coefficient */
//double u;                /* limb darkening coefficient */  
//double albedo;           /* 1 - (X-ray photons reprocessing efficiency) */	 
//double T_star;

/* neutron star parameter */
//double Lx_noniso;               /* erg/s */ 
//
//double h;                /* semi-thickness */ 
//double R;                /* radius */
//double y_tilt; 
//double z_tilt;  
//
//double y_tilt2; 
//double z_tilt2;  
//
//int picture;             /* print 3D picture */
//  
///* observer */
//double inclination;
//
//int lc_num; /* number of light curve points */
//
///* star and disk surface partition */
//int star_tiles;
//int disk_tiles;
//
//int threads;             /* number of OpenMP threads */
//
///* temperature of star and disk and spectral band */
//double T_disk;
//
//double a_cm;
//
//double PSI_pr;
//double kappa;
//
//double Lx_disk;
//double Lx_disk_2;
//double Lx_iso;
//
//int spot_disk;
//double T_spot;
//double spot_beg;
//double spot_end;
//
//double ns_theta;
//
//double spot_rho_in;
//double spot_rho_out;
//
//double drd_phi;
//double drd_theta;
//
//double rho_in;
//
//double A;
//
//double uniform_disk;
//
//double disk_flux_B;
//double disk_flux_V;
//
//double h_warp;
//
//char filter;
//
//double lc_start;
//double lc_end;


/**/



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



struct parameters {

  double q; /* mass ratio q = m_x/m_opt */
  double mu; /* roche lobe filling */
  double omega; /* dimentionless potential */ 
  double beta;    
  double u;   
  double X_albedo;  /* X-ray albedo */	 
  double T_star_polar;

  double a; /* distance between center of masses */
  double inclination;

  double Lx; /* NS X-ray luminosity */

  double NS_phi;
  double NS_kappa;
  double NS_theta;

  double h_out; /* semi-thickness of the outer edge of the disk */ 
  double r_out; /* outer radius of the disk */
  double r_in; /* inner radius of the disk */
  double gamma; /* h-profile exponent of the disk */

  double theta_out; 
  double phi_out;  

  double theta_in; 
  double phi_in;  

  double epsilon_0_out; /* tau scale */
  double epsilon_0_in;

  bool do_lc; /* If true calculate light curve. If false returns xyz-coordinates of points and corresponding temperatures */
  
  int N_lc; /* number of light curve points */

  int N_theta;
  int N_r;

  int OMP_threads; /* OpenMP threads */

  char filter[5];

  bool do_corona; /* calculate X-ray flux from corona */
  int N_corona; /* number of randomly distributed points in the corona */
  double h_corona; /* height of the corona */
  int rs_corona; /* random seed */

};

typedef struct parameters parameters;




struct disk {
  vec3 n;
  double h;
  double R;
  double r_out;
  double r_in;
  double theta_out;
  double theta_in;
  double phi_out;
  double phi_in;
  double gamma;
};

typedef struct disk disk;


struct ion {
  char * name[8];
  int * Z;
  int * N;
  float * E_th;
  float * E_max;
  float * E_0;
  float * sigma_0;
  float * y_a;
  float * P;
  float * y_w;
  float * y_0;
  float * y_1;
};
typedef struct ion ion;





bool ray_star(double omega, double q, vec3 o, vec3 p);
bool ray_disk(disk disk, vec3 o, vec3 p);
bool disk_shadow(vec3 p, disk disk);
double * rand_p(parameters parameters);

double len(vec3 p);
double dot(vec3 a, vec3 b);

vec3 sum(vec3 a, vec3 b);
vec3 scale(vec3 a, double s);


vec3 sp2dec(sp a);
sp dec2sp(vec3 a);

vec3 rotate(vec3 a, double Y, double Z);
vec3 R_x(vec3 a, double X);
vec3 R_y(vec3 a, double Y);
vec3 R_z(vec3 a, double Z);

vec3 axrot(vec3 a, vec3 u, double theta);

double F_filter(double T, char * filter);
double kappa_mm(double energy);
double x_ray_reflected_fraction(double zeta, double energy);
double normalized_x_ray_spectrum(double energy);
double x_ray_flux_integrated(double zeta);

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

double * phi_func(int steps_phi, int threads);
double * theta_func(int steps_theta, int threads);
double * shape_r(int steps_phi, int steps_theta, double * phi_array, double * theta_array, double q, double omega, int threads);
double * shape_g_abs(int steps_phi, int steps_theta, double * phi_array, double * theta_array, double q, double omega, int threads);

double * star_geometry(parameters parameters);


double * phi_func_disk(int steps_phi_disk);
double * theta_func_disk(int steps_theta_disk);

double eclipse_by_star(double omega, double q, vec3 o, vec3 p);
double radius_disk(disk disk, double phi, double theta);

double distance_to_disk(vec3 p, disk disk);
double distance_to_disk_inside(vec3 p, disk disk);

double eclipse_by_disk(disk disk, vec3 o, vec3 p);
double disk_h_profile(disk disk, double r, double gamma);
double disk_h_diff_profile(disk disk, double r, double gamma);
//double eclipse_by_disk_inside(disk disk, vec3 o, vec3 p);

double * disk_geometry(disk disk, parameters parameters);

double x_ray_corona(parameters parameters, disk disk, vec3 observer, double * corona_elements);

double star_F(double * star_elements, parameters parameters, disk disk, vec3 observer, vec3 neutron_star, double * Ix_dd);

double star_X(double * star_elements, parameters parameters, disk disk, vec3 observer, vec3 neutron_star, double * Ix_dd, double E);

double flux_star(vec3 o, double q, double omega, double beta, double u, disk disk, double Lx_noniso, double Lx_disk, double Lx_iso, double Lx_disk_2, double albedo, int star_tiles, double T_star, double a, vec3 neutron_star, double PSI_pr, int picture, sp disk_reflection_diagr, double * r_array, double * g_array, double * phi_array, double * theta_array, double * Ix_dd, double y_tilt, double y_tilt2, double z_tilt, double z_tilt2, double phi_orb, char * filter);

double B(disk disk, double A, double rho_in, double T);

double disk_F(double * disk_elements, parameters parameters, disk disk, vec3 observer);

double * flux_disk(vec3 o, disk disk, double rho_in, double A, double uniform_disk, double y_tilt, double z_tilt, double omega, double q, int disk_tiles, double phi_orb, double T_disk, double a, int picture, int spot_disk, double T_spot, double spot_beg, double spot_end, double spot_rho_in, double spot_rho_out, double h_warp, double * Ix_dd, double Lx_iso, double Lx_noniso, vec3 neutron_star, char * filter);

void ion_init();
double sigma_phot(double E, int i);


double integrate(double * f, int N);
double * Legendre_2(int N);
double * h_0(int N);
double * phi_scatter(double x_0, double x_2, double a, int N);
double I_h(double * h, double * phi, double mu, int N);
double * H(double x_0, double x_2, double a, int N);
double * rho(double x_0, double x_2, double a, int N);
double * I(double flux, double x_0, double x_2, double a, int N);
double * albedo_mu(double x_0, double x_2, double a, int N);



sp arcgen(double theta, double N_phi, int i);
double * x_ray_direction_diagram(double PSI_pr);



#endif //DEF_H
