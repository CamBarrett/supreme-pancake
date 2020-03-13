// header file for utility routines 
#include <iostream>
#include "stdlib.h"
#include  <cmath>
#include  <cstring>
//#include  <strings.h> // (GN Nov 15 2001) for SUN Blade 100
#include  <fstream>
using namespace std;

double wd_luminosity(double, double);

double wd_radius_Nauenberg(double);
double dlnm_dlnr_Nauenberg(double);

double wd_radius_ZS(double);
double dlnm_dlnr_ZS(double);

double wd_radius_PPE(double);
double dlnm_dlnr_PPE(double);

double wd_radius_semi_Savonije(double);
double dlnm_dlnr_Savonije(double);

double wd_radius_semi_TF(double);
double dlnm_dlnr_TF(double);

double wd_radius(double);
double dlnm_dlnr(double);

double wd_radius_semi(double);
double dlnm_dlnr_semi(double); 

double roche_lobe(double);
double zeta_L(double);

double get_rmin(double);

double eddington_limit(double, double, double, double);

double BP_SFR_disc(double, double);
//double BP_SFR_disc_int(double, double);
double BP_SFR_max(double);
double BP_SFR_int(double);
double BP_total_SFR(double);
double BP_total_SFR_disc(double);
double BP_total_SFR_bulge(double);

double birth_rate_coefficient(bool, int, double, int, double, double);
double next_birth_interval(double, double, double, double);
double next_birth_interval_BP(double, double, double, double*);
double next_birth_interval_constantSFR(double, double, double*);
//double OLD_next_birth_interval_BP(double, double, double, double*);

double BC_black_body(double);
double get_cross_section(double);
double X_ray_absorption(double, double);
double BB(double, double, double);
double Bnu(double, double, double);
void BB_trapezium(double, double, double, double*, int, double);
double BB_integrate(double, double, double, double);
double BB_int_simpson(double, double, double, double);
double BB_frac(double, double, double, double);

double get_A_V(double, double, double);
double EB_V(double, double);
double A_V_correction(double);
double get_N_H(double);

void random_position(double*, double*, double*);
double random_R_BP_disc(double);
void random_position_BP(double, double*, double*, double*);
void random_position_BP_GR(double, double*, double*, double*);
void random_position_halo(double, double*, double*, double*);
void random_position_around_Sun(double*, double*, double*);


double my_ran2(int&);
double random_number(double, double);
double gaussrand(double, double);

double bolometric_correction(double); 
double get_Mv(double, double);

