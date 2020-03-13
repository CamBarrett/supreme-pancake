// Utility functions for WDMS systems
// header file for utility routines 

#include "util.h"
#include "constants.h"


static double SFR_disc[16][13] = {
  {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
  {2.1400E+01, 2.1400E+01, 1.7490E+00, 6.6570E-01, 2.8459E-01, 1.3149E-01, 6.4415E-02, 3.3110E-02, 9.8793E-03, 3.3465E-03, 1.2592E-03, 4.0308E-04, 1.3619E-04},
  {5.2741E+01, 5.2741E+01, 3.2835E+00, 1.1407E+00, 4.5782E-01, 2.0215E-01, 9.5746E-02, 4.7957E-02, 1.3793E-02, 4.5833E-03, 1.7248E-03, 5.5220E-04, 1.8658E-04 },
  {1.8566E+02, 1.8566E+02, 1.0151E+01, 3.2820E+00, 1.2275E+00, 5.0583E-01, 2.3167E-01, 1.1051E-01, 3.0054E-02, 9.6826E-03, 3.7240E-03, 1.1927E-03, 4.0310E-04 },
  {4.1939E+02, 4.1939E+02, 3.4220E+01, 1.1431E+01, 4.2246E+00, 1.7242E+00, 7.6345E-01, 3.5981E-01, 9.4973E-02, 2.9762E-02, 1.1339E-02, 3.6367E-03, 1.2301E-03 },
  {4.0519E+02, 4.0519E+02, 6.4758E+01, 2.3491E+01, 9.4016E+00, 3.9277E+00, 1.7747E+00, 8.3396E-01, 2.2657E-01, 6.9737E-02, 2.6553E-02, 8.5404E-03, 2.8935E-03 },
  {3.7161E+02, 3.7161E+02, 7.8669E+01, 3.5347E+01, 1.5174E+01, 6.7169E+00, 3.2271E+00, 1.5787E+00, 4.2989E-01, 1.3695E-01, 5.2380E-02, 1.6724E-02, 5.6818E-03 },
  {2.9863E+02, 2.9863E+02, 8.5461E+01, 4.1794E+01, 2.0025E+01, 9.8152E+00, 4.8369E+00, 2.3978E+00, 6.7913E-01, 2.2340E-01, 8.5101E-02, 2.8139E-02, 9.5990E-03 },
  {2.2049E+02, 2.2049E+02, 8.3445E+01, 4.4238E+01, 2.2686E+01, 1.1554E+01, 6.0188E+00, 3.2392E+00, 9.3604E-01, 3.1813E-01, 1.2445E-01, 3.9695E-02, 1.3975E-02 },
  {1.2739E+02, 1.2739E+02, 6.7115E+01, 4.1712E+01, 2.4551E+01, 1.4001E+01, 7.7231E+00, 4.3227E+00, 1.4238E+00, 4.9025E-01, 1.9797E-01, 6.5754E-02, 2.2870E-02 },
  {7.0114E+01, 7.0114E+01, 5.2664E+01, 3.6050E+01, 2.3105E+01, 1.4220E+01, 8.5370E+00, 5.0712E+00, 1.7691E+00, 6.4403E-01, 2.6232E-01, 8.9170E-02, 3.2166E-02 },
  {3.4378E+01, 3.4378E+01, 3.5360E+01, 2.7245E+01, 1.9268E+01, 1.2973E+01, 8.4415E+00, 5.3780E+00, 2.1122E+00, 8.3154E-01, 3.5131E-01, 1.2299E-01, 4.4713E-02 },
  {1.7350E+01, 1.7350E+01, 2.0766E+01, 1.8165E+01, 1.4210E+01, 1.0456E+01, 7.3844E+00, 5.1058E+00, 2.2749E+00, 9.8425E-01, 4.3667E-01, 1.5936E-01, 5.9058E-02 },
  {1.0294E+01, 1.0294E+01, 1.4335E+01, 1.3521E+01, 1.1237E+01, 8.7075E+00, 6.4460E+00, 4.6259E+00, 2.2496E+00, 1.0344E+00, 4.7478E-01, 1.7860E-01, 6.7319E-02 },
  {8.0903E+00, 8.0903E+00, 9.1227E+00, 9.2050E+00, 8.1741E+00, 6.7149E+00, 5.2418E+00, 3.9979E+00, 2.1089E+00, 1.0432E+00, 4.9982E-01, 1.9573E-01, 7.5184E-02 },
  {6.4642E+00, 6.4642E+00, 5.7883E+00, 6.0321E+00, 5.6678E+00, 4.9212E+00, 4.0442E+00, 3.2010E+00, 1.8702E+00, 9.9303E-01, 4.9952E-01, 2.0530E-01, 8.1734E-02 }
};

static double SFR_R_disc[13] = {0., 1.3,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0, 10.0, 12.0, 14.0, 16.5, 19.0};

static double SFR_t[16] = {0., .012, .032, .101, .301, .604, 1.024, 1.520, 2.087, 3.045, 4.126, 5.612, 7.504, 9.125, 11.016, 13.583};

static double SFR_tot_disc[16] = {0, .377422, 0.874993, 2.99774, 7.21658, 8.35363, 8.93403, 8.67373, 7.91223, 6.58928, 5.46237, 4.30698, 3.2699, 2.69487, 2.16498, 1.69729};

static double SFR_tot_bulge[16] = {0., 0.34571, 0.838459, 2.93298, 6.77436, 6.96333, 6.64068, 5.63008, 4.42544, 2.80634, 1.75212, 0.984976, 0.535645, 0.344391, 0.242374, 0.17385 };

static double SFR_tot_int[16] = {0., 4.339e+06, 2.870e+07, 2.924e+08, 2.285e+09, 6.725e+09, 1.321e+10, 2.062e+10, 2.817e+10, 3.858e+10, 4.756e+10, 5.686e+10, 6.546e+10, 7.101e+10, 7.616e+10, 8.165e+10};

static double SFR_max[16] = { 0, 2.870e-08, 7.048e-08, 2.478e-07, 5.624e-07, 5.529e-07, 5.147e-07, 4.249e-07, 3.281e-07, 2.122e-07, 1.581e-07, 1.098e-07, 7.301e-08, 5.619e-08, 4.100e-08, 2.953e-08};

static int z_sign = 1;


//New cooling using Bergeron's models
//double wd_luminosity(double m, double t) {
//
//  // (GN Apr  8 1999) Driebe cooling
//  // //  double fit_mass = min(0.5, max(0.2, m));
//  //  //  double l_max_fit = 3. - 7.* (fit_mass - 0.2);
//  //  double fit_mass = min(0.6, max(0.2, m));
//  //  double l_max_fit = 3.83 - 4.77* fit_mass;
//  //  double coeff = 1.4;
//
//  //  // low mass guys loose envelope in thermal flashes if short period
//  //  if (m < 0.3) l_max_fit = 1.6;
//
//  // `Hansen' cooling
//  double l_max_fit = 1. - (0.9 - m);
//  if (m < 0.5) l_max_fit = 1.4 - 1.33 * (0.45 - m);
//  double coeff = 1.3333;
//
//  double l_max = 2.;
//  if (m < 0.5 ) l_max = -0.5;
//    
//  double l = min(l_max, l_max_fit - coeff * (log10(t) - 6.));
//
//  return l;
//}


//Nauenberg 1972 white dwarf radii
double wd_radius_Nauenberg(double m) {

  double mu = 2.;
  double m_rel23 = pow(m/CHANDRASEKHAR_MASS, TWO_THIRD);

  return (0.0225/mu) * sqrt(1./m_rel23 - m_rel23);

}

double dlnm_dlnr_Nauenberg(double m) {

  double m_rel23  = pow(m/CHANDRASEKHAR_MASS, -TWO_THIRD); 
  double beta   = ONE_THIRD - TWO_THIRD * (m_rel23/(m_rel23 - 1./m_rel23));

  return beta;
}


// Zapolsky & Salpeter for T=0 helium degenerates
//  return 8.76e8*pow(m,-0.3333)/R_SUN;
// Fit to eq. (3) of Warner (1995), ApSS, 225, 249
// which comes from Rappaport & Joss 1984, ApJ, 283, 232
// see also wdradius_RJ84.xmgr and wdradius.mws
double wd_radius_ZS(double m) {
  
  double r = 0.0106 - 0.0064 * log(m) + 0.0015 * pow(m, 2.);
  return r;

}

double dlnm_dlnr_ZS(double m) {

  double beta   = (-0.0064 + 0.003 * pow(m, 2.))/wd_radius_ZS(m);

  return beta;

}



// Eggleton 1986 fit to Nauenberg for high m and
// ZS for low m. From Rappaport & Verbunt 1988
double wd_radius_PPE(double m) {

  if (m >= 1.44) return 0.001;

  double fac1 = pow(m/CHANDRASEKHAR_MASS, 2./3.);
  double fac2 = 0.00057/m;
  double a = 3.5;
  double b = 1.;
  double two_third = 2./3.;

  double r = 0.0114 * sqrt(1./fac1 - fac1)
           * pow(1. + a*pow(fac2, 2./3.) + b*fac2, -two_third);

  return r;
      
}

double dlnm_dlnr_PPE( double m) {

  // From Tom Marsh, after Eggleton in Rappaport & Verbunt 1988 
 
  float fac1 = pow(m/CHANDRASEKHAR_MASS, 4./3.);
  double fac2 = 0.00057/m;
  double two_third = 2./3.;
  float fac3 = pow(fac2, two_third);
  float a = 3.5;
  float b = 1.;

  return -(1. + fac1) / (1. - fac1) / 3.+
    2.*(2.*a*fac3/3.+2*b*fac2)/(1.+a*fac3+b*fac2)/3.;

}



// Savonije et al 1986 for semi degenerate Helium stars
double wd_radius_semi_Savonije(double m) {

  return 2.03e9 * pow(m, -0.19)/SOLAR_RADIUS;
}

double dlnm_dlnr_Savonije(double m) {

  return -0.19;

}



// Fit to model 1.1 of Tutukov & Fedorova 1989
double wd_radius_semi_TF(double m) {

  return 2.99e9 * pow(m, -0.062)/SOLAR_RADIUS;
}

double dlnm_dlnr_TF(double m) {

  return -0.062;

}


double wd_radius(double m) {

  return wd_radius_PPE(m);

}


double wd_radius_semi(double m) {

  return wd_radius_semi_TF(m);

}




double dlnm_dlnr(double m) {

  return dlnm_dlnr_PPE(m);

}

double dlnm_dlnr_semi(double m) {

  return dlnm_dlnr_PPE(m);

}



//R_L(m) for q = m/M
double roche_lobe(double q) { 

  double q1_3 = pow(q, ONE_THIRD);
  double q2_3 = q1_3*q1_3;

  return 0.49*q2_3 / (0.6*q2_3 + log(1 + q1_3));

}

// zeta_L for R_L(m) if q = m/M
double zeta_L(double q) { 

  double q1_3 = pow(q, ONE_THIRD);
  double zeta = (1 + q)/3;
  zeta *= (2.*log(1 + q1_3) - q1_3/(1 + q1_3));
  zeta /= (0.6*q1_3*q1_3 + log(1 + q1_3));

  return zeta;

}

// Fit to Lubow & Shu 1975
double get_rmin(double q) {

  //A[0] = 0.0494792559
  //A[1] = 0.0381519871
  //A[2] = 0.0475218925
  //A[3] = 0.00697279913
  //R square = 0.9999623
  //Avg Y    = 0.1063


  double rmin = 0.04948 + 0.03815 * log10(q) 
          + 0.04752 * pow(log10(q), 2.) 
          + 0.006973 * pow(log10(q), 3.);

  return rmin;

}



// Eddington limit for star1
// From Han & Webbink 1999, A&A 349, L20
// solar units
double eddington_limit(double R1sun, double M1sun, double m2sun, double asun) {

  double M1   = M1sun * SOLAR_MASS;
  double m2   = m2sun * SOLAR_MASS;
  double R1   = R1sun * SOLAR_RADIUS;
  double a    = asun  * SOLAR_RADIUS;

  double q    = m2/M1;
  double q1_3 = pow(q, ONE_THIRD);
  double x    = a * (0.696 * q1_3 - 0.189 * q1_3*q1_3)
              / (1 + 0.014 * q);

  double mu   = M1/(M1 + m2);

  double phi_L1 = G * (-m2/x - M1/(a - x) 
	        - (M1 + m2) * pow((x - mu * a), 2.) / (2. * pow(a, 3.)));
		
  double phi_R1 = G * (-m2/a - M1/R1 
	        - (M1 + m2) * (2. * pow(R1, 2.)/3. + pow((a - mu * a), 2.)) 
	        / (2. * pow(a, 3.)));
		
  double X      = 0.; // Hydrogen mass fraction
  double kappa  = 0.2 * (1 + X);
			
  double g      = G * (M1/ pow(R1, 2.)  
		       - 2. * R1 * (M1 + m2) / (3. * pow(a, 3.)));

  double Md_Edd = 4. * PI * pow(R1, 2.) * C * g
                / (kappa * (phi_L1 - phi_R1));

  return Md_Edd * 3.1557e+7 / SOLAR_MASS;

}


// (GN Oct 17 2002) SFR(t, R) according to
// Boissier & Prantzos
double BP_SFR_disc(double t, double R) {


  t /= 1.e9;
  if (t > SFR_t[15] || R > SFR_R_disc[12]) {

    cerr << "ERROR in BP_SFR t or R too large" << endl;
    PRC(t);PRL(R);
    return 0.;
  }    

  // SFR_disc[t=16][R=13]

  int i = 0, j = 0;

  for (i = 15; i >= 0; i--) {
    
    //PRC(i);PRC(t);PRL(SFR_t[i]);
    if (t >= SFR_t[i]) break;

  }

  for (j = 12; j >= 0; j--) {
    
    //PRC(j);PRC(R);PRL(SFR_R_disc[j]);
    if (R >= SFR_R_disc[j]) break;

  }

  //PRC(i);PRL(j);

  double SFR_11 = SFR_disc[i][j];
  double SFR_22 = SFR_disc[i+1][j+1];
  double SFR_12 = SFR_disc[i][j+1];
  double SFR_21 = SFR_disc[i+1][j];

  //  PRC(SFR_11);PRL(SFR_21);
  //  PRC(SFR_12);PRL(SFR_22);

    
  double SFR_1 = SFR_11 
               + (t - SFR_t[i])/(SFR_t[i+1] - SFR_t[i]) * (SFR_21 - SFR_11);

  double SFR_2 = SFR_12 
               + (t - SFR_t[i])/(SFR_t[i+1] - SFR_t[i]) * (SFR_22 - SFR_12);

  //  PRL((t - SFR_t[i])/(SFR_t[i+1] - SFR_t[i]));
  //  PRC(SFR_1);PRL(SFR_2);

  //  PRL((R - SFR_R_disc[j])/(SFR_R_disc[j+1] - SFR_R_disc[j]));
  
  double SFR_Gyr = SFR_1 
                 + (R - SFR_R_disc[j])/(SFR_R_disc[j+1] - SFR_R_disc[j]) 
                 * (SFR_2 - SFR_1);

  return SFR_Gyr/1.e9; // in Msun/ pc^2  yr

}


double BP_SFR_int(double t) {

  t /= 1.e9;
  if (t > SFR_tot_int[15]) {

    cerr << "ERROR in BP_SFR_int t too large" << endl;
    PRL(t);
    return 0.;
  }    

  int i = 0;

  for (i = 15; i >= 0; i--) {
    
    //    PRC(i);PRC(t);PRL(SFR_t[i]);
    if (t > SFR_t[i]) break;

  }


  double SFR_int = SFR_tot_int[i] 
                   + (t - SFR_t[i])/(SFR_t[i+1] - SFR_t[i]) 
                   * (SFR_tot_int[i+1] - SFR_tot_int[i]);

  //PRC(SFR_t[i+1] - SFR_t[i]);

  //PRC(i);PRC(t-SFR_t[i]);PRC(SFR_tot_int[i]);PRL(SFR_int);
  return SFR_int;

}


double BP_SFR_max(double t) {

  t /= 1.e9;
  if (t > SFR_t[15]) {

    cerr << "ERROR in BP_SFR_int t too large" << endl;
    PRL(t);
    return 0.;
  }    

  int i = 0;

  for (i = 15; i >= 0; i--) {
    
    //    PRC(i);PRC(t);PRL(SFR_t[i]);
    if (t > SFR_t[i]) break;

  }

  double SFR_m = SFR_max[i] 
                   + (t - SFR_t[i])/(SFR_t[i+1] - SFR_t[i]) 
                   * (SFR_max[i+1] - SFR_max[i]);

  return SFR_m;

}


double BP_total_SFR_disc(double t) {


  t /= 1.e9;
  if (t > SFR_t[15]) {

    cerr << "ERROR in BP_total_SFR_disc t too large" << endl;
    PRL(t);
    return 0.;
  }    

  int i = 0;

  for (i = 15; i >= 0; i--) {
    
    //    PRC(i);PRC(t);PRL(SFR_t[i]);
    if (t > SFR_t[i]) break;

  }

  double SFR_total = SFR_tot_disc[i] 
                   + (t - SFR_t[i])/(SFR_t[i+1] - SFR_t[i]) 
                   * (SFR_tot_disc[i+1] - SFR_tot_disc[i]);

  return SFR_total;

}

double BP_total_SFR_bulge(double t) {


  t /= 1.e9;
  if (t > SFR_t[15]) {

    cerr << "ERROR in BP_total_SFR_bulge t too large" << endl;
    PRL(t);
    return 0.;
  }    

  int i = 0;

  for (i = 15; i >= 0; i--) {
    
    //    PRC(i);PRC(t);PRL(SFR_t[i]);
    if (t > SFR_t[i]) break;

  }

  double SFR_total = SFR_tot_bulge[i] 
                   + (t - SFR_t[i])/(SFR_t[i+1] - SFR_t[i]) 
                   * (SFR_tot_bulge[i+1] - SFR_tot_bulge[i]);

  return SFR_total;

}


double BP_total_SFR(double t) {


  return BP_total_SFR_disc(t) + BP_total_SFR_bulge(t);


}

double birth_rate_coefficient(bool massive, int function, double power, 
			      int number_of_points, double binary_fraction, double triple_fraction) {
  double f = binary_fraction;
  double f_t = triple_fraction;
  double f_p = (1. - f - f_t);
  double fraction_in_interval = 0, mean_mass = 0;

  switch (function) {

    case 1: // Power law

        if (power == -2.35) {
            cerr <<"ERROR: in powerlaw " << endl;
            cerr <<"please extend birth_rate_coefficient(..) " << endl;
            exit(-1);
            
            //values need to be checked
            mean_mass = f*0.58 + (1. - f)*0.35;
            if (massive) { fraction_in_interval = 0.0026; }
            else {  fraction_in_interval = 0.045; }
            
        } else if (power == -2.5) {
            //values need to be checked	 
            mean_mass = f*0.49 + (1. - f)*0.29;
            if (massive) { fraction_in_interval = 0.0014; }
            else {  fraction_in_interval = 0.033; }
	 
       } else {
            cerr <<"ERROR: power in powerlaw not -2.35 or -2.5! " << endl;
            cerr <<"please extend birth_rate_coefficient(..) " << endl;
       }
       break;

    case 2:  // Miller-Scalo 
        cerr <<"ERROR: in Miller-Scalo " << endl;
        cerr <<"please extend birth_rate_coefficient(..) " << endl;
        exit(-1);
        
        //values need to be checked 
        mean_mass = f*1.0 + (1. - f)*0.66;
        if (massive) { fraction_in_interval = 0.0063; }
        else {fraction_in_interval = 0.14;}
        break;

    case 3: // Scalo
        cerr <<"ERROR: Scalo IMF not yet implemented " << endl;
        cerr <<"please extend birth_rate_coefficient(..) " << endl;
        break;

    case 4: // Kroupa, Gilmore, Tout
//new normalization ./SeBa -f 4 -m 0.95 -M 10 -q 0
        mean_mass = f*0.74 + (1. - f)*0.49;
        if (massive) { fraction_in_interval = 0.0027; } //value need to be checked
        else {fraction_in_interval = 0.096;}
        break;
    
    case 5: //for it to suit triples
        double M_trip = 1.271;
        double M_bin = 0.895;
        double M_p = 0.593;
        double N_t = number_of_points;
        fraction_in_interval = 0.005466;
        double M_total = (N_t/fraction_in_interval) * (M_trip + (M_bin*f/(1.-f)) + (M_p*f_p/(f_t + f)));
        mean_mass = (f_t*M_trip) + (f*M_bin) + (f_p*M_p);
        
//old normalization ./SeBa -f 4 -m 0.96 -M 11 -q -1
//        mean_mass = f*0.79 + (1. - f)*0.49;
//        if (massive) { fraction_in_interval = 0.0027; }
//        else {fraction_in_interval = 0.095;}
        break;
  }

  double c = fraction_in_interval/(mean_mass * number_of_points);

  return c;

}



double rtsec_next_t(double time, double t_end, double coeff, double rand) {


  double eps = 0.01;
  double SFRint = BP_SFR_int(t_end - time);
  double t1 = time;
  double t2 = min(time + 10./(coeff * BP_total_SFR(t_end - time)), t_end-1.);

  double SFR1 = rand;
  double SFR2 = coeff * (BP_SFR_int(t_end - t2) - SFRint) + rand;
  
  //PRC(SFRint);PRC(rand);PRC(t2);PRL(SFR2);


  if (SFR2 > 0) return t_end; // last system

  double dx, swap, xl, rts;

  //PRC(time);PRC(t_end);PRC(t2);PRC(SFR1);PRL(SFR2);
  if (abs(SFR1) < abs(SFR2)) {
    
    rts = t1;
    xl = t2;
    swap = SFR1;
    SFR1 = SFR2;
    SFR2 = swap;

  } else {
   
    xl = t1;
    rts = t2;
  }

  for (int j = 1; j < 100 ; j++) {
    
    dx = (xl - rts) * SFR2 / (SFR2 - SFR1);
    xl = rts;
    SFR1 = SFR2;
    rts += dx;
    //PRC(dx);PRL(rts);
    if (rts >= t_end) rts = t_end - 1.;
    if (rts < time) rts = time;
    SFR2 = coeff * (BP_SFR_int(t_end - rts) - SFRint) + rand ;
    //PRC(rts);PRC(dx);PRC(coeff * (BP_SFR_int(t_end - rts) - SFRint));PRL(SFR2);
    if (abs(dx) < eps || SFR2 == 0 ) return rts;
  }
  return rts;
}


double next_birth_interval_BP(double time, double t_end, 
			      double coeff, double* t_random) {

  double rand = random_number(0.,1.);

  *t_random = rtsec_next_t(time, t_end, coeff, rand);

  double t = rtsec_next_t(time, t_end, coeff, 1.);

  if (*t_random == t_end) *t_random = 0.5 * (time + t_end);
    //cerr<< rand << "  " << t << "  " << *t_random<< " " << time << " " << t_end<< endl;
  return t;

}

double next_birth_interval_constantSFR(double time, double coeff, double* t_random) {

  //assuming a 1Msun/yr SFR
  //if a SFR of xMsun/yr is desired, multiply coeff with x

  double rand = random_number(0.,1.);
  *t_random = time + rand/coeff;  
  double t = time + 1./coeff;
  return t;
  
}


double next_birth_interval(double time, double t_end,
			   double coeff, double tau) {

  return tau * log(exp(t_end/tau)/(coeff*tau) + exp(time/tau));

}


double BC_black_body(double temp) {


  double bc = 0;
  temp *= 1000; // temp now in Kelvin


  // see /home/nelemans/Home/Work/SeBa/Fit/BC.xmgr
  if (temp < 7500) {
//y = -22.916 + 0.01553 x - 4.0937e-06 x^2 + 4.8969e-10 x^3 - 2.2114e-14 x^4


    bc =  -22.916 + temp * (0.01553
         	  + temp * (-4.0937e-06
		  + temp * (4.8969e-10
		  + temp * -2.2114e-14)));
      
  } else if (temp < 20000) {
	
    //Regression coefficient (SLOPE)		 = -0.0001309091
    //Regression constant (INTERCEPT)		 = 1.007818
    
    bc = 1.007818 - 0.000130909 * temp; // temp in K
    
  } else {

    double Emin = 2.066; // \lambda = 6000 \AA in eV
    double Emax = 2.478; // \lambda = 5000 \AA in eV
    double N_H = 0.;

    // 0.127516 is BB_frac for T_sun = 5775 to make BC(sun) = 0
    bc = -2.5 * log10(0.127516/BB_frac(temp, Emin, Emax, N_H)); 
    
  }
  
  return bc;
}

// B.C. for ms (Eggleton et al 1989)
// NOTE: temp in kK!!!!
double bolometric_correction(double temp) {

      double bc;

      if (temp<4.452)
         bc = 2.5*log10((6.859e-6*pow(temp,8.) + 9.316e-3)
                       /(1. + 5.975e-10*pow(temp,14.)));
      else if (temp>=4.452 && temp <=10.84)
         bc = 2.5*log10((3.407e-2*pow(temp,2.))
                       /(1. + 1.043e-4*pow(temp, 4.5)));
      else
         bc = 2.5*log10((2728/pow(temp,3.5) + 1.878e-2*temp)
                       /(1. + 5.362e-5*pow(temp,3.5)));
      return bc;
}  


double get_Mv(double lum, double R) {

  // temperature in kK !!
  double temperature = pow(1130. * lum/(R*R), 0.25);
  double BC = bolometric_correction(temperature);
  //  cerr << " t (1000 K) , bc " << t << " " << BC << endl; 
  double Mv = 4.74 - 2.5 * log10(lum) - BC;


  //  cerr << "temperature " << temperature;
  // (GN May 28 1999) test temperature limit
  //  if (temperature < 15.) Mv = 1000.;

  return Mv;
}



//interstellar photoelectric absorption cross-section.
//  Morrison and McCammon (1983, ApJ, 270, 119)
double get_cross_section(double energy) { // in keV


  if (energy > 10 || energy < 0.03 ) return 0.;



  static double emax[14] = { .1,.284,.4,.532,
			   .707,.867,1.303,1.84,2.471,
			   3.21,4.038,7.111,8.331,10. };
  static double c0[14] = { 17.3,34.6,78.1,71.4,
			 95.5,308.9,120.6,141.3,202.7,
			 342.7,352.2,433.9,629.,701.2 };
  static double c1[14] = { 608.1,267.9,18.8,66.8,
			 145.8,-380.6,169.3,146.8,104.7,
			 18.7,18.7,-2.4,30.9,25.2 };
  static double c2[14] = { -2150.,-476.1,4.3,
			   -51.4,-61.1,294.,-47.7,-31.5,
			   -17.,0.,0.,.75,0.,0. };

  int i = 1;
  for (i = 1; i <= 14; ++i) {
    
    if (energy < emax[i-1]) break;

  }

  i -= 1;
  return 1.e-24 * (c0[i] + energy*(c1[i] + energy*c2[i]))/pow(energy, 3);

}

double X_ray_absorption(double E, double N_H) { // N_H in cm^{-2}

  double cross_section = get_cross_section(E);
  PRL(cross_section);

  return exp(-N_H*cross_section);

}


double UV_absorption(double nu, double N_H) {

  double x = 1e-4*nu/C; //1/lambda in micron^-1
  

  if (x > 24)  return 1.;

  // linear fit to
  // Rumph et al. 1994 AJ 107 2108
  double cross_section = pow(10.,(24-x)/14.);
  double abs = exp(-N_H*cross_section);


  if (x > 8) return abs;


  double FMx = 0.;

  if (x >= 2.27) { // above B band

    //CCM
    double x0 = 4.595;
    double gamma = 1.051;
    double c1 = -0.384;
    double c2 = 0.739;
    double c3 = 3.961;
    double c4 = 0.265;
    
    double Fx = 0.;
    if (x > 5.9) Fx = 0.592*pow((x - 5.9),2.) + 0.0564*pow((x - 5.9),3.);
    
    double Dx = x*x/((x*x - x0*x0)*(x*x - x0*x0) + gamma*gamma*x*x);
  
    FMx = c1 + c2*x + c3*Dx + c4*Fx;
    PRC(x);PRL(FMx);

  }  else if (x < 1.82) { // below V band

    FMx = -3.0885 + 0.077775*x + 1.6061 *x*x - 0.39101*x*x*x;
    PRC(x);PRL(FMx);

  } else { // linear between V and B

    FMx = (x-1.82)/0.45;
    PRC(x);PRL(FMx);
  }

  double A_V = 5.59e-22 * N_H; // N_H = 0.179 A_V in 1e22 units
  double R_v = 3.1;

  double Alamb = A_V * (FMx/R_v + 1);

  PRC(x);PRL(abs);

  abs = pow(10, -0.4*Alamb);

  return abs;

}


double NIR_UV_absorption(double nu, double N_H) {


  double x = 1e-4*nu/C; //1/lambda in micron^-1
  PRL(x);

  if (x > 24.2 || x < 0.3)  return 1.;


  // from Ly break to 0.03 keV
  // linear fit to
  // Rumph et al. 1994 AJ 107 2108
  double cross_section = 1.245e-18*pow(10.,(24.2-x)/14.);
  double abs = exp(-N_H*cross_section);

  PRC(cross_section);PRL(abs);

  if (x > 11) return abs;

  // below Ly break

  double ax = 0.574*pow(x,1.61);
  double bx = -0.527*pow(x,1.61);
  
  if (x > 8) {

    double y = (x-8);
    ax = -1.073 - 0.628 * y + 0.137 * y*y - 0.070 * y*y*y;
    bx = 13.670 + 4.257 * y - 0.420 * y*y + 0.374 * y*y*y;

  } else if (x > 3.3) {

    double Fax = 0.;
    double Fbx = 0.;

    if (x > 5.9) {
      
      double y = x - 5.9;

      Fax = -0.04473 * y*y - 0.009779 * y*y*y;
      Fbx = 0.2130 * y*y + 0.1207 * y*y*y;
    }

    ax = 1.752 - 0.316*x -0.104/((x-4.67)*(x-4.67) + 0.341) + Fax;

    bx = -3.090 + 1.825*x + 1.206/((x-4.62)*(x-4.62) + 0.263) + Fbx;
    
  } else if (x > 1.1) {

    double y = (x - 1.82);

    ax = 1 + 0.1799*y - 0.50447 * y*y - 0.02427 * y*y*y + 0.72085 * pow(y, 4.)
      + 0.01979 * pow(y, 5.) - 0.77530 * pow(y, 6.) + 0.32999 * pow(y, 7.);

    bx = 1.41338 * y + 2.28305 * y*y + 1.07233 * y*y*y - 5.38434 * pow(y, 4.)
      - 0.62251 * pow(y, 5.) + 5.30260 * pow(y, 6.) - 2.09002 * pow(y, 7.);

  } 

  double A_V = 5.59e-22 * N_H; // N_H = 0.179 A_V in 1e22 units
  double R_v = 3.1;
  
  double Alamb = A_V * (ax + bx/R_v);

  abs = pow(10, -0.4*Alamb);

  
  return abs;
}


double BB(double Temp, double nu, double N_H) {


  // Fro consistency with old version

  return PI * Bnu(Temp, nu, N_H);

}


double Bnu(double Temp, double nu, double N_H) {

  double h = 6.6260755e-27;
  double k = 1.380658e-16;

  double E = h*nu/1.6021772e-9; // in keV
  double abs1 = X_ray_absorption(E, N_H);
  double abs2 = NIR_UV_absorption(nu, N_H);

  PRC(nu);PRC(abs1);PRL(abs2);
  
  double abs = min(abs1,abs2);

  double Bnu =  2*h/(C*C) * nu*nu*nu / (exp(h*nu/(k*Temp)) - 1.);

  return abs * Bnu;

}


// From Num Rec trapzd
void BB_trapezium(double Temp, double nu_min, double nu_max, 
			double* s, int n, double N_H) {

  if (n == 1) {



    //    PRC(Temp);PRL(nu_min);PRL(BB(Temp, nu_min, N_H));
    //    PRC(Temp);PRL(nu_max);PRL(BB(Temp, nu_max, N_H));

    *s = 0.5 * (nu_max - nu_min) 
             * (BB(Temp, nu_min, N_H) + BB(Temp, nu_max, N_H));

  } else {

    double two = 2.;
    int it = pow(two, n-1);
    double tnm = it;
    double del = (nu_max - nu_min)/tnm;
    double x = nu_min + 0.5*del;
    double sum = 0;

    for (int i = 1; i <= it; i++) {
      
      sum += BB(Temp, x, N_H);
      x += del;

    }

    *s = 0.5 * (*s + (nu_max - nu_min)*sum/tnm);

  }
  
}
  


// return int_{\nu_min}^{\nu_max} B(\nu) d \nu
double BB_integrate(double Temp, double nu_min, double nu_max, double N_H) {

  int JMAX = 30;
  double EPS = 1e-2;

  double olds = -1.e30;
  double s = 0;

  for (int j = 1; j <= JMAX; j++) {

    //PRC(j);
    BB_trapezium(Temp, nu_min, nu_max, &s, j, N_H);

    //PRC(olds);PRL(s);

    if (abs(s - olds) < EPS*abs(olds)) return s;
    if (s == 0 && olds == 0 & j > 6) return s;

    olds = s;

  }

  cerr << "ERROR: to many steps in BB_integrate!!!" << endl;

}

// return int_{\nu_min}^{\nu_max} B(\nu) d \nu
double BB_int_simpson(double Temp, double nu_min, double nu_max, double N_H) {

  int JMAX = 30;
  double EPS = 1e-2;

  double olds = -1.e30;
  double oldst = -1.e30;
  double st = 0;

  for (int j = 1; j <= JMAX; j++) {

    //PRC(j);
    BB_trapezium(Temp, nu_min, nu_max, &st, j, N_H);


    double s = (4.*st - oldst)/3.;

    //PRC(olds);PRL(s);

    //PRC(s - olds);PRC(abs(s - olds));PRC(olds);PRL(abs(olds));
    if (abs(s - olds) < EPS*abs(olds)) return s;
    if (s == 0 && olds == 0 & j > 6) return s;

    olds = s;
    oldst = st;

  }

  cerr << "ERROR: to many steps in BB_integrate!!!" << endl;

}


double BB_frac(double Temp, double Emin, double Emax, double N_H) {


  double sigma = 0.0000567051; //Stefan-Boltzmann constant 
  double eV_to_Hz = 2.41798814396244e14;
  double nu_min = eV_to_Hz*Emin;
  double nu_max = eV_to_Hz*Emax;

  double BB_in_band = BB_int_simpson(Temp, nu_min, nu_max, N_H);
  //double BB_in_band = BB_integrate(Temp, nu_min, nu_max, N_H);

  double BB_total = sigma*pow(Temp, 4.);

  double frac = BB_in_band/BB_total;

  //PRC(BB_in_band);PRC(BB_total);PRL(frac);

  return frac;

}




double get_A_V(double l, double b, double dist) {

  // (GN Dec 15 2006) 0.3  // in kpc see Nelemans et al 2001, A&A 365, 491 
  // (GN Mar 26 2007) dust is more concentrated!
  double h = 0.12; 

  double sinb = sin(abs(b)/rad_to_deg);

  // From Sandage (1972, ApJ, 178, 1)  
  double tanb = tan(abs(b)/rad_to_deg);
  double tan50 = 1.19175359; // tan 50 deg

  if (abs(b/rad_to_deg) > 0.87266463) return 0.;

  double A_V_inf = 0.165 * (tan50 - tanb) / sinb;
 
  // (GN Mar 11 2003) quick fix to increase absorption
  // towards the bulge, loosely based on SFK
  // increase A_V witing 1 deg
  //if (abs(l) < 1.) A_V_inf *= (2. - abs(l));

#if 0
  // From Slegel maps
  double A_V_inf = 3.315 * EB_V(l, b);
#endif 

  double h_max = min(h,sinb*23.5); //for very low b sample full dust layer

  return A_V_inf * tanh(dist*sinb/h_max);

}


// E(B-v) from SFK 1998 model
//double EB_V(double l, double b) {

//  float myl = (float)l;
//  float myb = (float)b;

//  float * EBV = lambert_getval("/home/nelemans/Compute/C++/AMCVn/SFD_dust_4096_ngp.fits",
//			       "/home/nelemans/Compute/C++/AMCVn/SFD_dust_4096_sgp.fits", 
//			       1, &myl, &myb, 0, 0, 0);

//  return *EBV;


//}


// from 6th order fit to fractional differences between
// Sandage and SFK 1998 model averaged over b
double A_V_correction(double l) {

  double x = 57.29577951*l; // in degrees

  double A1 = 3.9413;
  double A2 =-0.2204;
  double A3 = 0.004346;
  double A4 =-3.7895e-05;
  double A5 = 1.6213e-07;
  double A6 =-3.3668e-10;
  double A7 = 2.7322e-13;

  double correction = A1 + x*(A2 + x*(A3 + x*(A4 + x*(A5 + x*(A6 + x*A7)))));

  //  PRC(x);PRL(correction);

  return correction;

}

// calculate N_H from A_V Predehl & Schmitt, 1995, A&A 293, 889
double get_N_H(double A_V) {


  return 0.179 * A_V * 1.e22; // in cm^{-2}

}

// random position in the Galaxy with distribution:
// exp(-R/HR) sech^2(z/hz) see Nelemans et al. 2001a,b
void random_position(double* l, double* b, double* dist) {

  const double  HR    = 2.5;   // in kpc 
  const double  R_sun = 8.5;  // in kpc
  const double  hz    = 0.2;   // in kpc

  //    real R = -HR * log(random_number(0.,1.));
  // (GN Oct 23 2002) 
  // appox P(R) = R exp(-R/H)
  double R = pow(-13. * log(random_number(0.,1.)), 20./31.);
  double phi = random_number(0., 2 * PI);
  double X   = random_number(0.,1.);
  double z   = hz * log(-sqrt(-(X - 1) * (X + 1)) / (X - 1));
  double x   = R * cos(phi);
  double y   = R * sin(phi);
  
  double d_proj = sqrt(pow(x, 2) + pow(y-R_sun, 2));

  *dist  = sqrt(pow(x, 2) + pow(y-R_sun, 2) + pow(z, 2));    
  *l     = rad_to_deg * atan(x/(y - R_sun));
  *b     = rad_to_deg * atan(z/d_proj);
  
  if (y > R_sun)  *l -= (x/abs(x)) * 180.;

}


// gausrand: Return a random number distributed according to a Gaussian
//	     distribution with specified mean and standard deviation.
// (GN Oct 31 2002) from Starlab

#define BIG 10

double  gaussrand(double mean, double sdev)
{
    // Select a number with zero mean and unit variance using the
    // rejection method (see Numerical Recipes).

    for(;;) {
	double x = random_number(-BIG, BIG);
	if (random_number(0,1) < exp(-0.5*x*x)) return mean + x*sdev;
    }
}
#undef BIG


// (GN Jun  1 2007) new try: rejection method...
double random_R_BP_disc(double t_gal) {

  double R_max = SFR_R_disc[12];
  double SFR_max = 1.1 * BP_SFR_max(t_gal);
  double R = 0.;
  double SFR = 0.;

  do {
    R = random_number(0.,R_max);
    SFR = random_number(0.,SFR_max);

  } while (SFR > R * BP_SFR_disc(t_gal, R));

  return R;

}
  

// random position in the Galaxy with distribution:
// accroding to SFR of Boissier & Prantzos with
// P(z) =  sech^2(z/hz) 
void random_position_BP(double t_gal, double* l, double* b, double* dist) {

  const double  HR    = 2.5;   // in kpc 
  const double  R_sun = 8.5;  // in kpc
  double  hz    = 0.3;   // in kpc, e.g. Sparke & Gallagher, p.64


  //(GN Sept 2011) EXPERIMENT larger scale heigt for old systems!!!
  // experiment with thick disk. Are all old stars in the thick disk?  
  //  if (t_gal < 6.5e9) hz = 1.2;
  // !!!!!!!!!!!!!!

  double f = random_number(0.,1.);
  double x = 0, y = 0, z = 0;
  
  //PRL(t_gal);
  //PRC(f);PRC(f * BP_total_SFR(t_gal));PRL(BP_total_SFR_bulge(t_gal));


  if (f * BP_total_SFR(t_gal) < BP_total_SFR_bulge(t_gal)) { // system in bulge
    

    // (GN Aug 26 2010) BUG!!!! forgot the r^2 volume factor. 
    // Solution: 3 gaussion dists in x, y and z... 
    //double r = gaussrand(0., 0.5);
    // double phi = random_number(0., 2 * PI);
    //double theta = acos(random_number(-1., 1.));
    //
    //z = r * cos(theta);
    //double R = r * sin(theta);
    //
    //x = R * cos(phi);
    //y   = R * sin(phi);

    x = gaussrand(0., 0.5);
    y = gaussrand(0., 0.5);
    z = gaussrand(0., 0.5);


    //cerr << "Bulge: R,z = " << R << " " << z << endl; 

  } else { // system in disc

    double R   = random_R_BP_disc(t_gal);
    double phi = random_number(0., 2 * PI);
    double X   = random_number(0.,1.);
    
    //hz = 0.165 + 0.21*(R/R_sun - 5./8.); //Kent et al 1991, from Gal Astr. p611

    z   = z_sign * hz * log(-sqrt(-(X - 1) * (X + 1)) / (X - 1));
    x   = R * cos(phi);
    y   = R * sin(phi);
  
    z_sign = -z_sign;

    //cerr << "Disc: R,z = " << R << " " << z << endl; 

  }

  double d_proj = sqrt(pow(x, 2) + pow(y-R_sun, 2));

  *dist  = sqrt(pow(x, 2) + pow(y-R_sun, 2) + pow(z, 2));    
  *l     = rad_to_deg * atan(x/(y - R_sun));
  *b     = rad_to_deg * atan(z/d_proj);
  
  if (y > R_sun)  *l -= (x/abs(x)) * 180.;


}

// random position in the Galaxy with distribution:
// accroding to SFR of Boissier & Prantzos with
// P(z) =  sech^2(z/hz) 
void random_position_halo(double t_gal, double* l, double* b, double* dist) {

  const double  HR    = 2.5;   // in kpc 
  const double  R_sun = 8.5;  // in kpc
  double  hz    = 0.2;   // in kpc

  double f = random_number(0.,1.);
  double x = 0, y = 0, z = 0;
  
  //  PRL(t_gal);
  //  PRC(f);PRC(f * BP_total_SFR(t_gal));PRL(BP_total_SFR_bulge(t_gal));


    double r = gaussrand(0., 7.);
    double phi = random_number(0., 2 * PI);
    double theta = acos(random_number(-1., 1.));
    
    z = r * cos(theta);
    double R = r * sin(theta);

    x = R * cos(phi);
    y   = R * sin(phi);

    //    cerr << "Bulge: R,z = " << R << " " << z << endl; 

  double d_proj = sqrt(pow(x, 2) + pow(y-R_sun, 2));

  *dist  = sqrt(pow(x, 2) + pow(y-R_sun, 2) + pow(z, 2));    
  *l     = rad_to_deg * atan(x/(y - R_sun));
  *b     = rad_to_deg * atan(z/d_proj);
  
  if (y > R_sun)  *l -= (x/abs(x)) * 180.;

}

// random position in the Galaxy 
// in x pc around the Sun
void random_position_around_Sun(double* l, double* b, double* dist) {

    // x, y, z centred on the sun
    double x = 0, y = 0, z = 0;
    const double  R_sun = 8.5;  // in kpc
    const double r_sun = 0.05; // in kpc, max_distance_from_Sun 
      
    x = random_number(-r_sun, r_sun);
    y = random_number(-r_sun, r_sun);
    z = random_number(-r_sun, r_sun);
    
//    PRC(x);PRC(y);PRL(z);

    double d_proj = sqrt(pow(x, 2) + pow(y, 2));    
    *dist  = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));    
    *l     = rad_to_deg * atan(x/y);
    *b     = rad_to_deg * atan(z/d_proj);
  
  if (y > 0)  *l -= (x/abs(x)) * 180.;
}



// ran2 from Numerical Recipes, taken from Starlab
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.e-14
#define RNMX (1.0-EPS)

double my_ran2(int& idum)
{
    int j;
    int k;
    static int idum2 = 123456789;
    static int iy = 0;
    static int iv[NTAB];
    double temp;

    if (idum <= 0) {		// idum < 0 ==> initialize
	if (-(idum) < 1)
	    idum = 1;
	else
	    idum = -(idum);
	idum2 = (idum);

	for (j = NTAB+7; j >= 0; j--) {
	    k = (idum)/IQ1;
	    idum = IA1*(idum-k*IQ1) - k*IR1;
	    if (idum < 0) idum += IM1;
	    if (j < NTAB) iv[j] = idum;
	}
	iy = iv[0];
    }
    k = (idum)/IQ1;
    idum = IA1*(idum-k*IQ1) - k*IR1;
    if (idum < 0) idum += IM1;

    k = idum2/IQ2;
    idum2 = IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;

    j = iy/NDIV;
    iy = iv[j] - idum2;
    iv[j] = idum;
    if (iy < 1) iy += IMM1;

    if ((temp = AM*iy) > RNMX)
	return RNMX;
    else
	return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

double random_number(double a, double b) {

  //    static int randx = -281072;
  static int randx = time(0);

    return a + (b-a)*my_ran2(randx);

}

// (GN Oct 30 2006) for Gijs Roelofs SDSS article:
// only systems with phi = [
// random position in the Galaxy with distribution:
// accroding to SFR of Boissier & Prantzos with
// P(z) =  sech^2(z/hz) 
void random_position_BP_GR(double t_gal, double* l, double* b, double* dist) {

  const double  HR    = 2.5;   // in kpc 
  const double  R_sun = 8.5;  // in kpc
  double  hz    = 0.3;   // in kpc, e.g. Sparke & Gallagher, p.64

  double f = random_number(0.,1.);
  double x = 0, y = 0, z = 0;
  
  //PRL(t_gal);
  //PRC(f);PRC(f * BP_total_SFR(t_gal));PRL(BP_total_SFR_bulge(t_gal));


  if (f * BP_total_SFR(t_gal) < BP_total_SFR_bulge(t_gal)) { // system in bulge
    
    double r = gaussrand(0., 0.5);
    double phi = random_number(0., 2 * PI);
    double theta = acos(random_number(-1., 1.));
    
    z = r * cos(theta);
    double R = r * sin(theta);

    x = R * cos(phi);
    y   = R * sin(phi);

    //cerr << "x, y, phi =" << x << " " << y << " " << phi << endl;
    //cerr << "Bulge: R,z = " << R << " " << z << endl; 

  } else { // system in disc

    double R   = random_R_BP_disc(t_gal);
    double phi = random_number(0.5*PI-0.23, 0.5*PI+0.23);
    double X   = random_number(0.,1.);
    
    //hz = 0.165 + 0.21*(R/R_sun - 5./8.); //Kent et al 1991, from Gal Astr. p611

    z   = z_sign * hz * log(-sqrt(-(X - 1) * (X + 1)) / (X - 1));
    //z   = z_sign * random_number(0, 2.);
    x   = R * cos(phi);
    y   = R * sin(phi);
  
    z_sign = -z_sign;

    //cerr << "x, y, phi =" << x << " " << y << " " << phi << endl;
    //cerr << "Disc: R,z = " << R << " " << z << endl; 
    //cerr  << R << " " << z << " "  << x << " " << y << " " << phi << endl;


  }

  double d_proj = sqrt(pow(x, 2) + pow(y-R_sun, 2));

  *dist  = sqrt(pow(x, 2) + pow(y-R_sun, 2) + pow(z, 2));    
  *l     = rad_to_deg * atan(x/(y - R_sun));
  *b     = rad_to_deg * atan(z/d_proj);
  
  if (y > R_sun)  *l -= (x/abs(x)) * 180.;

}
