#include "my_main_sequence.h"
#include "constants.h"
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;

double base_main_sequence_radius(const double mass, const double z) { 
    double mx = pow(mass, 0.5);
    double teller = (smc.c(8,z)*pow(mass,2) + smc.c(9,z)*pow(mass,6))*mx + smc.c(10,z)*pow(mass,11) +(smc.c(11,z) +                               smc.c(12,z)*mx)*pow(mass,19);
    double noemer = smc.c(13,z) + smc.c(14,z)*pow(mass,2) + (smc.c(15,z)*pow(mass,8) + pow(mass,18) +                                             smc.c(16,z)*pow(mass,19))*mx;
    return teller/noemer;
}


double base_giant_branch_time(const double mass,
					 const double z) {

    double t_bgb;
    double pow_mass_7 = pow(mass, 7);
    double teller = smc.a(1, z) +smc.a(2, z)*pow(mass, 4) 
                      +smc.a(3, z)*pow(mass, 5.5)  
                      +       pow_mass_7;
    double noemer =  smc.a(4, z)*pow(mass, 2) +smc.a(5, z)*pow_mass_7; 
    t_bgb = teller/noemer;
  
    return t_bgb;
}


//Eq. 5
double main_sequence_time(const double mass, const double z) {

  double t_ms = base_giant_branch_time(mass, z)
            * max(stars_with_main_sequence_hook(mass, z),
                  stars_without_main_sequence_hook(mass, z)); 
 
  return t_ms;
}


double main_sequence_hook_time(const double mass,
					    const double z) { 

  double mu = stars_with_main_sequence_hook(mass, z);
  double t_bgb = base_giant_branch_time(mass, z);
  double t_hook = mu * t_bgb;

  return t_hook;
}


double main_sequence_hook_mass(const double z) {

  double zeta = log10(z/SOLAR_METALICITY);
  double m_hook = 1.0185 + zeta * (zeta * 0.0892 + 0.16015);

  return m_hook;
}


// Eq.6 identified as 'x' by Hurley
double stars_without_main_sequence_hook(const double mass, 
						    const double z) {
  
  double zeta = log10(z/SOLAR_METALICITY);
  double x = max(0.95, min(0.99,
                         0.95 - 0.03*(zeta + 0.30103)));
  
  return x;
} 

// Eq.7 identified as 'mu' by Hurley
double stars_with_main_sequence_hook(const double mass, 
						const double z) {
  
  double mu = max(0.5, 
		1.0 - 0.01*max(smc.a(6, z)/pow(mass, smc.a(7, z)),
			       smc.a(8, z) + smc.a(9, z)
			       /pow(mass, smc.a(10, z))));

  return mu;
}

//Eq.21a Eq.21b
double alpha_r_coefficient(const double mass, 
                                      const double z) {
    
    double alpha_r;
    if (mass<=smc.a(67, z) && mass>=smc.a(66, z)) {
        alpha_r = smc.a(58, z)*pow(mass, smc.a(60, z))
        / (smc.a(59, z) + pow(mass, smc.a(61, z)));
    }
    else if (mass>smc.a(67, z)) {
        
        double Cnstnt = alpha_r_coefficient(smc.a(67, z), z);
        alpha_r = (Cnstnt + smc.a(65, z)*(mass - smc.a(67, z)));
    }
    else if(mass<smc.a(66, z) && mass >=smc.a(68, z)) {
        
        double B = alpha_r_coefficient(smc.a(66, z), z);
        alpha_r = smc.a(64, z) + (B - smc.a(64, z))*(mass - smc.a(68, z))
        / (smc.a(66, z) - smc.a(68, z));
    }
    else if(mass<smc.a(68, z) && mass>=0.65) 
        alpha_r = smc.a(63, z) + (smc.a(64, z) - smc.a(63, z))*(mass - 0.65)
        / (smc.a(68, z) - 0.65);
    else if(mass<0.65 && mass>=0.50){
        alpha_r = smc.a(62, z) + (smc.a(63, z) 
                                  - smc.a(62, z))*(mass - 0.50) / 0.15;
    }
    else if(mass<0.50)
        alpha_r = smc.a(62, z);
    else {
        cerr << "WARNING: ill defined double "
        << "alpha_r_coefficient(const double mass, const double z) " << endl;
        PRC(mass);PRL(z);
        cerr << flush;
        return EXIT_FAILURE;
        //exit(1);
    }
    
    return alpha_r;
}

//Eq.22a Eq.22b
double beta_r_coefficient(const double mass, 
                                     const double z) {
    
    double beta_r;
    if (mass>=2 && mass<=16) {
        beta_r = smc.a(69, z)*pow(mass, 3.5)
        / (smc.a(70, z) + pow(mass, smc.a(71, z)));
    }
    else if (mass>16) {
        double Cnstnt = beta_r_coefficient(16, z)+1.0;      
        beta_r = Cnstnt + smc.a(73, z)*(mass - 16);
    }
    else if (mass<2 && mass>=smc.a(74, z)) { 
        double B = beta_r_coefficient(2, z)+1.0;
        beta_r = smc.a(72, z) + (B - smc.a(72, z))
        * (mass - smc.a(74, z))/(2 - smc.a(74, z));
    }
    else if (mass<smc.a(74, z) && mass>1)
        beta_r = 1.06 + (smc.a(72, z) - 1.06)*(mass - 1.0)/(smc.a(74, z) - 1.0);
    else if (mass<=1)
        beta_r = 1.06;
    else {
        cerr << "WARNING: ill defined double "
        << "beta_r_coefficient(const double mass, const double z) " << endl;
        cerr << flush;
        return EXIT_FAILURE;
        //exit(1);
    }
    return beta_r - 1;
}

//Eq.23 
double gamma_r_coefficient(const double mass, 
                                      const double z) {
    double gamma_r;
    if (mass<=1) {
        
        gamma_r = smc.a(76, z) + smc.a(77, z)
        * pow(mass - smc.a(78, z), smc.a(79, z));
    }
    else if (mass<=smc.a(75, z)) {
        
        double B = gamma_r_coefficient(1, z);
        gamma_r = B + (smc.a(80, z) - B)
        * pow( (mass - 1)/(smc.a(75, z) - 1), smc.a(81, z));
    }
    else if (mass<smc.a(75, z) + 0.1) {
        double Cnstnt;
        if (smc.a(75, z)<=1)
            Cnstnt = gamma_r_coefficient(1, z);
        else
            Cnstnt = smc.a(80, z);
        
        gamma_r = Cnstnt - 10*(mass - smc.a(75, z))*C;
    }
    else {
        gamma_r = 0;
    }
    
    if (mass>smc.a(75, z) + 0.1)
        gamma_r = 0;
        
    return max(0., gamma_r);
}



double ap(double zeta, double a, double b, double c, double d, double e) {

  double ai = a + zeta*(b + zeta*(c + zeta*(d + zeta*e)));

  return ai;
}

double stellar_model_constants::a(int index, double z) {

  double zeta = log10(z/SOLAR_METALICITY);

  double a = 0;
  switch(index) {
     case 1: a = ap(zeta, 1593.890, 2053.038, 1231.226, 232.7785);
             break;
     case 2: a = ap(zeta, 2.706708E+3, 1.483131E+3, 5.772723E+2, 7.411230E+1);  
             break;
     case 3: a = ap(zeta, 1.466143E+2, -1.048442E+2, -6.795374E+1, -1.391127E+1);
             break;
     case 4: a = ap(zeta, 4.141960E-2, 4.564888E-2, 2.958542E-2, 5.571483E-3);
             break;
     case 5: a = ap(zeta, 3.426349E-1);
             break;
     case 6: a = ap(zeta, 1.949814E+1, 1.758178E+0, -6.008212E+0, -4.470533E+0); 
             break;
     case 7: a = ap(zeta, 4.903830E+0);
             break;
     case 8: a = ap(zeta, 5.212154E-2, 3.166411E-2, -2.750074E-3, -2.271549E-3); 
             break;
     case 9: a = ap(zeta, 1.312179E+0, -3.294936E-1, 9.231860E-2, 2.610989E-2);
             break;
     case 10: a = ap(zeta, 8.073972E-1);
              
             break;
     case 17: {double s = log10(z);
               double log_a17 = max(0.097 - 0.1072*(s+3),
				 max(0.097, 
				     min(0.1461, 0.1461 + 0.1237*(s+2))));
              a = pow(10., log_a17);
     }
             break; 
     case 26: a = ap(zeta, 5.502535E+0, -6.601663E-2, 9.968707E-2, 3.599801E-2);
             break;
     case 27: a = ap(zeta, 9.511033E+1, 6.819618E+1, -1.045625E+1, -1.474939E+1);
             break;
     case 28: a = ap(zeta, 3.113458E+1, 1.012033E+1, -4.650511E+0, -2.463185E+0);
             break;
     case 29: {double a29 = ap(zeta, 1.413057E+0, 4.578814E-1, -6.850581E-2, 
		                  -5.588658E-2);
              a = pow(a29, smc.a(32, z)); 
     }
             break;
     case 30: a = ap(zeta, 3.910862E+1, 5.196646E+1, 2.264970E+1, 2.873680E+0);
             break;
     case 31: a = ap(zeta, 4.597479E+0, -2.855179E-1, 2.709724E-1);
             break;
     case 32: a = ap(zeta, 6.682518E+0, 2.827718E-1, -7.294429E-2);
             break;
     case 33: a = max(0.6355 - 0.4192*zeta, 
		      max(1.25, 
			  min(1.4, 1.5135 + 0.3769*zeta)));
             break;
     case 34: a = ap(zeta, 1.910302E-1, 1.158624E-1, 3.348990E-2, 2.599706E-3);
             break;
     case 35: a = ap(zeta, 3.931056E-1, 7.277637E-2, -1.366593E-1, -4.508946E-2);
             break;
     case 36: a = ap(zeta, 3.267776E-1, 1.204424E-1, 9.988332E-2, 2.455361E-2);
             break;
     case 37: a = ap(zeta, 5.990212E-1, 5.570264E-2, 6.207626E-2, 1.777283E-2);
             break;
     case 50: {double a50 = ap(zeta, 2.4000E-1, 1.8000E-1, 5.9500E-1);
              a = min(a50, 0.306 + 0.053*zeta);
     }
             break;
     case 51: {double a51 = ap(zeta, 3.3000E-1, 1.3200E-1, 2.1800E-1);
              a = min(a51, 0.3625 + 0.062*zeta);
     }
             break;
     case 52: {double a52 = ap(zeta, 1.1064E+0, 4.1500E-1, 1.8000E-1);
       			a = max(a52, 0.9);
                if(z>0.01)
					a = min(a, 1.0);
	 }
             break;
     case 53: {double a53 = ap(zeta, 1.1900E+0, 3.7700E-1, 1.7600E-1);
				a = max(a53, 1.0);
	            if(z>0.01)
					a = min(a, 1.1);
     }
             break;
     case 58: a = ap(zeta, 4.907546E-1, -1.683928E-1, -3.108742E-1, -7.202918E-2);
             break;
     case 59: a = ap(zeta, 4.537070E+0, -4.465455E+0, -1.612690E+0, -1.623246E+0);
             break;
     case 60: a = ap(zeta, 1.796220E+0, 2.814020E-1, 1.423325E+0, 3.421036E-1);
             break;
     case 61: a = ap(zeta, 2.256216E+0, 3.773400E-1, 1.537867E+0, 4.396373E-1);
             break;
     case 62: {double a62 = ap(zeta, 8.4300E-2, -4.7500E-2, -3.5200E-2);
              a = max(0.065, a62);
     }
             break;
     case 63: {a = ap(zeta, 7.3600E-2, 7.4900E-2, 4.4260E-2);
              if(z<0.004)
                  a = min(0.055, a);
    }
             break;
     case 64: {double a64 = ap(zeta, 1.3600E-1, 3.5200E-2);
              a = max(0.091, min(0.121, a64));
              if (smc.a(68, z) >= smc.a(66, z)) {
                  a =  smc.a(58, z)*pow(smc.a(66, z), smc.a(60, z))
                  / (smc.a(59, z) + pow(smc.a(66, z), smc.a(61, z)));
                                   
              }
     }
             break;
     case 65: {a = ap(zeta, 1.564231E-3, 1.653042E-3, -4.439786E-3, 
		     -4.951011E-3, -1.216530E-03);
     }
             break;
     case 66: {double a66 = ap(zeta, 1.4770E+0, 2.9600E-1);
             a = max(0.8, 
		     min(0.8 - 2.0*zeta, 
			 max(a66, 
			     min(1.6, -0.308 - 1.046*zeta))));
     }
             break;
     case 67: a = ap(zeta, 5.210157E+0, -4.143695E+0, -2.120870E+0);
             break;
     case 68: {double a68 = ap(zeta, 1.1160E+0, 1.6600E-1);
             a68 = max(0.9, min(a68, 1.0));
             a = min(a68, smc.a(66, z));
     }
             break;

     case 69: a = ap(zeta, 1.071489E+0, -1.164852E-1, -8.623831E-2, 
		     -1.582349E-2);
             break;
     case 70: a = ap(zeta, 7.108492E-1, 7.935927E-1, 3.926983E-1, 3.622146E-2);
             break;
     case 71: a = ap(zeta, 3.478514E+0, -2.585474E-2, -1.512955E-2, 
		     -2.833691E-3);
             break;
     case 72: {double a72 = ap(zeta, 9.132108E-1, -1.653695E-1, 0.0,  3.636784E-2);
             if(z>0.01)
	       a = max(a72, 0.95);
	     else
	       a = a72;
     }
             break;
     case 73: a = ap(zeta, 3.969331E-3, 4.539076E-3, 1.720906E-3, 
		     1.897857E-4);
             break;
     case 74: {double a74 = ap(zeta, 1.600E+0, 7.640E-1, 3.322E-1);
             a = max(1.4, min(a74, 1.6));
     }
             break;
     case 75: {double a75 = ap(zeta, 8.109E-1, -6.282E-1);
             a = max(max(1.0, 
			 min(a75, 1.27)), 0.6355 - 0.4192*zeta);
     }
             break;
     case 76: {double a76 = ap(zeta, 1.192334E-2, 1.083057E-2, 
			     1.230969E+0, 1.551656E+0);
             a = max(a76, 
		     -0.1015564 - 0.2161264*zeta - 0.05182516*pow(zeta, 2));
     }
             break;
     case 77: {double a77 = ap(zeta, -1.668868E-1, 5.818123E-1, -1.105027E+1, -1.668070E+1);
             a = max(-0.3868776 - 0.5457078*zeta - 0.1463472*pow(zeta, 2),
		     min(0.0, a77));
     }
             break;
     case 78: {double a78 = ap(zeta, 7.615495E-1, 1.068243E-1, 
			          -2.011333E-1, -9.371415E-2);
             a = max(0.0, min(a78, 7.454 + 9.046*zeta));
     }
             break;
     case 79: {double a79 = ap(zeta, 9.409838E+0, 1.522928E+0);
             a = min(a79, max(2.0, - 13.3 - 18.6*zeta));
     }
             break;
     case 80: {double a80 = ap(zeta, -2.7110E-1, -5.7560E-1, -8.3800E-2);
             a = max(0.0585542, a80);
     }
             break;
     case 81: {double a81 = ap(zeta, 2.4930E+0, 1.1475E+0);
             a = min(1.5, max(0.4, a81));
     }
             break;

     default:  cerr << "Not a valid index for "
		    << "constant.stellar_model_constant::a(int i= " 
		    << index << ")" << endl;
  }

  return a;
}


double bp(double zeta, double a, double b, double c, double d, double e) {

  // According to HPT2000 ap == bp (see Appendix A)
  return ap(zeta, a, b, c, d, e);
}

double stellar_model_constants::b(int index, double z) {

  double zeta = log10(z/SOLAR_METALICITY);

  double b = 0;
  switch(index) {
     case 1: {double b1 = bp(zeta, 3.9700E-1, 2.8826E-1, 5.2930E-1);
              b = min(b1, 0.54);}
             break;
     case 2: {double s = log10(z);
             double b2 = pow(10, -4.6739 - 0.9394*s);
	     b = min(max(b2, -0.04167+55.67*z), 0.4771 -9329.21*pow(z, 2.94));
     }
	     break;
     case 3: {double s = log10(z);
             double b3 = max(-0.1451, -2.2794 - s*(1.5175 + s*0.254));
             b = pow(10, b3);
	     if(z>0.004) 
	       b = max(b, 0.7307 + 14265.1*pow(z, 3.395));
             }
	     break;
     case 4: {double b4 = bp(zeta, 9.960283E-1, 8.164393E-1, 
			  2.383830E+0, 2.223436E+0, 8.638115E-1);
             b = b4 + 0.1231572*pow(zeta, 5);
     }
             break;
     case 5: b = bp(zeta, 2.561062E-1, 7.072646E-2, -5.444596E-2, -5.798167E-2,
		    -1.349129E-2);
             break;
     case 6: {double b6 = bp(zeta, 1.157338E+0, 1.467883E+0, 4.299661E+0, 
			  3.130500E+0, 6.992080E-1);
             b = b6 + 0.01640687*pow(zeta, 5);
     }
             break;
     case 7: b = bp(zeta, 4.022765E-1, 3.050010E-1, 9.962137E-1, 7.914079E-1,
		    1.728098E-1);
             break;

     case 8: cerr << "Unknown index for b8"<<endl;
          break;
          
     case 9: b = bp(zeta, 2.751631E+3, 3.557098E+2);
             break;
     case 10: b = bp(zeta, -3.820831E-2, 5.872664E-2);
              break;
     case 11: {double b11 = bp(zeta, 1.071738E+2, -8.970339E+1, -3.949739E+1);
              b = pow(b11, 2);
     }
             break;
     case 12: b = bp(zeta, 7.348793E+2, -1.531020E+2, -3.793700E+1);
             break;
     case 13: {double b13 = bp(zeta, 9.219293E+0, -2.005865E+0, -5.561309E-1);
              b = pow(b13, 2);
     }
             break;
     case 14: {double b14 = bp(zeta, 2.917412E+0, 1.575290E+0, 5.751814E-1);
              b = pow(b14, smc.b(15,z));
     }
             break;
     case 15: b = bp(zeta, 3.629118E+0, -9.112722E-1, 1.042291E+0);
             break;
     case 16: {double b16 = bp(zeta, 4.916389E+0, 2.862149E+0, 7.844850E-1);
              b = pow(b16, smc.b(15,z));
     }
	      break;
     case 17: b = 1;
              if (zeta>-1) 
		b = 1 - 0.3880523*pow(zeta+1, 2.862149);
	     break;
     case 18: b = bp(zeta, 5.496045E+1, -1.289968E+1, 6.385758E+0);
             break;
     case 19: b = bp(zeta, 1.832694E+0, -5.766608E-2, 5.696128E-2);
             break;
     case 20: b = bp(zeta, 1.211104E+2);
             break;
     case 21: b = bp(zeta, 2.214088E+2, 2.187113E+2, 1.170177E+1, -2.635340E+1);
             break;
     case 22: b = bp(zeta, 2.063983E+0, 7.363827E-1, 2.654323E-1, -6.140719E-2);
             break;
     case 23: b = bp(zeta, 2.003160E+0, 9.388871E-1, 9.656450E-1, 2.362266E-1);
             break;
     case 24: {double b24 = bp(zeta, 1.609901E+1, 7.391573E+0, 
			    2.277010E+1, 8.334227E+0);
              b = pow(b24, smc.b(28, z));
     }
             break;
     case 25: b = bp(zeta, 1.747500E-1, 6.271202E-2, -2.324229E-2, -1.844559E-2);
             break;
     case 26: b = 5 - 0.09138012*pow(z, -0.3671407);
             break;
     case 27: {double b27 = bp(zeta, 2.752869E+0, 2.729201E-2, 
			    4.996927E-1, 2.496551E-1);
              b = pow(b27, 2*smc.b(28, z));
     }
             break;
     case 28: b = bp(zeta, 3.518506E+0, 1.112440E+0, -4.556216E-1, -2.179426E-1);
             break;


     case 29: b = bp(zeta, 1.626062E+2, -1.168838E+1, -5.498343E+0);
             break;
     case 30: b = bp(zeta, 3.336833E-1, -1.458043E-1, -2.011751E-2);
             break;
     case 31: {double b31 = bp(zeta, 7.425137E+1, 1.790236E+1, 
			    3.033910E+1, 1.018259E+1);
              b = pow(b31, smc.b(33, z));
     }
             break;
     case 32: b = bp(zeta, 9.268325E+2, -9.739859E+1, -7.702152E+1, 
		     -3.158268E+1);
             break;
     case 33: b = bp(zeta, 2.474401E+0, 3.892972E-1);
             break;
     case 34: {double b34 = bp(zeta, 1.127018E+1, 1.622158E+0, 
			    -1.443664E+0, -9.474699E-1);
              b = pow(b34, smc.b(33, z));
     }
             break;
     case 35: cerr << "Unknown index for b35"<<endl;
	     break;
     case 36: {double b36 = bp(zeta, 1.445216E-1, -6.180219E-2, 
			    3.093878E-2, 1.567090E-2);
              b = pow(b36, 4);
     }
             break;
     case 37: {double b37 = bp(zeta, 1.304129E+0, 1.395919E-1, 
			    4.142455E-3, -9.732503E-3);
              b = 4*b37;
     }
             break;
     case 38: {double b38 = bp(zeta, 5.114149E-1, -1.160850E-2);
              b = pow(b38, 4);
     }
             break;
     case 39: b = bp(zeta, 1.314955E+2, 2.009258E+1, -5.143082E-1, -1.379140E+0);
             break;
     case 40: {double b40 = bp(zeta, 1.823973E+1, -3.074559E+0, -4.307878E+0);
             b = max(b40, 1.);
     }
             break;
     case 41: {double b41 = bp(zeta, 2.327037E+0, 2.403445E+0, 
			    1.208407E+0, 2.087263E-1);
             b = pow(b41, smc.b(42, z));
     }
             break;
     case 42: b = bp(zeta, 1.997378E+0, -8.126205E-1);
             break;
     case 43: b = bp(zeta, 1.079113E-1, 1.762409E-2, 1.096601E-2, 3.058818E-3);
             break;
     case 44: {double b44 = bp(zeta, 2.327409E+0, 6.901582E-1, 
			    -2.158431E-1, -1.084117E-1);
              b = pow(b44, 5);
     }
             break;
     case 45: {double p = zeta+1;
              if(p<=0)
		b = 1;
	      else
		b = 1 - (2.47162*p - 5.401682*pow(p, 2) + 3.247361*pow(p, 3));
     }
	      break;
     case 46: {b = bp(zeta, 2.214315E+0, -1.975747E+0);
     }
             break;
     case 47: {double p = zeta+1;
              b = 1.127733*p + 0.2344416*pow(p, 2) - 0.3793726*pow(p, 3);
     }
	      break;
     case 48: b = bp(zeta, 5.072525E+0, 1.146189E+1, 6.961724E+0, 1.316965E+0);
             break;
     case 49: b = bp(zeta, 5.139740E+0);
             break;
     case 51: {double b51 = bp(zeta, 1.125124E+0, 1.306486E+0, 
			    3.622359E+0, 2.601976E+0, 3.031270E-1);
             b = b51 - 0.1343798*pow(zeta, 5);
     }
             break;
     case 52: b = bp(zeta, 3.349489E-1, 4.531269E-3, 1.131793E-1, 2.300156E-1,
		     7.632745E-2);
             break;
     case 53: {double b53 = bp(zeta, 1.467794E+0, 2.798142E+0, 9.455580E+0, 
			    8.963904E+0, 3.339719E+0);
             b = b53 + 0.4426929*pow(zeta, 5);
     }
             break;
     case 54: b = bp(zeta, 4.658512E-1, 2.597451E-1, 9.048179E-1, 7.394505E-1,
		     1.607092E-1);
             break;
     case 55: {double b55 = bp(zeta, 1.0422E+0, 1.3156E-1, 4.5000E-2);
             b = min(b55, 0.99164 - 743.123*pow(z, 2.83));
     }
             break;
     case 56: {double b56 = bp(zeta, 1.110866E+0, 9.623856E-1, 2.735487E+0, 
			    2.445602E+0, 8.826352E-1);
              b = b56 + 0.1140142*pow(zeta, 5);
     }
             break;
     case 57: {double b57 = bp(zeta, -1.584333E-1, -1.728865E-1, -4.461431E-1, 
		     -3.925259E-1, -1.276203E-1);
              b = b57 - 0.01308728*pow(zeta, 5);
     }
             break;
     default:  cerr << "Not a valid index for "
		    << "constant.stellar_model_constant::b(int i= " 
		    << index << ")" << endl;
     }

return b;
}




double cp(double zeta, double a, double b, double c, double d, double e) {
  // According to Tout96 cp == ap (== bp) of HPT 2000 
  return ap(zeta, a, b, c, d, e);
}

double stellar_model_constants::c(int index, double z) {

  double zeta = log10(z/SOLAR_METALICITY);

  double c = 0;
  switch(index) {
     case 8: c = cp(zeta, 1.715359E+00, 6.2246212E-01, -9.2557761E-01, 
			-1.16996966E+00, -3.0631491E-01); 
             break;
     case 9: c = cp(zeta, 6.597788E+00, -4.2450044E-01,-1.213339427E+01,
			-1.073509484E+01, -2.51487077E+00);
             break;
     case 10: c = cp(zeta,1.008855000E+01, -7.11727086E+00,-3.167119479E+01, 
			-2.424848322E+01,-5.33608972E+00 );
             break;
     case 11: c = cp(zeta, 1.012495E+00, 3.2699690E-01, -9.23418E-03, 
			-3.876858E-02, -4.12750E-03);
             break;
     case 12: c = cp(zeta,7.490166E-02, 2.410413E-02, 7.233664E-02, 
			3.040467E-02, 1.97741E-03 );
 	     break;
     case 13: c = cp(zeta,1.077422E-02 );
             break;
     case 14: c = cp(zeta, 3.082234E+00, 9.447205E-01, -2.15200882E+00, 
			-2.49219496E+00, -6.3848738E-01);
             break;
     case 15: c = cp(zeta,1.784778E+01, -7.4534569E+00,-4.896066856E+01,
			-4.005386135E+01, -9.09331816E+00 ); 
             break;
     case 16: c = cp(zeta,2.2582E-04, -1.86899E-03, 3.88783E-03, 
			1.42402E-03,-7.671E-05 ); 
             break;
     default:  cerr << "Not a valid index for "
		    << "constant.stellar_model_constant::c(int i= " 
		    << index << ")" << endl;



  }
  return c;
}     





//              Standard lineair interpolation routine.
double lineair_interpolation(const double x,
			   const double x1, const double x2,
			   const double y1, const double y2) {

        double a = (y2-y1)/(x2-x1);
	double b = y1 - a*x1;

	double y = a*x + b;
	return y;
}
