//double base_main_sequence_luminosity(const double mass, const double z);
double base_main_sequence_radius(const double mass, const double z);
//double terminal_main_sequence_luminosity(const double mass, const double z);  
//double terminal_main_sequence_radius(const double mass, const double z);  
//double main_sequence_luminosity(const double time, const double mass, const double z);  
double main_sequence_radius(const double time, const double mass, const double z);  
//double zams_luminosity_correction(const double time, const double mass, const double z);  
double zams_radius_correction(const double time, const double mass, const double z);  

double base_giant_branch_time(const double mass, const double z);
double main_sequence_time(const double mass, const double z);
double main_sequence_hook_time(const double mass, const double z);
double main_sequence_hook_mass(const double z);
double stars_without_main_sequence_hook(const double mass, const double z); 
double stars_with_main_sequence_hook(const double mass, const double z); 

//double alpha_l_coefficient(const double mass, const double z);  
//double beta_l_coefficient(const double mass, const double z);  
double alpha_r_coefficient(const double mass, const double z);  
double beta_r_coefficient(const double mass, const double z);   
double gamma_r_coefficient(const double mass, const double z);   


static
class stellar_model_constants {  // Easy to have a name for compiling.

  public:
  double a(int, double);
  double b(int, double);
  double c(int, double);
} smc;

double ap(double zeta, double a, double b=0, double c=0, double d=0, double e=0);
double bp(double zeta, double a, double b=0, double c=0, double d=0, double e=0);
double cp(double zeta, double a, double b=0, double c=0, double d=0, double e=0);


//              Interpolation function
       double lineair_interpolation(const double, const double, const double,
				  const double, const double);
