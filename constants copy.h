// Constants 
#include  <cmath>


#define SOLAR_RADIUS 6.96e+10
#define SOLAR_MASS 1.989e+33
#define SOLAR_LUMINOSITY 3.862e+33
#define CHANDRASEKHAR_MASS 1.44
#define G 6.67e-8
#define C 2.9979e+10
#define YEAR 3.156e+7
#define M_year 3.156e+13
#define ONE_THIRD 0.33333333333333333333
#define TWO_THIRD 0.66666666666666666667
#define PI 3.14159265358979323846
#define kpc 3.0857e+21
#define rad_to_deg 57.295779513
#define SOLAR_METALICITY 0.02

inline double	max(double x, double y)		{return (x > y) ? x : y;}
inline int	max(int x, int y)            	{return (x > y) ? x : y;}
inline long int	max(long int x, long int y)  	{return (x > y) ? x : y;}

inline double	min(double x, double y)       	{return (x < y) ? x : y;}
inline int	min(int x, int y)            	{return (x < y) ? x : y;}
inline long int	min(long int x, long int y)  	{return (x < y) ? x : y;}

//inline long int  abs(long int x)            	{return (x < 0) ? -x : x;}
//inline double    abs(double x)                	{return (x < 0) ? -x : x;}

//inline long int  abs(long int x)            	{return labs(x);}
//inline double    abs(double x)                	{return fabs(x);}




#define PR(x)  cerr << #x << " = " << x << " "
#define PRC(x) cerr << #x << " = " << x << ",  "
#define PRL(x) cerr << #x << " = " << x << endl
