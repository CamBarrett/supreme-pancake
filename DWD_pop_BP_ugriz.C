#include "util.h"
#include "constants.h"
#include "my_main_sequence.h"
#include "mb.h"
#include <cstdlib>
#include <vector>
#include <iostream>
#include <cstdio>
#include <math.h>

#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
using namespace std;

double m_wd_array2[] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2};
static vector<double> m_wd_array(m_wd_array2, m_wd_array2 + sizeof(m_wd_array2)/sizeof(double));
static vector< vector<double> > table_Mwd_time;
static vector< vector<double> > table_Mwd_temp;
static vector< vector<double> > table_Mwd_u;
static vector< vector<double> > table_Mwd_g;
static vector< vector<double> > table_Mwd_r;
static vector< vector<double> > table_Mwd_i;
static vector< vector<double> > table_Mwd_z;
static vector< vector<double> > table_Mwd_Mbol;
static vector< vector<double> > table_G_;
static vector< vector<double> > table_G_BP;
static vector< vector<double> > table_G_RP;
static vector< vector<double> > table_Mms;

void read_global_table(string filename, vector< vector<double> > &table){

    ifstream ifs(filename.c_str());
    if (!ifs) cerr << "error: couldn't open file "<< filename <<endl;
    
    while (!ifs.eof()){
        string line;
        getline(ifs, line);
        //   cout<<line<<endl;
        istringstream iss(line);
        vector<double> row;
        copy(istream_iterator<double>(iss), istream_iterator<double>(), back_inserter(row));
        table.push_back(row);
    };    
    ifs.close();
}


//i_wd table of lower limit mass
//i_wd+1 table of upper limit mass
//j_wd lower limit mass, lower limit time
//j_wd+1 lower limit mass, upper limit time
//k_wd upper limit mass, lower limit time
//k_wd+1 upper limit mass, upper limit time
void find_nearest_neighbours_3d(int& ii_1, int& ii_2, int& jj_1, int& jj_2, int& kk_1, int& kk_2, 
                                    const double mass, const double time, 
                                    const vector<double>& m_array, const vector< vector<double> > &table_time){
     
//    cerr<< "find_nearest_neighbours_3d" <<endl;                                   
    int i_1, j_1, k_1; 
    int i_2, j_2, k_2; 
    if (ii_1 < 0){
        // first iteration  
        
        int length_mass = m_array.size();
        if (mass <= m_array[0]){
             i_1 = 0;       
             i_2 = 0;       
        }
        else if (mass >= m_array[length_mass-1]){
             i_1 = length_mass-1;
             i_2 = length_mass-1;             
        }
        else{        
            for (i_2 = 0; i_2 < length_mass; ++i_2){
                if (mass < m_array[i_2]) break;    
            }
            i_1 = i_2-1;
        }

        int length = table_time[i_1].size();
        if (time <= table_time[i_1][0]){
            j_1 = 0;
            j_2 = 0;
        }
        else if (time >= table_time[i_1][length-1]){
            j_1 = length-1;
            j_2 = length-1;
        }
        else{
            for (j_2 = 0; j_2 < length; ++j_2){
                if (time < table_time[i_1][j_2]) break;                 
            }        
            j_1 = j_2-1;
        }
        
        int length2 = table_time[i_2].size(); 
        if (time <= table_time[i_2][0]){
             k_1 = 0;
             k_2 = 0;
        }
        else if (time >= table_time[i_2][length2-1]){ 
            k_1 = length2 -1;
            k_2 = length2 -1;
        }
        else{
            for (k_2 = 0; k_2 < length2; ++k_2){
                if (time < table_time[i_2][k_2]) break;                 
            }     
            k_1 = k_2-1;
        } 
    }
    else{// after first iteration
        i_1 = ii_1;     
        i_2 = ii_2;     
        int length = table_time[i_1].size();        
        if (time <= table_time[i_1][0]){
             j_1 = 0;
             j_2 = 0;
        }
        else if(time >= table_time[i_1][length-1]){
             j_1 = length-1;
             j_2 = length-1;
        }
        else{
            for (j_2 = jj_1; j_2 < length; ++j_2){
                if (time < table_time[i_1][j_2]) break;                 
            }     
            j_1 = j_2-1;
        }
        
        int length2 = table_time[i_2].size();                
        if (time <= table_time[i_2][0]){
             k_1 = 0;
             k_2 = 0;             
        }
        else if(time >= table_time[i_2][length2-1]){
             k_1 = length2-1;
             k_2 = length2-1;
        }
        else{
            for (k_2 = kk_1; k_2 < length2; ++k_2){
                if (time < table_time[i_2][k_2]) break;                 
            }     
            k_1 = k_2-1;           
        }
    }
    
    ii_1 = i_1;
    ii_2 = i_2;
    jj_1 = j_1;
    jj_2 = j_2;
    kk_1 = k_1;
    kk_2 = k_2;
}

void find_nearest_neighbours_2d(int& ii_1, int& ii_2, const double x, const vector<double>& x_array){
      
    int length = x_array.size();
    if (x <= x_array[0]){
         ii_1 = 0;
         ii_2 = 0;
    }
    else if (x >= x_array[length-1]){
         ii_1 = length -1;
         ii_2 = length -1;
    }
    else{
        for (ii_2 = 0; ii_2 < length; ++ii_2){
            if (x < x_array[ii_2]) break;
        }
        ii_1 = ii_2-1;
    }
}

double interpol_irregular_grid(const int ii_1, const int ii_2, const int jj_1, 
                                const int jj_2, const int kk_1, const int kk_2, 
                                const double mass, const double time, 
                                const vector<double>& m_array, vector< vector<double> > &table_time, 
                                vector< vector<double> > &table_interpol){

//    cerr<< "interpol_irregular_grid" <<endl;                                   
    double grid_1, grid_2;  
    double grid_11 = table_interpol[ii_1][jj_1];
    double time_11 = table_time[ii_1][jj_1];  
    if (jj_1 != jj_2){     
        double grid_12 = table_interpol[ii_1][jj_2];  
        double time_12 = table_time[ii_1][jj_2];   
        grid_1 = grid_11 + (time - time_11)/(time_12 - time_11) * (grid_12 - grid_11);
    }
    else grid_1 = grid_11;
    
    if (ii_1 == ii_2){
        return grid_1;
    }

    double grid_21 = table_interpol[ii_2][kk_1];  
    double time_21 = table_time[ii_2][kk_1];   
    if (kk_1 != kk_2){
        double grid_22 = table_interpol[ii_2][kk_2];  
        double time_22 = table_time[ii_2][kk_2];  
        grid_2 = grid_21 + (time - time_21)/(time_22 - time_21) * (grid_22 - grid_21);   
    }
    else grid_2 = grid_21;
    
    double grid_interpol = grid_1 + (mass - m_array[ii_1])/(m_array[ii_2] - m_array[ii_1]) * (grid_2 - grid_1);   
    return grid_interpol; 
    
}

double interpol_2d(const double x, const double x1, const double x2,
               const double y1, const double y2){
    
    if (x1 == x2){
        return y1;
    }

    double a = (y2-y1)/(x2-x1);
    double b = y1 - a*x1;
    double y = a*x + b;
    return y;
        
}


//Rappaport, Verbunt, Joss 1983
// P**y = Pi**y + y*xi*(t-ti)
// t = ti + (Pmax**y - Pi**y)/y/xi
double P_next(double Mstar_1, double Mstar_2, double Pi, double ti, double t, double z) {
    // ti = t_form, t = t
    double P = 0;
    Pi *= 8.64e4; // in sec
    double yr_to_sec = 3.1536E7;
    
    double xi_gw = -3.68e-6 * Mstar_1 * Mstar_2 / pow(Mstar_1 + Mstar_2, ONE_THIRD);
    double y_gw = 8./3.; 
    P = pow(pow(Pi,y_gw) + (t-ti)*yr_to_sec*y_gw*xi_gw, 1/y_gw);
    
    P /= 8.64e4;
    return P;

}

double start_time(double Mstar_1, double Mstar_2, double Pmax, double Pi, double ti, double z) {

    double t = 0;
    Pi *= 8.64e4;
    Pmax *= 8.64e4;
    double yr_to_sec = 3.1536E7;
    
    double xi_gw = -3.68e-6 * Mstar_1 * Mstar_2 / pow(Mstar_1 + Mstar_2, ONE_THIRD);
    double y_gw = 8./3.; 
    t = ti + (pow(Pmax,y_gw)-pow(Pi,y_gw)) /yr_to_sec/y_gw/xi_gw; 
    
    return t;
}


double period_RLOF(double Mstar_1, double Mstar_2, double Rstar_1, double Rstar_2) {

    double q1_3   = pow(Mstar_2/Mstar_1, ONE_THIRD);
    double q2_3   = q1_3*q1_3;  
    double ap_RLOF = Rstar_1*(0.6*(1./q2_3) + log(1 + (1./q1_3)))/(0.49*(1./q2_3));
    double as_RLOF = Rstar_2*(0.6*q2_3 + log(1 + q1_3))/(0.49*q2_3);
    double a_RLOF  = max(ap_RLOF,as_RLOF);
    double P = 0.116 * sqrt(pow(a_RLOF, 3.)/(Mstar_1 + Mstar_2)); // in days
        
    return P;

}


int main(int argc, char* argv[]) {    
    //assuming file contains 250000 binaries, and a binary fraction of 50%, otherwise change coeff    
    //assuming a SeBa run with: ./SeBa -f 4 -m 0.95 -M 10 -> needed for birth_rate_coefficient in myutil.C 
    //char filename[] = "/Users/silviato/Development/process_data/data_run_SeBa/data_r105/SeBa_r105_ag_wdwd.data";
    char filename[] = "/Users/Cam/Desktop/Uni_work/Year_4/Project/Original/stable.txt";
//    double coeff= birth_rate_coefficient(false, 4, 0, 500000, 0.5);    
    
    // if M_wd<0.2, the colors of a 0.2Msun are taken, similarly to a WD mass more massive then 1.2
    double g_lim = 23; 
    double P_max = 100; // in days
    double z = SOLAR_METALICITY; 
    int set_resolution = 0; // desired resolution
        // 0: default
        // 1: 10*MW resolution, P_max = 100
        // 2: 50*MW resolution, all periods, distance < 200 pc
    
    double max_ms_mass=100;
    int set_ms_mass = 0;
    // if M_MS<0.088, the colors of a 0.088Msun are taken, similarly to a WD mass out of the range 0.2-1.2Msun.
    // for now only for stars with M < 3.8Msun, although be careful with M>2.9, I don't trust them

    for (int i = 1; i < argc; ++i) { 
        if (i + 1 != argc){
            if ( strcmp(argv[i], "-f") == 0) {
                strcpy(filename, argv[i + 1]);
                i++;
            } else if ( strcmp(argv[i], "-g") == 0) {
                g_lim = atof(argv[i + 1]);
                i++;
            } else if ( strcmp(argv[i], "-P") == 0) {
                P_max = atof(argv[i + 1]);
                i++;
            } else if (strcmp(argv[i], "-R") == 0) {
                set_resolution = atoi(argv[i + 1]);
                i++;
            } else if (strcmp(argv[i], "-M") == 0) {
                max_ms_mass = atof(argv[i + 1]);
                i++;
            } else if (strcmp(argv[i], "-m") == 0) {
                set_ms_mass = atoi(argv[i + 1]);
                i++;
            } else {
                cout << "Not enough or invalid arguments, please try again.\n";
                exit(0);
            }
        }
    }
    /*
    read_global_table("Bergeron_time.txt", table_Mwd_time);
    read_global_table("Bergeron_temp.txt", table_Mwd_temp);
    read_global_table("Bergeron_u.txt", table_Mwd_u);
    read_global_table("Bergeron_g.txt", table_Mwd_g);
    read_global_table("Bergeron_r.txt", table_Mwd_r);
    read_global_table("Bergeron_i.txt", table_Mwd_i);
    read_global_table("Bergeron_z.txt", table_Mwd_z);
    read_global_table("Bergeron_Mbol.txt", table_Mwd_Mbol);
    read_global_table("G_.txt", table_G_);
    */
    read_global_table("Tables/time.txt", table_Mwd_time);
    read_global_table("Tables/temp.txt", table_Mwd_temp);
    read_global_table("Tables/u.txt", table_Mwd_u);
    read_global_table("Tables/g.txt", table_Mwd_g);
    read_global_table("Tables/r.txt", table_Mwd_r);
    read_global_table("Tables/i.txt", table_Mwd_i);
    read_global_table("Tables/z.txt", table_Mwd_z);
    read_global_table("Tables/Mbol.txt", table_Mwd_Mbol);
    read_global_table("Tables/G_.txt", table_G_);
    read_global_table("Tables/G_BP.txt", table_G_BP);
    read_global_table("Tables/G_RP.txt", table_G_RP);
    
    if (set_ms_mass == 0){
        read_global_table("Tables/Kraus.txt", table_Mms);
        max_ms_mass = min(max_ms_mass, 3.8);
    } else if (set_ms_mass == 1){
        read_global_table("Tables/Kraus_rebassa_mansergas.txt", table_Mms);
        max_ms_mass = min(max_ms_mass, 0.472);
    } else {
        cout << "Conversion from spectral type to MS mass not defined, please try again.\n";
        exit(0);
    }
    
    double coeff;
    //cerr << 'test'<<endl;
    if (set_resolution == 0){
        coeff= birth_rate_coefficient(false, 5, 0, 1*50007, 0.4, 0.1); //coeff= birth_rate_coefficient(false, 5, 0, 250000, 0.5);
        // This is MW resolution so is ideal to keep this?
    } else if (set_resolution == 1) {//high resolution
        coeff= birth_rate_coefficient(false, 5, 0, 0.1*50007, 0.4, 0.1);
    } else if (set_resolution == 2) {//high resolution & d<0.2
        coeff= birth_rate_coefficient(false, 5, 0, 0.02*50007, 0.4, 0.1);
        P_max = 1e10;
    } else {
        cout << "Resolution not defined, please try again.\n";
        exit(0);        
    }
    
    //cerr<<coeff<<endl;
    
    cerr<<"# The g band magnitude limit is: "<<g_lim<<endl;
    cerr<<"# The maximum period is: "<<P_max<<"d"<<endl;
//    cerr<<"# The maximum donor mass is: "<<max_ms_mass<<"solar mass"<<endl;
//    cerr<<"# SpT to MS mass from: "<<set_ms_mass<<endl;
    cerr<<"# Data from: " << filename <<endl;
    cerr<<"# Resolution is: " << set_resolution << endl;
        
    double T_end = 1.35e10; //in yrs
    double inner_sep =0, ecc = 0, formation_time_star_1 =0, formation_time_star_2 =0, formation_time_star_3 =0, evolution_time =0;
    double l, b, d = 0.;    

    double star_1_mass = 0, star_2_mass = 0, star_3_mass = 0;
    double inclination =0, inner_ecc =0, outer_sep =0, outer_ecc =0;
    double star_1_id = 0, star_2_id = 0, star_3_id = 0;
    int bin_id = 0, id =0, line_id =0;
    int star_1_ii1 = 0, star_1_jj1 = 0, star_1_kk1 = 0;
    int star_1_ii2 = 0, star_1_jj2 = 0, star_1_kk2 = 0;
    int star_2_ii1 = 0, star_2_jj1 = 0, star_2_kk1 = 0;
    int star_2_ii2 = 0, star_2_jj2 = 0, star_2_kk2 = 0;
    int star_3_ii1 = 0, star_3_jj1 = 0, star_3_kk1 = 0;
    int star_3_ii2 = 0, star_3_jj2 = 0, star_3_kk2 = 0;

    int num = 0;   
    double dt;//to read in data
    
    ifstream is(filename, ios::in);
    if (!is) cerr << "error: couldn't open file "<< filename <<endl;
    
    do {
        
        is >> line_id >> dt >> inclination >> inner_sep >> inner_ecc >> dt >> dt >> outer_sep >> outer_ecc >> dt >> dt >> star_1_mass >> star_2_mass >> star_3_mass >> dt >> dt >> dt >> dt >> dt >> dt >> dt >> dt >> dt >> star_1_id >> star_2_id >> star_3_id >> dt >> dt >> dt >> evolution_time >> dt >> dt >> formation_time_star_1 >> formation_time_star_2;


        double inner_mass = star_1_mass + star_2_mass;
        //double formation_time_star_3 = evolution_time;
        
        
        
        star_1_ii1 = -1, star_1_jj1 = -1, star_1_kk1 = -1;
        star_1_ii2 = -1, star_1_jj2 = -1, star_1_kk2 = -1;
        star_2_ii1 = -1, star_2_jj1 = -1, star_2_kk1 = -1;
        star_2_ii2 = -1, star_2_jj2 = -1, star_2_kk2 = -1;
        star_3_ii1 = -1, star_3_jj1 = -1, star_3_kk1 = -1;
        star_3_ii2 = -1, star_3_jj2 = -1, star_3_kk2 = -1;
        
        if ( (((2<star_1_id && star_1_id<10) || star_1_id > 12) && ((2<star_2_id && star_2_id<10) || star_2_id > 12)) || (((2<star_1_id && star_1_id<10) || star_1_id > 12) && ((2<star_3_id && star_3_id<10) || star_3_id > 12)) || (((2<star_2_id && star_2_id<10) || star_2_id > 12) && ((2<star_3_id && star_3_id<10) || star_3_id > 12)) ) {
            cerr << "ERROR something wrong!!! Using Cam version of id cut" << endl;
            PRC(star_1_id);PRL(star_1_mass);
            PRC(star_2_id);PRL(star_2_mass);
            PRC(star_3_id);PRL(star_3_mass);
            double inner_sep = outer_sep;
            PRC(formation_time_star_1);PRC(formation_time_star_2);PRC(formation_time_star_3);PRL(inner_sep); //inner_sep?
            exit(1);
            
        }
        
        
        //cout << bin_id << endl;
        double P = 0.116 * sqrt(pow(inner_sep, 3.)/(star_1_mass + star_2_mass)); // in days
        double P_out = 0.116 * sqrt(pow(outer_sep, 3.)/(inner_mass + star_3_mass)); // in days
        //cout << "Starting Inner " << P << " Starting outer " << P_out << " " << endl;
        double Pi = P;
        double Po = P_out;
        
        double star_1_rad =0;
        double star_2_rad =0;
        double star_3_rad =0;
        double t_ms =0;
        double form_time_star_1 =0;
        double form_time_star_2 =0;
        double form_time_star_3 =0;
        double t_form =0;
        
        if (star_1_id < 2) {
            // Star 1 = MS
            star_1_rad = base_main_sequence_radius(star_1_mass, z);
            star_2_rad = wd_radius_PPE(star_2_mass);
            star_3_rad = wd_radius_PPE(star_3_mass);
            
            t_ms = main_sequence_time(star_1_mass, z) * 1e6; // time on ms stage
            
            form_time_star_1 = 0;
            form_time_star_2 = formation_time_star_1 * 1e6;
            form_time_star_3 = formation_time_star_2 * 1e6;
            
            t_form = formation_time_star_2 * 1e6;
        }
        
        if (star_2_id < 2) {
            // Star 2 = MS
            star_1_rad = wd_radius_PPE(star_1_mass);
            star_2_rad = base_main_sequence_radius(star_2_mass, z);
            star_3_rad = wd_radius_PPE(star_3_mass);
            
            t_ms = main_sequence_time(star_2_mass, z) * 1e6;
            
            form_time_star_1 = formation_time_star_1 * 1e6;
            form_time_star_2 = 0;
            form_time_star_3 = formation_time_star_2 * 1e6;
            
            t_form = formation_time_star_2 * 1e6;
        }
        
        if (star_3_id < 2) {
            // Star 3 = MS
            star_1_rad = wd_radius_PPE(star_1_mass);
            star_2_rad = wd_radius_PPE(star_2_mass);
            star_3_rad = base_main_sequence_radius(star_3_mass, z);
            
            t_ms = main_sequence_time(star_3_mass, z) * 1e6;
            
            form_time_star_1 = formation_time_star_1 * 1e6;
            form_time_star_2 = formation_time_star_2 * 1e6;
            form_time_star_3 = 0;
            
            t_form = formation_time_star_2 * 1e6;
        }
        
        if ((star_1_id==10 || star_1_id==11 || star_1_id==12) && (star_2_id==10 || star_2_id==11 || star_2_id==12) && (star_3_id==10 || star_3_id==11 || star_3_id==12)) {
            // If all the stars are WDs
            star_1_rad = wd_radius_PPE(star_1_mass);
            star_2_rad = wd_radius_PPE(star_2_mass);
            star_3_rad = wd_radius_PPE(star_3_mass);
            
            form_time_star_1 = formation_time_star_1 * 1e6;
            form_time_star_2 = formation_time_star_2 * 1e6;
            form_time_star_3 = evolution_time * 1e6;
            t_form = form_time_star_3;
        }
        
        
        double t_form_star_1 = form_time_star_1; // formation of first white dwarf --now in this case the inner binary??
        double t_form_star_2 = form_time_star_2;
        double t_form_star_3 = form_time_star_3;
        double t = t_form; // formation time of current state of the system?
        //cout << "The time for this system is " << t << " " << endl;
        if (Pi > P_max) t = start_time(star_1_mass, star_2_mass, P_max, Pi, t_form, z);
        
        
        do {
            
            double next_random_time;
            //cout << "NEXT RANDOM TIME " << next_random_time << " " << endl;
            double t_next = next_birth_interval_BP(t, 
            		     T_end, coeff, &next_random_time);
            num++;
            //cout << t_next;
            double t_gal = T_end - next_random_time;
            
            // position in the galaxy  
            random_position_BP(t_gal, &l, &b, &d);
            double A_V = get_A_V(l, b, d);
            double logd = log10(d);           

            if (set_resolution < 3){

            	if (next_random_time <= T_end && next_random_time <= t_ms) {
      
                    P = P_next(star_1_mass, star_2_mass, Pi, t_form, next_random_time, z);
                    P_out = P_next(inner_mass, star_3_mass, Po, t_form, next_random_time, z);
                    
                    if (P < P_max) { 
    
                        // RLOF
                        double P_RLOF = period_RLOF(star_1_mass, star_2_mass, star_1_rad, star_2_rad);
                        if (P <= P_RLOF) break; 

//                    double t_gal = T_end - next_random_time;
//                    
//                    // position in the galaxy  
//                    random_position_BP(t_gal, &l, &b, &d);
//                    double A_V = get_A_V(l, b, d);
//                    double logd = log10(d);           
                                            
                        //  colors
                        
                        double u_star_1, g_star_1, r_star_1, i_star_1, z_star_1, T_star_1, Mbol_star_1, G_1, G_BP_1, G_RP_1;
                        if (star_1_id < 2) {
                            //Must assign MS colours to the star
                            int ii_ms1, ii_ms2;
                            find_nearest_neighbours_2d(ii_ms1, ii_ms2, star_1_mass, table_Mms[10]);
                            T_star_1 = interpol_2d(star_1_mass, table_Mms[10][ii_ms1], table_Mms[10][ii_ms2],
                                                        table_Mms[9][ii_ms1], table_Mms[9][ii_ms2]);
                            u_star_1 = interpol_2d(star_1_mass, table_Mms[10][ii_ms1], table_Mms[10][ii_ms2],
                                                        table_Mms[0][ii_ms1], table_Mms[0][ii_ms2]);
                            g_star_1 = interpol_2d(star_1_mass, table_Mms[10][ii_ms1], table_Mms[10][ii_ms2],
                                                        table_Mms[1][ii_ms1], table_Mms[1][ii_ms2]);
                            r_star_1 = interpol_2d(star_1_mass, table_Mms[10][ii_ms1], table_Mms[10][ii_ms2],
                                                        table_Mms[2][ii_ms1], table_Mms[2][ii_ms2]);
                            i_star_1 = interpol_2d(star_1_mass, table_Mms[10][ii_ms1], table_Mms[10][ii_ms2],
                                                        table_Mms[3][ii_ms1], table_Mms[3][ii_ms2]);
                            z_star_1 = interpol_2d(star_1_mass, table_Mms[10][ii_ms1], table_Mms[10][ii_ms2],
                                                        table_Mms[4][ii_ms1], table_Mms[4][ii_ms2]);
                            double C2 = g_star_1 - i_star_1;
                            double C3 = g_star_1 - r_star_1;
                            // Colours transformed into Gaia colours according to Jordi et al. 2010
                            G_1 = g_star_1 + (-0.1005) + (-0.5358*C2) + (-0.1207*pow(C2, 2)) + (0.0082*pow(C2, 3)) + (-0.0272*C3) + (0.1270*pow(C3, 2)) + (-0.0205*pow(C3, 3)) + (-0.0176*C2*C3);
                            G_BP_1 = g_star_1 - (-0.0428) - (0.0158*C2) - (0.0122*pow(C2, 2)) - (0.0005*pow(C2, 3)) - (0.2382*C3) - (0.1855*pow(C3, 2)) - (-0.0096*pow(C3, 3)) - (-0.0493*C2*C3);
                            G_RP_1 = g_star_1 - (0.3731) - (1.3335*C2) - (-0.0242*pow(C2, 2)) - (-0.0047*pow(C2, 3)) - (-0.4086*C3) - (-0.2729*pow(C3, 2)) - (0.0230*pow(C3, 3)) - (0.2166*C2*C3);
                            Mbol_star_1 = 0;
                        }
                        else {
                            find_nearest_neighbours_3d(star_1_ii1, star_1_ii2, star_1_jj1, star_1_jj2, star_1_kk1, star_1_kk2,
                                    star_1_mass, next_random_time - t_form_star_1, m_wd_array, table_Mwd_time);
                            u_star_1 = interpol_irregular_grid(star_1_ii1, star_1_ii2, star_1_jj1, star_1_jj2, star_1_kk1, star_1_kk2,
                                    star_1_mass, next_random_time-t_form_star_1, m_wd_array, table_Mwd_time, table_Mwd_u);
                            g_star_1 = interpol_irregular_grid(star_1_ii1, star_1_ii2, star_1_jj1, star_1_jj2, star_1_kk1, star_1_kk2,
                                    star_1_mass, next_random_time-t_form_star_1, m_wd_array, table_Mwd_time, table_Mwd_g);
                            r_star_1 = interpol_irregular_grid(star_1_ii1, star_1_ii2, star_1_jj1, star_1_jj2, star_1_kk1, star_1_kk2,
                                    star_1_mass, next_random_time-t_form_star_1, m_wd_array, table_Mwd_time, table_Mwd_r);
                            i_star_1 = interpol_irregular_grid(star_1_ii1, star_1_ii2, star_1_jj1, star_1_jj2, star_1_kk1, star_1_kk2,
                                    star_1_mass, next_random_time-t_form_star_1, m_wd_array, table_Mwd_time, table_Mwd_i);
                            z_star_1 = interpol_irregular_grid(star_1_ii1, star_1_ii2, star_1_jj1, star_1_jj2, star_1_kk1, star_1_kk2,
                                    star_1_mass, next_random_time-t_form_star_1, m_wd_array, table_Mwd_time, table_Mwd_z);
                            T_star_1 = interpol_irregular_grid(star_1_ii1, star_1_ii2, star_1_jj1, star_1_jj2, star_1_kk1, star_1_kk2,
                                    star_1_mass, next_random_time-t_form_star_1, m_wd_array, table_Mwd_time, table_Mwd_temp);
                            Mbol_star_1 = interpol_irregular_grid(star_1_ii1, star_1_ii2, star_1_jj1, star_1_jj2, star_1_kk1, star_1_kk2,
                                   star_1_mass, next_random_time-t_form_star_1, m_wd_array, table_Mwd_time, table_Mwd_Mbol);
                            G_1 = interpol_irregular_grid(star_1_ii1, star_1_ii2, star_1_jj1, star_1_jj2, star_1_kk1, star_1_kk2,
                                    star_1_mass, next_random_time-t_form_star_1, m_wd_array, table_Mwd_time, table_G_);
                            G_BP_1 = interpol_irregular_grid(star_1_ii1, star_1_ii2, star_1_jj1, star_1_jj2, star_1_kk1, star_1_kk2,
                                    star_1_mass, next_random_time-t_form_star_1, m_wd_array, table_Mwd_time, table_G_BP);
                            G_RP_1 = interpol_irregular_grid(star_1_ii1, star_1_ii2, star_1_jj1, star_1_jj2, star_1_kk1, star_1_kk2,
                                    star_1_mass, next_random_time-t_form_star_1, m_wd_array, table_Mwd_time, table_G_RP);
                        }
                            //   d in kpc
                        u_star_1 += 10. + 5*logd + 1.58*A_V;
                        g_star_1 += 10. + 5*logd + 1.16*A_V;
                        r_star_1 += 10. + 5*logd + 0.84*A_V;
                        i_star_1 += 10. + 5*logd + 0.64*A_V;
                        z_star_1 += 10. + 5*logd + 0.45*A_V;
                        G_1 += 10. + 5*logd + 0.789*A_V;
                        G_BP_1 += 10. + 5*logd + 1.002*A_V;
                        G_RP_1 += 10. + 5*logd + 0.589*A_V;
                        
                        
                        double u_star_2, g_star_2, r_star_2, i_star_2, z_star_2, T_star_2, Mbol_star_2, G_2, G_BP_2, G_RP_2;
                        if (star_2_id<2) {
                            //Must assign MS colours to the star
                            int ii_ms1, ii_ms2;
                            find_nearest_neighbours_2d(ii_ms1, ii_ms2, star_2_mass, table_Mms[10]);
                            T_star_2 = interpol_2d(star_2_mass, table_Mms[10][ii_ms1], table_Mms[10][ii_ms2],
                                                        table_Mms[9][ii_ms1], table_Mms[9][ii_ms2]);
                            u_star_2 = interpol_2d(star_2_mass, table_Mms[10][ii_ms1], table_Mms[10][ii_ms2],
                                                        table_Mms[0][ii_ms1], table_Mms[0][ii_ms2]);
                            g_star_2 = interpol_2d(star_2_mass, table_Mms[10][ii_ms1], table_Mms[10][ii_ms2],
                                                        table_Mms[1][ii_ms1], table_Mms[1][ii_ms2]);
                            r_star_2 = interpol_2d(star_2_mass, table_Mms[10][ii_ms1], table_Mms[10][ii_ms2],
                                                        table_Mms[2][ii_ms1], table_Mms[2][ii_ms2]);
                            i_star_2 = interpol_2d(star_2_mass, table_Mms[10][ii_ms1], table_Mms[10][ii_ms2],
                                                        table_Mms[3][ii_ms1], table_Mms[3][ii_ms2]);
                            z_star_2 = interpol_2d(star_2_mass, table_Mms[10][ii_ms1], table_Mms[10][ii_ms2],
                                                        table_Mms[4][ii_ms1], table_Mms[4][ii_ms2]);
                            
                            double C2 = g_star_2 - i_star_2;
                            double C3 = g_star_2 - r_star_2;
                            
                            G_2 = g_star_2 + (-0.1005) + (-0.5358*C2) + (-0.1207*pow(C2, 2)) + (0.0082*pow(C2, 3)) + (-0.0272*C3) + (0.1270*pow(C3, 2)) + (-0.0205*pow(C3, 3)) + (-0.0176*C2*C3);
                            G_BP_2 = g_star_2 - (-0.0428) - (0.0158*C2) - (0.0122*pow(C2, 2)) - (0.0005*pow(C2, 3)) - (0.2382*C3) - (0.1855*pow(C3, 2)) - (-0.0096*pow(C3, 3)) - (-0.0493*C2*C3);
                            G_RP_2 = g_star_2 - (0.3731) - (1.3335*C2) - (-0.0242*pow(C2, 2)) - (-0.0047*pow(C2, 3)) - (-0.4086*C3) - (-0.2729*pow(C3, 2)) - (0.0230*pow(C3, 3)) - (0.2166*C2*C3);
                            Mbol_star_2 = 0;
                        }
                        else {
                            find_nearest_neighbours_3d(star_2_ii1, star_2_ii2, star_2_jj1, star_2_jj2, star_2_kk1, star_2_kk2,
                                    star_2_mass, next_random_time - t_form_star_2, m_wd_array, table_Mwd_time);
                            u_star_2 = interpol_irregular_grid(star_2_ii1, star_2_ii2, star_2_jj1, star_2_jj2, star_2_kk1, star_2_kk2,
                                    star_2_mass, next_random_time-t_form_star_2, m_wd_array, table_Mwd_time, table_Mwd_u);
                            g_star_2 = interpol_irregular_grid(star_2_ii1, star_2_ii2, star_2_jj1, star_2_jj2, star_2_kk1, star_2_kk2,
                                    star_2_mass, next_random_time-t_form_star_2, m_wd_array, table_Mwd_time, table_Mwd_g);
                            r_star_2 = interpol_irregular_grid(star_2_ii1, star_2_ii2, star_2_jj1, star_2_jj2, star_2_kk1, star_2_kk2,
                                    star_2_mass, next_random_time-t_form_star_2, m_wd_array, table_Mwd_time, table_Mwd_r);
                            i_star_2 = interpol_irregular_grid(star_2_ii1, star_2_ii2, star_2_jj1, star_2_jj2, star_2_kk1, star_2_kk2,
                                    star_2_mass, next_random_time-t_form_star_2, m_wd_array, table_Mwd_time, table_Mwd_i);
                            z_star_2 = interpol_irregular_grid(star_2_ii1, star_2_ii2, star_2_jj1, star_2_jj2, star_2_kk1, star_2_kk2,
                                    star_2_mass, next_random_time-t_form_star_2, m_wd_array, table_Mwd_time, table_Mwd_z);
                            T_star_2 = interpol_irregular_grid(star_2_ii1, star_2_ii2, star_2_jj1, star_2_jj2, star_2_kk1, star_2_kk2,
                                    star_2_mass, next_random_time-t_form_star_2, m_wd_array, table_Mwd_time, table_Mwd_temp);
                            Mbol_star_2 = interpol_irregular_grid(star_2_ii1, star_2_ii2, star_2_jj1, star_2_jj2, star_2_kk1, star_2_kk2,
                                    star_2_mass, next_random_time-t_form_star_2, m_wd_array, table_Mwd_time, table_Mwd_Mbol);
                            G_2 = interpol_irregular_grid(star_2_ii1, star_2_ii2, star_2_jj1, star_2_jj2, star_2_kk1, star_2_kk2,
                                    star_2_mass, next_random_time-t_form_star_2, m_wd_array, table_Mwd_time, table_G_);
                            G_BP_2 = interpol_irregular_grid(star_2_ii1, star_2_ii2, star_2_jj1, star_2_jj2, star_2_kk1, star_2_kk2,
                                    star_2_mass, next_random_time-t_form_star_2, m_wd_array, table_Mwd_time, table_G_BP);
                            G_RP_2 = interpol_irregular_grid(star_2_ii1, star_2_ii2, star_2_jj1, star_2_jj2, star_2_kk1, star_2_kk2,
                                    star_2_mass, next_random_time-t_form_star_2, m_wd_array, table_Mwd_time, table_G_RP);
                        }
                        //   d in kpc
                        u_star_2 += 10. + 5*logd + 1.58*A_V;
                        g_star_2 += 10. + 5*logd + 1.16*A_V;
                        r_star_2 += 10. + 5*logd + 0.84*A_V;
                        i_star_2 += 10. + 5*logd + 0.64*A_V;
                        z_star_2 += 10. + 5*logd + 0.45*A_V;
                        G_2 += 10. + 5*logd + 0.789*A_V;
                        G_BP_2 += 10. + 5*logd + 1.002*A_V;
                        G_RP_2 += 10. + 5*logd + 0.589*A_V;
                        
                        
                        double u_star_3, g_star_3, r_star_3, i_star_3, z_star_3, T_star_3, Mbol_star_3, G_3, G_BP_3, G_RP_3;
                        if (star_3_id<2) {
                            //Must assign MS colours to the star
                            int ii_ms1, ii_ms2;
                            find_nearest_neighbours_2d(ii_ms1, ii_ms2, star_3_mass, table_Mms[10]);
                            T_star_3 = interpol_2d(star_3_mass, table_Mms[10][ii_ms1], table_Mms[10][ii_ms2],
                                                        table_Mms[9][ii_ms1], table_Mms[9][ii_ms2]);
                            u_star_3 = interpol_2d(star_3_mass, table_Mms[10][ii_ms1], table_Mms[10][ii_ms2],
                                                        table_Mms[0][ii_ms1], table_Mms[0][ii_ms2]);
                            g_star_3 = interpol_2d(star_3_mass, table_Mms[10][ii_ms1], table_Mms[10][ii_ms2],
                                                        table_Mms[1][ii_ms1], table_Mms[1][ii_ms2]);
                            r_star_3 = interpol_2d(star_3_mass, table_Mms[10][ii_ms1], table_Mms[10][ii_ms2],
                                                        table_Mms[2][ii_ms1], table_Mms[2][ii_ms2]);
                            i_star_3 = interpol_2d(star_3_mass, table_Mms[10][ii_ms1], table_Mms[10][ii_ms2],
                                                        table_Mms[3][ii_ms1], table_Mms[3][ii_ms2]);
                            z_star_3 = interpol_2d(star_3_mass, table_Mms[10][ii_ms1], table_Mms[10][ii_ms2],
                                                        table_Mms[4][ii_ms1], table_Mms[4][ii_ms2]);
                            
                            double C2 = g_star_3 - i_star_3;
                            double C3 = g_star_3 - r_star_3;
                            
                            G_3 = g_star_3 + (-0.1005) + (-0.5358*C2) + (-0.1207*pow(C2, 2)) + (0.0082*pow(C2, 3)) + (-0.0272*C3) + (0.1270*pow(C3, 2)) + (-0.0205*pow(C3, 3)) + (-0.0176*C2*C3);
                            G_BP_3 = g_star_3 - (-0.0428) - (0.0158*C2) - (0.0122*pow(C2, 2)) - (0.0005*pow(C2, 3)) - (0.2382*C3) - (0.1855*pow(C3, 2)) - (-0.0096*pow(C3, 3)) - (-0.0493*C2*C3);
                            G_RP_3 = g_star_3 - (0.3731) - (1.3335*C2) - (-0.0242*pow(C2, 2)) - (-0.0047*pow(C2, 3)) - (-0.4086*C3) - (-0.2729*pow(C3, 2)) - (0.0230*pow(C3, 3)) - (0.2166*C2*C3);
                            Mbol_star_3 = 0;
                        }
                        else {
                            find_nearest_neighbours_3d(star_3_ii1, star_3_ii2, star_3_jj1, star_3_jj2, star_3_kk1, star_3_kk2,
                                    star_3_mass, next_random_time - t_form_star_3, m_wd_array, table_Mwd_time);
                            u_star_3 = interpol_irregular_grid(star_3_ii1, star_3_ii2, star_3_jj1, star_3_jj2, star_3_kk1, star_3_kk2,
                                    star_3_mass, next_random_time-t_form_star_3, m_wd_array, table_Mwd_time, table_Mwd_u);
                            g_star_3 = interpol_irregular_grid(star_3_ii1, star_3_ii2, star_3_jj1, star_3_jj2, star_3_kk1, star_3_kk2,
                                    star_3_mass, next_random_time-t_form_star_3, m_wd_array, table_Mwd_time, table_Mwd_g);
                            r_star_3 = interpol_irregular_grid(star_3_ii1, star_3_ii2, star_3_jj1, star_3_jj2, star_3_kk1, star_3_kk2,
                                    star_3_mass, next_random_time-t_form_star_3, m_wd_array, table_Mwd_time, table_Mwd_r);
                            i_star_3 = interpol_irregular_grid(star_3_ii1, star_3_ii2, star_3_jj1, star_3_jj2, star_3_kk1, star_3_kk2,
                                    star_3_mass, next_random_time-t_form_star_3, m_wd_array, table_Mwd_time, table_Mwd_i);
                            z_star_3 = interpol_irregular_grid(star_3_ii1, star_3_ii2, star_3_jj1, star_3_jj2, star_3_kk1, star_3_kk2,
                                    star_3_mass, next_random_time-t_form_star_3, m_wd_array, table_Mwd_time, table_Mwd_z);
                            T_star_3 = interpol_irregular_grid(star_3_ii1, star_3_ii2, star_3_jj1, star_3_jj2, star_3_kk1, star_3_kk2,
                                    star_3_mass, next_random_time-t_form_star_3, m_wd_array, table_Mwd_time, table_Mwd_temp);
                            Mbol_star_3 = interpol_irregular_grid(star_3_ii1, star_3_ii2, star_3_jj1, star_3_jj2, star_3_kk1, star_3_kk2,
                                    star_3_mass, next_random_time-t_form_star_3, m_wd_array, table_Mwd_time, table_Mwd_Mbol);
                            G_3 = interpol_irregular_grid(star_3_ii1, star_3_ii2, star_3_jj1, star_3_jj2, star_3_kk1, star_3_kk2,
                                    star_3_mass, next_random_time-t_form_star_3, m_wd_array, table_Mwd_time, table_G_);
                            G_BP_3 = interpol_irregular_grid(star_3_ii1, star_3_ii2, star_3_jj1, star_3_jj2, star_3_kk1, star_3_kk2,
                                    star_3_mass, next_random_time-t_form_star_3, m_wd_array, table_Mwd_time, table_G_BP);
                            G_RP_3 = interpol_irregular_grid(star_3_ii1, star_3_ii2, star_3_jj1, star_3_jj2, star_3_kk1, star_3_kk2,
                                    star_3_mass, next_random_time-t_form_star_3, m_wd_array, table_Mwd_time, table_G_RP);
                        }
                        
                        //   d in kpc
                        u_star_3 += 10. + 5*logd + 1.58*A_V;
                        g_star_3 += 10. + 5*logd + 1.16*A_V;
                        r_star_3 += 10. + 5*logd + 0.84*A_V;
                        i_star_3 += 10. + 5*logd + 0.64*A_V;
                        z_star_3 += 10. + 5*logd + 0.45*A_V;
                        G_3 += 10. + 5*logd + 0.789*A_V;
                        G_BP_3 += 10. + 5*logd + 1.002*A_V;
                        G_RP_3 += 10. + 5*logd + 0.589*A_V;
                        
//                        double r_in = outer_sep * ( star_3_mass / (inner_mass + star_3_mass) ) ;
//                        double r_out = outer_sep * ( inner_mass / (inner_mass + star_3_mass) ) ;
//
//                        if (r_in<0) {r_in *= -1;}
//                        if (r_out<0) {r_out *= -1;}
//
//                        double inn_1 = inner_sep * ( star_2_mass / inner_mass ) ;
//                        double inn_2 = inner_sep * ( star_3_mass / inner_mass ) ;
//
//                        if (inn_1<0) {inn_1 *= -1;}
//                        if (inn_2<0) {inn_2 *= -1;}
//
//                        double cm_1 = 0;
//                        double cm_2 = 0;
//
//                        if (inclination = 0) {
//                            cm_1 = r_in + inn_1; // solar radii
//                            cm_2 = r_in - inn_2;
//                        }
//                        else {
//                            cm_1 = sqrt( pow(inn_1, 2) + pow(r_in, 2)    ); // keeps it in solar radii
//                            cm_2 = ;
//                        }
//
//                        double d1 = sqrt( ( pow(d*3.1e19 , 2) ) + ( pow( (( r_in + inn_1 )*6.96e8) , 2) ) ) / (3.1e19) ; // in kpc to match d
//                        double d2 = sqrt( ( pow(d*3.1e19 , 2) ) + ( pow( (( r_in - inn_2 )*6.96e8) , 2) ) ) / (3.1e19) ;
//                        double d3 = sqrt( ( pow(d*3.1e19 , 2) ) + ( pow( r_out*6.96e8, 2) ) ) / (3.1e19) ;
//
//                        double common_motion_1 = atan((r_in+inn_1)*(6.96e8) / (d*3.1e19)) * 3600 * (180/M_PI) * 4 / (P_out*3.1536E7); //in arcsecs per year
//                        double common_motion_2 = atan((r_in-inn_1)*(6.96e8) / (d*3.1e19)) * 3600 * (180/M_PI) * 4 / (P_out*3.1536E7);
//                        double common_motion_3 = atan(r_out*(6.96e8) / (d*3.1e19)) * 3600 * (180/M_PI) * 4 / (P_out*3.1536E7);
                        //3.1536E7 seconds in a year, so divide by this to go from seconds to year
                        
                        //cout << d << " " << d1 << " " << d2 << " " << d3 << " " << endl;
                        
                        //cout << id << " " << star_1_mass << " " << star_2_mass << " " << star_1_id << " " << star_2_id << " " << star_3_id << endl
                        
                        
                        // depending on whether or not I need them later on I could just eliminate things from tge output, add a header for each column type, and
                        // then read in python using pandas
                        if (G_1 < g_lim || G_2 < g_lim || G_3 < g_lim) {
                            cout << bin_id << " " << P <<  " " << P_out << " "
                            << l << " " << b << " " << d << " "
                            << A_V << " "
                            << star_1_mass << " "  << star_1_rad << " " << T_star_1 << " "
                            << star_2_mass << " "  << star_2_rad << " " << T_star_2 << " "
                            << star_3_mass << " "  << star_3_rad << " " << T_star_3 << " "
                            << u_star_1 << " " << g_star_1 << " " << r_star_1 << " "
                            << i_star_1 << " " << z_star_1 << " "
                            << u_star_2 << " " << g_star_2 << " " << r_star_2 << " "
                            << i_star_2 << " " << z_star_2 << " "
                            << u_star_3 << " " << g_star_3 << " " << r_star_3 << " "
                            << i_star_3 << " " << z_star_3 << " "
                            << G_1 << " " << G_BP_1 << " " << G_RP_1 << " " << G_2 << " " << G_BP_2 << " " << G_RP_2 << " " << G_3 << " " << G_BP_3 << " " << G_RP_3 << " " << " " << star_1_id << " " << star_2_id << " " << star_3_id << " " << t << " " << line_id << " " << inner_sep << " " << outer_sep << " " << inner_ecc << " " << outer_ecc << " " << inclination << " " <<
                             endl;
                        } // cout
                    } // P < P_max
                } // next_random_time <= T_end
            } // resolution = 0,1 or (resolution = 2 & d<200 pc)  
            // evolve to end of "next birth interval"
            t = t_next;
            //cout << "Time is now " << t << " " << endl;
            //cout << t << " " << endl;
            P = P_next(star_1_mass, star_2_mass, P, t_form, t, z); // change to Pi =
            P_out = P_next(inner_mass, star_3_mass, P_out, t_form, t, z); // change to Po =
            //cout << "Inner is " << P << " The outer is " << P_out << " " << endl;
        } while (t < T_end && t < t_ms); //&& t < t_ms); // do loop

    } while (!is.eof()); // end of file
    is.close();    
    cout << "# Total = " << num << endl;
}
