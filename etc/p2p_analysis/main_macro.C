#include <iostream>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include <math.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TClonesArray.h"
#include "TChain.h"
#include "find_z.C"
#include "find_isotopes.C"
#include "all_params.C"
#include "find_psi_par.C"
using namespace std;

char fname[500];
bool sortcol( const vector<double>& v1,const vector<double>& v2 ) {
                                    return v1[1] > v2[1];
                                    }



void main_macro(char const count_i[50]){
cout << "thsi is count_i " << count_i << endl;
cout << "and this is count_i[50] " << count_i[50] << endl;
stringstream ss;
ss << count_i;
string my_runs;
ss >> my_runs;
all_params(my_runs);
char f_out_name[500];
char hist_name[500];
char f_parameters[500];
//DEFINITIONS OF THE HISTOGRAMS:

TH1F* h1_z_all_iso;
TH1F* h1_a_q_z5;
TH1F* h1_a_q_z6;

TH2F* h2_theta_out_vs_mw3_10b;
TH2F* h2_theta_out_vs_mw3_11b;
TH2F* h2_theta_out_vs_mw3_11c;
TH2F* h2_theta_out_vs_mw3_12c;

TH2F* h2_z_vs_a_q;
TH1F* h1_gamma_energyE_10B;
TH1F* h1_gamma_energyE_11B;
TH1F* h1_gamma_energyE_11B_no_border;
TH1F* h1_gamma_energyE_max_val_11B;
TH2F* h2_gamma_fst_vs_snd_11B;
TH1F* h1_theta_1_plus_theta_2_CALIFA_11b;
TH1F* h1_theta_1_plus_theta_2_CALIFA_10b;
TH2F* h2_E1_vs_E2_CALIFA_10b;
TH2F* h2_E1_vs_E2_CALIFA_11b;

TH1F* h1_E1_plus_E2_CALIFA_10b;
TH1F* h1_E1_plus_E2_CALIFA_11b;

TH2F* h2_long_mom_p1p2_long_mom11b;
TH2F* h2_long_mom_p1p2_long_mom10b;
TH1F* h1_transv_mom_difference_p1_p2_10b;
TH1F* h1_transv_mom_difference_p1_p2_11b;
TH2F* h2_tof_vs_aq_fix_g_11b;
TH2F* h2_tof_vs_aq_fix_g_10b;
TH2F* h2_tof_vs_aq_fix_g_11c;
TH2F* h2_tof_vs_aq_fix_g_12c;
TH2F* h2_E1_plus_E2_CALIFA_vs_full_mom11b;

TH2F* h2_E1_vs_E2_CALIFA_11b_cut;
TH2F* h2_theta_1_vs_theta_2_CALIFA_11b_cut;
TH1F* h1_binding_energy_11b;
TH1F* h1_cos_gamma_11b_p_i;
TH1F* h1_time_diff_gamma_11b;

TH2F* h2_theta1_theta2_11b;
TH2F* h2_psi1_psi2_11b;
TH2F* h2_psi1_psi2_11b_2pi;
TH1F* h1_tof_detnr_strange_12c;
TH1F* h1_tof_detnr_strange_12c_large_aq;
TH1F* h1_mw3_fy_12c;
TH1F* h1_mw3_fy_strange_12c;
TH1F* h1_mw3_fy_strange_12c_large_aq;
TH1F* h1_tof_detnr_12c;
TH2F* h2_z_vs_a_q_nocut;
TH1F* h1_theta_1_CALIFA_11b;
TH1F* h1_theta_2_CALIFA_11b;

TH2F* h2_gamma_energyE_psi_11B;
TH2F* h2_long_mom_p1p2_long_mom11b_no_border;
TH2F* h2_long_mom_p1p2_long_mom11b_high;
TH2F* h2_long_mom_p1p2_long_mom11b_low;
TH2F* h2_long_mom_p1p2_long_mom11b_no_border_high;
TH2F* h2_long_mom_p1p2_long_mom11b_no_border_low;
TH2F* h2_long_mom_p1p2_long_mom11b_small_opening;
TH2F* h2_long_mom_p1p2_long_mom11b_theta_phi_constr;
TH2F* h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom;
TH2F* h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom200;
TH2F* h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom150;
TH1F* h1_m_mom_250_gamma_spec;
TH1F* h1_m_mom_200_gamma_spec;
TH1F* h1_m_mom_150_gamma_spec;
TH1F* h1_cluster_gamma_11b;
TH1F* h1_cluster_proton_11b;
TH2F* h2_long_mom_p1p2_long_mom11b_large_opening;
TH1F* h1_cos_gamma_11b_p_i_cms;
TH1F* h1_cos_gamma_11b_p_i_cms_optimized;
TH1F* h1_psi_in_11B_cut_triangle;
TH1F* h1_psi_in_11B;
TH1F* h1_gamma_energyE_max_val_11B_angle_cut;
TH1F* h1_gamma_spec_escape_angle_cut;
TH2F* h2_gamma_energy_vs_theta_11B;
TH1F* h1_abs_phi_diff_11b;

TH2F* h2_tof_path_strange12c;
TH2F* h2_tof_path_strange_12c_large_aq;
TH2F* h2_tof_path_12c;
TH2F* h2_radius_vs_beta_12c;

TH2F* h2_radius_path_strange12c;
TH2F* h2_radius_path_strange_12c_large_aq;
TH2F* h2_radius_path_12c;
TH2F* h2_mw2x_vs_tof_12c;
TH2F* h2_mw3x_vs_tof_12c;
TH2F* h2_charge_vs_tof_12c;
TH2F* h2_psi_in_vs_tof_12c;
TH2F* h2_detnr_vs_tof_12c;
TH2F* h2_long_mom_p1p2_long_mom11b_large_gamma;
TH2F* h2_radius_vs_tof_12c;

//further analyis histograms for mw2,3 and psi_in..
TH1F* h1_mw3_fy_11b_cut_triangle;
TH1F* h1_mw2_fy_11b_cut_triangle;
TH1F* h1_mw3_fy_11b;
TH1F* h1_mw2_fy_11b;
TH1F* h1_psi_in_11b_cut_triangle;
TH1F* h1_psi_in_11b;

TH1F* h1_angl_diff_y;
TH1F* h1_angl_diff_y_cut_triangle;
//-------------------------------------
TH2F* h2_long_mom_p1p2_long_mom11b_low_optimized;
TH2F* h2_long_mom_p1p2_long_mom11b_high_optimized;
TH2F* h2_long_mom_p1p2_long_mom11b_low_optimized_08;
TH1F* h1_p_missing_11B;
TH1F* h1_cos_low_gamma;
TH2F* h2_theta_sum_vs_phi_diff_11b;

TH1F* h1_tofwall_posy_slow[NbDets];
TH1F* h1_tofwall_posy_medium[NbDets];
TH1F* h1_tofwall_posy_fast[NbDets];

//--------------------------------------------------

//IMPLEMENTATION OF THE HISTOGRAMS:-----------------


//checking tof in sofia TOFW
for (Int_t i = 0; i < NbDets; i++){

    sprintf(hist_name, "TofWall Detnr.%i,pos in y direction in ns for A/q > 1.997, 12C", i + 1);
    h1_tofwall_posy_slow[i] = new TH1F(hist_name, hist_name, 40000, -20, 20);
    h1_tofwall_posy_slow[i]->GetXaxis()->SetTitle("Raw position in y direction [ns]");
    h1_tofwall_posy_slow[i]->GetYaxis()->SetTitle("Counts per bin");
    h1_tofwall_posy_slow[i]->GetXaxis()->CenterTitle(true);
    h1_tofwall_posy_slow[i]->GetYaxis()->CenterTitle(true);
    h1_tofwall_posy_slow[i]->GetYaxis()->SetLabelSize(0.045);
    h1_tofwall_posy_slow[i]->GetYaxis()->SetTitleSize(0.045);
	
	sprintf(hist_name, "TofWall Detnr.%i,pos in y direction in ns for 35.35 < tof < 35.61 ns, 12C", i + 1);
    h1_tofwall_posy_medium[i] = new TH1F(hist_name, hist_name, 40000, -20, 20);
    h1_tofwall_posy_medium[i]->GetXaxis()->SetTitle("Raw position in y direction [ns]");
    h1_tofwall_posy_medium[i]->GetYaxis()->SetTitle("Counts per bin");
    h1_tofwall_posy_medium[i]->GetXaxis()->CenterTitle(true);
    h1_tofwall_posy_medium[i]->GetYaxis()->CenterTitle(true);
    h1_tofwall_posy_medium[i]->GetYaxis()->SetLabelSize(0.045);
    h1_tofwall_posy_medium[i]->GetYaxis()->SetTitleSize(0.045);
	
	sprintf(hist_name, "TofWall Detnr.%i,pos in y direction in ns for tof < 35.13 ns, 12C", i + 1);
    h1_tofwall_posy_fast[i] = new TH1F(hist_name, hist_name, 40000, -20, 20);
    h1_tofwall_posy_fast[i]->GetXaxis()->SetTitle("Raw position in y direction [ns]");
    h1_tofwall_posy_fast[i]->GetYaxis()->SetTitle("Counts per bin");
    h1_tofwall_posy_fast[i]->GetXaxis()->CenterTitle(true);
    h1_tofwall_posy_fast[i]->GetYaxis()->CenterTitle(true);
    h1_tofwall_posy_fast[i]->GetYaxis()->SetLabelSize(0.045);
    h1_tofwall_posy_fast[i]->GetYaxis()->SetTitleSize(0.045);
}



//

sprintf(hist_name, "Z for all isotopes");
h1_z_all_iso = new TH1F(hist_name,hist_name,1000,1,10);
h1_z_all_iso->GetXaxis()->SetTitle("Z (charge)");
h1_z_all_iso->GetYaxis()->SetTitle("# Events ");
h1_z_all_iso->GetXaxis()->CenterTitle(true);
h1_z_all_iso->GetYaxis()->CenterTitle(true);
h1_z_all_iso->GetYaxis()->SetLabelSize(0.045);
h1_z_all_iso->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "A/q for Z = 5");
h1_a_q_z5 = new TH1F(hist_name,hist_name,400,1,3);
h1_a_q_z5->GetXaxis()->SetTitle("A/q for Z = 5");
h1_a_q_z5->GetYaxis()->SetTitle(" # Events");
h1_a_q_z5->GetXaxis()->CenterTitle(true);
h1_a_q_z5->GetYaxis()->CenterTitle(true);
h1_a_q_z5->GetYaxis()->SetLabelSize(0.045);
h1_a_q_z5->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "A/q for Z = 6");
h1_a_q_z6 = new TH1F(hist_name,hist_name,400,1,3);
h1_a_q_z6->GetXaxis()->SetTitle("A/q for Z = 6");
h1_a_q_z6->GetYaxis()->SetTitle(" # Events");
h1_a_q_z6->GetXaxis()->CenterTitle(true);
h1_a_q_z6->GetYaxis()->CenterTitle(true);
h1_a_q_z6->GetYaxis()->SetLabelSize(0.045);
h1_a_q_z6->GetYaxis()->SetTitleSize(0.045);

//- Plots for the third step to get psi_in parameters
sprintf(hist_name, "Theta_out versus MWPC3.fX for 12C");
h2_theta_out_vs_mw3_12c = new TH2F(hist_name,hist_name,800,-400,400,1200,0,0.6);
h2_theta_out_vs_mw3_12c->GetXaxis()->SetTitle("MWPC3 f.X [mm]");
h2_theta_out_vs_mw3_12c->GetYaxis()->SetTitle("Theta_out [degrees]");
h2_theta_out_vs_mw3_12c->GetXaxis()->CenterTitle(true);
h2_theta_out_vs_mw3_12c->GetYaxis()->CenterTitle(true);
h2_theta_out_vs_mw3_12c->GetYaxis()->SetLabelSize(0.045);
h2_theta_out_vs_mw3_12c->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Theta_out versus MWPC3.fX for 11C");
h2_theta_out_vs_mw3_11c = new TH2F(hist_name,hist_name,800,-400,400,1200,0,0.6);
h2_theta_out_vs_mw3_11c->GetXaxis()->SetTitle("MWPC3 f.X [mm]");
h2_theta_out_vs_mw3_11c->GetYaxis()->SetTitle("Theta_out [degrees]");
h2_theta_out_vs_mw3_11c->GetXaxis()->CenterTitle(true);
h2_theta_out_vs_mw3_11c->GetYaxis()->CenterTitle(true);
h2_theta_out_vs_mw3_11c->GetYaxis()->SetLabelSize(0.045);
h2_theta_out_vs_mw3_11c->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Theta_out versus MWPC3.fX for 11B");
h2_theta_out_vs_mw3_11b = new TH2F(hist_name,hist_name,800,-400,400,1200,0,0.6);
h2_theta_out_vs_mw3_11b->GetXaxis()->SetTitle("MWPC3 f.X [mm]");
h2_theta_out_vs_mw3_11b->GetYaxis()->SetTitle("Theta_out [degrees]");
h2_theta_out_vs_mw3_11b->GetXaxis()->CenterTitle(true);
h2_theta_out_vs_mw3_11b->GetYaxis()->CenterTitle(true);
h2_theta_out_vs_mw3_11b->GetYaxis()->SetLabelSize(0.045);
h2_theta_out_vs_mw3_11b->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Theta_out versus MWPC3.fX for 10B");
h2_theta_out_vs_mw3_10b = new TH2F(hist_name,hist_name,800,-400,400,1200,0,0.6);
h2_theta_out_vs_mw3_10b->GetXaxis()->SetTitle("MWPC3 f.X [mm]");
h2_theta_out_vs_mw3_10b->GetYaxis()->SetTitle("Theta_out [degrees]");
h2_theta_out_vs_mw3_10b->GetXaxis()->CenterTitle(true);
h2_theta_out_vs_mw3_10b->GetYaxis()->CenterTitle(true);
h2_theta_out_vs_mw3_10b->GetYaxis()->SetLabelSize(0.045);
h2_theta_out_vs_mw3_10b->GetYaxis()->SetTitleSize(0.045);
//
//first serious try:
sprintf(hist_name, "Z versus A/q");
h2_z_vs_a_q = new TH2F(hist_name,hist_name,200,1,3,100,0,10);
h2_z_vs_a_q->GetXaxis()->SetTitle("A/q");
h2_z_vs_a_q->GetYaxis()->SetTitle("Z (charge) ");
h2_z_vs_a_q->GetXaxis()->CenterTitle(true);
h2_z_vs_a_q->GetYaxis()->CenterTitle(true);
h2_z_vs_a_q->GetYaxis()->SetLabelSize(0.045);
h2_z_vs_a_q->GetYaxis()->SetTitleSize(0.045);

//gamma analysis from CALIFA

sprintf(hist_name, "Gamma Energies E for 10B");
h1_gamma_energyE_10B = new TH1F(hist_name,hist_name,200,0,10);
h1_gamma_energyE_10B->GetXaxis()->SetTitle("Energy E [MeV]");
h1_gamma_energyE_10B->GetYaxis()->SetTitle("Number of entries");
h1_gamma_energyE_10B->GetXaxis()->CenterTitle(true);
h1_gamma_energyE_10B->GetYaxis()->CenterTitle(true);
h1_gamma_energyE_10B->GetYaxis()->SetLabelSize(0.045);
h1_gamma_energyE_10B->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Gamma Energies E for 11B");
h1_gamma_energyE_11B = new TH1F(hist_name,hist_name,200,0,10);
h1_gamma_energyE_11B->GetXaxis()->SetTitle("Energy E [MeV]");
h1_gamma_energyE_11B->GetYaxis()->SetTitle("Number of entries");
h1_gamma_energyE_11B->GetXaxis()->CenterTitle(true);
h1_gamma_energyE_11B->GetYaxis()->CenterTitle(true);
h1_gamma_energyE_11B->GetYaxis()->SetLabelSize(0.045);
h1_gamma_energyE_11B->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Highest Gamma Energy E for 11B");
h1_gamma_energyE_max_val_11B = new TH1F(hist_name,hist_name,200,0,10);
h1_gamma_energyE_max_val_11B->GetXaxis()->SetTitle("Energy E [MeV]");
h1_gamma_energyE_max_val_11B->GetYaxis()->SetTitle("Number of entries");
h1_gamma_energyE_max_val_11B->GetXaxis()->CenterTitle(true);
h1_gamma_energyE_max_val_11B->GetYaxis()->CenterTitle(true);
h1_gamma_energyE_max_val_11B->GetYaxis()->SetLabelSize(0.045);
h1_gamma_energyE_max_val_11B->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Highest Gamma Energy E  vs second highest for 11B");
h2_gamma_fst_vs_snd_11B = new TH2F(hist_name,hist_name,200,0,10,200,0,10);
h2_gamma_fst_vs_snd_11B->GetXaxis()->SetTitle("Gamma 2nd highestEnergy E [MeV]");
h2_gamma_fst_vs_snd_11B->GetYaxis()->SetTitle("Gamma highest Energy E [MeV]");
h2_gamma_fst_vs_snd_11B->GetXaxis()->CenterTitle(true);
h2_gamma_fst_vs_snd_11B->GetYaxis()->CenterTitle(true);
h2_gamma_fst_vs_snd_11B->GetYaxis()->SetLabelSize(0.045);
h2_gamma_fst_vs_snd_11B->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Highest Gamma Energy E for 11B with theta and phi cut");
h1_gamma_energyE_max_val_11B_angle_cut = new TH1F(hist_name,hist_name,200,0,10);
h1_gamma_energyE_max_val_11B_angle_cut->GetXaxis()->SetTitle("Energy E [MeV]");
h1_gamma_energyE_max_val_11B_angle_cut->GetYaxis()->SetTitle("Number of entries");
h1_gamma_energyE_max_val_11B_angle_cut->GetXaxis()->CenterTitle(true);
h1_gamma_energyE_max_val_11B_angle_cut->GetYaxis()->CenterTitle(true);
h1_gamma_energyE_max_val_11B_angle_cut->GetYaxis()->SetLabelSize(0.045);
h1_gamma_energyE_max_val_11B_angle_cut->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Gamma spectrum with single and double escape for  11B with theta and phi cut");
h1_gamma_spec_escape_angle_cut = new TH1F(hist_name,hist_name,200,0,10);
h1_gamma_spec_escape_angle_cut->GetXaxis()->SetTitle("Energy E [MeV]");
h1_gamma_spec_escape_angle_cut->GetYaxis()->SetTitle("Number of entries");
h1_gamma_spec_escape_angle_cut->GetXaxis()->CenterTitle(true);
h1_gamma_spec_escape_angle_cut->GetYaxis()->CenterTitle(true);
h1_gamma_spec_escape_angle_cut->GetYaxis()->SetLabelSize(0.045);
h1_gamma_spec_escape_angle_cut->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Gamma Energies versus Theta for 11B with theta and phi cut");
h2_gamma_energy_vs_theta_11B = new TH2F(hist_name,hist_name,200,0,10,12,0,120);
h2_gamma_energy_vs_theta_11B->GetXaxis()->SetTitle("Energy E [MeV]");
h2_gamma_energy_vs_theta_11B->GetYaxis()->SetTitle("Theta [degrees]");
h2_gamma_energy_vs_theta_11B->GetXaxis()->CenterTitle(true);
h2_gamma_energy_vs_theta_11B->GetYaxis()->CenterTitle(true);
h2_gamma_energy_vs_theta_11B->GetYaxis()->SetLabelSize(0.045);
h2_gamma_energy_vs_theta_11B->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Gamma Energies E for 11B without border");
h1_gamma_energyE_11B_no_border = new TH1F(hist_name,hist_name,200,0,10);
h1_gamma_energyE_11B_no_border->GetXaxis()->SetTitle("Energy E [MeV]");
h1_gamma_energyE_11B_no_border->GetYaxis()->SetTitle("Number of entries");
h1_gamma_energyE_11B_no_border->GetXaxis()->CenterTitle(true);
h1_gamma_energyE_11B_no_border->GetYaxis()->CenterTitle(true);
h1_gamma_energyE_11B_no_border->GetYaxis()->SetLabelSize(0.045);
h1_gamma_energyE_11B_no_border->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Gamma Energies versus Psi for  11B");
h2_gamma_energyE_psi_11B = new TH2F(hist_name,hist_name,200,0,10,60,0,360);
h2_gamma_energyE_psi_11B->GetXaxis()->SetTitle("Energy E [MeV]");
h2_gamma_energyE_psi_11B->GetYaxis()->SetTitle("Psi [degrees]");
h2_gamma_energyE_psi_11B->GetXaxis()->CenterTitle(true);
h2_gamma_energyE_psi_11B->GetYaxis()->CenterTitle(true);
h2_gamma_energyE_psi_11B->GetYaxis()->SetLabelSize(0.045);
h2_gamma_energyE_psi_11B->GetYaxis()->SetTitleSize(0.045);

//check psi_in for all 11B and for the triangle...
sprintf(hist_name, "Theta_in for 11B");
h1_psi_in_11B = new TH1F(hist_name,hist_name,200,-0.1,0.1);
h1_psi_in_11B->GetXaxis()->SetTitle("Theta_in [rad]");
h1_psi_in_11B->GetYaxis()->SetTitle("Number of entries");
h1_psi_in_11B->GetXaxis()->CenterTitle(true);
h1_psi_in_11B->GetYaxis()->CenterTitle(true);
h1_psi_in_11B->GetYaxis()->SetLabelSize(0.045);
h1_psi_in_11B->GetYaxis()->SetTitleSize(0.045);

//cut triangle:
//check psi_in for all 11B and for the triangle...
sprintf(hist_name, "Theta_in for 11B in cut triangle");
h1_psi_in_11B_cut_triangle = new TH1F(hist_name,hist_name,200,-0.1,0.1);
h1_psi_in_11B_cut_triangle->GetXaxis()->SetTitle("Theta_in [rad]");
h1_psi_in_11B_cut_triangle->GetYaxis()->SetTitle("Number of entries");
h1_psi_in_11B_cut_triangle->GetXaxis()->CenterTitle(true);
h1_psi_in_11B_cut_triangle->GetYaxis()->CenterTitle(true);
h1_psi_in_11B_cut_triangle->GetYaxis()->SetLabelSize(0.045);
h1_psi_in_11B_cut_triangle->GetYaxis()->SetTitleSize(0.045);

//here I do some Mw3 and mw2 and psi_in analysis for the cut triangle and the rest and compare it, to see if the strange cut triangle comes from straggling in Y direction...

sprintf(hist_name, "MW3.fY for the cut triangle , 11B");
h1_mw3_fy_11b_cut_triangle = new TH1F(hist_name,hist_name,400,-400,400);
h1_mw3_fy_11b_cut_triangle->GetXaxis()->SetTitle("MW3.fY [mm]");
h1_mw3_fy_11b_cut_triangle->GetYaxis()->SetTitle("No. of Events");
h1_mw3_fy_11b_cut_triangle->GetXaxis()->CenterTitle(true);
h1_mw3_fy_11b_cut_triangle->GetYaxis()->CenterTitle(true);
h1_mw3_fy_11b_cut_triangle->GetYaxis()->SetLabelSize(0.045);
h1_mw3_fy_11b_cut_triangle->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "MW2.fY for the cut triangle , 11B");
h1_mw2_fy_11b_cut_triangle = new TH1F(hist_name,hist_name,400,-400,400);
h1_mw2_fy_11b_cut_triangle->GetXaxis()->SetTitle("MW2.fY [mm]");
h1_mw2_fy_11b_cut_triangle->GetYaxis()->SetTitle("No. of Events");
h1_mw2_fy_11b_cut_triangle->GetXaxis()->CenterTitle(true);
h1_mw2_fy_11b_cut_triangle->GetYaxis()->CenterTitle(true);
h1_mw2_fy_11b_cut_triangle->GetYaxis()->SetLabelSize(0.045);
h1_mw2_fy_11b_cut_triangle->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "MW3.fY for anti-correlated line , 11B");
h1_mw3_fy_11b = new TH1F(hist_name,hist_name,400,-400,400);
h1_mw3_fy_11b->GetXaxis()->SetTitle("MW3.fY [mm]");
h1_mw3_fy_11b->GetYaxis()->SetTitle("No. of Events");
h1_mw3_fy_11b->GetXaxis()->CenterTitle(true);
h1_mw3_fy_11b->GetYaxis()->CenterTitle(true);
h1_mw3_fy_11b->GetYaxis()->SetLabelSize(0.045);
h1_mw3_fy_11b->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "MW2.fY for anti-correlated line , 11B");
h1_mw2_fy_11b = new TH1F(hist_name,hist_name,400,-400,400);
h1_mw2_fy_11b->GetXaxis()->SetTitle("MW2.fY [mm]");
h1_mw2_fy_11b->GetYaxis()->SetTitle("No. of Events");
h1_mw2_fy_11b->GetXaxis()->CenterTitle(true);
h1_mw2_fy_11b->GetYaxis()->CenterTitle(true);
h1_mw2_fy_11b->GetYaxis()->SetLabelSize(0.045);
h1_mw2_fy_11b->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Psi_in for the cut triangle , 11B");
h1_psi_in_11b_cut_triangle = new TH1F(hist_name,hist_name,2000,-0.1,0.1);
h1_psi_in_11b_cut_triangle->GetXaxis()->SetTitle("Psi_in  [degrees]");
h1_psi_in_11b_cut_triangle->GetYaxis()->SetTitle("No. of Events");
h1_psi_in_11b_cut_triangle->GetXaxis()->CenterTitle(true);
h1_psi_in_11b_cut_triangle->GetYaxis()->CenterTitle(true);
h1_psi_in_11b_cut_triangle->GetYaxis()->SetLabelSize(0.045);
h1_psi_in_11b_cut_triangle->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Psi_in for anti-correlated line, 11B");
h1_psi_in_11b = new TH1F(hist_name,hist_name,2000,-0.1,0.1);
h1_psi_in_11b->GetXaxis()->SetTitle("Psi_in [degrees]");
h1_psi_in_11b->GetYaxis()->SetTitle("No. of Events");
h1_psi_in_11b->GetXaxis()->CenterTitle(true);
h1_psi_in_11b->GetYaxis()->CenterTitle(true);
h1_psi_in_11b->GetYaxis()->SetLabelSize(0.045);
h1_psi_in_11b->GetYaxis()->SetTitleSize(0.045);


//-------------------------------------------------

//more CALIFA
sprintf(hist_name, "Theta1 plus Theta2 for CALIFA 11B");
h1_theta_1_plus_theta_2_CALIFA_11b = new TH1F(hist_name,hist_name,52,22.15,152.15);
h1_theta_1_plus_theta_2_CALIFA_11b->GetXaxis()->SetTitle("Theta1 plus Theta2 [degrees]");
h1_theta_1_plus_theta_2_CALIFA_11b->GetYaxis()->SetTitle("No. of Events");
h1_theta_1_plus_theta_2_CALIFA_11b->GetXaxis()->CenterTitle(true);
h1_theta_1_plus_theta_2_CALIFA_11b->GetYaxis()->CenterTitle(true);
h1_theta_1_plus_theta_2_CALIFA_11b->GetYaxis()->SetLabelSize(0.045);
h1_theta_1_plus_theta_2_CALIFA_11b->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Theta1 plus Theta2 versus phi difference for 11B");
h2_theta_sum_vs_phi_diff_11b = new TH2F(hist_name,hist_name,52,22.15,152.15,60,0,360);
h2_theta_sum_vs_phi_diff_11b->GetXaxis()->SetTitle("Theta sum [degrees]");
h2_theta_sum_vs_phi_diff_11b->GetYaxis()->SetTitle(" Phi difference [degrees]");
h2_theta_sum_vs_phi_diff_11b->GetXaxis()->CenterTitle(true);
h2_theta_sum_vs_phi_diff_11b->GetYaxis()->CenterTitle(true);
h2_theta_sum_vs_phi_diff_11b->GetYaxis()->SetLabelSize(0.045);
h2_theta_sum_vs_phi_diff_11b->GetYaxis()->SetTitleSize(0.045);



//here to check binning of angles:
sprintf(hist_name, "Theta1 for  CALIFA 11B");
h1_theta_1_CALIFA_11b = new TH1F(hist_name,hist_name,1500,0,150);
h1_theta_1_CALIFA_11b->GetXaxis()->SetTitle("Theta1 [degrees]");
h1_theta_1_CALIFA_11b->GetYaxis()->SetTitle("No. of Events");
h1_theta_1_CALIFA_11b->GetXaxis()->CenterTitle(true);
h1_theta_1_CALIFA_11b->GetYaxis()->CenterTitle(true);
h1_theta_1_CALIFA_11b->GetYaxis()->SetLabelSize(0.045);
h1_theta_1_CALIFA_11b->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Theta2 for  CALIFA 11B");
h1_theta_2_CALIFA_11b = new TH1F(hist_name,hist_name,1500,0,150);
h1_theta_2_CALIFA_11b->GetXaxis()->SetTitle("Theta2 [degrees]");
h1_theta_2_CALIFA_11b->GetYaxis()->SetTitle("No. of Events");
h1_theta_2_CALIFA_11b->GetXaxis()->CenterTitle(true);
h1_theta_2_CALIFA_11b->GetYaxis()->CenterTitle(true);
h1_theta_2_CALIFA_11b->GetYaxis()->SetLabelSize(0.045);
h1_theta_2_CALIFA_11b->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Abs(phi_difference) for identified 11B");
h1_abs_phi_diff_11b = new TH1F(hist_name,hist_name,60,0,360);
h1_abs_phi_diff_11b->GetXaxis()->SetTitle("Absolute phi difference of the two protons [degrees]");
h1_abs_phi_diff_11b->GetYaxis()->SetTitle("No. of Events");
h1_abs_phi_diff_11b->GetXaxis()->CenterTitle(true);
h1_abs_phi_diff_11b->GetYaxis()->CenterTitle(true);
h1_abs_phi_diff_11b->GetYaxis()->SetLabelSize(0.045);
h1_abs_phi_diff_11b->GetYaxis()->SetTitleSize(0.045);


//--------------------
sprintf(hist_name, "Theta1 plus Theta2 for CALIFA 10B");
h1_theta_1_plus_theta_2_CALIFA_10b = new TH1F(hist_name,hist_name,52,22.15,152.15);
h1_theta_1_plus_theta_2_CALIFA_10b->GetXaxis()->SetTitle("Theta1 plus Theta2 [degrees]");
h1_theta_1_plus_theta_2_CALIFA_10b->GetYaxis()->SetTitle("No. of Events");
h1_theta_1_plus_theta_2_CALIFA_10b->GetXaxis()->CenterTitle(true);
h1_theta_1_plus_theta_2_CALIFA_10b->GetYaxis()->CenterTitle(true);
h1_theta_1_plus_theta_2_CALIFA_10b->GetYaxis()->SetLabelSize(0.045);
h1_theta_1_plus_theta_2_CALIFA_10b->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "E1 versus E2 for CALIFA for 11B");
h2_E1_vs_E2_CALIFA_11b = new TH2F(hist_name,hist_name,50,0,500,50,0,500);
h2_E1_vs_E2_CALIFA_11b->GetXaxis()->SetTitle("E1 [MeV]");
h2_E1_vs_E2_CALIFA_11b->GetYaxis()->SetTitle("E2 [MeV]");
h2_E1_vs_E2_CALIFA_11b->GetXaxis()->CenterTitle(true);
h2_E1_vs_E2_CALIFA_11b->GetYaxis()->CenterTitle(true);
h2_E1_vs_E2_CALIFA_11b->GetYaxis()->SetLabelSize(0.045);
h2_E1_vs_E2_CALIFA_11b->GetYaxis()->SetTitleSize(0.045);



sprintf(hist_name, "E1 versus E2 for CALIFA for 10B");
h2_E1_vs_E2_CALIFA_10b = new TH2F(hist_name,hist_name,50,0,500,50,0,500);
h2_E1_vs_E2_CALIFA_10b->GetXaxis()->SetTitle("E1 [MeV]");
h2_E1_vs_E2_CALIFA_10b->GetYaxis()->SetTitle("E2 [MeV]");
h2_E1_vs_E2_CALIFA_10b->GetXaxis()->CenterTitle(true);
h2_E1_vs_E2_CALIFA_10b->GetYaxis()->CenterTitle(true);
h2_E1_vs_E2_CALIFA_10b->GetYaxis()->SetLabelSize(0.045);
h2_E1_vs_E2_CALIFA_10b->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "E1 plus E2 for CALIFA 10B");
h1_E1_plus_E2_CALIFA_10b = new TH1F(hist_name,hist_name,100,0,1000);
h1_E1_plus_E2_CALIFA_10b->GetXaxis()->SetTitle("E1 plus E2 [MeV]");
h1_E1_plus_E2_CALIFA_10b->GetYaxis()->SetTitle("No. of Events");
h1_E1_plus_E2_CALIFA_10b->GetXaxis()->CenterTitle(true);
h1_E1_plus_E2_CALIFA_10b->GetYaxis()->CenterTitle(true);
h1_E1_plus_E2_CALIFA_10b->GetYaxis()->SetLabelSize(0.045);
h1_E1_plus_E2_CALIFA_10b->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "E1 plus E2 for CALIFA 11B");
h1_E1_plus_E2_CALIFA_11b = new TH1F(hist_name,hist_name,100,0,1000);
h1_E1_plus_E2_CALIFA_11b->GetXaxis()->SetTitle("E1 plus E2 [MeV]");
h1_E1_plus_E2_CALIFA_11b->GetYaxis()->SetTitle("No. of Events");
h1_E1_plus_E2_CALIFA_11b->GetXaxis()->CenterTitle(true);
h1_E1_plus_E2_CALIFA_11b->GetYaxis()->CenterTitle(true);
h1_E1_plus_E2_CALIFA_11b->GetYaxis()->SetLabelSize(0.045);
h1_E1_plus_E2_CALIFA_11b->GetYaxis()->SetTitleSize(0.045);

if (gamma_given == 1.42942){

//these histos are for checking strange behavior of tof vs A/q at 12C:
sprintf(hist_name, "TOF Wall Detector Number for tof < 35.13 ns , 12C");
h1_tof_detnr_strange_12c = new TH1F(hist_name,hist_name,30,0,30);
h1_tof_detnr_strange_12c->GetXaxis()->SetTitle("DetectorNumber TOFW");
h1_tof_detnr_strange_12c->GetYaxis()->SetTitle("No. of Events");
h1_tof_detnr_strange_12c->GetXaxis()->CenterTitle(true);
h1_tof_detnr_strange_12c->GetYaxis()->CenterTitle(true);
h1_tof_detnr_strange_12c->GetYaxis()->SetLabelSize(0.045);
h1_tof_detnr_strange_12c->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "TOF Wall Detector Number for 35.35 < tof < 35.61 ns , 12C");
h1_tof_detnr_12c = new TH1F(hist_name,hist_name,30,0,30);
h1_tof_detnr_12c->GetXaxis()->SetTitle("DetectorNumber TOFW");
h1_tof_detnr_12c->GetYaxis()->SetTitle("No. of Events");
h1_tof_detnr_12c->GetXaxis()->CenterTitle(true);
h1_tof_detnr_12c->GetYaxis()->CenterTitle(true);
h1_tof_detnr_12c->GetYaxis()->SetLabelSize(0.045);
h1_tof_detnr_12c->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "TOF Wall Detector Number for A/q > 1.997, 12C");
h1_tof_detnr_strange_12c_large_aq = new TH1F(hist_name,hist_name,30,0,30);
h1_tof_detnr_strange_12c_large_aq->GetXaxis()->SetTitle("DetectorNumber TOFW");
h1_tof_detnr_strange_12c_large_aq->GetYaxis()->SetTitle("No. of Events");
h1_tof_detnr_strange_12c_large_aq->GetXaxis()->CenterTitle(true);
h1_tof_detnr_strange_12c_large_aq->GetYaxis()->CenterTitle(true);
h1_tof_detnr_strange_12c_large_aq->GetYaxis()->SetLabelSize(0.045);
h1_tof_detnr_strange_12c_large_aq->GetYaxis()->SetTitleSize(0.045);

//Y inspection in MW3
sprintf(hist_name, "MW3.fY for tof < 35.13 ns , 12C");
h1_mw3_fy_strange_12c = new TH1F(hist_name,hist_name,400,-400,400);
h1_mw3_fy_strange_12c->GetXaxis()->SetTitle("MW3.fY [mm]");
h1_mw3_fy_strange_12c->GetYaxis()->SetTitle("No. of Events");
h1_mw3_fy_strange_12c->GetXaxis()->CenterTitle(true);
h1_mw3_fy_strange_12c->GetYaxis()->CenterTitle(true);
h1_mw3_fy_strange_12c->GetYaxis()->SetLabelSize(0.045);
h1_mw3_fy_strange_12c->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "MW3.fY for 35.35 < tof < 35.61 ns , 12C");
h1_mw3_fy_12c = new TH1F(hist_name,hist_name,400,-400,400);
h1_mw3_fy_12c->GetXaxis()->SetTitle("MW3.fY [mm]");
h1_mw3_fy_12c->GetYaxis()->SetTitle("No. of Events");
h1_mw3_fy_12c->GetXaxis()->CenterTitle(true);
h1_mw3_fy_12c->GetYaxis()->CenterTitle(true);
h1_mw3_fy_12c->GetYaxis()->SetLabelSize(0.045);
h1_mw3_fy_12c->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "MW3.fY for  A/q > 1.997, 12C");
h1_mw3_fy_strange_12c_large_aq = new TH1F(hist_name,hist_name,400,-400,400);
h1_mw3_fy_strange_12c_large_aq->GetXaxis()->SetTitle("MW3.fY [mm]");
h1_mw3_fy_strange_12c_large_aq->GetYaxis()->SetTitle("No. of Events");
h1_mw3_fy_strange_12c_large_aq->GetXaxis()->CenterTitle(true);
h1_mw3_fy_strange_12c_large_aq->GetYaxis()->CenterTitle(true);
h1_mw3_fy_strange_12c_large_aq->GetYaxis()->SetLabelSize(0.045);
h1_mw3_fy_strange_12c_large_aq->GetYaxis()->SetTitleSize(0.045);

//check tof vs pathlength...
sprintf(hist_name, "Time of flight vs pathlength, for tof < 35.13 ns , 12C");
h2_tof_path_strange12c = new TH2F(hist_name,hist_name,400,7300,7700,1000,30,40);
h2_tof_path_strange12c->GetXaxis()->SetTitle("Pathlength [mm]");
h2_tof_path_strange12c->GetYaxis()->SetTitle("ToF [ns]");
h2_tof_path_strange12c->GetXaxis()->CenterTitle(true);
h2_tof_path_strange12c->GetYaxis()->CenterTitle(true);
h2_tof_path_strange12c->GetYaxis()->SetLabelSize(0.045);
h2_tof_path_strange12c->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Time of flight vs pathlength, for A/q > 1.997, 12C");
h2_tof_path_strange_12c_large_aq = new TH2F(hist_name,hist_name,400,7300,7700,1000,30,40);
h2_tof_path_strange_12c_large_aq->GetXaxis()->SetTitle("Pathlength [mm]");
h2_tof_path_strange_12c_large_aq->GetYaxis()->SetTitle("ToF [ns]");
h2_tof_path_strange_12c_large_aq->GetXaxis()->CenterTitle(true);
h2_tof_path_strange_12c_large_aq->GetYaxis()->CenterTitle(true);
h2_tof_path_strange_12c_large_aq->GetYaxis()->SetLabelSize(0.045);
h2_tof_path_strange_12c_large_aq->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Time of flight vs pathlength, for 35.35 < tof < 35.61 ns , 12C");
h2_tof_path_12c = new TH2F(hist_name,hist_name,400,7300,7700,1000,30,40);
h2_tof_path_12c->GetXaxis()->SetTitle("Pathlength [mm]");
h2_tof_path_12c->GetYaxis()->SetTitle("ToF [ns]");
h2_tof_path_12c->GetXaxis()->CenterTitle(true);
h2_tof_path_12c->GetYaxis()->CenterTitle(true);
h2_tof_path_12c->GetYaxis()->SetLabelSize(0.045);
h2_tof_path_12c->GetYaxis()->SetTitleSize(0.045);

//check radius vs time of flight....
sprintf(hist_name, "Radius vs time of flight, for tof < 35.13 ns , 12C");
h2_radius_path_strange12c = new TH2F(hist_name,hist_name,400,6300,6700,1000,30,40);
h2_radius_path_strange12c->GetXaxis()->SetTitle("Radius [mm]");
h2_radius_path_strange12c->GetYaxis()->SetTitle("ToF [ns]");
h2_radius_path_strange12c->GetXaxis()->CenterTitle(true);
h2_radius_path_strange12c->GetYaxis()->CenterTitle(true);
h2_radius_path_strange12c->GetYaxis()->SetLabelSize(0.045);
h2_radius_path_strange12c->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Radius vs time of flight, for A/q > 1.997, 12C");
h2_radius_path_strange_12c_large_aq = new TH2F(hist_name,hist_name,400,6300,6700,1000,30,40);
h2_radius_path_strange_12c_large_aq->GetXaxis()->SetTitle("Radius [mm]");
h2_radius_path_strange_12c_large_aq->GetYaxis()->SetTitle("ToF [ns]");
h2_radius_path_strange_12c_large_aq->GetXaxis()->CenterTitle(true);
h2_radius_path_strange_12c_large_aq->GetYaxis()->CenterTitle(true);
h2_radius_path_strange_12c_large_aq->GetYaxis()->SetLabelSize(0.045);
h2_radius_path_strange_12c_large_aq->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Radius vs time of flight, for 35.35 < tof < 35.61 ns , 12C");
h2_radius_path_12c = new TH2F(hist_name,hist_name,400,6300,6700,1000,30,40);
h2_radius_path_12c->GetXaxis()->SetTitle("Radius [mm]");
h2_radius_path_12c->GetYaxis()->SetTitle("ToF [ns]");
h2_radius_path_12c->GetXaxis()->CenterTitle(true);
h2_radius_path_12c->GetYaxis()->CenterTitle(true);
h2_radius_path_12c->GetYaxis()->SetLabelSize(0.045);
h2_radius_path_12c->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Radius vs beta, for 12C");
h2_radius_vs_beta_12c = new TH2F(hist_name,hist_name,1000,0,1,400,6300,6700);
h2_radius_vs_beta_12c->GetXaxis()->SetTitle("beta");
h2_radius_vs_beta_12c->GetYaxis()->SetTitle("Radius [mm]");
h2_radius_vs_beta_12c->GetXaxis()->CenterTitle(true);
h2_radius_vs_beta_12c->GetYaxis()->CenterTitle(true);
h2_radius_vs_beta_12c->GetYaxis()->SetLabelSize(0.045);
h2_radius_vs_beta_12c->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Radius vs time of flight for  12C");
h2_radius_vs_tof_12c = new TH2F(hist_name,hist_name,400,6300,6700,1000,30,40);
h2_radius_vs_tof_12c->GetXaxis()->SetTitle("Radius [mm]");
h2_radius_vs_tof_12c->GetYaxis()->SetTitle("ToF [ns]");
h2_radius_vs_tof_12c->GetXaxis()->CenterTitle(true);
h2_radius_vs_tof_12c->GetYaxis()->CenterTitle(true);
h2_radius_vs_tof_12c->GetYaxis()->SetLabelSize(0.045);
h2_radius_vs_tof_12c->GetYaxis()->SetTitleSize(0.045);

//check mw2.x
sprintf(hist_name, "MW2.X vs time of flight for 12C");
h2_mw2x_vs_tof_12c = new TH2F(hist_name,hist_name,400,-400,400,1000,30,40);
h2_mw2x_vs_tof_12c->GetXaxis()->SetTitle("MW2.fX [mm]");
h2_mw2x_vs_tof_12c->GetYaxis()->SetTitle("ToF [ns]");
h2_mw2x_vs_tof_12c->GetXaxis()->CenterTitle(true);
h2_mw2x_vs_tof_12c->GetYaxis()->CenterTitle(true);
h2_mw2x_vs_tof_12c->GetYaxis()->SetLabelSize(0.045);
h2_mw2x_vs_tof_12c->GetYaxis()->SetTitleSize(0.045);

//check mw3.x
sprintf(hist_name, "MW3.X vs time of flight for 12C");
h2_mw3x_vs_tof_12c = new TH2F(hist_name,hist_name,400,-400,400,1000,30,40);
h2_mw3x_vs_tof_12c->GetXaxis()->SetTitle("MW3.fX [mm]");
h2_mw3x_vs_tof_12c->GetYaxis()->SetTitle("ToF [ns]");
h2_mw3x_vs_tof_12c->GetXaxis()->CenterTitle(true);
h2_mw3x_vs_tof_12c->GetYaxis()->CenterTitle(true);
h2_mw3x_vs_tof_12c->GetYaxis()->SetLabelSize(0.045);
h2_mw3x_vs_tof_12c->GetYaxis()->SetTitleSize(0.045);

//check charge...
sprintf(hist_name, "Charge vs time of flight for 12C");
h2_charge_vs_tof_12c = new TH2F(hist_name,hist_name,100,0,10,1000,30,40);
h2_charge_vs_tof_12c->GetXaxis()->SetTitle("Charge");
h2_charge_vs_tof_12c->GetYaxis()->SetTitle("ToF [ns]");
h2_charge_vs_tof_12c->GetXaxis()->CenterTitle(true);
h2_charge_vs_tof_12c->GetYaxis()->CenterTitle(true);
h2_charge_vs_tof_12c->GetYaxis()->SetLabelSize(0.045);
h2_charge_vs_tof_12c->GetYaxis()->SetTitleSize(0.045);

//check psi_in
sprintf(hist_name, "Psi_in vs time of flight for 12C");
h2_psi_in_vs_tof_12c = new TH2F(hist_name,hist_name,200,-0.1,0.1,1000,30,40);
h2_psi_in_vs_tof_12c->GetXaxis()->SetTitle("Psi_in [degrees]");
h2_psi_in_vs_tof_12c->GetYaxis()->SetTitle("ToF [ns]");
h2_psi_in_vs_tof_12c->GetXaxis()->CenterTitle(true);
h2_psi_in_vs_tof_12c->GetYaxis()->CenterTitle(true);
h2_psi_in_vs_tof_12c->GetYaxis()->SetLabelSize(0.045);
h2_psi_in_vs_tof_12c->GetYaxis()->SetTitleSize(0.045);


//check detnr...
sprintf(hist_name, "Detnr TOFW vs time of flight for 12C");
h2_detnr_vs_tof_12c = new TH2F(hist_name,hist_name,30,0,30,1000,30,40);
h2_detnr_vs_tof_12c->GetXaxis()->SetTitle("Detnr TOFW");
h2_detnr_vs_tof_12c->GetYaxis()->SetTitle("ToF [ns]");
h2_detnr_vs_tof_12c->GetXaxis()->CenterTitle(true);
h2_detnr_vs_tof_12c->GetYaxis()->CenterTitle(true);
h2_detnr_vs_tof_12c->GetYaxis()->SetLabelSize(0.045);
h2_detnr_vs_tof_12c->GetYaxis()->SetTitleSize(0.045);



//

sprintf(hist_name, "Z versus A/q all isotopes");
h2_z_vs_a_q_nocut = new TH2F(hist_name,hist_name,200,1,3,100,0,10);
h2_z_vs_a_q_nocut->GetXaxis()->SetTitle("A/q");
h2_z_vs_a_q_nocut->GetYaxis()->SetTitle("Z (charge) ");
h2_z_vs_a_q_nocut->GetXaxis()->CenterTitle(true);
h2_z_vs_a_q_nocut->GetYaxis()->CenterTitle(true);
h2_z_vs_a_q_nocut->GetYaxis()->SetLabelSize(0.045);
h2_z_vs_a_q_nocut->GetYaxis()->SetTitleSize(0.045);


//

sprintf(hist_name, "Theta1 vs Theta2 in CALIFA");
h2_theta1_theta2_11b = new TH2F(hist_name,hist_name,25,22.15,84.65,25,22.15,84.65);
h2_theta1_theta2_11b->GetXaxis()->SetTitle("Theta1 [degrees]");
h2_theta1_theta2_11b->GetYaxis()->SetTitle("Theta2 [degrees]");
h2_theta1_theta2_11b->GetXaxis()->CenterTitle(true);
h2_theta1_theta2_11b->GetYaxis()->CenterTitle(true);
h2_theta1_theta2_11b->GetYaxis()->SetLabelSize(0.045);
h2_theta1_theta2_11b->GetYaxis()->SetTitleSize(0.045);

//


sprintf(hist_name, "Phi1 vs Phi2 in CALIFA");
h2_psi1_psi2_11b = new TH2F(hist_name,hist_name,60,0,360,60,0,360);
h2_psi1_psi2_11b->GetXaxis()->SetTitle("Psi1 [degrees]");
h2_psi1_psi2_11b->GetYaxis()->SetTitle("Psi2 [degrees]");
h2_psi1_psi2_11b->GetXaxis()->CenterTitle(true);
h2_psi1_psi2_11b->GetYaxis()->CenterTitle(true);
h2_psi1_psi2_11b->GetYaxis()->SetLabelSize(0.045);
h2_psi1_psi2_11b->GetYaxis()->SetTitleSize(0.045);


//2pi test for phi...

sprintf(hist_name, "Phi1 vs Phi2 in CALIFA with 2Pi");
h2_psi1_psi2_11b_2pi = new TH2F(hist_name,hist_name,60,0,360,60,0,360);
h2_psi1_psi2_11b_2pi->GetXaxis()->SetTitle("Psi1 [degrees]");
h2_psi1_psi2_11b_2pi->GetYaxis()->SetTitle("Psi2 [degrees]");
h2_psi1_psi2_11b_2pi->GetXaxis()->CenterTitle(true);
h2_psi1_psi2_11b_2pi->GetYaxis()->CenterTitle(true);
h2_psi1_psi2_11b_2pi->GetYaxis()->SetLabelSize(0.045);
h2_psi1_psi2_11b_2pi->GetYaxis()->SetTitleSize(0.045);
//


sprintf(hist_name, "cos(gamma) in the z-x plane for 11B and p_i");
h1_cos_gamma_11b_p_i = new TH1F(hist_name,hist_name,40,-2,2);
h1_cos_gamma_11b_p_i->GetXaxis()->SetTitle("cos(gamma)");
h1_cos_gamma_11b_p_i->GetYaxis()->SetTitle("No. of Events");
h1_cos_gamma_11b_p_i->GetXaxis()->CenterTitle(true);
h1_cos_gamma_11b_p_i->GetYaxis()->CenterTitle(true);
h1_cos_gamma_11b_p_i->GetYaxis()->SetLabelSize(0.045);
h1_cos_gamma_11b_p_i->GetYaxis()->SetTitleSize(0.045);


//angle between p_in and p_11b in 12C frame...
sprintf(hist_name, "cos(gamma) in the z-x plane for 11B and p_i in 12C rest frame");
h1_cos_gamma_11b_p_i_cms = new TH1F(hist_name,hist_name,50,-1,1);
h1_cos_gamma_11b_p_i_cms->GetXaxis()->SetTitle("cos(gamma)");
h1_cos_gamma_11b_p_i_cms->GetYaxis()->SetTitle("No. of Events");
h1_cos_gamma_11b_p_i_cms->GetXaxis()->CenterTitle(true);
h1_cos_gamma_11b_p_i_cms->GetYaxis()->CenterTitle(true);
h1_cos_gamma_11b_p_i_cms->GetYaxis()->SetLabelSize(0.045);
h1_cos_gamma_11b_p_i_cms->GetYaxis()->SetTitleSize(0.045);

//here the optimized version of h1_cos_gamma_11b_p_i_cms sweeping both phi1 and ph2 +- 2.8 degrees...
sprintf(hist_name, "cos(gamma) in the z-x plane for 11B and p_i in 12C rest frame  optimized version");
h1_cos_gamma_11b_p_i_cms_optimized = new TH1F(hist_name,hist_name,50,-1,1);
h1_cos_gamma_11b_p_i_cms_optimized->GetXaxis()->SetTitle("cos(gamma)");
h1_cos_gamma_11b_p_i_cms_optimized->GetYaxis()->SetTitle("No. of Events");
h1_cos_gamma_11b_p_i_cms_optimized->GetXaxis()->CenterTitle(true);
h1_cos_gamma_11b_p_i_cms_optimized->GetYaxis()->CenterTitle(true);
h1_cos_gamma_11b_p_i_cms_optimized->GetYaxis()->SetLabelSize(0.045);
h1_cos_gamma_11b_p_i_cms_optimized->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Time difference between START and CALIFA for the gammas");
h1_time_diff_gamma_11b = new TH1F(hist_name,hist_name,200,-100,100);
h1_time_diff_gamma_11b->GetXaxis()->SetTitle("Time difference [ns]");
h1_time_diff_gamma_11b->GetYaxis()->SetTitle("No. of Events");
h1_time_diff_gamma_11b->GetXaxis()->CenterTitle(true);
h1_time_diff_gamma_11b->GetYaxis()->CenterTitle(true);
h1_time_diff_gamma_11b->GetYaxis()->SetLabelSize(0.045);
h1_time_diff_gamma_11b->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Bindungsenergie B_p = T_0 - (T_1+T_2+T_A-1)");
h1_binding_energy_11b = new TH1F(hist_name,hist_name,1000,-500,500);
h1_binding_energy_11b->GetXaxis()->SetTitle("Calc. Binding energy [MeV]");
h1_binding_energy_11b->GetYaxis()->SetTitle("No. of Events");
h1_binding_energy_11b->GetXaxis()->CenterTitle(true);
h1_binding_energy_11b->GetYaxis()->CenterTitle(true);
h1_binding_energy_11b->GetYaxis()->SetLabelSize(0.045);
h1_binding_energy_11b->GetYaxis()->SetTitleSize(0.045);
//here I look at the triangle cut, see if these events have special angles, Energies,...


sprintf(hist_name, "E1 versus E2 for CALIFA for 11B in the cut triangle");
h2_E1_vs_E2_CALIFA_11b_cut = new TH2F(hist_name,hist_name,50,0,500,50,0,500);
h2_E1_vs_E2_CALIFA_11b_cut->GetXaxis()->SetTitle("E1 [MeV]");
h2_E1_vs_E2_CALIFA_11b_cut->GetYaxis()->SetTitle("E2 [MeV]");
h2_E1_vs_E2_CALIFA_11b_cut->GetXaxis()->CenterTitle(true);
h2_E1_vs_E2_CALIFA_11b_cut->GetYaxis()->CenterTitle(true);
h2_E1_vs_E2_CALIFA_11b_cut->GetYaxis()->SetLabelSize(0.045);
h2_E1_vs_E2_CALIFA_11b_cut->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Theta1 vs Theta2 for CALIFA 11B in the cut triangle");
h2_theta_1_vs_theta_2_CALIFA_11b_cut = new TH2F(hist_name,hist_name,75,0,150,75,0,150);
h2_theta_1_vs_theta_2_CALIFA_11b_cut->GetXaxis()->SetTitle("Theta1 [degrees]");
h2_theta_1_vs_theta_2_CALIFA_11b_cut->GetYaxis()->SetTitle("Theta2 [degrees]");
h2_theta_1_vs_theta_2_CALIFA_11b_cut->GetXaxis()->CenterTitle(true);
h2_theta_1_vs_theta_2_CALIFA_11b_cut->GetYaxis()->CenterTitle(true);
h2_theta_1_vs_theta_2_CALIFA_11b_cut->GetYaxis()->SetLabelSize(0.045);
h2_theta_1_vs_theta_2_CALIFA_11b_cut->GetYaxis()->SetTitleSize(0.045);

//-----------------------------------------------------------------------------------------

sprintf(hist_name, "E1+E2 versus full momentum 11B");
h2_E1_plus_E2_CALIFA_vs_full_mom11b = new TH2F(hist_name,hist_name,100,0,1000,200,9500,11500);
h2_E1_plus_E2_CALIFA_vs_full_mom11b->GetXaxis()->SetTitle("E1+E2 [MeV]");
h2_E1_plus_E2_CALIFA_vs_full_mom11b->GetYaxis()->SetTitle("Full Momentum 11B");
h2_E1_plus_E2_CALIFA_vs_full_mom11b->GetXaxis()->CenterTitle(true);
h2_E1_plus_E2_CALIFA_vs_full_mom11b->GetYaxis()->CenterTitle(true);
h2_E1_plus_E2_CALIFA_vs_full_mom11b->GetYaxis()->SetLabelSize(0.045);
h2_E1_plus_E2_CALIFA_vs_full_mom11b->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Long.  mom. p1p2 vs long. mom. 11B");
h2_long_mom_p1p2_long_mom11b = new TH2F(hist_name,hist_name,200,0,2000,200,9500,11500);
h2_long_mom_p1p2_long_mom11b->GetXaxis()->SetTitle("Longitudinal mom. p1+p2");
h2_long_mom_p1p2_long_mom11b->GetYaxis()->SetTitle("Longitudinal mom. 11B");
h2_long_mom_p1p2_long_mom11b->GetXaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b->GetYaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b->GetYaxis()->SetLabelSize(0.045);
h2_long_mom_p1p2_long_mom11b->GetYaxis()->SetTitleSize(0.045);

//here a cut on opening angle > < 84 degreees..
sprintf(hist_name, "Long.  mom. p1p2 vs long. mom. 11B for opening angles smaller than 84 degrees");
h2_long_mom_p1p2_long_mom11b_small_opening = new TH2F(hist_name,hist_name,200,0,2000,200,9500,11500);
h2_long_mom_p1p2_long_mom11b_small_opening->GetXaxis()->SetTitle("Longitudinal mom. p1+p2");
h2_long_mom_p1p2_long_mom11b_small_opening->GetYaxis()->SetTitle("Longitudinal mom. 11B");
h2_long_mom_p1p2_long_mom11b_small_opening->GetXaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_small_opening->GetYaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_small_opening->GetYaxis()->SetLabelSize(0.045);
h2_long_mom_p1p2_long_mom11b_small_opening->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Long.  mom. p1p2 vs long. mom. 11B for opening angles larger than 84 degrees");
h2_long_mom_p1p2_long_mom11b_large_opening = new TH2F(hist_name,hist_name,200,0,2000,200,9500,11500);
h2_long_mom_p1p2_long_mom11b_large_opening->GetXaxis()->SetTitle("Longitudinal mom. p1+p2");
h2_long_mom_p1p2_long_mom11b_large_opening->GetYaxis()->SetTitle("Longitudinal mom. 11B");
h2_long_mom_p1p2_long_mom11b_large_opening->GetXaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_large_opening->GetYaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_large_opening->GetYaxis()->SetLabelSize(0.045);
h2_long_mom_p1p2_long_mom11b_large_opening->GetYaxis()->SetTitleSize(0.045);

//constraints both on theta and phi...
sprintf(hist_name, "Long.  mom. p1p2 vs long. mom. 11B for opening angles smaller than 84 degrees and phi-diff +-30degrees");
h2_long_mom_p1p2_long_mom11b_theta_phi_constr = new TH2F(hist_name,hist_name,200,0,2000,200,9500,11500);
h2_long_mom_p1p2_long_mom11b_theta_phi_constr->GetXaxis()->SetTitle("Longitudinal mom. p1+p2");
h2_long_mom_p1p2_long_mom11b_theta_phi_constr->GetYaxis()->SetTitle("Longitudinal mom. 11B");
h2_long_mom_p1p2_long_mom11b_theta_phi_constr->GetXaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_theta_phi_constr->GetYaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_theta_phi_constr->GetYaxis()->SetLabelSize(0.045);
h2_long_mom_p1p2_long_mom11b_theta_phi_constr->GetYaxis()->SetTitleSize(0.045);

//making also constraint on missing mass..
sprintf(hist_name, "Long.  mom. p1p2 vs long. mom. 11B for opening angles smaller than 84 degrees and phi-diff +-30degrees and missing mom < 250");
h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom = new TH2F(hist_name,hist_name,200,0,2000,200,9500,11500);
h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom->GetXaxis()->SetTitle("Longitudinal mom. p1+p2");
h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom->GetYaxis()->SetTitle("Longitudinal mom. 11B");
h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom->GetXaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom->GetYaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom->GetYaxis()->SetLabelSize(0.045);
h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Long.  mom. p1p2 vs long. mom. 11B for opening angles smaller than 84 degrees and phi-diff +-30degrees and missing mom < 200");
h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom200 = new TH2F(hist_name,hist_name,200,0,2000,200,9500,11500);
h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom200->GetXaxis()->SetTitle("Longitudinal mom. p1+p2");
h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom200->GetYaxis()->SetTitle("Longitudinal mom. 11B");
h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom200->GetXaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom200->GetYaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom200->GetYaxis()->SetLabelSize(0.045);
h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom200->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Long.  mom. p1p2 vs long. mom. 11B for opening angles smaller than 84 degrees and phi-diff +-30degrees and missing mom < 150");
h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom150 = new TH2F(hist_name,hist_name,200,0,2000,200,9500,11500);
h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom150->GetXaxis()->SetTitle("Longitudinal mom. p1+p2");
h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom150->GetYaxis()->SetTitle("Longitudinal mom. 11B");
h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom150->GetXaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom150->GetYaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom150->GetYaxis()->SetLabelSize(0.045);
h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom150->GetYaxis()->SetTitleSize(0.045);

//gamma spectrum for according missing mass...


sprintf(hist_name, "Highest Gamma Energy E for 11B with missing mom < 250");
h1_m_mom_250_gamma_spec = new TH1F(hist_name,hist_name,200,0,10);
h1_m_mom_250_gamma_spec->GetXaxis()->SetTitle("Energy E [MeV]");
h1_m_mom_250_gamma_spec->GetYaxis()->SetTitle("Number of entries");
h1_m_mom_250_gamma_spec->GetXaxis()->CenterTitle(true);
h1_m_mom_250_gamma_spec->GetYaxis()->CenterTitle(true);
h1_m_mom_250_gamma_spec->GetYaxis()->SetLabelSize(0.045);
h1_m_mom_250_gamma_spec->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Highest Gamma Energy E for 11B with missing mom < 200");
h1_m_mom_200_gamma_spec = new TH1F(hist_name,hist_name,200,0,10);
h1_m_mom_200_gamma_spec->GetXaxis()->SetTitle("Energy E [MeV]");
h1_m_mom_200_gamma_spec->GetYaxis()->SetTitle("Number of entries");
h1_m_mom_200_gamma_spec->GetXaxis()->CenterTitle(true);
h1_m_mom_200_gamma_spec->GetYaxis()->CenterTitle(true);
h1_m_mom_200_gamma_spec->GetYaxis()->SetLabelSize(0.045);
h1_m_mom_200_gamma_spec->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Highest Gamma Energy E for 11B with missing mom < 150");
h1_m_mom_150_gamma_spec = new TH1F(hist_name,hist_name,200,0,10);
h1_m_mom_150_gamma_spec->GetXaxis()->SetTitle("Energy E [MeV]");
h1_m_mom_150_gamma_spec->GetYaxis()->SetTitle("Number of entries");
h1_m_mom_150_gamma_spec->GetXaxis()->CenterTitle(true);
h1_m_mom_150_gamma_spec->GetYaxis()->CenterTitle(true);
h1_m_mom_150_gamma_spec->GetYaxis()->SetLabelSize(0.045);
h1_m_mom_150_gamma_spec->GetYaxis()->SetTitleSize(0.045);

//cluster checking...
sprintf(hist_name, "Number of gamma clusters for 11B with missing mom < 250");
h1_cluster_gamma_11b = new TH1F(hist_name,hist_name,10,0,10);
h1_cluster_gamma_11b->GetXaxis()->SetTitle("Number of Clusters");
h1_cluster_gamma_11b->GetYaxis()->SetTitle("Counts");
h1_cluster_gamma_11b->GetXaxis()->CenterTitle(true);
h1_cluster_gamma_11b->GetYaxis()->CenterTitle(true);
h1_cluster_gamma_11b->GetYaxis()->SetLabelSize(0.045);
h1_cluster_gamma_11b->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Number of proton clusters for 11B with missing mom < 250");
h1_cluster_proton_11b = new TH1F(hist_name,hist_name,10,0,10);
h1_cluster_proton_11b->GetXaxis()->SetTitle("Number of Clusters");
h1_cluster_proton_11b->GetYaxis()->SetTitle("Counts");
h1_cluster_proton_11b->GetXaxis()->CenterTitle(true);
h1_cluster_proton_11b->GetYaxis()->CenterTitle(true);
h1_cluster_proton_11b->GetYaxis()->SetLabelSize(0.045);
h1_cluster_proton_11b->GetYaxis()->SetTitleSize(0.045);


//now just constraining on gamma energy larger 1MeV...
sprintf(hist_name, "Long.  mom. p1p2 vs long. mom. 11B for E_gamma > 1MeV");
h2_long_mom_p1p2_long_mom11b_large_gamma = new TH2F(hist_name,hist_name,200,0,2000,200,9500,11500);
h2_long_mom_p1p2_long_mom11b_large_gamma->GetXaxis()->SetTitle("Longitudinal mom. p1+p2");
h2_long_mom_p1p2_long_mom11b_large_gamma->GetYaxis()->SetTitle("Longitudinal mom. 11B");
h2_long_mom_p1p2_long_mom11b_large_gamma->GetXaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_large_gamma->GetYaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_large_gamma->GetYaxis()->SetLabelSize(0.045);
h2_long_mom_p1p2_long_mom11b_large_gamma->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Long.  mom. p1p2 vs long. mom. 11B cos_high");
h2_long_mom_p1p2_long_mom11b_high = new TH2F(hist_name,hist_name,200,0,2000,200,9500,11500);
h2_long_mom_p1p2_long_mom11b_high->GetXaxis()->SetTitle("Longitudinal mom. p1+p2");
h2_long_mom_p1p2_long_mom11b_high->GetYaxis()->SetTitle("Longitudinal mom. 11B");
h2_long_mom_p1p2_long_mom11b_high->GetXaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_high->GetYaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_high->GetYaxis()->SetLabelSize(0.045);
h2_long_mom_p1p2_long_mom11b_high->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Long.  mom. p1p2 vs long. mom. 11B cos_low");
h2_long_mom_p1p2_long_mom11b_low = new TH2F(hist_name,hist_name,200,0,2000,200,9500,11500);
h2_long_mom_p1p2_long_mom11b_low->GetXaxis()->SetTitle("Longitudinal mom. p1+p2");
h2_long_mom_p1p2_long_mom11b_low->GetYaxis()->SetTitle("Longitudinal mom. 11B");
h2_long_mom_p1p2_long_mom11b_low->GetXaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_low->GetYaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_low->GetYaxis()->SetLabelSize(0.045);
h2_long_mom_p1p2_long_mom11b_low->GetYaxis()->SetTitleSize(0.045);

//do p1+p2 long.mom. vs 11b long. mom. plots with optimizations.... and restrictions on phi and theta....
sprintf(hist_name, "Long.  mom. p1p2 vs long. mom. 11B cos_low with optimized");
h2_long_mom_p1p2_long_mom11b_low_optimized = new TH2F(hist_name,hist_name,200,0,2000,200,9500,11500);
h2_long_mom_p1p2_long_mom11b_low_optimized->GetXaxis()->SetTitle("Longitudinal mom. p1+p2");
h2_long_mom_p1p2_long_mom11b_low_optimized->GetYaxis()->SetTitle("Longitudinal mom. 11B");
h2_long_mom_p1p2_long_mom11b_low_optimized->GetXaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_low_optimized->GetYaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_low_optimized->GetYaxis()->SetLabelSize(0.045);
h2_long_mom_p1p2_long_mom11b_low_optimized->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Long.  mom. p1p2 vs long. mom. 11B with cos < -0.8");
h2_long_mom_p1p2_long_mom11b_low_optimized_08 = new TH2F(hist_name,hist_name,200,0,2000,200,9500,11500);
h2_long_mom_p1p2_long_mom11b_low_optimized_08->GetXaxis()->SetTitle("Longitudinal mom. p1+p2");
h2_long_mom_p1p2_long_mom11b_low_optimized_08->GetYaxis()->SetTitle("Longitudinal mom. 11B");
h2_long_mom_p1p2_long_mom11b_low_optimized_08->GetXaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_low_optimized_08->GetYaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_low_optimized_08->GetYaxis()->SetLabelSize(0.045);
h2_long_mom_p1p2_long_mom11b_low_optimized_08->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Long.  mom. p1p2 vs long. mom. 11B cos_high with optimized");
h2_long_mom_p1p2_long_mom11b_high_optimized = new TH2F(hist_name,hist_name,200,0,2000,200,9500,11500);
h2_long_mom_p1p2_long_mom11b_high_optimized->GetXaxis()->SetTitle("Longitudinal mom. p1+p2");
h2_long_mom_p1p2_long_mom11b_high_optimized->GetYaxis()->SetTitle("Longitudinal mom. 11B");
h2_long_mom_p1p2_long_mom11b_high_optimized->GetXaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_high_optimized->GetYaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_high_optimized->GetYaxis()->SetLabelSize(0.045);
h2_long_mom_p1p2_long_mom11b_high_optimized->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "missing four momentum p_missing 11B");
h1_p_missing_11B = new TH1F(hist_name,hist_name,750,0,1500);
h1_p_missing_11B->GetXaxis()->SetTitle("p_missing [MeV/c]");
h1_p_missing_11B->GetYaxis()->SetTitle("No. of Events");
h1_p_missing_11B->GetXaxis()->CenterTitle(true);
h1_p_missing_11B->GetYaxis()->CenterTitle(true);
h1_p_missing_11B->GetYaxis()->SetLabelSize(0.045);
h1_p_missing_11B->GetYaxis()->SetTitleSize(0.045);


//------------------------------------------------------------------------------------------------

sprintf(hist_name, "Long.  mom. p1p2 vs long. mom. 11B no borders");
h2_long_mom_p1p2_long_mom11b_no_border = new TH2F(hist_name,hist_name,200,0,2000,200,9500,11500);
h2_long_mom_p1p2_long_mom11b_no_border->GetXaxis()->SetTitle("Longitudinal mom. p1+p2");
h2_long_mom_p1p2_long_mom11b_no_border->GetYaxis()->SetTitle("Longitudinal mom. 11B");
h2_long_mom_p1p2_long_mom11b_no_border->GetXaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_no_border->GetYaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_no_border->GetYaxis()->SetLabelSize(0.045);
h2_long_mom_p1p2_long_mom11b_no_border->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Long.  mom. p1p2 vs long. mom. 11B cos_high no borders");
h2_long_mom_p1p2_long_mom11b_no_border_high = new TH2F(hist_name,hist_name,200,0,2000,200,9500,11500);
h2_long_mom_p1p2_long_mom11b_no_border_high->GetXaxis()->SetTitle("Longitudinal mom. p1+p2");
h2_long_mom_p1p2_long_mom11b_no_border_high->GetYaxis()->SetTitle("Longitudinal mom. 11B");
h2_long_mom_p1p2_long_mom11b_no_border_high->GetXaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_no_border_high->GetYaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_no_border_high->GetYaxis()->SetLabelSize(0.045);
h2_long_mom_p1p2_long_mom11b_no_border_high->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Long.  mom. p1p2 vs long. mom. 11B cos_low no borders");
h2_long_mom_p1p2_long_mom11b_no_border_low = new TH2F(hist_name,hist_name,200,0,2000,200,9500,11500);
h2_long_mom_p1p2_long_mom11b_no_border_low->GetXaxis()->SetTitle("Longitudinal mom. p1+p2");
h2_long_mom_p1p2_long_mom11b_no_border_low->GetYaxis()->SetTitle("Longitudinal mom. 11B");
h2_long_mom_p1p2_long_mom11b_no_border_low->GetXaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_no_border_low->GetYaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b_no_border_low->GetYaxis()->SetLabelSize(0.045);
h2_long_mom_p1p2_long_mom11b_no_border_low->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Long.  mom. p1p2 vs long. mom. 10B");
h2_long_mom_p1p2_long_mom10b = new TH2F(hist_name,hist_name,200,0,2000,200,8500,10500);
h2_long_mom_p1p2_long_mom10b->GetXaxis()->SetTitle("Longitudinal mom. p1+p2");
h2_long_mom_p1p2_long_mom10b->GetYaxis()->SetTitle("Longitudinal mom. 10B");
h2_long_mom_p1p2_long_mom10b->GetXaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom10b->GetYaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom10b->GetYaxis()->SetLabelSize(0.045);
h2_long_mom_p1p2_long_mom10b->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "B*rho/(beta*fixed_gamma) versus Time of Flight 11B");
h2_tof_vs_aq_fix_g_11b = new TH2F(hist_name,hist_name,900,1,3,1000,30,40);
h2_tof_vs_aq_fix_g_11b->GetXaxis()->SetTitle("B*rho/(beta*fixed_gamma)");
h2_tof_vs_aq_fix_g_11b->GetYaxis()->SetTitle("Time of Flight [ns]");
h2_tof_vs_aq_fix_g_11b->GetXaxis()->CenterTitle(true);
h2_tof_vs_aq_fix_g_11b->GetYaxis()->CenterTitle(true);
h2_tof_vs_aq_fix_g_11b->GetYaxis()->SetLabelSize(0.045);
h2_tof_vs_aq_fix_g_11b->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "B*rho/(beta*fixed_gamma) versus Time of Flight 10B");
h2_tof_vs_aq_fix_g_10b = new TH2F(hist_name,hist_name,900,1,3,1000,30,40);
h2_tof_vs_aq_fix_g_10b->GetXaxis()->SetTitle("B*rho/(beta*fixed_gamma)");
h2_tof_vs_aq_fix_g_10b->GetYaxis()->SetTitle("Time of Flight [ns]");
h2_tof_vs_aq_fix_g_10b->GetXaxis()->CenterTitle(true);
h2_tof_vs_aq_fix_g_10b->GetYaxis()->CenterTitle(true);
h2_tof_vs_aq_fix_g_10b->GetYaxis()->SetLabelSize(0.045);
h2_tof_vs_aq_fix_g_10b->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "B*rho/(beta*fixed_gamma) versus Time of Flight 11C");
h2_tof_vs_aq_fix_g_11c = new TH2F(hist_name,hist_name,900,1,3,1000,30,40);
h2_tof_vs_aq_fix_g_11c->GetXaxis()->SetTitle("B*rho/(beta*fixed_gamma)");
h2_tof_vs_aq_fix_g_11c->GetYaxis()->SetTitle("Time of Flight [ns]");
h2_tof_vs_aq_fix_g_11c->GetXaxis()->CenterTitle(true);
h2_tof_vs_aq_fix_g_11c->GetYaxis()->CenterTitle(true);
h2_tof_vs_aq_fix_g_11c->GetYaxis()->SetLabelSize(0.045);
h2_tof_vs_aq_fix_g_11c->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "B*rho/(beta*fixed_gamma) versus Time of Flight 12C");
h2_tof_vs_aq_fix_g_12c = new TH2F(hist_name,hist_name,900,1,3,1000,30,40);
h2_tof_vs_aq_fix_g_12c->GetXaxis()->SetTitle("B*rho/(beta*fixed_gamma)");
h2_tof_vs_aq_fix_g_12c->GetYaxis()->SetTitle("Time of Flight [ns]");
h2_tof_vs_aq_fix_g_12c->GetXaxis()->CenterTitle(true);
h2_tof_vs_aq_fix_g_12c->GetYaxis()->CenterTitle(true);
h2_tof_vs_aq_fix_g_12c->GetYaxis()->SetLabelSize(0.045);
h2_tof_vs_aq_fix_g_12c->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Abs. Angle Difference before and after GLAD in y-Plane without cut triangle");
h1_angl_diff_y = new TH1F(hist_name,hist_name,100,0,25);
h1_angl_diff_y->GetXaxis()->SetTitle("Angle differnce (abs.) [degrees]");
h1_angl_diff_y->GetYaxis()->SetTitle("No. of Events");
h1_angl_diff_y->GetXaxis()->CenterTitle(true);
h1_angl_diff_y->GetYaxis()->CenterTitle(true);
h1_angl_diff_y->GetYaxis()->SetLabelSize(0.045);
h1_angl_diff_y->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Abs. Angle Difference before and after GLAD in y-Plane for cut triangle");
h1_angl_diff_y_cut_triangle = new TH1F(hist_name,hist_name,100,0,25);
h1_angl_diff_y_cut_triangle->GetXaxis()->SetTitle("Angle differnce (abs.) [degrees]");
h1_angl_diff_y_cut_triangle->GetYaxis()->SetTitle("No. of Events");
h1_angl_diff_y_cut_triangle->GetXaxis()->CenterTitle(true);
h1_angl_diff_y_cut_triangle->GetYaxis()->CenterTitle(true);
h1_angl_diff_y_cut_triangle->GetYaxis()->SetLabelSize(0.045);
h1_angl_diff_y_cut_triangle->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Gamma spectrum for 12C(p,2p)11B reaction  for cos < -0.6");
h1_cos_low_gamma = new TH1F(hist_name,hist_name,200,0,10);
h1_cos_low_gamma->GetXaxis()->SetTitle("Energy [MeV]");
h1_cos_low_gamma->GetYaxis()->SetTitle("No. of Events");
h1_cos_low_gamma->GetXaxis()->CenterTitle(true);
h1_cos_low_gamma->GetYaxis()->CenterTitle(true);
h1_cos_low_gamma->GetYaxis()->SetLabelSize(0.045);
h1_cos_low_gamma->GetYaxis()->SetTitleSize(0.045);

}
else {

sprintf(hist_name, "E1+E2 versus full momentum 11B");
h2_E1_plus_E2_CALIFA_vs_full_mom11b = new TH2F(hist_name,hist_name,100,0,1000,300,15000,18000);
h2_E1_plus_E2_CALIFA_vs_full_mom11b->GetXaxis()->SetTitle("E1+E2 [MeV]");
h2_E1_plus_E2_CALIFA_vs_full_mom11b->GetYaxis()->SetTitle("Full Momentum 11B");
h2_E1_plus_E2_CALIFA_vs_full_mom11b->GetXaxis()->CenterTitle(true);
h2_E1_plus_E2_CALIFA_vs_full_mom11b->GetYaxis()->CenterTitle(true);
h2_E1_plus_E2_CALIFA_vs_full_mom11b->GetYaxis()->SetLabelSize(0.045);
h2_E1_plus_E2_CALIFA_vs_full_mom11b->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Long.  mom. p1p2 vs long. mom. 11B");
h2_long_mom_p1p2_long_mom11b = new TH2F(hist_name,hist_name,200,0,2000,300,15000,18000);
h2_long_mom_p1p2_long_mom11b->GetXaxis()->SetTitle("Longitudinal mom. p1+p2");
h2_long_mom_p1p2_long_mom11b->GetYaxis()->SetTitle("Longitudinal mom. 11B");
h2_long_mom_p1p2_long_mom11b->GetXaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b->GetYaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom11b->GetYaxis()->SetLabelSize(0.045);
h2_long_mom_p1p2_long_mom11b->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Long.  mom. p1p2 vs long. mom. 10B");
h2_long_mom_p1p2_long_mom10b = new TH2F(hist_name,hist_name,200,0,2000,300,13500,16500);
h2_long_mom_p1p2_long_mom10b->GetXaxis()->SetTitle("Longitudinal mom. p1+p2");
h2_long_mom_p1p2_long_mom10b->GetYaxis()->SetTitle("Longitudinal mom. 10B");
h2_long_mom_p1p2_long_mom10b->GetXaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom10b->GetYaxis()->CenterTitle(true);
h2_long_mom_p1p2_long_mom10b->GetYaxis()->SetLabelSize(0.045);
h2_long_mom_p1p2_long_mom10b->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "B*rho/(beta*fixed_gamma) versus Time of Flight 11B");
h2_tof_vs_aq_fix_g_11b = new TH2F(hist_name,hist_name,900,1,3,1000,23,33);
h2_tof_vs_aq_fix_g_11b->GetXaxis()->SetTitle("B*rho/(beta*fixed_gamma)");
h2_tof_vs_aq_fix_g_11b->GetYaxis()->SetTitle("Time of Flight [ns]");
h2_tof_vs_aq_fix_g_11b->GetXaxis()->CenterTitle(true);
h2_tof_vs_aq_fix_g_11b->GetYaxis()->CenterTitle(true);
h2_tof_vs_aq_fix_g_11b->GetYaxis()->SetLabelSize(0.045);
h2_tof_vs_aq_fix_g_11b->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "B*rho/(beta*fixed_gamma) versus Time of Flight 10B");
h2_tof_vs_aq_fix_g_10b = new TH2F(hist_name,hist_name,900,1,3,1000,23,33);
h2_tof_vs_aq_fix_g_10b->GetXaxis()->SetTitle("B*rho/(beta*fixed_gamma)");
h2_tof_vs_aq_fix_g_10b->GetYaxis()->SetTitle("Time of Flight [ns]");
h2_tof_vs_aq_fix_g_10b->GetXaxis()->CenterTitle(true);
h2_tof_vs_aq_fix_g_10b->GetYaxis()->CenterTitle(true);
h2_tof_vs_aq_fix_g_10b->GetYaxis()->SetLabelSize(0.045);
h2_tof_vs_aq_fix_g_10b->GetYaxis()->SetTitleSize(0.045);

//sprintf(hist_name, "B*rho/(beta*fixed_gamma) versus Time of Flight 11C");
//h2_tof_vs_aq_fix_g_11c = new TH2F(hist_name,hist_name,900,1,3,1000,23,33);
//h2_tof_vs_aq_fix_g_11c->GetXaxis()->SetTitle("B*rho/(beta*fixed_gamma)");
//h2_tof_vs_aq_fix_g_11c->GetYaxis()->SetTitle("Time of Flight [ns]");
//h2_tof_vs_aq_fix_g_11c->GetXaxis()->CenterTitle(true);
//h2_tof_vs_aq_fix_g_11c->GetYaxis()->CenterTitle(true);
//h2_tof_vs_aq_fix_g_11c->GetYaxis()->SetLabelSize(0.045);
//h2_tof_vs_aq_fix_g_11c->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "B*rho/(beta*fixed_gamma) versus Time of Flight 12C");
h2_tof_vs_aq_fix_g_12c = new TH2F(hist_name,hist_name,900,1,3,1000,23,33);
h2_tof_vs_aq_fix_g_12c->GetXaxis()->SetTitle("B*rho/(beta*fixed_gamma)");
h2_tof_vs_aq_fix_g_12c->GetYaxis()->SetTitle("Time of Flight [ns]");
h2_tof_vs_aq_fix_g_12c->GetXaxis()->CenterTitle(true);
h2_tof_vs_aq_fix_g_12c->GetYaxis()->CenterTitle(true);
h2_tof_vs_aq_fix_g_12c->GetYaxis()->SetLabelSize(0.045);
h2_tof_vs_aq_fix_g_12c->GetYaxis()->SetTitleSize(0.045);

}

sprintf(hist_name, "Transv. mom. difference p1 p2 for 10B");
h1_transv_mom_difference_p1_p2_10b = new TH1F(hist_name,hist_name,100,0,500);
h1_transv_mom_difference_p1_p2_10b->GetXaxis()->SetTitle("Transv. mom. difference p1p2 [MeV/c]");
h1_transv_mom_difference_p1_p2_10b->GetYaxis()->SetTitle("Number of entries");
h1_transv_mom_difference_p1_p2_10b->GetXaxis()->CenterTitle(true);
h1_transv_mom_difference_p1_p2_10b->GetYaxis()->CenterTitle(true);
h1_transv_mom_difference_p1_p2_10b->GetYaxis()->SetLabelSize(0.045);
h1_transv_mom_difference_p1_p2_10b->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Transv. mom. difference p1 p2 for 11B");
h1_transv_mom_difference_p1_p2_11b = new TH1F(hist_name,hist_name,100,0,500);
h1_transv_mom_difference_p1_p2_11b->GetXaxis()->SetTitle("Transv. mom. difference p1p2 [MeV/c]");
h1_transv_mom_difference_p1_p2_11b->GetYaxis()->SetTitle("Number of entries");
h1_transv_mom_difference_p1_p2_11b->GetXaxis()->CenterTitle(true);
h1_transv_mom_difference_p1_p2_11b->GetYaxis()->CenterTitle(true);
h1_transv_mom_difference_p1_p2_11b->GetYaxis()->SetLabelSize(0.045);
h1_transv_mom_difference_p1_p2_11b->GetYaxis()->SetTitleSize(0.045);

//END OF IMPLEMENTATION OF HISTOS--------------------------------------------------


sprintf(f_out_name,"/home/ge37liw/plots/histos/zero_with_miss/fast_macro/main_macro_phi_test_automatic_src_califa_%s.root", count_i);
TFile * f = new TFile(f_out_name,"RECREATE");
//Path to input file
sprintf(fname,"/scratch8/ge37liw/workingspace/data/root_files/all_ts_unpack/sweep_target/ts_cone_cluster_%s_zero_cm.root",count_i);
sprintf(f_parameters,"/home/ge37liw/plots/histos/zero_with_miss/par_files_dir/par_run_%s.txt",count_i);
double t = find_z(f_parameters,fname);
cout << "this is the return value\t" << t << endl;
vector<double> s = find_isotopes(f_parameters,fname,t);
cout << "the isotope cuts are:\t" <<endl;\
cout << "10B:\t " << s[0] << endl;\
cout << "11B:\t" << s[1] << endl;\
cout << "11C:\t" << s[2] << endl;\
cout << "12C:\t" << s[3] << endl;

vector<double> u = find_psi_par(f_parameters,fname,t,s);
cout << "job done, now we can proceede with analysis" << endl;

//Finally the precise calculation with the parameters from Step1 -> Step3
cout << "finally precise calculation!---------------------------------" << endl;

TChain* chain = new TChain("evt");
chain->Reset();
chain->Add(fname);

//TClonesArrays containing the TObjects R3BSofToFWTcalData,... which allow access to data over function calls
//
TClonesArray* SofMwpc3HitData = new TClonesArray("R3BSofMwpcHitData",5);
R3BSofMwpcHitData** sofmwpc3hitdata;
TBranch *branchSofMwpc3HitData = chain->GetBranch("Mwpc3HitData");
branchSofMwpc3HitData->SetAddress(&SofMwpc3HitData);
//
TClonesArray* SofSciTcalData = new TClonesArray("R3BSofSciTcalData",2);
R3BSofSciTcalData** sofscitcaldata;
TBranch *branchSofSciTcalData = chain->GetBranch("SofSciTcalData");
branchSofSciTcalData->SetAddress(&SofSciTcalData);
//
TClonesArray* SofToFWTcalData = new TClonesArray("R3BSofToFWTcalData",2);
R3BSofToFWTcalData** softofwtcaldata;
TBranch *branchSofToFWTcalData = chain->GetBranch("SofToFWTcalData");
branchSofToFWTcalData->SetAddress(&SofToFWTcalData);
//
TClonesArray* SofMwpc0HitData = new TClonesArray("R3BSofMwpcHitData",5);
R3BSofMwpcHitData** sofmwpc0hitdata;
TBranch *branchSofMwpc0HitData = chain->GetBranch("Mwpc0HitData");
branchSofMwpc0HitData->SetAddress(&SofMwpc0HitData);
//
TClonesArray* SofMwpc1HitData = new TClonesArray("R3BSofMwpcHitData",5);
R3BSofMwpcHitData** sofmwpc1hitdata;
TBranch *branchSofMwpc1HitData = chain->GetBranch("Mwpc1HitData");
branchSofMwpc1HitData->SetAddress(&SofMwpc1HitData);
//
TClonesArray* SofMwpc2HitData = new TClonesArray("R3BSofMwpcHitData",5);
R3BSofMwpcHitData** sofmwpc2hitdata;
TBranch *branchSofMwpc2HitData = chain->GetBranch("Mwpc2HitData");
branchSofMwpc2HitData->SetAddress(&SofMwpc2HitData);
//
TClonesArray* SofToFWMappedData = new TClonesArray("R3BSofToFWMappedData",2);
R3BSofToFWMappedData** softofwmappeddata;
TBranch *branchSofToFWMappedData = chain->GetBranch("SofToFWMappedData");
branchSofToFWMappedData->SetAddress(&SofToFWMappedData);
//
TClonesArray* SofTwimHitData = new TClonesArray("R3BSofTwimHitData",2);
R3BSofTwimHitData** softwimhitdata;
TBranch *branchSofTwimHitData = chain->GetBranch("TwimHitData");
branchSofTwimHitData->SetAddress(&SofTwimHitData);
//
TClonesArray* CalifaHitData = new TClonesArray("R3BCalifaHitData",3);
R3BCalifaHitData** califahitdata;
TBranch *branchCalifaHitData = chain->GetBranch("CalifaHitData");
branchCalifaHitData->SetAddress(&CalifaHitData);
//
Long64_t nevents = chain->GetEntries(); //number of events in file with name "fname"
for(Long64_t i=0;i< nevents;i++){
    Long64_t evtnr = i;

    if (i%100000==0)
        cout<<"Processing event "<<i<<endl;

	chain->GetEvent(i); //event_i_call
	entries_mw0 = SofMwpc0HitData->GetEntries();
    entries_mw1 = SofMwpc1HitData->GetEntries();
    entries_mw2 = SofMwpc2HitData->GetEntries();
    entries_mw3 = SofMwpc3HitData->GetEntries();
    entries_start = SofSciTcalData->GetEntriesFast();
    entries_tofw = SofToFWTcalData->GetEntriesFast();
    entries_califa = CalifaHitData->GetEntries();
    softwimhitdata = new R3BSofTwimHitData*[1];
    softwimhitdata[0] = (R3BSofTwimHitData*)SofTwimHitData->At(0);
    if (entries_start ==2 && entries_tofw == 2 && entries_mw0 == 1 && entries_mw1 == 1 && entries_mw2 == 1 && entries_mw3 == 1 && softwimhitdata[0]){
        charge_val = softwimhitdata[0]->GetZcharge();
        sofmwpc0hitdata = new R3BSofMwpcHitData*[1];
        sofmwpc1hitdata = new R3BSofMwpcHitData*[1];
        sofmwpc2hitdata = new R3BSofMwpcHitData*[1];
        sofmwpc3hitdata = new R3BSofMwpcHitData*[1];
        softofwtcaldata = new R3BSofToFWTcalData*[1];
        softofwmappeddata = new R3BSofToFWMappedData*[2];
        califahitdata = new R3BCalifaHitData*[entries_califa];
        sofmwpc0hitdata[0] = (R3BSofMwpcHitData*)SofMwpc0HitData->At(0);
        sofmwpc1hitdata[0] = (R3BSofMwpcHitData*)SofMwpc1HitData->At(0);
        sofmwpc2hitdata[0] = (R3BSofMwpcHitData*)SofMwpc2HitData->At(0);
        sofmwpc3hitdata[0] = (R3BSofMwpcHitData*)SofMwpc3HitData->At(0);
        softofwtcaldata[0] = (R3BSofToFWTcalData*)SofToFWTcalData->At(0);
        Double_t xMW0 = sofmwpc0hitdata[0]->GetX();
        Double_t xMW1 = sofmwpc1hitdata[0]->GetX();
        Double_t xMW2 = sofmwpc2hitdata[0]->GetX();
        Double_t xMW3 = sofmwpc3hitdata[0]->GetX();
        Double_t yMW1 = sofmwpc1hitdata[0]->GetY();
        Double_t yMW2 = sofmwpc2hitdata[0]->GetY();
        Double_t yMW3 = sofmwpc3hitdata[0]->GetY();
if (xMW1 != -1000 && xMW2 != -1000 && xMW3 != -1000  && xMW0 != -1000){
        Double_t psi_in = atan(((xMW2+xMW2_shift)-(xMW1+xMW1_shift))/(zM2-zM1));
        Double_t T_to_M2 = (zM2-zT)/(cos(psi_in));
        Double_t b1_cutoff = (xMW1+xMW1_shift)-tan(psi_in)*zM1;
        Double_t z_labMW3 = middle_zM3+cos(2*PI/5.)*(xMW3-xMW3_shift);
        Double_t x_labMW3 = middle_xM3+sin(2*PI/5.)*(xMW3-xMW3_shift);
        Double_t zB = (bGLAD_cutoff-b1_cutoff)/(tan(psi_in)-tan(PI/2.-alpha_G));
        Double_t xB = tan(psi_in)*zB+b1_cutoff;
        for (Int_t rec = 0; rec < 50; rec++){
            zC = (zGm_cutoff-b1_cutoff)/(tan(psi_in)-tan(PI/2.-alpha_G)) + 2*rec;
            xC = tan(psi_in)*zC+b1_cutoff;
            slope = (x_labMW3-xC)/(z_labMW3-zC);
            offset_slope = x_labMW3-slope*z_labMW3;
            psi_out_rec[rec] = -atan(slope);
            z_D = (D_cutoff-offset_slope)/(slope-tan(PI/2.-alpha_G));
            x_D = tan(PI/2.-alpha_G)*z_D +D_cutoff;
            l_diff[rec] = (sqrt((xC-xB)*(xC-xB)+(zC-zB)*(zC-zB))-sqrt((x_D-xC)*(x_D-xC)+(z_D-zC)*(z_D-zC)));
            z_pos_shift[rec] = zC;
            x_pos_shift[rec] = xC;
        }
        Int_t n = 50;
        TGraph *gr1 = new TGraph(n,psi_out_rec,l_diff);
        gr1->Fit("pol1","Q");
        TF1 *f3 =gr1->GetFunction("pol1");
        Double_t f3_slope = f3->GetParameter(1);
        Double_t f3_offset = f3->GetParameter(0);
        Double_t psi_out = -f3_offset/f3_slope;
        TGraph *gr2 = new TGraph(n,z_pos_shift,l_diff);
        gr2->Fit("pol1", "Q");
        TF1 *f4 = gr2->GetFunction("pol1");
        Double_t f4_slope = f4->GetParameter(1);
        Double_t f4_offset = f4->GetParameter(0);
        Double_t new_z_C = -f4_offset/f4_slope;

        TGraph *gr3 = new TGraph(n,x_pos_shift,l_diff);
        gr3->Fit("pol1","Q");
        TF1 *f5 = gr3->GetFunction("pol1");
        Double_t f5_slope = f5->GetParameter(1);
        Double_t f5_offset = f5->GetParameter(0);
        Double_t new_x_C = -f5_offset/f5_slope;
        Double_t z_diff = new_z_C-z_pos_shift[0];
        Double_t offset_new = new_x_C+tan(psi_out)*new_z_C;
        z_D = (D_cutoff-offset_new)/(-tan(psi_out)-tan(PI/2.-alpha_G));
        x_D = -tan(psi_out)*z_D + offset_new;
        Double_t zToFW = (cutoff_ToFW-offset_new)/(-tan(psi_out)-tan(2*PI/5.));
        Double_t xToFW = -tan(psi_out)*zToFW+offset_new;
        Double_t path_D_to_TOFW = (zToFW-z_D)/(cos(psi_out));
        Double_t BD = sqrt((x_D-xB)*(x_D-xB)+(z_D-zB)*(z_D-zB));

        Double_t m = (cos(psi_out)-cos(psi_in))/(sin(psi_out)+sin(psi_in));
        Double_t rho =  (Leff/(2*sin((psi_in+psi_out)/2)))*sqrt(pow(((m+tan(alpha_G))/(1-m*tan(alpha_G))),2)+1);
        Double_t w = 2*abs(asin(BD/(2*rho)));


        tot_length = start_to_target+ abs((zB-zT)/(cos(psi_in)))+rho*w+abs((zToFW-z_D)/cos(psi_out));
        delete gr1;
        delete gr2;
        delete gr3;

        if (SofToFWTcalData && SofToFWTcalData->GetEntriesFast())
            {
            for (UShort_t i = 0; i < NbDets; i++)
                {
                for (UShort_t j = 0; j < NbChs; j++)
                    {
                    mult[i * NbChs + j] = 0;
                    }
                }
            for (Int_t ihit = 0; ihit < entries_tofw; ihit++)
                {
                R3BSofToFWMappedData* hitmapped = (R3BSofToFWMappedData*)SofToFWMappedData->At(ihit);
                if (!hitmapped)
                    continue;
                iDet = hitmapped->GetDetector() - 1;
                iCh = hitmapped->GetPmt() - 1;
                mult[iDet * NbChs + iCh]++;
                }
            for (Int_t ihit = 0; ihit < entries_tofw; ihit++)
                {
                R3BSofToFWTcalData* hittcal = (R3BSofToFWTcalData*)SofToFWTcalData->At(ihit);
                if (!hittcal)
                    continue;
                if (hittcal->GetPmt() == 3)
                    continue;
                iDet = hittcal->GetDetector() - 1;
                iCh = hittcal->GetPmt() - 1;
                iRawTimeNs[iDet * 2 + iCh] = hittcal->GetRawTimeNs();
                }
            for (Int_t s_hit = 0 ; s_hit < entries_start; s_hit++)
                {
                R3BSofSciTcalData* start_hits = (R3BSofSciTcalData*)SofSciTcalData->At(s_hit);
                if(!start_hits)
                    continue;
                iRawTimeStartNs[s_hit] = start_hits->GetRawTimeNs();

                }
                sofmwpc3hitdata[0]=(R3BSofMwpcHitData*)SofMwpc3HitData->At(0);
            if ((sofmwpc3hitdata[0]->GetY() > -500.) && (sofmwpc3hitdata[0]->GetX() > -500.))  //cut on the mwpc3
            {
            Double_t raw_mwpc3 = sofmwpc3hitdata[0]->GetY();
            Double_t tofpos_ns = -1000;
            Double_t raw_tofpos_ns = -1000;
            Double_t raw_t_start = -1000000.;
            Double_t raw_time_of_flight = 0.;
            for (UShort_t i = 0; i < NbDets; i++)
                {
                if ((mult[i * NbChs] == 1) && (mult[i * NbChs + 1] == 1))
                    {
                    tofpos_ns = 0.5*(iRawTimeNs[i * NbChs + 1] - iRawTimeNs[i * NbChs]-offsets[i+1]);
                    raw_tofpos_ns = 0.5*(iRawTimeNs[i * NbChs + 1] - iRawTimeNs[i * NbChs]);
                    raw_t_start = 0.5*(iRawTimeStartNs[0]-iRawTimeStartNs[1]);
                    if (raw_t_start !=-1000000. && raw_tofpos_ns != -1000 && tofpos_ns != -1000)
                        {
                        Double_t tof_diff_up_down = 0.5*(iRawTimeNs[i * NbChs + 1] - iRawTimeNs[i * NbChs]);
                        raw_time_of_flight=(0.5*(iRawTimeNs[i*NbChs+1]+iRawTimeNs[i*NbChs]))-0.5*(iRawTimeStartNs[0]+iRawTimeStartNs[1]);
                        Double_t time_target_tof = raw_time_of_flight+tof_offs[i] -time_start_target;
                        Double_t path_from_target = abs((zB-zT)/(cos(psi_in)))+rho*w+abs((zToFW-z_D)/cos(psi_out));
                        Double_t beta = ((path_from_target/time_target_tof)*pow(10,6))/light_c;
                        Double_t mag_field = current*0.0006527728074785267;
                        Double_t a_q = (pow(10,-3)*((((mag_field*rho)/(beta))*1.602176634*pow(10,-19))/(1.66053906660*pow(10,-27)*light_c)))/gamma_given;
                        if (charge_val< t && charge_val> t-1.) //Z = 5 before it was charge_val< charge_cut && charge_val> charge_cut-0.8, too straight ....
                        {
                        if (a_q > s[0] && a_q < (s[0]+2*s[3]) && entries_califa >= 2){  //11B
                        //psi_in = mean_psi_out_11b -slope_11b*xMW3 - angle_offs_11b;
						psi_in = u[3] - u[4]*xMW3 - u[5];
                        m = (cos(psi_out)-cos(psi_in))/(sin(psi_out)+sin(psi_in));
                        rho =  (Leff/(2*sin((psi_in+psi_out)/2)))*sqrt(pow(((m+tan(alpha_G))/(1-m*tan(alpha_G))),2)+1);
                        w = 2*abs(asin(BD/(2*rho)));
                        Double_t pathlength_precise = abs((zB-zT)/(cos(psi_in)))+rho*w+abs((zToFW-z_D)/cos(psi_out));
                        Double_t beta_precise = ((pathlength_precise/time_target_tof)*pow(10,6))/light_c;
                        Double_t gamma_precise = 1/(sqrt(1-beta_precise*beta_precise));
                        Double_t a_q_precise = (pow(10,-3)*((((mag_field*rho)/(beta_precise))*1.602176634*pow(10,-19))/(1.66053906660*pow(10,-27)*light_c)))/gamma_given;
                        Double_t mom_11b = 11*931.5*(beta_precise/(sqrt(1-(beta_precise*beta_precise))));
                        h2_z_vs_a_q->Fill(a_q_precise,charge_val);
                        h2_tof_vs_aq_fix_g_11b->Fill(a_q_precise,time_target_tof);
						//CALIFA handling
						                        Double_t psi_CALIFA;
                        Double_t theta_CALIFA;
                        Double_t time_CALIFA;
                        Double_t E_Califa;
                        vector<double> vE1;
                        vector<double> vtheta1CALIFA;
                        vector<double> vpsi1CALIFA;
                        vector<double> vE2;
                        vector<double> vtheta2CALIFA;
                        vector<double> vpsi2CALIFA;
                        vector<double> vCALIFA_time_1;
                        vector<double> vCALIFA_time_2;
                        vector<double> v_E_gamma;
                        vector<double> vpsi;
                        vector<double> vtheta;
                        vector<double> vE;
                        vector<double> v_gamma_sorted;
                        vector<double> vtime;
                        vector<double> v_origin_E;
                        Double_t ind_high_energy;
                        Double_t ind_second_high_energy;
                        Double_t ind_first_gamma;
                        Double_t ind_second_gamma;
                        Double_t theta_gamma;
                        Double_t psi1;
                        Double_t psi2;
                        Double_t psi1_2pi;
                        Double_t psi2_2pi;
                        Double_t psi1_2pi_degr;
                        Double_t psi2_2pi_degr;
                        Double_t theta1_degr;
                        Double_t theta2_degr;
                        Double_t complete_gamma_deposited = 0.;
                        for (Int_t j = 0.;j<entries_califa;j++){
                            califahitdata[j] = (R3BCalifaHitData*)CalifaHitData->At(j);
                            psi_CALIFA = califahitdata[j]->GetPhi();
                            theta_CALIFA = califahitdata[j]->GetTheta();
                            E_Califa = califahitdata[j]->GetEnergy();
                            time_CALIFA = califahitdata[j]->GetTime();
                            vpsi.push_back(psi_CALIFA);
                            vtheta.push_back(theta_CALIFA);
                            vE.push_back(E_Califa);
                            v_origin_E.push_back(E_Califa);
                            vtime.push_back(time_CALIFA);
                            }
                        sort(vE.begin(),vE.end());
                        if (vE[vE.size()-1] >30000 && vE[vE.size()-2] > 30000)
                        {
                        vector<vector<double> > energy_theta_v;
                        for (Int_t k =0;k< vE.size();k++)
                        {
                            if (v_origin_E[k] < 10000){
                            vector<double> temp_vec;
                            temp_vec.push_back(vtheta[k]);                                                                                      //column 0 for angle theta
                            temp_vec.push_back((v_origin_E[k]/1000)*(1/(sqrt(1-beta_precise*beta_precise)))*(1-beta_precise*cos(vtheta[k])));   //column 1 for energy
                            energy_theta_v.push_back(temp_vec);
                            temp_vec.clear();
                            }
                            if (v_origin_E[k] == vE[vE.size()-1]) {
                            ind_high_energy = k;
                            }
                            else if (v_origin_E[k] == vE[vE.size()-2]){
                            ind_second_high_energy = k;
                            }

                            else if(v_origin_E[k] < 10000)
                            {
                            theta_gamma = vtheta[k];
                            Double_t energy_gamma = (v_origin_E[k]/1000)*(1/(sqrt(1-beta_precise*beta_precise)))*(1-beta_precise*cos(theta_gamma));
                            v_E_gamma.push_back(energy_gamma);
                            h1_gamma_energyE_11B->Fill(energy_gamma);
                            complete_gamma_deposited +=v_origin_E[k]/1000;
                            if (vpsi[k]< 0.){
                            h2_gamma_energyE_psi_11B->Fill(energy_gamma,((2*PI+vpsi[k])/PI)*180);
                            }
                            if (vpsi[k] > 0.){
                            h2_gamma_energyE_psi_11B->Fill(energy_gamma,(vpsi[k]/PI)*180);

                            }
                            }



                        }
                        Double_t E1 = v_origin_E[ind_high_energy]/1000.;
                        Double_t E2 = v_origin_E[ind_second_high_energy]/1000.;
                        Double_t theta1 = vtheta[ind_high_energy];
                        Double_t theta2 = vtheta[ind_second_high_energy];
                        cout << "------------------------------------------------------" << endl;
                        cout << "------------------------------------------------------" << endl;
                        cout << "phi1 = \t" << (vpsi[ind_high_energy]/PI)*180. << endl;
                        cout << "phi2 = \t" << (vpsi[ind_second_high_energy]/PI)*180. << endl;
                        cout << "------------------------------------------------------" << endl;
                        if (vpsi[ind_high_energy] < 0.){
                        psi1 = PI - vpsi[ind_high_energy];
                        psi1_2pi = 2*PI+vpsi[ind_high_energy];
                        }
                        if (vpsi[ind_high_energy] > 0.){
                        psi1_2pi = vpsi[ind_high_energy];
                        }
                        if (vpsi[ind_second_high_energy] < 0.){
                        psi2 = PI - vpsi[ind_second_high_energy];
                        psi2_2pi = 2*PI+vpsi[ind_second_high_energy];
                        }
                        if (vpsi[ind_second_high_energy] > 0.){
                        psi2_2pi = vpsi[ind_second_high_energy];
                        }

                        psi1_2pi_degr = (psi1_2pi/PI)*180.;
                        psi2_2pi_degr = (psi2_2pi/PI)*180.;
                        theta1_degr = (theta1/PI)*180.;
                        theta2_degr = (theta2/PI)*180.;

                        h1_abs_phi_diff_11b->Fill(abs(psi1_2pi_degr-psi2_2pi_degr));


                        Double_t mom_p1 = sqrt(pow((E1+938.272),2)-pow(938.272,2));
                        Double_t mom_p2 = sqrt(pow((E2+938.272),2)-pow(938.272,2));
                        Double_t beta_1 = mom_p1/(sqrt(mom_p1*mom_p1+938.272*938.272));
                        Double_t gamma_1 = 1/(sqrt(1-beta_1*beta_1));
                        Double_t beta_2 = mom_p2/(sqrt(mom_p2*mom_p2+938.272*938.272));
                        Double_t gamma_2 = 1/(sqrt(1-beta_2*beta_2));
                        Double_t angle_p1_11b = abs(atan(tan(theta1)*cos(psi1_2pi) )-psi_in);
                        Double_t angle_p2_11b = abs(atan(tan(theta2)*cos(psi2_2pi) )-psi_in);
                        Double_t angle_p2_p1 = abs(atan(tan(theta1)*cos(psi1_2pi) )-atan(tan(theta2)*cos(psi2_2pi)));
                        Double_t exc_energy = sqrt(2*938.272*938.272+(11*931.5)*(11*931.5)+2*gamma_1*gamma_2*938.272*938.272*(1-beta_1*beta_2*cos(angle_p2_p1))+gamma_1*gamma_precise*(11*931.5)*938.272*(1-beta_precise*beta_1*cos(angle_p1_11b))+gamma_2*gamma_precise*(11*931.5)*938.272*(1-beta_precise*beta_2*cos(angle_p2_11b)) + 938.272*0.714549*1.42942*938.272*0.714549*1.42942)+complete_gamma_deposited-(12*931.5)-938.272-1.42942*938;
                        cout << "exec energy\t "<< exc_energy << endl;
                        h1_psi_in_11B->Fill(psi_in);
                        TVector2 k_0 =   TVector2(11414.5,0);
                        TVector2 k_1 =   TVector2(mom_p1*cos(theta1),mom_p1*sin(theta1));
                        TVector2 k_2 =   TVector2(mom_p2*cos(theta2),mom_p2*sin(theta2));
                        TVector2 k_a_1 = TVector2(mom_11b*cos(psi_in),mom_11b*sin(psi_in));
                        Double_t cos_distr = ((k_0 - k_1 - k_2)*k_a_1)/(sqrt((k_0 - k_1 - k_2)*(k_0 - k_1 - k_2))*sqrt(k_a_1*k_a_1));
                        //CALIFA histos..
                        if (v_E_gamma.size()){
						cout << "Event with gamma from p2p:\t" << evtnr << endl;
						cout << "p2p gamma energy:\t" << *max_element(v_E_gamma.begin(),v_E_gamma.end()) << endl;
                        h1_gamma_energyE_max_val_11B->Fill(*max_element(v_E_gamma.begin(),v_E_gamma.end()));
                        if(v_E_gamma.size() > 1){
                        sort(v_E_gamma.begin(),v_E_gamma.end());
                        h2_gamma_fst_vs_snd_11B->Fill(v_E_gamma[v_E_gamma.size()-2],v_E_gamma[v_E_gamma.size()-1]);
                        }
                        }
                        h1_cos_gamma_11b_p_i->Fill(cos_distr);
                        h1_theta_1_plus_theta_2_CALIFA_11b->Fill(((theta1+theta2)/PI)*180.);
                        h2_E1_vs_E2_CALIFA_11b->Fill(E1,E2);
                        h1_E1_plus_E2_CALIFA_11b->Fill(E1+E2);
                        h2_long_mom_p1p2_long_mom11b->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                        if (v_E_gamma.size() &&  *max_element(v_E_gamma.begin(),v_E_gamma.end()) > 1.){
                        h2_long_mom_p1p2_long_mom11b_large_gamma->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                        }
                        h1_transv_mom_difference_p1_p2_11b->Fill(abs(mom_p1*sin(theta1)-mom_p2*sin(theta2)));
                        h2_E1_plus_E2_CALIFA_vs_full_mom11b->Fill(E1+E2,mom_11b);
                        h1_binding_energy_11b->Fill(11*400-(E1+E2+sqrt(pow(mom_11b,2)+pow(11*938.272,2))-(11*938.272)));
                        h2_theta1_theta2_11b->Fill((theta1/PI)*180.,(theta2/PI)*180.);
                        h2_psi1_psi2_11b->Fill((psi1/PI)*180.,(psi2/PI)*180.);
                        h2_psi1_psi2_11b_2pi->Fill((psi1_2pi/PI)*180.,(psi2_2pi/PI)*180.);
                        h1_theta_1_CALIFA_11b->Fill((theta1/PI)*180.);
                        h1_theta_2_CALIFA_11b->Fill((theta2/PI)*180.);
                        h2_theta_sum_vs_phi_diff_11b->Fill(theta1_degr+theta2_degr,abs(psi1_2pi_degr-psi2_2pi_degr));
                        if (-0.9973*(mom_p1*cos(theta1)+mom_p2*cos(theta2))+10913 > (mom_11b*cos(psi_in)) && (mom_11b*cos(psi_in)) > 9820 && (mom_p1*cos(theta1)+mom_p2*cos(theta2)) > 750) {
                        h2_theta_1_vs_theta_2_CALIFA_11b_cut->Fill((theta1/PI)*180.,(theta2/PI)*180.);
                        h2_E1_vs_E2_CALIFA_11b_cut->Fill(E1,E2);
                        h1_psi_in_11B_cut_triangle->Fill(psi_in);

                        }
                        else {
						}
                        TVector3 p_11b (11*931.5*gamma_precise*beta_precise*sin(psi_in),0.,-gamma_given*0.714549*11*931.5*gamma_precise+11*931.5*gamma_precise*beta_precise*cos(psi_in)*gamma_given);
                        TVector3 p_p1(938.272*gamma_1*beta_1*sin(theta1)*cos(psi1_2pi),0.,-gamma_given*0.714549*938.272*gamma_1+gamma_given*938.272*gamma_1*beta_1*cos(theta1));
                        TVector3 p_p2(938.272*gamma_2*beta_2*sin(theta2)*cos(psi2_2pi),0.,-gamma_given*0.714549*938.272*gamma_2+gamma_given*938.272*gamma_2*beta_2*cos(theta2));
                        TVector3 p_tar(0.,0.,-gamma_given*0.714549*938.272);
                        TVector3 p_initial(p_p1+p_p2-p_tar);
                        cout << "this is momentum of initial projectile\t" << p_initial.Mag() << endl;
                        h1_p_missing_11B->Fill(p_initial.Mag());
                        Double_t my_sweep = 0.048869219056;
                        vector<double> cos_sweep;
                        for (Int_t sweep_1 = -1; sweep_1 < 2; sweep_1++){
                            for (Int_t sweep_2 = -1; sweep_2 < 2; sweep_2++){

                        TVector3 p_11b_sweep (11*931.5*gamma_precise*beta_precise*sin(psi_in),0.,-gamma_given*0.714549*11*931.5*gamma_precise+11*931.5*gamma_precise*beta_precise*cos(psi_in)*gamma_given);
                        TVector3 p_p1_sweep(938.272*gamma_1*beta_1*sin(theta1)*cos(psi1_2pi+sweep_1*my_sweep),0.,-gamma_given*0.714549*938.272*gamma_1+gamma_given*938.272*gamma_1*beta_1*cos(theta1));
                        TVector3 p_p2_sweep(938.272*gamma_2*beta_2*sin(theta2)*cos(psi2_2pi+sweep_2*my_sweep),0.,-gamma_given*0.714549*938.272*gamma_2+gamma_given*938.272*gamma_2*beta_2*cos(theta2));
                        TVector3 p_tar_sweep(0.,0.,-gamma_given*0.714549*938.272);
                        TVector3 p_initial_sweep(p_p1_sweep+p_p2_sweep-p_tar_sweep);
                        Double_t angle_cms_sweep = p_initial_sweep.Angle(p_11b_sweep);
                        cos_sweep.push_back(cos(angle_cms_sweep));
                            }


                        }
                        Double_t min_cos_angle_cms = *min_element(cos_sweep.begin(),cos_sweep.end());
                        cos_sweep.clear();
                        Double_t angle_cms = p_initial.Angle(p_11b);
                        cout << "cos(gamma) questioning....." << endl;
                        if (cos(angle_cms) < -0.6){
                        h2_long_mom_p1p2_long_mom11b_low->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                        }
                        else if (cos(angle_cms)>-0.6){
                        h2_long_mom_p1p2_long_mom11b_high->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                        }
                        if ((theta1_degr+theta2_degr) < 84 ){
                        h2_long_mom_p1p2_long_mom11b_small_opening->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                        }
                        else if ((theta1_degr+theta2_degr) > 84 ){
                        h2_long_mom_p1p2_long_mom11b_large_opening->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                        }
                        cout << "asking for angular restrictions...." << endl;
                        if ((theta1_degr+theta2_degr) < 84 && abs(180-abs(psi1_2pi_degr-psi2_2pi_degr)) < 30){
                            if(energy_theta_v.size() > 1){
                                cout << "AFTER asking!" << endl;
                                sort(energy_theta_v.begin(),energy_theta_v.end(),sortcol);
                                cout << "you sorted it riht..."<< endl;

                                Double_t mev_lab_2_1 = 2.124/((1/(sqrt(1-beta_precise*beta_precise)))*(1-beta_precise*cos(energy_theta_v[0][0])));
                                Double_t mev_lab_4_4 = 4.445/((1/(sqrt(1-beta_precise*beta_precise)))*(1-beta_precise*cos(energy_theta_v[0][0])));
                                Double_t mev_lab_5_0 = 5.020/((1/(sqrt(1-beta_precise*beta_precise)))*(1-beta_precise*cos(energy_theta_v[0][0])));
                                Double_t escape_peak = energy_theta_v[1][1]*(sqrt(1-beta_precise*beta_precise))/(1-beta_precise*cos(energy_theta_v[1][0]));
                                Double_t full_spec_energy = energy_theta_v[0][1];
                                cout << "I reach the lab stuff..." << endl;
                                if ((energy_theta_v[0][1] > (mev_lab_2_1-0.6) && energy_theta_v[0][1] < (mev_lab_2_1-0.4) && escape_peak < 0.6)||( energy_theta_v[0][1] > (mev_lab_2_1-1.1) && energy_theta_v[0][1] < (mev_lab_2_1-0.9) && escape_peak < 1.1)\
                                || (energy_theta_v[0][1] > (mev_lab_4_4-0.6) && energy_theta_v[0][1] < (mev_lab_4_4-0.4) && escape_peak < 0.6)|| (energy_theta_v[0][1] > (mev_lab_4_4-1.1) && energy_theta_v[0][1] < (mev_lab_4_4-0.9) && escape_peak < 1.1)\
                                || (energy_theta_v[0][1] > (mev_lab_5_0-0.6) && energy_theta_v[0][1] < (mev_lab_5_0-0.4) && escape_peak < 0.6)|| (energy_theta_v[0][1] > (mev_lab_5_0-1.1) && energy_theta_v[0][1] < (mev_lab_5_0-0.9) && escape_peak < 1.1)){
                                cout << "this is spec_energy:\t" << full_spec_energy <<endl;
                                for (Int_t spec = 1; spec < energy_theta_v.size(); spec++){
                                    full_spec_energy += (energy_theta_v[spec][1]*(sqrt(1-beta_precise*beta_precise))/(1-beta_precise*cos(energy_theta_v[spec][0])))*((1/(sqrt(1-beta_precise*beta_precise)))*(1-beta_precise*cos(energy_theta_v[0][0]))); //here the energy of 2nd gamma is first transformed back to laboratory then back to the particle system taking the angle of the gamma with the largest engergy
                                    cout << " in the spec loop" <<endl;
                                }
                                cout << "fillig spec histo...." << endl;
                                h1_gamma_spec_escape_angle_cut->Fill(full_spec_energy);
                                h2_gamma_energy_vs_theta_11B->Fill(full_spec_energy,((energy_theta_v[0][0])/PI)*180);
                                }
                                else {
                                    h1_gamma_spec_escape_angle_cut->Fill(full_spec_energy);
                                    h2_gamma_energy_vs_theta_11B->Fill(full_spec_energy,((energy_theta_v[0][0])/PI)*180);
                                }
                            }
                            else if (energy_theta_v.size() == 1){
                                h1_gamma_spec_escape_angle_cut->Fill(energy_theta_v[0][1]);
                                h2_gamma_energy_vs_theta_11B->Fill(energy_theta_v[0][1],((energy_theta_v[0][0])/PI)*180);
                            }
                            h2_long_mom_p1p2_long_mom11b_theta_phi_constr->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                            if (p_initial.Mag() < 250){
                                h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                                if (v_E_gamma.size()){
                                cout << "after questioning vE_gamma_size..." << endl;
                                cout << "this is max element:\t" << *max_element(v_E_gamma.begin(),v_E_gamma.end()) << endl;
                                Double_t max_gamma_test = *max_element(v_E_gamma.begin(),v_E_gamma.end());
                                h1_m_mom_250_gamma_spec->Fill(max_gamma_test);
                                cout << "after filling h1_m_mom_250_gamma_spec ... " << endl;
                                cout << "after h1_m_mom_250_gamma_spec->Fill .." << endl;
                                Int_t gamma_cluster = 0;
                                Int_t proton_cluster = 0;
                                cout << "right before 3rd hello...." << endl;
                                if (vE.size()){
                                    cout << "3rdHELLO" << endl;
                                    for (Int_t cl = 0;cl <vE.size();cl++){
                                        if (vE[cl] > 30000){
                                        cout << "TJ, hello from inside proton cluster" << endl;
                                        proton_cluster+= 1;
                                        }
                                        else{
                                        gamma_cluster+=1;
                                        cout << "TJ, hello from inside gamma cluster" << endl;
                                        }

                                    }
                                }
                                h1_cluster_gamma_11b->Fill(gamma_cluster);
                                h1_cluster_proton_11b->Fill(proton_cluster);
                                }
                                }
                            if (p_initial.Mag() < 200 && v_E_gamma.size()){
                                h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom200->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                                h1_m_mom_200_gamma_spec->Fill(*max_element(v_E_gamma.begin(),v_E_gamma.end()));
                                }
                            if (p_initial.Mag() < 150 &&  v_E_gamma.size()){
                                h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom150->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                                h1_m_mom_150_gamma_spec->Fill(*max_element(v_E_gamma.begin(),v_E_gamma.end()));
                                }
                            if (min_cos_angle_cms < -0.6){
                                h2_long_mom_p1p2_long_mom11b_low_optimized->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                                if (v_E_gamma.size()){
                                    h1_cos_low_gamma->Fill(*max_element(v_E_gamma.begin(),v_E_gamma.end()));
                                }
                            }
                            if (min_cos_angle_cms > 0.6){
                                h2_long_mom_p1p2_long_mom11b_high_optimized->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                            }

                            if (min_cos_angle_cms < -0.8){
                                h2_long_mom_p1p2_long_mom11b_low_optimized_08->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                            }
                            cout << "before triangle cut" <<endl;
                            if (-0.9973*(mom_p1*cos(theta1)+mom_p2*cos(theta2))+10913 > (mom_11b*cos(psi_in)) && (mom_11b*cos(psi_in)) > 9820 && (mom_p1*cos(theta1)+mom_p2*cos(theta2)) > 750) {
                            h1_mw3_fy_11b_cut_triangle->Fill(raw_mwpc3);
                            h1_mw2_fy_11b_cut_triangle->Fill(yMW2);
                            h1_psi_in_11b_cut_triangle->Fill(psi_in);

                            if (v_E_gamma.size()){
                            h1_gamma_energyE_max_val_11B_angle_cut->Fill(v_E_gamma[v_E_gamma.size()-1]);
                            }

                            }
                            else {
                            h1_mw3_fy_11b->Fill(raw_mwpc3);
                            h1_mw2_fy_11b->Fill(yMW2);
                            h1_psi_in_11b->Fill(psi_in);
                            if (v_E_gamma.size()){
                            h1_gamma_energyE_max_val_11B_angle_cut->Fill(v_E_gamma[v_E_gamma.size()-1]);
                            }
                            }
                          if (v_E_gamma.size()){
                          v_E_gamma.clear();
                          }
                        }
                        cout << "before border constraints...." << endl;
                        cout << "psi1 in degrees:\t" << psi1_2pi_degr<< endl;
                        cout << "psi2 in degrees:\t" << psi2_2pi_degr << endl;
                        cout << "theta1 in degrees:\t" << theta1_degr << endl;
                        cout << "theta2 in degrees:\t" << theta2_degr << endl;
                         if( (psi1_2pi_degr > 84 && psi1_2pi_degr < 96) || (psi1_2pi_degr > 264 && psi1_2pi_degr < 276) || (psi2_2pi_degr > 84 && psi2_2pi_degr < 96) || (psi2_2pi_degr > 264 && psi2_2pi_degr)){
                            cout << "hello in first condition" << endl;
                            continue;
                            }
                        if((psi1_2pi_degr> 42 &&  psi1_2pi_degr < 160.5 && theta1_degr < 46) || (psi1_2pi_degr> 199.5 &&  psi1_2pi_degr < 318 && theta1_degr < 46)){
                            cout << "hello in second contidtion" << endl;
                            continue;
                            }
                        if((psi2_2pi_degr> 42 &&  psi2_2pi_degr < 160.5 && theta2_degr < 46) || (psi2_2pi_degr> 199.5 &&  psi2_2pi_degr < 318 && theta2_degr < 46)){
                            cout << "hello in third condition " << endl;
                            continue;
                            }
                        if (theta1_degr < 25 || theta2_degr < 25){
                            cout << "hello in forth condition" << endl;
                            continue;
                            }
                        cout << "After border constraints ...." << endl;
                        h2_long_mom_p1p2_long_mom11b_no_border->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                        h1_cos_gamma_11b_p_i_cms->Fill(cos(angle_cms));
                        h1_cos_gamma_11b_p_i_cms_optimized->Fill(min_cos_angle_cms);
                        cout << "after filling some histograms..." << endl;
                        if (cos(angle_cms) < -0.6){
                        h2_long_mom_p1p2_long_mom11b_no_border_low->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                        }
                        else if (cos(angle_cms)>0.6){
                        h2_long_mom_p1p2_long_mom11b_no_border_high->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                        }
                        cout << "After asking cos(gamma) <> -0.6...."<< endl;
                        cout << "last calculation....." << endl;
                        for (Int_t k =0;k< vE.size();k++)
                        {
                            cout << "inside the last for loop in 11B" << endl;
                            if(v_origin_E[k] < 10000)
                            {
                            cout << "inside the last if questioning..." << endl;
                            cout << "this is v_origin_E[k]" << v_origin_E[k] << endl;
                            cout << " I shouldn't see this line if v_origin_E[k] is empty..." << endl;
                            Double_t theta_gamma_noborder = vtheta[k];
                            Double_t energy_gamma = (v_origin_E[k]/1000)*(1/(sqrt(1-beta*beta)))*(1-beta*cos(theta_gamma_noborder));
                            h1_gamma_energyE_11B_no_border->Fill(energy_gamma);

                            }



                        }
					}
                        }
                        if (a_q < s[0] && a_q > (s[0]-0.15) && entries_califa >= 2){  //10B
                        //psi_in = mean_psi_out_10b-slope_10b*xMW3 - angle_offs_10b;
						psi_in = u[0]-u[1]*xMW3 -u[2];
                        m = (cos(psi_out)-cos(psi_in))/(sin(psi_out)+sin(psi_in));
                        rho =  (Leff/(2*sin((psi_in+psi_out)/2)))*sqrt(pow(((m+tan(alpha_G))/(1-m*tan(alpha_G))),2)+1);
                        w = 2*abs(asin(BD/(2*rho)));
                        Double_t pathlength_precise = abs((zB-zT)/(cos(psi_in)))+rho*w+abs((zToFW-z_D)/cos(psi_out));
                        Double_t beta_precise = ((pathlength_precise/time_target_tof)*pow(10,6))/light_c;
                        Double_t a_q_precise = (pow(10,-3)*((((mag_field*rho)/(beta_precise))*1.602176634*pow(10,-19))/(1.66053906660*pow(10,-27)*light_c)))/gamma_given;
                        Double_t mom_10b = a_q_precise*5*(beta/(sqrt(1-(beta*beta))))*938.272;

                        Double_t psi_CALIFA;
                        Double_t theta_CALIFA;
                        Double_t time_CALIFA;
                        Double_t E_Califa;
                        vector<double> vE1;
                        vector<double> vtheta1CALIFA;
                        vector<double> vpsi1CALIFA;
                        vector<double> vE2;
                        vector<double> vtheta2CALIFA;
                        vector<double> vpsi2CALIFA;
                        vector<double> vpsi;
                        vector<double> vtheta;
                        vector<double> vE;
                        vector<double> vtime;
                        vector<double> v_origin_E;
                        Double_t ind_high_energy;
                        Double_t ind_second_high_energy;
                        Double_t theta_gamma;
                        Double_t psi1;
                        Double_t psi2;
                        for (Int_t j = 0.;j<entries_califa;j++){
                            califahitdata[j] = (R3BCalifaHitData*)CalifaHitData->At(j);
                            psi_CALIFA = califahitdata[j]->GetPhi();
                            theta_CALIFA = califahitdata[j]->GetTheta();
                            E_Califa = califahitdata[j]->GetEnergy();
                            time_CALIFA = califahitdata[j]->GetTime();
                            vpsi.push_back(psi_CALIFA);
                            vtheta.push_back(theta_CALIFA);
                            vE.push_back(E_Califa);
                            v_origin_E.push_back(E_Califa);
                            vtime.push_back(time_CALIFA);
                            }
                        sort(vE.begin(),vE.end());
                        if (vE[vE.size()-1] >30000 && vE[vE.size()-2] > 30000)
                        {
                        for (Int_t k =0;k<vE.size();k++)
                        {
                            if (v_origin_E[k] == vE[vE.size()-1]) {
                            ind_high_energy = k;
                            }
                            else if (v_origin_E[k] == vE[vE.size()-2]){
                            ind_second_high_energy = k;
                            }
                            else if(v_origin_E[k] < 10000)
                            {
                            theta_gamma = vtheta[k];
                            Double_t energy_gamma = (v_origin_E[k]/1000)*(1/(sqrt(1-beta*beta)))*(1-beta*cos(theta_gamma));
                            h1_gamma_energyE_10B->Fill(energy_gamma);

                            }

                        }
                        Double_t E1 = v_origin_E[ind_high_energy]/1000.;
                        Double_t E2 = v_origin_E[ind_second_high_energy]/1000.;
                        Double_t theta1 = vtheta[ind_high_energy];
                        Double_t theta2 = vtheta[ind_second_high_energy];

                        if (vpsi[ind_high_energy] < 0.){
                        psi1 = 2*PI + vpsi[ind_high_energy];
                        }
                        if (vpsi[ind_high_energy] > 0.){
                        psi1 = vpsi[ind_high_energy];
                        }
                        if (vpsi[ind_second_high_energy] < 0.){
                        psi2 = 2*PI + vpsi[ind_second_high_energy];
                        }
                        if (vpsi[ind_second_high_energy] > 0.){
                        psi2 = vpsi[ind_second_high_energy];
                        }

                        Double_t mom_p1 = sqrt(pow((E1+938.272),2)-pow(938.272,2));
                        Double_t mom_p2 = sqrt(pow((E2+938.272),2)-pow(938.272,2));
                        h1_theta_1_plus_theta_2_CALIFA_10b->Fill(((theta1+theta2)/PI)*180.);
                        h2_E1_vs_E2_CALIFA_10b->Fill(E1,E2);
                        h1_E1_plus_E2_CALIFA_10b->Fill(E1+E2);
                        h2_long_mom_p1p2_long_mom10b->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_10b*cos(psi_in));
                        h1_transv_mom_difference_p1_p2_10b->Fill(abs(mom_p1*sin(theta1)-mom_p2*sin(theta2)));
                        }
                        h2_z_vs_a_q->Fill(a_q_precise,charge_val);
                        h2_tof_vs_aq_fix_g_10b->Fill(a_q_precise,time_target_tof);
                        }
                        }
                        if (charge_val< t+0.8 && charge_val> t) //Z = 6
                        {
                        if (a_q > s[1] && a_q < (s[1] +2*s[5])){  //12C
                        //psi_in = mean_psi_out_12c-slope_12c*xMW3 - angle_offs_12c;
						psi_in = u[9] -u[10]*xMW3 - u[11];
                        m = (cos(psi_out)-cos(psi_in))/(sin(psi_out)+sin(psi_in));
                        rho =  (Leff/(2*sin((psi_in+psi_out)/2)))*sqrt(pow(((m+tan(alpha_G))/(1-m*tan(alpha_G))),2)+1);
                        w = 2*abs(asin(BD/(2*rho)));
                        Double_t pathlength_precise = abs((zB-zT)/(cos(psi_in)))+rho*w+abs((zToFW-z_D)/cos(psi_out));
                        Double_t beta_precise = ((pathlength_precise/time_target_tof)*pow(10,6))/light_c;
                        Double_t a_q_precise = (pow(10,-3)*((((mag_field*rho)/(beta_precise))*1.602176634*pow(10,-19))/(1.66053906660*pow(10,-27)*light_c)))/gamma_given;
                        if (time_target_tof < 35.13){
                        h1_tof_detnr_strange_12c->Fill(i+1);
                        h1_mw3_fy_strange_12c->Fill(raw_mwpc3);
                        h2_tof_path_strange12c->Fill(pathlength_precise,time_target_tof);
                        h2_radius_path_strange12c->Fill(rho,time_target_tof);
                        h1_tofwall_posy_fast[i]->Fill(tof_diff_up_down);
                        }
                        if (time_target_tof> 35.35 && time_target_tof<35.61){
						//cout << "HEEEEEEE-----------------:\t" << i+1 << endl;
                        h1_tof_detnr_12c->Fill(i+1);
						//cout << "AFTER FILLING!!!!!!!!!!" << endl;
                        h1_mw3_fy_12c->Fill(raw_mwpc3);
						//cout << "after h1_mw3_fy_12c->Fill(raw_mwpc3)" << endl;
                        h2_tof_path_12c->Fill(pathlength_precise,time_target_tof);
						//cout << "after h2_tof_path_12c->Fill(pathlength_precise,time_target_tof)" << endl;
                        h2_radius_path_12c->Fill(rho,time_target_tof);
						//cout << "after h2_radius_path_12c->Fill(rho,time_target_tof)" << endl;
                        h1_tofwall_posy_medium[i]->Fill(tof_diff_up_down);
						//cout << "at end of if condition..." << endl;
                        }
                        if (a_q_precise > 1.997){
                        h1_tof_detnr_strange_12c_large_aq->Fill(i+1);
                        h1_mw3_fy_strange_12c_large_aq->Fill(raw_mwpc3);
                        h2_tof_path_strange_12c_large_aq->Fill(pathlength_precise,time_target_tof);
                        h2_radius_path_strange_12c_large_aq->Fill(rho,time_target_tof);
                        h1_tofwall_posy_slow[i]->Fill(tof_diff_up_down);
                        }
                        h2_z_vs_a_q->Fill(a_q_precise,charge_val);
                        h2_tof_vs_aq_fix_g_12c->Fill(a_q_precise,time_target_tof);
                        h2_mw2x_vs_tof_12c->Fill(xMW2,time_target_tof);
                        h2_mw3x_vs_tof_12c->Fill(xMW3,time_target_tof);
                        h2_charge_vs_tof_12c->Fill(charge_val,time_target_tof);
                        h2_psi_in_vs_tof_12c->Fill(psi_in,time_target_tof);
                        h2_detnr_vs_tof_12c->Fill(i+1,time_target_tof);
                        h2_radius_vs_beta_12c->Fill(beta_precise,rho);
                        h2_radius_vs_tof_12c->Fill(rho,time_target_tof);
                        }
                        if (a_q < s[1]-0.06 && a_q > (s[1]-1.65*s[4]) && entries_califa >= 2){  //11C before it was a_q < aq_cut_z_6 && a_q > (aq_cut_z_6-2*width_11c), too bad for 11c...
                        //psi_in = mean_psi_out_11c-slope_11c*xMW3 - angle_offs_11c;
						psi_in = u[6] - u[7]*xMW3 - u[8];
                        m = (cos(psi_out)-cos(psi_in))/(sin(psi_out)+sin(psi_in));
                        rho =  (Leff/(2*sin((psi_in+psi_out)/2)))*sqrt(pow(((m+tan(alpha_G))/(1-m*tan(alpha_G))),2)+1);
                        w = 2*abs(asin(BD/(2*rho)));
                        Double_t pathlength_precise = abs((zB-zT)/(cos(psi_in)))+rho*w+abs((zToFW-z_D)/cos(psi_out));
                        Double_t beta_precise = ((pathlength_precise/time_target_tof)*pow(10,6))/light_c;
                        Double_t a_q_precise = (pow(10,-3)*((((mag_field*rho)/(beta_precise))*1.602176634*pow(10,-19))/(1.66053906660*pow(10,-27)*light_c)))/gamma_given;
                        h2_z_vs_a_q->Fill(a_q_precise,charge_val);
                        h2_tof_vs_aq_fix_g_11c->Fill(a_q_precise,time_target_tof);
                        }

                        }
}
}
}
}
}
}
 delete [] sofmwpc0hitdata;
 delete [] sofmwpc1hitdata;
 delete [] sofmwpc2hitdata;
 delete [] sofmwpc3hitdata;
 delete [] softofwtcaldata;
 delete [] softofwmappeddata;
 delete [] califahitdata;
}
delete [] softwimhitdata;
}

cout << "before adding everything to TList..." << endl;
TList *l = new TList();
l->Add(h2_z_vs_a_q);
l->Add(h1_a_q_z6);
l->Add(h1_a_q_z5);
l->Add(h1_gamma_energyE_10B);
l->Add(h1_gamma_energyE_11B);
l->Add(h1_theta_1_plus_theta_2_CALIFA_10b);
l->Add(h1_theta_1_plus_theta_2_CALIFA_11b);
l->Add(h2_E1_vs_E2_CALIFA_10b);
l->Add(h2_E1_vs_E2_CALIFA_11b);
l->Add(h1_E1_plus_E2_CALIFA_10b);
l->Add(h1_E1_plus_E2_CALIFA_11b);
l->Add(h2_long_mom_p1p2_long_mom11b);
l->Add(h2_long_mom_p1p2_long_mom10b);
l->Add(h1_transv_mom_difference_p1_p2_10b);
l->Add(h1_transv_mom_difference_p1_p2_11b);
l->Add(h2_tof_vs_aq_fix_g_11b);
l->Add(h2_tof_vs_aq_fix_g_10b);
l->Add(h2_tof_vs_aq_fix_g_11c);
l->Add(h2_tof_vs_aq_fix_g_12c);
l->Add(h2_theta_out_vs_mw3_10b);
l->Add(h2_theta_out_vs_mw3_11b);
l->Add(h2_theta_out_vs_mw3_11c);
l->Add(h2_theta_out_vs_mw3_12c);
l->Add(h2_E1_plus_E2_CALIFA_vs_full_mom11b);
l->Add(h2_theta_1_vs_theta_2_CALIFA_11b_cut);
l->Add(h2_E1_vs_E2_CALIFA_11b_cut);
l->Add(h1_binding_energy_11b);
l->Add(h1_cos_gamma_11b_p_i);
l->Add(h1_time_diff_gamma_11b);
l->Add(h2_theta1_theta2_11b);
l->Add(h2_psi1_psi2_11b);
l->Add(h2_psi1_psi2_11b_2pi);
l->Add(h1_tof_detnr_strange_12c);
l->Add(h1_tof_detnr_strange_12c_large_aq);
l->Add(h1_tof_detnr_12c);
l->Add(h1_mw3_fy_strange_12c);
l->Add(h1_mw3_fy_12c);
l->Add(h1_mw3_fy_strange_12c_large_aq);
l->Add(h2_z_vs_a_q_nocut);
l->Add(h1_theta_1_CALIFA_11b);
l->Add(h1_theta_2_CALIFA_11b);
l->Add(h2_gamma_energyE_psi_11B);
l->Add(h2_long_mom_p1p2_long_mom11b_no_border);
l->Add(h1_cos_gamma_11b_p_i_cms);
l->Add(h1_cos_gamma_11b_p_i_cms_optimized);
l->Add(h1_gamma_energyE_11B_no_border);
l->Add(h1_psi_in_11B_cut_triangle);
l->Add(h1_psi_in_11B);
l->Add(h2_long_mom_p1p2_long_mom11b_no_border_high);
l->Add(h2_long_mom_p1p2_long_mom11b_no_border_low);
l->Add(h2_long_mom_p1p2_long_mom11b_high);
l->Add(h2_long_mom_p1p2_long_mom11b_low);
l->Add(h1_gamma_energyE_max_val_11B);
l->Add(h2_gamma_fst_vs_snd_11B);
l->Add(h2_long_mom_p1p2_long_mom11b_small_opening);
l->Add(h2_long_mom_p1p2_long_mom11b_theta_phi_constr);
l->Add(h2_long_mom_p1p2_long_mom11b_large_opening);
l->Add(h1_gamma_energyE_max_val_11B_angle_cut);
l->Add(h1_gamma_spec_escape_angle_cut);
l->Add(h2_gamma_energy_vs_theta_11B);
//---------------------------------------------------
l->Add(h1_mw3_fy_11b_cut_triangle);
l->Add(h1_mw2_fy_11b_cut_triangle);
l->Add(h1_psi_in_11b_cut_triangle);
l->Add(h1_mw3_fy_11b);
l->Add(h1_mw2_fy_11b);
l->Add(h1_psi_in_11b);
l->Add(h2_tof_path_strange12c);
l->Add(h2_tof_path_strange_12c_large_aq);
l->Add(h2_tof_path_12c);
l->Add(h2_radius_path_strange12c);
l->Add(h2_radius_path_strange_12c_large_aq);
l->Add(h2_radius_path_12c);
l->Add(h2_mw2x_vs_tof_12c);
l->Add(h2_mw3x_vs_tof_12c);
l->Add(h2_charge_vs_tof_12c);
l->Add(h2_psi_in_vs_tof_12c);
l->Add(h2_detnr_vs_tof_12c);
l->Add(h2_long_mom_p1p2_long_mom11b_low_optimized);
l->Add(h2_long_mom_p1p2_long_mom11b_high_optimized);
l->Add(h2_long_mom_p1p2_long_mom11b_low_optimized_08);
l->Add(h1_p_missing_11B);
l->Add(h1_cos_low_gamma);
l->Add(h2_theta_sum_vs_phi_diff_11b);

l->Add(h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom);
l->Add(h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom200);
l->Add(h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom150);
cout << "before actually adding ..spec to TList (now It doesn't work..)" << endl;
l->Add(h1_m_mom_250_gamma_spec);
cout << "shortly after..." << endl;
l->Add(h1_m_mom_200_gamma_spec);
l->Add(h1_m_mom_150_gamma_spec);
l->Add(h1_cluster_gamma_11b);
l->Add(h1_cluster_proton_11b);
l->Add(h2_radius_vs_beta_12c);
l->Add(h2_radius_vs_tof_12c);
l->Add(h1_abs_phi_diff_11b);
l->Add(h2_long_mom_p1p2_long_mom11b_large_gamma);
for (Int_t i = 0; i < NbDets; i++){
l->Add(h1_tofwall_posy_fast[i]);
l->Add(h1_tofwall_posy_medium[i]);
l->Add(h1_tofwall_posy_slow[i]);
}
cout << "shortly after for loop in TList..." << endl;

l->Write("histlist", TObject::kSingleKey);

cout << "after creating file f..." << endl;
cout << "after l->Write...." << endl;
cout << "before deleting TList..." << endl;
cout << "end of process, successful" <<endl;
delete SofMwpc3HitData;
delete SofSciTcalData;
delete SofToFWTcalData;
delete SofMwpc0HitData;
delete SofMwpc1HitData;
delete SofMwpc2HitData;
delete SofToFWMappedData;
delete SofTwimHitData;
delete CalifaHitData;
cout << "please end!!" << endl;
}



