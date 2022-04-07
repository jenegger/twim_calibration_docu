#include <iostream>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include <math.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TChain.h"
#include "find_z.C"
#include "find_isotopes.C"
#include "all_params.C"
#include "find_psi_par.C"
#include <string>
using namespace std;

char fname[500];
bool sortcol( const vector<double>& v1,const vector<double>& v2 ) {
                                    return v1[1] > v2[1];
                                    }
bool sortcol_pseudo( const vector<double>& v1,const vector<double>& v2 ) {
                                    return v1[0] > v2[0];
                                    }
bool sortcol_lorentz( const TLorentzVector& v1,const TLorentzVector& v2 ) {
                                    return v1[3] > v2[3];
                                    }



void to_test_macro(char const count_i[50]){
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
TH1F* h1_gamma_energyE_sum_11B;
TH1F* h1_gamma_energyE_sum_11B_cut;
TH1F* h1_gamma_energyE_max_val_11B_m_cut;
TH1F* h1_gamma_energyE_max_val_10B;
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
TH1F* h1_cos_low_gamma;
TH2F* h2_theta_sum_vs_phi_diff_11b;

TH1F* h1_tofwall_posy_slow[NbDets];
TH1F* h1_tofwall_posy_medium[NbDets];
TH1F* h1_tofwall_posy_fast[NbDets];


//here are the histograms I added after doing the revision of the code
TH1F* h1_missing_energy_11B_p2p;
TH1F* h1_missing_energy_11B_p2p_cut;
TH1F* h1_missing_mass_11B_p2p_cut;
TH1F* h1_missing_mass_10B;
TH1F* h1_missing_mass_11B;
TH1F* h1_missing_mass_10B_neutron_miss_cut;
TH1F* h1_missing_mass_10B_rec;
TH1F* h1_missing_mass_10B_rec_miss_cut;
TH2F* h2_missing_energy_vs_t1_t2_11B_p2p;
TH1F* h1_p_missing_11B;
TH1F* h1_p_missing_11B_cut;
TH1F* h1_p_missing_11B_x;
TH1F* h1_p_missing_11B_x_cut;
TH1F* h1_p_missing_11B_y_cut;
TH1F* h1_p_missing_11B_z;
TH1F* h1_p_missing_11B_z_cut;
TH1F* h1_p_missing_10B;
TH2F* h2_missing_energy_vs_t1_t2_10B;
TH1F* h1_angle_pi_11B;
TH1F* h1_angle_pi_pn_10B;
TH1F* h1_angle_pi_pn_10B_cos;
TH1F* h1_angle_src_vs_10B_cos;
TH1F* h1_angle_pi_11B_cos;
TH1F* h1_angle_pi_11B_cos_cut;
TH1F* h1_angle_pi_11B_cos_cut_fake_y;
TH1F* h1_angle_pi_11B_cos_cut_plane;
TH1F* h1_angle_pi_11B_cos_hcut_plane;
TH1F* h1_gamma_energyE_max_val_11B_lab;
TH2F* h2_gamma_energy_11B_vs_mult;
TH2F* h2_gamma_energy_11B_vs_mult_cut;
TH2F* h2_gamma_energy_11_vs_angle;
TH2F* h2_theta1_theta2_10b;
TH2F* h2_theta1_theta2_10b_m_cut;
TH1F* h1_p_missing_10B_m_cut;
TH1F* h1_p_missing_10B_p_i_cut;
TH2F* h2_theta1_theta2_10b_double_cut;
TH2F* h2_theta1_theta2_10b_miss_mom_cut;
TH1F* h1_gamma_energyE_max_val_11B_cut;
TH1F* h1_gamma_energyE_max_val_11Bframe_cut;
TH2F* h2_p_missing_vs_p_11B_cut;
TH2F* h2_p_missing_vs_p_11B_cut_fake_y;
TH1F* h1_excit_energy_11B_cut;
TH2F* h2_corr_y_proton_fragm;
TH2F* h2_corr_y_proton_fragm_tj;
TH2F* h2_corr_y_proton_fragm_cos_less;
TH2F* h2_corr_y_proton_fragm_cos_more;
TH2F* h2_corr_x_proton_fragm;
TH2F* h2_missing_e_2meth_11B_p2p_cut;
TH2F* h2_corr_y_mws;
TH2F* h2_corr_y_mws_minusmw3;
TH2F* h2_E_sum_vs_theta_sum;
TH1F* h1_mom11B_z_cut;
TH1F* h1_mom11B_x_cut;
TH1F* h1_mom11B_y_cut;
TH1F* h1_mom11B_cut;
TH2F* h2_theta1_theta2_10b_lower250;
TH2F* h2_theta1_theta2_10b_lower650;
TH2F* h2_theta1_theta2_10b_higher650;
TH2F* h2_n_f_n_s_10b;
TH2F* h2_n_f_n_s_10b_lower250;
TH2F* h2_n_f_n_s_10b_lower650;
TH2F* h2_n_f_n_s_10b_higher650;

TH1F* h1_angle_pi_p_n_xy_plane;
TH1F* h1_angle_pi_p_n_xy_plane_mom_cut;
TH1F* h1_angle_src_vs_10B_cos_xy_plane;
TH1F* h1_angle_src_vs_10B_cos_xy_plane_mom_cut;
TH1F* h1_wr_first_wr_master;
TH1F* h1_wr_second_wr_master;
//with NEULAND ANALYSIS
TH2F* h2_theta1_theta2_10b_neu;
TH2F* h2_theta1_theta2_10b_lower250_neu;
TH2F* h2_theta1_theta2_10b_lower650_neu;
TH2F* h2_theta1_theta2_10b_higher650_neu;
//--------------------------------------------------
TH2F* h2_reco_mass_angle_comb_cut;
TH1F* h1_reco_mass_11B_p2p_cut_two_iPhos;
TH1F* h1_reco_mass_11B_p2p_cut_iPhos_barrel;
TH1F* h1_reco_mass_11B_p2p_cut_two_barrel;

TH2F* h2_reco_mass_11B_iphos_angle;
TH2F* h2_reco_mass_11B_barrel_angle;

TH1F* h1_mom_abs_10B;
TH1F* h1_mom_x_10B;
TH1F* h1_mom_y_10B;
TH1F* h1_mom_z_10B;
TH2F* h2_sum_mom_vs_cos_p_i_p_n_xy;
TH1F* h1_angle_src_vs_10B_cos_3d;
TH1F* h1_angle_src_vs_10B_cos_3d_paper_cut;
TH2F* h2_e_miss_vs_theta_sum_p2p_11B;
TH2F* h2_e_miss_vs_theta_sum_p2p_10B;
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

sprintf(hist_name, "Highest Gamma Energy E for 11B with multiplicity cut < 3 and angular cut");
h1_gamma_energyE_max_val_11B_m_cut = new TH1F(hist_name,hist_name,200,0,10);
h1_gamma_energyE_max_val_11B_m_cut->GetXaxis()->SetTitle("Energy E [MeV]");
h1_gamma_energyE_max_val_11B_m_cut->GetYaxis()->SetTitle("Number of entries");
h1_gamma_energyE_max_val_11B_m_cut->GetXaxis()->CenterTitle(true);
h1_gamma_energyE_max_val_11B_m_cut->GetYaxis()->CenterTitle(true);
h1_gamma_energyE_max_val_11B_m_cut->GetYaxis()->SetLabelSize(0.045);
h1_gamma_energyE_max_val_11B_m_cut->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Highest Gamma Energy E for 11B with angular cuts applied");
h1_gamma_energyE_max_val_11B_cut = new TH1F(hist_name,hist_name,200,0,10);
h1_gamma_energyE_max_val_11B_cut->GetXaxis()->SetTitle("Energy E [MeV]");
h1_gamma_energyE_max_val_11B_cut->GetYaxis()->SetTitle("Number of entries");
h1_gamma_energyE_max_val_11B_cut->GetXaxis()->CenterTitle(true);
h1_gamma_energyE_max_val_11B_cut->GetYaxis()->CenterTitle(true);
h1_gamma_energyE_max_val_11B_cut->GetYaxis()->SetLabelSize(0.045);
h1_gamma_energyE_max_val_11B_cut->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Highest Gamma Energy E for 11B with angular cuts applied boosting to the 11B frame");
h1_gamma_energyE_max_val_11Bframe_cut = new TH1F(hist_name,hist_name,200,0,10);
h1_gamma_energyE_max_val_11Bframe_cut->GetXaxis()->SetTitle("Energy E [MeV]");
h1_gamma_energyE_max_val_11Bframe_cut->GetYaxis()->SetTitle("Number of entries");
h1_gamma_energyE_max_val_11Bframe_cut->GetXaxis()->CenterTitle(true);
h1_gamma_energyE_max_val_11Bframe_cut->GetYaxis()->CenterTitle(true);
h1_gamma_energyE_max_val_11Bframe_cut->GetYaxis()->SetLabelSize(0.045);
h1_gamma_energyE_max_val_11Bframe_cut->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "SUM Gamma Energy E for 11B ");
h1_gamma_energyE_sum_11B = new TH1F(hist_name,hist_name,200,0,10);
h1_gamma_energyE_sum_11B->GetXaxis()->SetTitle("SUM Energy E [MeV]");
h1_gamma_energyE_sum_11B->GetYaxis()->SetTitle("Number of entries");
h1_gamma_energyE_sum_11B->GetXaxis()->CenterTitle(true);
h1_gamma_energyE_sum_11B->GetYaxis()->CenterTitle(true);
h1_gamma_energyE_sum_11B->GetYaxis()->SetLabelSize(0.045);
h1_gamma_energyE_sum_11B->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "SUM Gamma Energy E for 11B with angular cuts applied");
h1_gamma_energyE_sum_11B_cut = new TH1F(hist_name,hist_name,200,0,10);
h1_gamma_energyE_sum_11B_cut->GetXaxis()->SetTitle("SUM Energy E [MeV]");
h1_gamma_energyE_sum_11B_cut->GetYaxis()->SetTitle("Number of entries");
h1_gamma_energyE_sum_11B_cut->GetXaxis()->CenterTitle(true);
h1_gamma_energyE_sum_11B_cut->GetYaxis()->CenterTitle(true);
h1_gamma_energyE_sum_11B_cut->GetYaxis()->SetLabelSize(0.045);
h1_gamma_energyE_sum_11B_cut->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Highest Gamma Energy E for 10B");
h1_gamma_energyE_max_val_10B = new TH1F(hist_name,hist_name,200,0,10);
h1_gamma_energyE_max_val_10B->GetXaxis()->SetTitle("Energy E [MeV]");
h1_gamma_energyE_max_val_10B->GetYaxis()->SetTitle("Number of entries");
h1_gamma_energyE_max_val_10B->GetXaxis()->CenterTitle(true);
h1_gamma_energyE_max_val_10B->GetYaxis()->CenterTitle(true);
h1_gamma_energyE_max_val_10B->GetYaxis()->SetLabelSize(0.045);
h1_gamma_energyE_max_val_10B->GetYaxis()->SetTitleSize(0.045);

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

if (gamma_given == 1.42942 || gamma_given == 1.41868){

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


//These are the histograms added after the code revision was done...
sprintf(hist_name, "Missing Energy calculated in the 12C rest frame 12C(p,2p)11B");
h1_missing_energy_11B_p2p = new TH1F(hist_name,hist_name,3000,-1500,1500);
h1_missing_energy_11B_p2p->GetXaxis()->SetTitle("E_miss [MeV]");
h1_missing_energy_11B_p2p->GetYaxis()->SetTitle("Counts");
h1_missing_energy_11B_p2p->GetXaxis()->CenterTitle(true);
h1_missing_energy_11B_p2p->GetYaxis()->CenterTitle(true);
h1_missing_energy_11B_p2p->GetYaxis()->SetLabelSize(0.045);
h1_missing_energy_11B_p2p->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "E_miss versus Theta_sum(of the two outgoing protons), 12C(p,2p)11B");
h2_e_miss_vs_theta_sum_p2p_11B = new TH2F(hist_name,hist_name,400,-1000,1000,52,22.15,152.15);
h2_e_miss_vs_theta_sum_p2p_11B->GetXaxis()->SetTitle("E_miss ( = m_p - e_miss)");
h2_e_miss_vs_theta_sum_p2p_11B->GetYaxis()->SetTitle("Theta1 + Theta2 [degr]");
h2_e_miss_vs_theta_sum_p2p_11B->GetXaxis()->CenterTitle(true);
h2_e_miss_vs_theta_sum_p2p_11B->GetYaxis()->CenterTitle(true);
h2_e_miss_vs_theta_sum_p2p_11B->GetYaxis()->SetLabelSize(0.045);
h2_e_miss_vs_theta_sum_p2p_11B->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "E_miss versus Theta_sum(of the two outgoing protons),cut deltaphi > 100 degr 10B");
h2_e_miss_vs_theta_sum_p2p_10B = new TH2F(hist_name,hist_name,400,-1000,1000,52,22.15,152.15);
h2_e_miss_vs_theta_sum_p2p_10B->GetXaxis()->SetTitle("E_miss ( = m_p - e_miss)");
h2_e_miss_vs_theta_sum_p2p_10B->GetYaxis()->SetTitle("Theta1 + Theta2 [degr]");
h2_e_miss_vs_theta_sum_p2p_10B->GetXaxis()->CenterTitle(true);
h2_e_miss_vs_theta_sum_p2p_10B->GetYaxis()->CenterTitle(true);
h2_e_miss_vs_theta_sum_p2p_10B->GetYaxis()->SetLabelSize(0.045);
h2_e_miss_vs_theta_sum_p2p_10B->GetYaxis()->SetTitleSize(0.045);



sprintf(hist_name, "Missing Energy calculated in the 12C rest frame 12C(p,2p)11B with angular cut applied");
h1_missing_energy_11B_p2p_cut = new TH1F(hist_name,hist_name,3000,-1500,1500);
h1_missing_energy_11B_p2p_cut->GetXaxis()->SetTitle("E_miss [MeV]");
h1_missing_energy_11B_p2p_cut->GetYaxis()->SetTitle("Counts");
h1_missing_energy_11B_p2p_cut->GetXaxis()->CenterTitle(true);
h1_missing_energy_11B_p2p_cut->GetYaxis()->CenterTitle(true);
h1_missing_energy_11B_p2p_cut->GetYaxis()->SetLabelSize(0.045);
h1_missing_energy_11B_p2p_cut->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Are the missing energy and Sep. Energy the same 12C(p,2p)11B with angular cut");
h2_missing_e_2meth_11B_p2p_cut = new TH2F(hist_name,hist_name,3000,-1500,1500,3000,-1500,1500);
h2_missing_e_2meth_11B_p2p_cut->GetXaxis()->SetTitle("E_miss [MeV]");
h2_missing_e_2meth_11B_p2p_cut->GetYaxis()->SetTitle("Separation Energy [MeV/c]");
h2_missing_e_2meth_11B_p2p_cut->GetXaxis()->CenterTitle(true);
h2_missing_e_2meth_11B_p2p_cut->GetYaxis()->CenterTitle(true);
h2_missing_e_2meth_11B_p2p_cut->GetYaxis()->SetLabelSize(0.045);
h2_missing_e_2meth_11B_p2p_cut->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Missing Mass = proton inital mass for 12C(p,2p)11B reaction with angular cut applied");
h1_missing_mass_11B_p2p_cut = new TH1F(hist_name,hist_name,3000,-1500,1500);
h1_missing_mass_11B_p2p_cut->GetXaxis()->SetTitle("Reconstructed initial proton mass [MeV]");
h1_missing_mass_11B_p2p_cut->GetYaxis()->SetTitle("Counts");
h1_missing_mass_11B_p2p_cut->GetXaxis()->CenterTitle(true);
h1_missing_mass_11B_p2p_cut->GetYaxis()->CenterTitle(true);
h1_missing_mass_11B_p2p_cut->GetYaxis()->SetLabelSize(0.045);
h1_missing_mass_11B_p2p_cut->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Separation/Excitation Energy of 11B with angular cut applied");
h1_excit_energy_11B_cut = new TH1F(hist_name,hist_name,3000,-1500,1500);
h1_excit_energy_11B_cut->GetXaxis()->SetTitle("Excitation Energy [MeV]");
h1_excit_energy_11B_cut->GetYaxis()->SetTitle("Counts");
h1_excit_energy_11B_cut->GetXaxis()->CenterTitle(true);
h1_excit_energy_11B_cut->GetYaxis()->CenterTitle(true);
h1_excit_energy_11B_cut->GetYaxis()->SetLabelSize(0.045);
h1_excit_energy_11B_cut->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Missing Energy calculated in the 12C rest frame vs (theta1+theta2) 12C(p,2p)11B with angular cut applied");
h2_missing_energy_vs_t1_t2_11B_p2p = new TH2F(hist_name,hist_name,400,-1000,1000,52,22.15,152.15);
h2_missing_energy_vs_t1_t2_11B_p2p->GetXaxis()->SetTitle("E_miss [MeV]");
h2_missing_energy_vs_t1_t2_11B_p2p->GetYaxis()->SetTitle("theta1+theta2 [degr]");
h2_missing_energy_vs_t1_t2_11B_p2p->GetXaxis()->CenterTitle(true);
h2_missing_energy_vs_t1_t2_11B_p2p->GetYaxis()->CenterTitle(true);
h2_missing_energy_vs_t1_t2_11B_p2p->GetYaxis()->SetLabelSize(0.045);
h2_missing_energy_vs_t1_t2_11B_p2p->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "missing momentum p_missing 11B");
h1_p_missing_11B = new TH1F(hist_name,hist_name,750,0,1500);
h1_p_missing_11B->GetXaxis()->SetTitle("p_missing [MeV/c]");
h1_p_missing_11B->GetYaxis()->SetTitle("No. of Events");
h1_p_missing_11B->GetXaxis()->CenterTitle(true);
h1_p_missing_11B->GetYaxis()->CenterTitle(true);
h1_p_missing_11B->GetYaxis()->SetLabelSize(0.045);
h1_p_missing_11B->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "missing momentum p_missing 11B with angular cuts applied");
h1_p_missing_11B_cut = new TH1F(hist_name,hist_name,750,0,1500);
h1_p_missing_11B_cut->GetXaxis()->SetTitle("p_missing [MeV/c]");
h1_p_missing_11B_cut->GetYaxis()->SetTitle("No. of Events");
h1_p_missing_11B_cut->GetXaxis()->CenterTitle(true);
h1_p_missing_11B_cut->GetYaxis()->CenterTitle(true);
h1_p_missing_11B_cut->GetYaxis()->SetLabelSize(0.045);
h1_p_missing_11B_cut->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Missing momentum p_missing vs momentum 11B in 12C frame with angular cuts applied");
h2_p_missing_vs_p_11B_cut = new TH2F(hist_name,hist_name,750,0,1500,750,0,1500);
h2_p_missing_vs_p_11B_cut->GetXaxis()->SetTitle("p_missing [MeV/c]");
h2_p_missing_vs_p_11B_cut->GetYaxis()->SetTitle("momentum 11B in c.m. frame [MeV/c]");
h2_p_missing_vs_p_11B_cut->GetXaxis()->CenterTitle(true);
h2_p_missing_vs_p_11B_cut->GetYaxis()->CenterTitle(true);
h2_p_missing_vs_p_11B_cut->GetYaxis()->SetLabelSize(0.045);
h2_p_missing_vs_p_11B_cut->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Missing momentum p_missing vs momentum 11B in 12C frame with angular cuts applied with fake y_11b");
h2_p_missing_vs_p_11B_cut_fake_y = new TH2F(hist_name,hist_name,750,0,1500,750,0,1500);
h2_p_missing_vs_p_11B_cut_fake_y->GetXaxis()->SetTitle("p_missing [MeV/c]");
h2_p_missing_vs_p_11B_cut_fake_y->GetYaxis()->SetTitle("momentum 11B in c.m. frame [MeV/c]");
h2_p_missing_vs_p_11B_cut_fake_y->GetXaxis()->CenterTitle(true);
h2_p_missing_vs_p_11B_cut_fake_y->GetYaxis()->CenterTitle(true);
h2_p_missing_vs_p_11B_cut_fake_y->GetYaxis()->SetLabelSize(0.045);
h2_p_missing_vs_p_11B_cut_fake_y->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "missing momentum p_missing x component 11B");
h1_p_missing_11B_x = new TH1F(hist_name,hist_name,1000,-1000,1000);
h1_p_missing_11B_x->GetXaxis()->SetTitle("p_missing x component [MeV/c]");
h1_p_missing_11B_x->GetYaxis()->SetTitle("No. of Events");
h1_p_missing_11B_x->GetXaxis()->CenterTitle(true);
h1_p_missing_11B_x->GetYaxis()->CenterTitle(true);
h1_p_missing_11B_x->GetYaxis()->SetLabelSize(0.045);
h1_p_missing_11B_x->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "missing momentum p_missing x component 11B with angular cuts applied");
h1_p_missing_11B_x_cut = new TH1F(hist_name,hist_name,1000,-1000,1000);
h1_p_missing_11B_x_cut->GetXaxis()->SetTitle("p_missing x component [MeV/c]");
h1_p_missing_11B_x_cut->GetYaxis()->SetTitle("No. of Events");
h1_p_missing_11B_x_cut->GetXaxis()->CenterTitle(true);
h1_p_missing_11B_x_cut->GetYaxis()->CenterTitle(true);
h1_p_missing_11B_x_cut->GetYaxis()->SetLabelSize(0.045);
h1_p_missing_11B_x_cut->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "missing momentum p_missing y component 11B with angular cuts applied");
h1_p_missing_11B_y_cut = new TH1F(hist_name,hist_name,1000,-1000,1000);
h1_p_missing_11B_y_cut->GetXaxis()->SetTitle("p_missing y component [MeV/c]");
h1_p_missing_11B_y_cut->GetYaxis()->SetTitle("No. of Events");
h1_p_missing_11B_y_cut->GetXaxis()->CenterTitle(true);
h1_p_missing_11B_y_cut->GetYaxis()->CenterTitle(true);
h1_p_missing_11B_y_cut->GetYaxis()->SetLabelSize(0.045);
h1_p_missing_11B_y_cut->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "missing momentum p_missing z component 11B");
h1_p_missing_11B_z = new TH1F(hist_name,hist_name,1000,-1000,1000);
h1_p_missing_11B_z->GetXaxis()->SetTitle("p_missing z component [MeV/c]");
h1_p_missing_11B_z->GetYaxis()->SetTitle("No. of Events");
h1_p_missing_11B_z->GetXaxis()->CenterTitle(true);
h1_p_missing_11B_z->GetYaxis()->CenterTitle(true);
h1_p_missing_11B_z->GetYaxis()->SetLabelSize(0.045);
h1_p_missing_11B_z->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "missing momentum p_missing z component 11B with angular cuts applied");
h1_p_missing_11B_z_cut = new TH1F(hist_name,hist_name,1000,-1000,1000);
h1_p_missing_11B_z_cut->GetXaxis()->SetTitle("p_missing z component [MeV/c]");
h1_p_missing_11B_z_cut->GetYaxis()->SetTitle("No. of Events");
h1_p_missing_11B_z_cut->GetXaxis()->CenterTitle(true);
h1_p_missing_11B_z_cut->GetYaxis()->CenterTitle(true);
h1_p_missing_11B_z_cut->GetYaxis()->SetLabelSize(0.045);
h1_p_missing_11B_z_cut->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Momentum z-component of 11B with angular cuts applied");
h1_mom11B_z_cut = new TH1F(hist_name,hist_name,1000,-1000,1000);
h1_mom11B_z_cut->GetXaxis()->SetTitle("11B momentum z-component [MeV/c]");
h1_mom11B_z_cut->GetYaxis()->SetTitle("No. of Events");
h1_mom11B_z_cut->GetXaxis()->CenterTitle(true);
h1_mom11B_z_cut->GetYaxis()->CenterTitle(true);
h1_mom11B_z_cut->GetYaxis()->SetLabelSize(0.045);
h1_mom11B_z_cut->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Momentum x-component of 11B with angular cuts applied");
h1_mom11B_x_cut = new TH1F(hist_name,hist_name,1000,-1000,1000);
h1_mom11B_x_cut->GetXaxis()->SetTitle("11B momentum x-component [MeV/c]");
h1_mom11B_x_cut->GetYaxis()->SetTitle("No. of Events");
h1_mom11B_x_cut->GetXaxis()->CenterTitle(true);
h1_mom11B_x_cut->GetYaxis()->CenterTitle(true);
h1_mom11B_x_cut->GetYaxis()->SetLabelSize(0.045);
h1_mom11B_x_cut->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Momentum y-component of 11B with angular cuts applied");
h1_mom11B_y_cut = new TH1F(hist_name,hist_name,1000,-1000,1000);
h1_mom11B_y_cut->GetXaxis()->SetTitle("11B momentum y-component [MeV/c]");
h1_mom11B_y_cut->GetYaxis()->SetTitle("No. of Events");
h1_mom11B_y_cut->GetXaxis()->CenterTitle(true);
h1_mom11B_y_cut->GetYaxis()->CenterTitle(true);
h1_mom11B_y_cut->GetYaxis()->SetLabelSize(0.045);
h1_mom11B_y_cut->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Full Momentum 11B with angular cuts applied");
h1_mom11B_cut = new TH1F(hist_name,hist_name,1000,-1000,1000);
h1_mom11B_cut->GetXaxis()->SetTitle("Full 11B momentum [MeV/c]");
h1_mom11B_cut->GetYaxis()->SetTitle("No. of Events");
h1_mom11B_cut->GetXaxis()->CenterTitle(true);
h1_mom11B_cut->GetYaxis()->CenterTitle(true);
h1_mom11B_cut->GetYaxis()->SetLabelSize(0.045);
h1_mom11B_cut->GetYaxis()->SetTitleSize(0.045);



sprintf(hist_name, "missing momentum p_missing 10B");
h1_p_missing_10B = new TH1F(hist_name,hist_name,750,0,1500);
h1_p_missing_10B->GetXaxis()->SetTitle("p_missing [MeV/c]");
h1_p_missing_10B->GetYaxis()->SetTitle("No. of Events");
h1_p_missing_10B->GetXaxis()->CenterTitle(true);
h1_p_missing_10B->GetYaxis()->CenterTitle(true);
h1_p_missing_10B->GetYaxis()->SetLabelSize(0.045);
h1_p_missing_10B->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Missing Mass in 12C(p,2p)10B reaction");
h1_missing_mass_10B = new TH1F(hist_name,hist_name,2000,-2000,2000);
h1_missing_mass_10B->GetXaxis()->SetTitle("M_miss [MeV]");
h1_missing_mass_10B->GetYaxis()->SetTitle("Counts");
h1_missing_mass_10B->GetXaxis()->CenterTitle(true);
h1_missing_mass_10B->GetYaxis()->CenterTitle(true);
h1_missing_mass_10B->GetYaxis()->SetLabelSize(0.045);
h1_missing_mass_10B->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Missing Mass in 12C(p,2p)11B reaction (should have peak at around 0)");
h1_missing_mass_11B = new TH1F(hist_name,hist_name,2000,-2000,2000);
h1_missing_mass_11B->GetXaxis()->SetTitle("M_miss [MeV]");
h1_missing_mass_11B->GetYaxis()->SetTitle("Counts");
h1_missing_mass_11B->GetXaxis()->CenterTitle(true);
h1_missing_mass_11B->GetYaxis()->CenterTitle(true);
h1_missing_mass_11B->GetYaxis()->SetLabelSize(0.045);
h1_missing_mass_11B->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Reconstructed neutron mass  in 12C(p,2p)10B reaction for p_i momentum > 900");
h1_missing_mass_10B_neutron_miss_cut = new TH1F(hist_name,hist_name,2000,-2000,2000);
h1_missing_mass_10B_neutron_miss_cut->GetXaxis()->SetTitle("M_miss [MeV]");
h1_missing_mass_10B_neutron_miss_cut->GetYaxis()->SetTitle("Counts");
h1_missing_mass_10B_neutron_miss_cut->GetXaxis()->CenterTitle(true);
h1_missing_mass_10B_neutron_miss_cut->GetYaxis()->CenterTitle(true);
h1_missing_mass_10B_neutron_miss_cut->GetYaxis()->SetLabelSize(0.045);
h1_missing_mass_10B_neutron_miss_cut->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Mass in 12C(p,2p)10B reaction reconstructed from p_miss");
h1_missing_mass_10B_rec = new TH1F(hist_name,hist_name,2000,0,2000);
h1_missing_mass_10B_rec->GetXaxis()->SetTitle("M_miss [MeV]");
h1_missing_mass_10B_rec->GetYaxis()->SetTitle("Counts");
h1_missing_mass_10B_rec->GetXaxis()->CenterTitle(true);
h1_missing_mass_10B_rec->GetYaxis()->CenterTitle(true);
h1_missing_mass_10B_rec->GetYaxis()->SetLabelSize(0.045);
h1_missing_mass_10B_rec->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Mass in 12C(p,2p)10B reaction reconstructed from p_miss for momentum of p_miss > 900");
h1_missing_mass_10B_rec_miss_cut = new TH1F(hist_name,hist_name,2000,0,2000);
h1_missing_mass_10B_rec_miss_cut->GetXaxis()->SetTitle("M_miss [MeV]");
h1_missing_mass_10B_rec_miss_cut->GetYaxis()->SetTitle("Counts");
h1_missing_mass_10B_rec_miss_cut->GetXaxis()->CenterTitle(true);
h1_missing_mass_10B_rec_miss_cut->GetYaxis()->CenterTitle(true);
h1_missing_mass_10B_rec_miss_cut->GetYaxis()->SetLabelSize(0.045);
h1_missing_mass_10B_rec_miss_cut->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Missing Momentum for 10B calculated in the 12C rest frame vs (theta1+theta2)");
h2_missing_energy_vs_t1_t2_10B = new TH2F(hist_name,hist_name,750,0,1500,52,22.15,152.15);
h2_missing_energy_vs_t1_t2_10B->GetXaxis()->SetTitle("p_missing [MeV/c]");
h2_missing_energy_vs_t1_t2_10B->GetYaxis()->SetTitle("theta1+theta2 [degr]");
h2_missing_energy_vs_t1_t2_10B->GetXaxis()->CenterTitle(true);
h2_missing_energy_vs_t1_t2_10B->GetYaxis()->CenterTitle(true);
h2_missing_energy_vs_t1_t2_10B->GetYaxis()->SetLabelSize(0.045);
h2_missing_energy_vs_t1_t2_10B->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Angle in the CMS between 11B and p_i(projectile proton) in the xz-plane");
h1_angle_pi_11B = new TH1F(hist_name,hist_name,90,-180,180);
h1_angle_pi_11B->GetXaxis()->SetTitle("Angle [degr]");
h1_angle_pi_11B->GetYaxis()->SetTitle("Counts");
h1_angle_pi_11B->GetXaxis()->CenterTitle(true);
h1_angle_pi_11B->GetYaxis()->CenterTitle(true);
h1_angle_pi_11B->GetYaxis()->SetLabelSize(0.045);
h1_angle_pi_11B->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Cosine of the angle in the CMS between 11B and p_i(projectile proton) in the xz-plane");
h1_angle_pi_11B_cos = new TH1F(hist_name,hist_name,25,-1,1);
h1_angle_pi_11B_cos->GetXaxis()->SetTitle("Cosine");
h1_angle_pi_11B_cos->GetYaxis()->SetTitle("Counts");
h1_angle_pi_11B_cos->GetXaxis()->CenterTitle(true);
h1_angle_pi_11B_cos->GetYaxis()->CenterTitle(true);
h1_angle_pi_11B_cos->GetYaxis()->SetLabelSize(0.045);
h1_angle_pi_11B_cos->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Cosine of the angle in the CMS between 11B and p_i(projectile proton) in 3D with angular cut");
h1_angle_pi_11B_cos_cut = new TH1F(hist_name,hist_name,25,-1,1);
h1_angle_pi_11B_cos_cut->GetXaxis()->SetTitle("Cosine");
h1_angle_pi_11B_cos_cut->GetYaxis()->SetTitle("Counts");
h1_angle_pi_11B_cos_cut->GetXaxis()->CenterTitle(true);
h1_angle_pi_11B_cos_cut->GetYaxis()->CenterTitle(true);
h1_angle_pi_11B_cos_cut->GetYaxis()->SetLabelSize(0.045);
h1_angle_pi_11B_cos_cut->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Cosine of the angle in the CMS between 11B and p_i(projectile proton) in x-y plane with angular cut");
h1_angle_pi_11B_cos_cut_plane = new TH1F(hist_name,hist_name,25,-1,1);
h1_angle_pi_11B_cos_cut_plane->GetXaxis()->SetTitle("Cosine");
h1_angle_pi_11B_cos_cut_plane->GetYaxis()->SetTitle("Counts");
h1_angle_pi_11B_cos_cut_plane->GetXaxis()->CenterTitle(true);
h1_angle_pi_11B_cos_cut_plane->GetYaxis()->CenterTitle(true);
h1_angle_pi_11B_cos_cut_plane->GetYaxis()->SetLabelSize(0.045);
h1_angle_pi_11B_cos_cut_plane->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Cosine of the angle in the  12C CMS between the initial proton and the reconstructed neutron in the x-y plane");
h1_angle_pi_p_n_xy_plane = new TH1F(hist_name,hist_name,25,-1,1);
h1_angle_pi_p_n_xy_plane->GetXaxis()->SetTitle("Cosine");
h1_angle_pi_p_n_xy_plane->GetYaxis()->SetTitle("Counts");
h1_angle_pi_p_n_xy_plane->GetXaxis()->CenterTitle(true);
h1_angle_pi_p_n_xy_plane->GetYaxis()->CenterTitle(true);
h1_angle_pi_p_n_xy_plane->GetYaxis()->SetLabelSize(0.045);
h1_angle_pi_p_n_xy_plane->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Cosine of the angle in the  12C CMS between the initial proton and the reconstructed neutron in the x-y plane for p_i momentum > 650 MeV/c");
h1_angle_pi_p_n_xy_plane_mom_cut = new TH1F(hist_name,hist_name,25,-1,1);
h1_angle_pi_p_n_xy_plane_mom_cut->GetXaxis()->SetTitle("Cosine");
h1_angle_pi_p_n_xy_plane_mom_cut->GetYaxis()->SetTitle("Counts");
h1_angle_pi_p_n_xy_plane_mom_cut->GetXaxis()->CenterTitle(true);
h1_angle_pi_p_n_xy_plane_mom_cut->GetYaxis()->CenterTitle(true);
h1_angle_pi_p_n_xy_plane_mom_cut->GetYaxis()->SetLabelSize(0.045);
h1_angle_pi_p_n_xy_plane_mom_cut->GetYaxis()->SetTitleSize(0.045);



sprintf(hist_name, "Cosine of the angle in the CMS between 11B and p_i(projectile proton) in x-y plane with restrictive cuts");
h1_angle_pi_11B_cos_hcut_plane = new TH1F(hist_name,hist_name,25,-1,1);
h1_angle_pi_11B_cos_hcut_plane->GetXaxis()->SetTitle("Cosine");
h1_angle_pi_11B_cos_hcut_plane->GetYaxis()->SetTitle("Counts");
h1_angle_pi_11B_cos_hcut_plane->GetXaxis()->CenterTitle(true);
h1_angle_pi_11B_cos_hcut_plane->GetYaxis()->CenterTitle(true);
h1_angle_pi_11B_cos_hcut_plane->GetYaxis()->SetLabelSize(0.045);
h1_angle_pi_11B_cos_hcut_plane->GetYaxis()->SetTitleSize(0.045);



sprintf(hist_name, "Cosine of the angle in the CMS between 11B and p_i(projectile proton) in 3D with angular cut, and fake y");
h1_angle_pi_11B_cos_cut_fake_y = new TH1F(hist_name,hist_name,25,-1,1);
h1_angle_pi_11B_cos_cut_fake_y->GetXaxis()->SetTitle("Cosine");
h1_angle_pi_11B_cos_cut_fake_y->GetYaxis()->SetTitle("Counts");
h1_angle_pi_11B_cos_cut_fake_y->GetXaxis()->CenterTitle(true);
h1_angle_pi_11B_cos_cut_fake_y->GetYaxis()->CenterTitle(true);
h1_angle_pi_11B_cos_cut_fake_y->GetYaxis()->SetLabelSize(0.045);
h1_angle_pi_11B_cos_cut_fake_y->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Angle in the CMS between p_i(projectile proton) and p_n in the xz-plane");
h1_angle_pi_pn_10B = new TH1F(hist_name,hist_name,90,-180,180);
h1_angle_pi_pn_10B->GetXaxis()->SetTitle("Angle [degr]");
h1_angle_pi_pn_10B->GetYaxis()->SetTitle("Counts");
h1_angle_pi_pn_10B->GetXaxis()->CenterTitle(true);
h1_angle_pi_pn_10B->GetYaxis()->CenterTitle(true);
h1_angle_pi_pn_10B->GetYaxis()->SetLabelSize(0.045);
h1_angle_pi_pn_10B->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Cosine of the angle in the CMS between the momentum of SRC proton and the SRC neutron in the xz-plane");
h1_angle_pi_pn_10B_cos = new TH1F(hist_name,hist_name,25,-1,1);
h1_angle_pi_pn_10B_cos->GetXaxis()->SetTitle("Cosine");
h1_angle_pi_pn_10B_cos->GetYaxis()->SetTitle("Counts");
h1_angle_pi_pn_10B_cos->GetXaxis()->CenterTitle(true);
h1_angle_pi_pn_10B_cos->GetYaxis()->CenterTitle(true);
h1_angle_pi_pn_10B_cos->GetYaxis()->SetLabelSize(0.045);
h1_angle_pi_pn_10B_cos->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Cosine of the angle in the CMS between the SRC pair and the 10B fragment  in the xz-plane");
h1_angle_src_vs_10B_cos = new TH1F(hist_name,hist_name,25,-1,1);
h1_angle_src_vs_10B_cos->GetXaxis()->SetTitle("Cosine");
h1_angle_src_vs_10B_cos->GetYaxis()->SetTitle("Counts");
h1_angle_src_vs_10B_cos->GetXaxis()->CenterTitle(true);
h1_angle_src_vs_10B_cos->GetYaxis()->CenterTitle(true);
h1_angle_src_vs_10B_cos->GetYaxis()->SetLabelSize(0.045);
h1_angle_src_vs_10B_cos->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Cosine of the angle in the CMS between the SRC pair and the 10B fragment  in the xy-plane");
h1_angle_src_vs_10B_cos_xy_plane = new TH1F(hist_name,hist_name,25,-1,1);
h1_angle_src_vs_10B_cos_xy_plane->GetXaxis()->SetTitle("Cosine");
h1_angle_src_vs_10B_cos_xy_plane->GetYaxis()->SetTitle("Counts");
h1_angle_src_vs_10B_cos_xy_plane->GetXaxis()->CenterTitle(true);
h1_angle_src_vs_10B_cos_xy_plane->GetYaxis()->CenterTitle(true);
h1_angle_src_vs_10B_cos_xy_plane->GetYaxis()->SetLabelSize(0.045);
h1_angle_src_vs_10B_cos_xy_plane->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Cosine of the angle in the CMS between the SRC pair and the 10B fragment  in three dimensions");
h1_angle_src_vs_10B_cos_3d = new TH1F(hist_name,hist_name,25,-1,1);
h1_angle_src_vs_10B_cos_3d->GetXaxis()->SetTitle("Cosine");
h1_angle_src_vs_10B_cos_3d->GetYaxis()->SetTitle("Counts");
h1_angle_src_vs_10B_cos_3d->GetXaxis()->CenterTitle(true);
h1_angle_src_vs_10B_cos_3d->GetYaxis()->CenterTitle(true);
h1_angle_src_vs_10B_cos_3d->GetYaxis()->SetLabelSize(0.045);
h1_angle_src_vs_10B_cos_3d->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Cosine of the angle in the CMS between the SRC pair and the 10B fragment  in three dimensions paper cut");
h1_angle_src_vs_10B_cos_3d_paper_cut = new TH1F(hist_name,hist_name,25,-1,1);
h1_angle_src_vs_10B_cos_3d_paper_cut->GetXaxis()->SetTitle("Cosine");
h1_angle_src_vs_10B_cos_3d_paper_cut->GetYaxis()->SetTitle("Counts");
h1_angle_src_vs_10B_cos_3d_paper_cut->GetXaxis()->CenterTitle(true);
h1_angle_src_vs_10B_cos_3d_paper_cut->GetYaxis()->CenterTitle(true);
h1_angle_src_vs_10B_cos_3d_paper_cut->GetYaxis()->SetLabelSize(0.045);
h1_angle_src_vs_10B_cos_3d_paper_cut->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Cosine of the angle in the CMS between the SRC pair and the 10B fragment  in the xy-plane with p_i momentum > 650 MeV/c");
h1_angle_src_vs_10B_cos_xy_plane_mom_cut = new TH1F(hist_name,hist_name,25,-1,1);
h1_angle_src_vs_10B_cos_xy_plane_mom_cut->GetXaxis()->SetTitle("Cosine");
h1_angle_src_vs_10B_cos_xy_plane_mom_cut->GetYaxis()->SetTitle("Counts");
h1_angle_src_vs_10B_cos_xy_plane_mom_cut->GetXaxis()->CenterTitle(true);
h1_angle_src_vs_10B_cos_xy_plane_mom_cut->GetYaxis()->CenterTitle(true);
h1_angle_src_vs_10B_cos_xy_plane_mom_cut->GetYaxis()->SetLabelSize(0.045);
h1_angle_src_vs_10B_cos_xy_plane_mom_cut->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Highest Gamma Energy E for 11B in lab. system");
h1_gamma_energyE_max_val_11B_lab = new TH1F(hist_name,hist_name,200,0,10);
h1_gamma_energyE_max_val_11B_lab->GetXaxis()->SetTitle("Energy E [MeV]");
h1_gamma_energyE_max_val_11B_lab->GetYaxis()->SetTitle("Number of entries");
h1_gamma_energyE_max_val_11B_lab->GetXaxis()->CenterTitle(true);
h1_gamma_energyE_max_val_11B_lab->GetYaxis()->CenterTitle(true);
h1_gamma_energyE_max_val_11B_lab->GetYaxis()->SetLabelSize(0.045);
h1_gamma_energyE_max_val_11B_lab->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Highest Gamma Energy E for 11B in cm frame vs multiplicity for gammas E_lab < 10 MeV");
h2_gamma_energy_11B_vs_mult = new TH2F(hist_name,hist_name,200,0,10,20,0,20);
h2_gamma_energy_11B_vs_mult->GetXaxis()->SetTitle("Energy E [MeV]");
h2_gamma_energy_11B_vs_mult->GetYaxis()->SetTitle("theta1+theta2 [degr]");
h2_gamma_energy_11B_vs_mult->GetXaxis()->CenterTitle(true);
h2_gamma_energy_11B_vs_mult->GetYaxis()->CenterTitle(true);
h2_gamma_energy_11B_vs_mult->GetYaxis()->SetLabelSize(0.045);
h2_gamma_energy_11B_vs_mult->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Highest Gamma Energy E for 11B in cm frame vs multiplicity for gammas E_lab < 10 MeV with angular cut");
h2_gamma_energy_11B_vs_mult_cut = new TH2F(hist_name,hist_name,200,0,10,20,0,20);
h2_gamma_energy_11B_vs_mult_cut->GetXaxis()->SetTitle("Energy E [MeV]");
h2_gamma_energy_11B_vs_mult_cut->GetYaxis()->SetTitle("theta1+theta2 [degr]");
h2_gamma_energy_11B_vs_mult_cut->GetXaxis()->CenterTitle(true);
h2_gamma_energy_11B_vs_mult_cut->GetYaxis()->CenterTitle(true);
h2_gamma_energy_11B_vs_mult_cut->GetYaxis()->SetLabelSize(0.045);
h2_gamma_energy_11B_vs_mult_cut->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Highest Gamma Energy E for 11B in cm frame vs opening angle");
h2_gamma_energy_11_vs_angle = new TH2F(hist_name,hist_name,200,0,10,25,22.15,84.65);
h2_gamma_energy_11_vs_angle->GetXaxis()->SetTitle("Energy E [MeV]");
h2_gamma_energy_11_vs_angle->GetYaxis()->SetTitle("Escaping angle [degrees]");
h2_gamma_energy_11_vs_angle->GetXaxis()->CenterTitle(true);
h2_gamma_energy_11_vs_angle->GetYaxis()->CenterTitle(true);
h2_gamma_energy_11_vs_angle->GetYaxis()->SetLabelSize(0.045);
h2_gamma_energy_11_vs_angle->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Theta1 vs Theta2 for the 10B reaction");
h2_theta1_theta2_10b = new TH2F(hist_name,hist_name,25,22.15,84.65,25,22.15,84.65);
h2_theta1_theta2_10b->GetXaxis()->SetTitle("Theta1 [degrees]");
h2_theta1_theta2_10b->GetYaxis()->SetTitle("Theta2 [degrees]");
h2_theta1_theta2_10b->GetXaxis()->CenterTitle(true);
h2_theta1_theta2_10b->GetYaxis()->CenterTitle(true);
h2_theta1_theta2_10b->GetYaxis()->SetLabelSize(0.045);
h2_theta1_theta2_10b->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Theta1 vs Theta2 for the 10B reaction for p_missing < 250 MeV/c");
h2_theta1_theta2_10b_lower250 = new TH2F(hist_name,hist_name,25,22.15,84.65,25,22.15,84.65);
h2_theta1_theta2_10b_lower250->GetXaxis()->SetTitle("Theta1 [degrees]");
h2_theta1_theta2_10b_lower250->GetYaxis()->SetTitle("Theta2 [degrees]");
h2_theta1_theta2_10b_lower250->GetXaxis()->CenterTitle(true);
h2_theta1_theta2_10b_lower250->GetYaxis()->CenterTitle(true);
h2_theta1_theta2_10b_lower250->GetYaxis()->SetLabelSize(0.045);
h2_theta1_theta2_10b_lower250->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Theta1 vs Theta2 for the 10B reaction for 250 <  p_missing < 650 MeV/c");
h2_theta1_theta2_10b_lower650 = new TH2F(hist_name,hist_name,25,22.15,84.65,25,22.15,84.65);
h2_theta1_theta2_10b_lower650->GetXaxis()->SetTitle("Theta1 [degrees]");
h2_theta1_theta2_10b_lower650->GetYaxis()->SetTitle("Theta2 [degrees]");
h2_theta1_theta2_10b_lower650->GetXaxis()->CenterTitle(true);
h2_theta1_theta2_10b_lower650->GetYaxis()->CenterTitle(true);
h2_theta1_theta2_10b_lower650->GetYaxis()->SetLabelSize(0.045);
h2_theta1_theta2_10b_lower650->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Theta1 vs Theta2 for the 10B reaction for p_missing > 650 MeV/c");
h2_theta1_theta2_10b_higher650 = new TH2F(hist_name,hist_name,25,22.15,84.65,25,22.15,84.65);
h2_theta1_theta2_10b_higher650->GetXaxis()->SetTitle("Theta1 [degrees]");
h2_theta1_theta2_10b_higher650->GetYaxis()->SetTitle("Theta2 [degrees]");
h2_theta1_theta2_10b_higher650->GetXaxis()->CenterTitle(true);
h2_theta1_theta2_10b_higher650->GetYaxis()->CenterTitle(true);
h2_theta1_theta2_10b_higher650->GetYaxis()->SetLabelSize(0.045);
h2_theta1_theta2_10b_higher650->GetYaxis()->SetTitleSize(0.045);

//analyse if hit in NEULAND seeen:

sprintf(hist_name, "Theta1 vs Theta2 for the 10B reactioni with entry in NEULAND");
h2_theta1_theta2_10b_neu = new TH2F(hist_name,hist_name,25,22.15,84.65,25,22.15,84.65);
h2_theta1_theta2_10b_neu->GetXaxis()->SetTitle("Theta1 [degrees]");
h2_theta1_theta2_10b_neu->GetYaxis()->SetTitle("Theta2 [degrees]");
h2_theta1_theta2_10b_neu->GetXaxis()->CenterTitle(true);
h2_theta1_theta2_10b_neu->GetYaxis()->CenterTitle(true);
h2_theta1_theta2_10b_neu->GetYaxis()->SetLabelSize(0.045);
h2_theta1_theta2_10b_neu->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Theta1 vs Theta2 for the 10B reaction for p_missing < 250 MeV/c with entry in NEULAND");
h2_theta1_theta2_10b_lower250_neu = new TH2F(hist_name,hist_name,25,22.15,84.65,25,22.15,84.65);
h2_theta1_theta2_10b_lower250_neu->GetXaxis()->SetTitle("Theta1 [degrees]");
h2_theta1_theta2_10b_lower250_neu->GetYaxis()->SetTitle("Theta2 [degrees]");
h2_theta1_theta2_10b_lower250_neu->GetXaxis()->CenterTitle(true);
h2_theta1_theta2_10b_lower250_neu->GetYaxis()->CenterTitle(true);
h2_theta1_theta2_10b_lower250_neu->GetYaxis()->SetLabelSize(0.045);
h2_theta1_theta2_10b_lower250_neu->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Theta1 vs Theta2 for the 10B reaction for 250 <  p_missing < 650 MeV/c with entry in NEULAND");
h2_theta1_theta2_10b_lower650_neu = new TH2F(hist_name,hist_name,25,22.15,84.65,25,22.15,84.65);
h2_theta1_theta2_10b_lower650_neu->GetXaxis()->SetTitle("Theta1 [degrees]");
h2_theta1_theta2_10b_lower650_neu->GetYaxis()->SetTitle("Theta2 [degrees]");
h2_theta1_theta2_10b_lower650_neu->GetXaxis()->CenterTitle(true);
h2_theta1_theta2_10b_lower650_neu->GetYaxis()->CenterTitle(true);
h2_theta1_theta2_10b_lower650_neu->GetYaxis()->SetLabelSize(0.045);
h2_theta1_theta2_10b_lower650_neu->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Theta1 vs Theta2 for the 10B reaction for p_missing > 650 MeV/c with entry in NEULAND");
h2_theta1_theta2_10b_higher650_neu = new TH2F(hist_name,hist_name,25,22.15,84.65,25,22.15,84.65);
h2_theta1_theta2_10b_higher650_neu->GetXaxis()->SetTitle("Theta1 [degrees]");
h2_theta1_theta2_10b_higher650_neu->GetYaxis()->SetTitle("Theta2 [degrees]");
h2_theta1_theta2_10b_higher650_neu->GetXaxis()->CenterTitle(true);
h2_theta1_theta2_10b_higher650_neu->GetYaxis()->CenterTitle(true);
h2_theta1_theta2_10b_higher650_neu->GetYaxis()->SetLabelSize(0.045);
h2_theta1_theta2_10b_higher650_neu->GetYaxis()->SetTitleSize(0.045);
//------------------------------

sprintf(hist_name, "Theta1 vs Theta2 for the 10B reaction with reconstructed neutron mass cut 850 < n_M < 1100 ");
h2_theta1_theta2_10b_m_cut = new TH2F(hist_name,hist_name,25,22.15,84.65,25,22.15,84.65);
h2_theta1_theta2_10b_m_cut->GetXaxis()->SetTitle("Theta1 [degrees]");
h2_theta1_theta2_10b_m_cut->GetYaxis()->SetTitle("Theta2 [degrees]");
h2_theta1_theta2_10b_m_cut->GetXaxis()->CenterTitle(true);
h2_theta1_theta2_10b_m_cut->GetYaxis()->CenterTitle(true);
h2_theta1_theta2_10b_m_cut->GetYaxis()->SetLabelSize(0.045);
h2_theta1_theta2_10b_m_cut->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "missing momentum p_missing 10B reaction with reconstructed neutron mass cut 850 < n_M < 1100");
h1_p_missing_10B_m_cut = new TH1F(hist_name,hist_name,750,0,1500);
h1_p_missing_10B_m_cut->GetXaxis()->SetTitle("p_missing [MeV/c]");
h1_p_missing_10B_m_cut->GetYaxis()->SetTitle("No. of Events");
h1_p_missing_10B_m_cut->GetXaxis()->CenterTitle(true);
h1_p_missing_10B_m_cut->GetYaxis()->CenterTitle(true);
h1_p_missing_10B_m_cut->GetYaxis()->SetLabelSize(0.045);
h1_p_missing_10B_m_cut->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "missing momentum p_missing 10B reaction with reconstructed initial proton  830 < m_p_i < 970");
h1_p_missing_10B_p_i_cut = new TH1F(hist_name,hist_name,750,0,1500);
h1_p_missing_10B_p_i_cut->GetXaxis()->SetTitle("p_missing [MeV/c]");
h1_p_missing_10B_p_i_cut->GetYaxis()->SetTitle("No. of Events");
h1_p_missing_10B_p_i_cut->GetXaxis()->CenterTitle(true);
h1_p_missing_10B_p_i_cut->GetYaxis()->CenterTitle(true);
h1_p_missing_10B_p_i_cut->GetYaxis()->SetLabelSize(0.045);
h1_p_missing_10B_p_i_cut->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Theta1 vs Theta2 for the 10B reaction with reconstructed neutron mass cut 850 < n_M < 1100 and proton initial mass 844 < p_M < 960");
h2_theta1_theta2_10b_double_cut = new TH2F(hist_name,hist_name,25,22.15,84.65,25,22.15,84.65);
h2_theta1_theta2_10b_double_cut->GetXaxis()->SetTitle("Theta1 [degrees]");
h2_theta1_theta2_10b_double_cut->GetYaxis()->SetTitle("Theta2 [degrees]");
h2_theta1_theta2_10b_double_cut->GetXaxis()->CenterTitle(true);
h2_theta1_theta2_10b_double_cut->GetYaxis()->CenterTitle(true);
h2_theta1_theta2_10b_double_cut->GetYaxis()->SetLabelSize(0.045);
h2_theta1_theta2_10b_double_cut->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Theta1 vs Theta2 for the 10B reaction with cut on initial proton mom > 900 MeV/c");
h2_theta1_theta2_10b_miss_mom_cut = new TH2F(hist_name,hist_name,25,22.15,84.65,25,22.15,84.65);
h2_theta1_theta2_10b_miss_mom_cut->GetXaxis()->SetTitle("Theta1 [degrees]");
h2_theta1_theta2_10b_miss_mom_cut->GetYaxis()->SetTitle("Theta2 [degrees]");
h2_theta1_theta2_10b_miss_mom_cut->GetXaxis()->CenterTitle(true);
h2_theta1_theta2_10b_miss_mom_cut->GetYaxis()->CenterTitle(true);
h2_theta1_theta2_10b_miss_mom_cut->GetYaxis()->SetLabelSize(0.045);
h2_theta1_theta2_10b_miss_mom_cut->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Correlation of the cartesian components of the internal momentum of the initial proton and fragment with angular cut");
h2_corr_y_proton_fragm = new TH2F(hist_name,hist_name,1000,-500,500,1000,-500,500);
h2_corr_y_proton_fragm->GetXaxis()->SetTitle("P_y fragment [MeV/c]");
h2_corr_y_proton_fragm->GetYaxis()->SetTitle("P_y proton [MeV/c]");
h2_corr_y_proton_fragm->GetXaxis()->CenterTitle(true);
h2_corr_y_proton_fragm->GetYaxis()->CenterTitle(true);
h2_corr_y_proton_fragm->GetYaxis()->SetLabelSize(0.045);
h2_corr_y_proton_fragm->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Correlation of the cartesian components of the internal momentum of the initial proton and fragment with angular cut cos < 0");
h2_corr_y_proton_fragm_cos_less = new TH2F(hist_name,hist_name,1000,-500,500,1000,-500,500);
h2_corr_y_proton_fragm_cos_less->GetXaxis()->SetTitle("P_y fragment [MeV/c]");
h2_corr_y_proton_fragm_cos_less->GetYaxis()->SetTitle("P_y proton [MeV/c]");
h2_corr_y_proton_fragm_cos_less->GetXaxis()->CenterTitle(true);
h2_corr_y_proton_fragm_cos_less->GetYaxis()->CenterTitle(true);
h2_corr_y_proton_fragm_cos_less->GetYaxis()->SetLabelSize(0.045);
h2_corr_y_proton_fragm_cos_less->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Correlation of the cartesian components of the internal momentum of the initial proton and fragment with angular cut cos > 0"); 
h2_corr_y_proton_fragm_cos_more = new TH2F(hist_name,hist_name,1000,-500,500,1000,-500,500);
h2_corr_y_proton_fragm_cos_more->GetXaxis()->SetTitle("P_y fragment [MeV/c]");
h2_corr_y_proton_fragm_cos_more->GetYaxis()->SetTitle("P_y proton [MeV/c]");
h2_corr_y_proton_fragm_cos_more->GetXaxis()->CenterTitle(true);
h2_corr_y_proton_fragm_cos_more->GetYaxis()->CenterTitle(true);
h2_corr_y_proton_fragm_cos_more->GetYaxis()->SetLabelSize(0.045);
h2_corr_y_proton_fragm_cos_more->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Correlation of the cartesian components of the internal momentum of the initial proton and fragment with angular cut in x");
h2_corr_x_proton_fragm = new TH2F(hist_name,hist_name,1000,-500,500,1000,-500,500);
h2_corr_x_proton_fragm->GetXaxis()->SetTitle("P_x fragment [MeV/c]");
h2_corr_x_proton_fragm->GetYaxis()->SetTitle("P_x proton [MeV/c]");
h2_corr_x_proton_fragm->GetXaxis()->CenterTitle(true);
h2_corr_x_proton_fragm->GetYaxis()->CenterTitle(true);
h2_corr_x_proton_fragm->GetYaxis()->SetLabelSize(0.045);
h2_corr_x_proton_fragm->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Correlation of the cartesian components of the internal momentum of the initial proton and fragment with angular cut using own formula");
h2_corr_y_proton_fragm_tj = new TH2F(hist_name,hist_name,1000,-500,500,1000,-500,500);
h2_corr_y_proton_fragm_tj->GetXaxis()->SetTitle("P_y fragment [MeV/c]");
h2_corr_y_proton_fragm_tj->GetYaxis()->SetTitle("P_y proton [MeV/c]");
h2_corr_y_proton_fragm_tj->GetXaxis()->CenterTitle(true);
h2_corr_y_proton_fragm_tj->GetYaxis()->CenterTitle(true);
h2_corr_y_proton_fragm_tj->GetYaxis()->SetLabelSize(0.045);
h2_corr_y_proton_fragm_tj->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Correlation (MWPC2.Y-MWPC1.Y) vs (MWPC3.Y-MWPC1.Y)");
h2_corr_y_mws = new TH2F(hist_name,hist_name,1000,-500,500,1000,-500,500);
h2_corr_y_mws->GetXaxis()->SetTitle("(MWPC2.Y-MWPC1.Y)[mm]");
h2_corr_y_mws->GetYaxis()->SetTitle("(MWPC2.Y-MWPC1.Y)[mm]");
h2_corr_y_mws->GetXaxis()->CenterTitle(true);
h2_corr_y_mws->GetYaxis()->CenterTitle(true);
h2_corr_y_mws->GetYaxis()->SetLabelSize(0.045);
h2_corr_y_mws->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Correlation (MWPC2.Y-MWPC1.Y) vs (MWPC3.Y-MWPC1.Y) with reversed MWPC3.Y");
h2_corr_y_mws_minusmw3 = new TH2F(hist_name,hist_name,1000,-500,500,1000,-500,500);
h2_corr_y_mws_minusmw3->GetXaxis()->SetTitle("(MWPC2.Y-MWPC1.Y)[mm]");
h2_corr_y_mws_minusmw3->GetYaxis()->SetTitle("(MWPC2.Y-MWPC1.Y)[mm]");
h2_corr_y_mws_minusmw3->GetXaxis()->CenterTitle(true);
h2_corr_y_mws_minusmw3->GetYaxis()->CenterTitle(true);
h2_corr_y_mws_minusmw3->GetYaxis()->SetLabelSize(0.045);
h2_corr_y_mws_minusmw3->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "E1_kin + E2_kin (of proton 1 and proton 2 in 12cFrame) vs. Theta_sum (in lab frame)");
h2_E_sum_vs_theta_sum = new TH2F(hist_name,hist_name,1000,-1000,1000,1000,-500,500);
h2_E_sum_vs_theta_sum->GetXaxis()->SetTitle("E1 + E2 [MeV]");
h2_E_sum_vs_theta_sum->GetYaxis()->SetTitle("Theta_1 + Theta_2 [degrees]");
h2_E_sum_vs_theta_sum->GetXaxis()->CenterTitle(true);
h2_E_sum_vs_theta_sum->GetYaxis()->CenterTitle(true);
h2_E_sum_vs_theta_sum->GetYaxis()->SetLabelSize(0.045);
h2_E_sum_vs_theta_sum->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "WRCALIFA ts first entry minus WRMaster ts");
h1_wr_first_wr_master = new TH1F(hist_name,hist_name,1000,-5000,5000);
h1_wr_first_wr_master->GetXaxis()->SetTitle("Time Difference [ns]");
h1_wr_first_wr_master->GetYaxis()->SetTitle("Counts");
h1_wr_first_wr_master->GetXaxis()->CenterTitle(true);
h1_wr_first_wr_master->GetYaxis()->CenterTitle(true);
h1_wr_first_wr_master->GetYaxis()->SetLabelSize(0.045);
h1_wr_first_wr_master->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "WRCALIFA ts second entry minus WRMaster ts");
h1_wr_second_wr_master = new TH1F(hist_name,hist_name,1000,-5000,5000);
h1_wr_second_wr_master->GetXaxis()->SetTitle("Time Difference [ns]");
h1_wr_second_wr_master->GetYaxis()->SetTitle("Counts");
h1_wr_second_wr_master->GetXaxis()->CenterTitle(true);
h1_wr_second_wr_master->GetYaxis()->CenterTitle(true);
h1_wr_second_wr_master->GetYaxis()->SetLabelSize(0.045);
h1_wr_second_wr_master->GetYaxis()->SetTitleSize(0.045);

//N_f vs N_s analysis

sprintf(hist_name, "N_f vs N_s  for the 10B reaction (filling behaviour of both light particles)");
h2_n_f_n_s_10b = new TH2F(hist_name,hist_name,250,0,250,250,0,250);
h2_n_f_n_s_10b->GetXaxis()->SetTitle("N_f [MeV]");
h2_n_f_n_s_10b->GetYaxis()->SetTitle("N_s [MeV]");
h2_n_f_n_s_10b->GetXaxis()->CenterTitle(true);
h2_n_f_n_s_10b->GetYaxis()->CenterTitle(true);
h2_n_f_n_s_10b->GetYaxis()->SetLabelSize(0.045);
h2_n_f_n_s_10b->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "N_f vs N_s  for the 10B reaction (filling behaviour of both light particles) for p_missing < 250 MeV/c");
h2_n_f_n_s_10b_lower250 = new TH2F(hist_name,hist_name,250,0,250,250,0,250);
h2_n_f_n_s_10b_lower250->GetXaxis()->SetTitle("N_f [MeV]");
h2_n_f_n_s_10b_lower250->GetYaxis()->SetTitle("N_s [MeV]");
h2_n_f_n_s_10b_lower250->GetXaxis()->CenterTitle(true);
h2_n_f_n_s_10b_lower250->GetYaxis()->CenterTitle(true);
h2_n_f_n_s_10b_lower250->GetYaxis()->SetLabelSize(0.045);
h2_n_f_n_s_10b_lower250->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "N_f vs N_s  for the 10B reaction (filling behaviour of both light particles) for 250 <  p_missing < 650 MeV/c");
h2_n_f_n_s_10b_lower650 = new TH2F(hist_name,hist_name,250,0,250,250,0,250);
h2_n_f_n_s_10b_lower650->GetXaxis()->SetTitle("N_f [MeV]");
h2_n_f_n_s_10b_lower650->GetYaxis()->SetTitle("N_s [MeV]");
h2_n_f_n_s_10b_lower650->GetXaxis()->CenterTitle(true);
h2_n_f_n_s_10b_lower650->GetYaxis()->CenterTitle(true);
h2_n_f_n_s_10b_lower650->GetYaxis()->SetLabelSize(0.045);
h2_n_f_n_s_10b_lower650->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "N_f vs N_s  for the 10B reaction (filling behaviour of both light particles) for p_missing > 650 MeV/c");
h2_n_f_n_s_10b_higher650 = new TH2F(hist_name,hist_name,250,0,250,250,0,250);
h2_n_f_n_s_10b_higher650->GetXaxis()->SetTitle("N_f [MeV]");
h2_n_f_n_s_10b_higher650->GetYaxis()->SetTitle("N_s [MeV]");
h2_n_f_n_s_10b_higher650->GetXaxis()->CenterTitle(true);
h2_n_f_n_s_10b_higher650->GetYaxis()->CenterTitle(true);
h2_n_f_n_s_10b_higher650->GetYaxis()->SetLabelSize(0.045);
h2_n_f_n_s_10b_higher650->GetYaxis()->SetTitleSize(0.045);


//-----------------------
const Int_t comb_nr = 3;
const char *angle_comb[comb_nr] = {"iPhos + iPhos" , "iPhos + Barrel", "Barrel + Barrel"};
sprintf(hist_name, "Reconstructed mass of the 12C knocked out proton with angular cuts applied vs opening angles of the two outgoing protons");
h2_reco_mass_angle_comb_cut = new TH2F(hist_name,hist_name,3000,-1500,1500,3,0,3);
h2_reco_mass_angle_comb_cut->GetXaxis()->SetTitle("Reconstructed Proton Mass [MeV]");
h2_reco_mass_angle_comb_cut->GetYaxis()->SetTitle("Opening angle combination of two outgoing Protons");
h2_reco_mass_angle_comb_cut->GetXaxis()->CenterTitle(true);
h2_reco_mass_angle_comb_cut->GetYaxis()->CenterTitle(true);
h2_reco_mass_angle_comb_cut->GetYaxis()->SetLabelSize(0.045);
h2_reco_mass_angle_comb_cut->GetYaxis()->SetTitleSize(0.045);
//h2_reco_mass_angle_comb_cut->SetCanExtend(TH1::kAllAxes);
//h2_reco_mass_angle_comb_cut->LabelsOption("a");

h2_reco_mass_angle_comb_cut->GetYaxis()->SetBinLabel(1,angle_comb[0]);
h2_reco_mass_angle_comb_cut->GetYaxis()->SetBinLabel(2,angle_comb[1]);
h2_reco_mass_angle_comb_cut->GetYaxis()->SetBinLabel(3,angle_comb[2]);


sprintf(hist_name, "Reconstructed mass of the 12C knocked out proton with angular cuts for both protons in iPhos");
h1_reco_mass_11B_p2p_cut_two_iPhos = new TH1F(hist_name,hist_name,3000,-1500,1500);
h1_reco_mass_11B_p2p_cut_two_iPhos->GetXaxis()->SetTitle("Reconstructed initial proton mass [MeV]");
h1_reco_mass_11B_p2p_cut_two_iPhos->GetYaxis()->SetTitle("Counts");
h1_reco_mass_11B_p2p_cut_two_iPhos->GetXaxis()->CenterTitle(true);
h1_reco_mass_11B_p2p_cut_two_iPhos->GetYaxis()->CenterTitle(true);
h1_reco_mass_11B_p2p_cut_two_iPhos->GetYaxis()->SetLabelSize(0.045);
h1_reco_mass_11B_p2p_cut_two_iPhos->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Reconstructed mass of the 12C knocked out proton with angular cuts for protons in iPhos and Barrel");
h1_reco_mass_11B_p2p_cut_iPhos_barrel = new TH1F(hist_name,hist_name,3000,-1500,1500);
h1_reco_mass_11B_p2p_cut_iPhos_barrel->GetXaxis()->SetTitle("Reconstructed initial proton mass [MeV]");
h1_reco_mass_11B_p2p_cut_iPhos_barrel->GetYaxis()->SetTitle("Counts");
h1_reco_mass_11B_p2p_cut_iPhos_barrel->GetXaxis()->CenterTitle(true);
h1_reco_mass_11B_p2p_cut_iPhos_barrel->GetYaxis()->CenterTitle(true);
h1_reco_mass_11B_p2p_cut_iPhos_barrel->GetYaxis()->SetLabelSize(0.045);
h1_reco_mass_11B_p2p_cut_iPhos_barrel->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "IPhos Angle vs reconstucted initial proton mass (for protons in iPhos and Barrel)");
h2_reco_mass_11B_iphos_angle = new TH2F(hist_name,hist_name,1000,500,1500,120,15,45);
h2_reco_mass_11B_iphos_angle->GetXaxis()->SetTitle("Reconstructed initial proton mass [MeV]");
h2_reco_mass_11B_iphos_angle->GetYaxis()->SetTitle("iPhos angle [degr]");
h2_reco_mass_11B_iphos_angle->GetXaxis()->CenterTitle(true);
h2_reco_mass_11B_iphos_angle->GetYaxis()->CenterTitle(true);
h2_reco_mass_11B_iphos_angle->GetYaxis()->SetLabelSize(0.045);
h2_reco_mass_11B_iphos_angle->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Barrel Angle  vs reconstucted initial proton mass (for protons in iPhos and Barrel)");
h2_reco_mass_11B_barrel_angle = new TH2F(hist_name,hist_name,1000,500,1500,240,40,100);
h2_reco_mass_11B_barrel_angle->GetXaxis()->SetTitle("Reconstructed initial proton mass [MeV]");
h2_reco_mass_11B_barrel_angle->GetYaxis()->SetTitle("Barrel angle [degr]");
h2_reco_mass_11B_barrel_angle->GetXaxis()->CenterTitle(true);
h2_reco_mass_11B_barrel_angle->GetYaxis()->CenterTitle(true);
h2_reco_mass_11B_barrel_angle->GetYaxis()->SetLabelSize(0.045);
h2_reco_mass_11B_barrel_angle->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Reconstructed mass of the 12C knocked out proton with angular cuts for both protons in Barrel");
h1_reco_mass_11B_p2p_cut_two_barrel = new TH1F(hist_name,hist_name,3000,-1500,1500);
h1_reco_mass_11B_p2p_cut_two_barrel->GetXaxis()->SetTitle("Reconstructed initial proton mass [MeV]");
h1_reco_mass_11B_p2p_cut_two_barrel->GetYaxis()->SetTitle("Counts");
h1_reco_mass_11B_p2p_cut_two_barrel->GetXaxis()->CenterTitle(true);
h1_reco_mass_11B_p2p_cut_two_barrel->GetYaxis()->CenterTitle(true);
h1_reco_mass_11B_p2p_cut_two_barrel->GetYaxis()->SetLabelSize(0.045);
h1_reco_mass_11B_p2p_cut_two_barrel->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Momentum of 10B (in 12C frame)");
h1_mom_abs_10B = new TH1F(hist_name,hist_name,750,0,1500);
h1_mom_abs_10B->GetXaxis()->SetTitle("abs(p_10B) [MeV/c]");
h1_mom_abs_10B->GetYaxis()->SetTitle("Counts");
h1_mom_abs_10B->GetXaxis()->CenterTitle(true);
h1_mom_abs_10B->GetYaxis()->CenterTitle(true);
h1_mom_abs_10B->GetYaxis()->SetLabelSize(0.045);
h1_mom_abs_10B->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "X component Momentum of 10B (in 12C frame)");
h1_mom_x_10B = new TH1F(hist_name,hist_name,1000,-1000,1000);
h1_mom_x_10B->GetXaxis()->SetTitle("p_10B x component [MeV/c]");
h1_mom_x_10B->GetYaxis()->SetTitle("Counts");
h1_mom_x_10B->GetXaxis()->CenterTitle(true);
h1_mom_x_10B->GetYaxis()->CenterTitle(true);
h1_mom_x_10B->GetYaxis()->SetLabelSize(0.045);
h1_mom_x_10B->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Y component Momentum of 10B (in 12C frame)");
h1_mom_y_10B = new TH1F(hist_name,hist_name,1000,-1000,1000);
h1_mom_y_10B->GetXaxis()->SetTitle("p_10B y component [MeV/c]");
h1_mom_y_10B->GetYaxis()->SetTitle("Counts");
h1_mom_y_10B->GetXaxis()->CenterTitle(true);
h1_mom_y_10B->GetYaxis()->CenterTitle(true);
h1_mom_y_10B->GetYaxis()->SetLabelSize(0.045);
h1_mom_y_10B->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Z component Momentum of 10B (in 12C frame)");
h1_mom_z_10B = new TH1F(hist_name,hist_name,1000,-1000,1000);
h1_mom_z_10B->GetXaxis()->SetTitle("p_10B z component [MeV/c]");
h1_mom_z_10B->GetYaxis()->SetTitle("Counts");
h1_mom_z_10B->GetXaxis()->CenterTitle(true);
h1_mom_z_10B->GetYaxis()->CenterTitle(true);
h1_mom_z_10B->GetYaxis()->SetLabelSize(0.045);
h1_mom_z_10B->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Summed up 3-momentum p_i + p_n + p_10B vs cos(p_i,p_n)");
h2_sum_mom_vs_cos_p_i_p_n_xy = new TH2F(hist_name,hist_name,25,-1,1,1000,0,1000);
h2_sum_mom_vs_cos_p_i_p_n_xy->GetXaxis()->SetTitle("Cos(p_i,p_n) in x-y plane");
h2_sum_mom_vs_cos_p_i_p_n_xy->GetYaxis()->SetTitle("Summed momentum (p_i+p_n+p_10B) [MeV/c]");
h2_sum_mom_vs_cos_p_i_p_n_xy->GetXaxis()->CenterTitle(true);
h2_sum_mom_vs_cos_p_i_p_n_xy->GetYaxis()->CenterTitle(true);
h2_sum_mom_vs_cos_p_i_p_n_xy->GetYaxis()->SetLabelSize(0.045);
h2_sum_mom_vs_cos_p_i_p_n_xy->GetYaxis()->SetTitleSize(0.045);


//END OF IMPLEMENTATION OF HISTOS--------------------------------------------------


sprintf(f_out_name,"/home/ge37liw/plots/histos/zero_with_miss/fast_macro/phi_test_automatic_src_califa_%s.root", count_i);
TFile * f = new TFile(f_out_name,"RECREATE");
//Path to input file
sprintf(fname,"/scratch8/ge37liw/workingspace/data/root_files/all_ts_unpack/sweep_target/ts_cone_cluster_%s_zero_cm_y_corr.root",count_i);
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
R3BEventHeader* DataCA = new R3BEventHeader();
TBranch* branchData = chain->GetBranch("R3BEventHeader");
branchData->SetAddress(&DataCA);
//
TClonesArray* WRMasterData = new TClonesArray("R3BWRMasterData",1);
R3BWRMasterData** wrmasterdata;
TBranch *branchWRMasterData = chain->GetBranch("WRMasterData");
branchWRMasterData->SetAddress(&WRMasterData);
//
TClonesArray* WRCalifaData = new TClonesArray("R3BWRCalifaData",3);
R3BWRCalifaData** wrcalifadata;
TBranch *branchWRCalifaData = chain->GetBranch("WRCalifaData");
branchWRCalifaData->SetAddress(&WRCalifaData);
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
		Double_t yMW0 = sofmwpc0hitdata[0]->GetY();
        Double_t yMW1 = sofmwpc1hitdata[0]->GetY();
        Double_t yMW2 = sofmwpc2hitdata[0]->GetY();
        Double_t yMW3 = sofmwpc3hitdata[0]->GetY();
if (xMW1 != -1000 && xMW2 != -1000 && xMW3 != -1000  && xMW0 != -1000 && yMW0 != -1000){
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

						//it has to beconsidered that the beam does not go parallel to the z-direction, but has deflections in x and y. Therefore also for later analysis and boosting the right trajectory 
						//before the 12C hits the target has to be found. I do this by using the info of the target point and the x-y info of the MW0
						TVector3 v_mw0(xMW0+xMW0_shift,yMW0+yMW0_shift,-2520); //vector of the MW0
						TVector3 v_target(0,0,-684.5); //vector of target position
						TVector3 v_diff_mw0_target;
						v_diff_mw0_target = v_mw0 - v_target;
						Double_t c12_arzimuth = v_diff_mw0_target.Phi(); // get arzimuthal angle of the incoming 12C beam (careful, it's negative..., or?)
						Double_t c12_theta = v_diff_mw0_target.Theta(); // get polar angle of the incoming 12C beam (careful, it's negative ..., or?)
						Double_t beta_x = abs(sin(c12_theta)*cos(c12_arzimuth));
						Double_t beta_y = abs(sin(c12_theta)*sin(c12_arzimuth));
						Double_t beta_z = abs(cos(c12_theta));
						TVector3 v_beta_12c (beta_x*beta_given,beta_y*beta_given,beta_z*beta_given);
						//compute position direction of fragment after target using info of mw1 and mw2
						TVector3 v_mw1(xMW1+xMW1_shift,yMW1+yMW1_shift,zM1);
						TVector3 v_mw2(xMW2+xMW2_shift,yMW2+yMW2_shift,zM2);
						TVector3 v_diff_mw2_mw1;
						v_diff_mw2_mw1 = v_mw2 - v_mw1;
						Double_t fragment_theta = v_diff_mw2_mw1.Theta();
						Double_t fragment_arzimuth = v_diff_mw2_mw1.Phi();
						Double_t fragment_x = sin(fragment_theta)*cos(fragment_arzimuth);
						Double_t fragment_y = sin(fragment_theta)*sin(fragment_arzimuth);
						Double_t fragment_z = cos(fragment_theta);


                        if (charge_val< t && charge_val> t-1.) //Z = 5 before it was charge_val< charge_cut && charge_val> charge_cut-0.8, too straight ....
                        {

//############################################################################################################
//------------BEGIN OF 11B ANALYSIS---------------------------------------------------------------------------
//############################################################################################################
                        if (a_q > s[0] && a_q < (s[0]+2*s[3]) && entries_califa >= 2){  //11B
						psi_in = u[3] - u[4]*xMW3 - u[5];
                        m = (cos(psi_out)-cos(psi_in))/(sin(psi_out)+sin(psi_in));
                        rho =  (Leff/(2*sin((psi_in+psi_out)/2)))*sqrt(pow(((m+tan(alpha_G))/(1-m*tan(alpha_G))),2)+1);
                        w = 2*abs(asin(BD/(2*rho)));
                        Double_t pathlength_precise = abs((zB-zT)/(cos(psi_in)))+rho*w+abs((zToFW-z_D)/cos(psi_out));
                        Double_t beta_precise = ((pathlength_precise/time_target_tof)*pow(10,6))/light_c;
                        Double_t gamma_precise = 1/(sqrt(1-beta_precise*beta_precise));
                        Double_t a_q_precise = (pow(10,-3)*((((mag_field*rho)/(beta_precise))*1.602176634*pow(10,-19))/(1.66053906660*pow(10,-27)*light_c)))/gamma_given;
						// checking MWPCs Y-orientation:
						h2_corr_y_mws->Fill((yMW2+yMW2_shift)-(yMW1+yMW1_shift),(yMW3+yMW3_shift)-(yMW1+yMW1_shift));
						h2_corr_y_mws_minusmw3->Fill((yMW2+yMW2_shift)-(yMW1+yMW1_shift),-(yMW3+yMW3_shift)-(yMW1+yMW1_shift));
						
					    
						//compute position direction of fragment after the target using the target pos. and the MW1/2 (x-pos) and MWPC3 (y-pos.)
//--> dangerous!! set fragment y to -y, TODO: change it back after usage
						fragment_y = -(yMW3+yMW3_shift)/(pathlength_precise-760.4); //substract path between MW3 and TOFW
						fragment_x =  v_diff_mw2_mw1.X()/v_diff_mw2_mw1.Mag();
						fragment_z = v_diff_mw2_mw1.Z()/v_diff_mw2_mw1.Mag();	
						//beginning after isotope has been selected:
						//TLorentzVector of 11B in the lab and cms 
						Double_t tot_mom11B = 10255.09731131*(beta_precise/(sqrt(1-(beta_precise*beta_precise)))); //here now right mass with mass defect
						//TODO: this energy calculation is wrong, be careful and recheck it..
						//Double_t energy_11B = sqrt(pow(tot_mom11B,2)+pow(11*938.272,2));
						Double_t energy_11B = (1./(sqrt(1-(beta_precise*beta_precise))))*10255.09731131;  //here now right mass with mass defect
						TLorentzVector v_11B_lab(tot_mom11B*fragment_x,tot_mom11B*fragment_y,tot_mom11B*fragment_z,energy_11B);
						TLorentzVector v_11B_cms = v_11B_lab;
						v_11B_cms.Boost(-v_beta_12c);
						//Lorentz vector of the 12C in rest 
						TLorentzVector p_12C_cms(0,0,0,11177.92337806);

						//TLorentzVector of the target in lab and cms frame
						TLorentzVector v_target_lab(0,0,0,938.272);
						TLorentzVector v_target_cms = v_target_lab;
						v_target_cms.Boost(-v_beta_12c);
						//fill entries of CALIFA in vectors, no matter if they are protons or gammas..
						//variales needed:
						vector<vector<double> > v_CALIFA_evts_lab_energy; //here I insert vectors v = (E_lab,theta_lab,phi_lab)
						vector<vector<double> > v_CALIFA_evts_cms_energy; //here I inserv vectors v = (E_cms,theta_lab,phi_lab)
						//for latest update I change to: v_CALIFA_evts_cms_energy(E_cms,theta_lab,phi_lab,Nf,Ns) 
						vector<TLorentzVector> v_lorentz_CALIFA_evts_lab; //here I insert Lorentz vector v = (p_x_lab,p_y_lab,p_z_lab,E_lab)
						vector<TLorentzVector> v_lorentz_CALIFA_evts_cms; //here I insert Lorentz vector v = (p_x_cms,p_y_cms,p_z_cms,E_cms)
						vector<TLorentzVector> v_lorentz_PROTONS_lab; //here I inserrt lorentz vector v = (p_x_lab,p_y_lab,p_z_lab,E_lab) where E_lab > 30MeV
						vector<TLorentzVector> v_lorentz_PROTONS_cms; //here I inserrt lorentz vector v = (p_x_cms,p_y_cms,p_z_cms,E_cms) 
						vector<TLorentzVector> v_lorentz_GAMMAS_lab; //here I inserrt lorentz vector v = (p_x_lab,p_y_lab,p_z_lab,E_lab) where E_lab < 10MeV for now... (change to 30MeV??)
						vector<TLorentzVector> v_lorentz_GAMMAS_cms; //here I inserrt lorentz vector v = (p_x_lab,p_y_lab,p_z_lab,E_cms) 
						
						for (Int_t j = 0.;j<entries_califa;j++){
							vector<double> v_temp(3);
							vector<double> v_temp_cms(5);
							califahitdata[j] = (R3BCalifaHitData*)CalifaHitData->At(j);
							v_temp[0] = (califahitdata[j]->GetEnergy())/1000.;
							//randomize angles to remove artifacts...
							//if (califahitdata[j]->GetTheta() > 0.45){
							//	Int_t rand_angl = rand() % 4500 + 1;
							//	v_temp[1] = califahitdata[j]->GetTheta() - 0.0225 + 0.00001*rand_angl;
							//	}
							//else {
							//	v_temp[1] = califahitdata[j]->GetTheta();
							//	}
							//
							v_temp[1] = califahitdata[j]->GetTheta();
							v_temp[2] = califahitdata[j]->GetPhi();
							v_temp_cms[0] = v_temp[0]*(1/(sqrt(1-beta_precise*beta_precise)))*(1-beta_precise*cos(v_temp[1])); //Doppler correction
							v_temp_cms[1] = v_temp[1];
						    v_temp_cms[2] = v_temp[2];
							v_temp_cms[3] = (califahitdata[j]->GetNf());
							v_temp_cms[4] = (califahitdata[j]->GetNs());
							v_CALIFA_evts_lab_energy.push_back(v_temp);
							v_CALIFA_evts_cms_energy.push_back(v_temp_cms);
							v_temp.clear();
							v_temp_cms.clear();
							
						}
						
						
						for (int CALIFA_events = 0; CALIFA_events < entries_califa; CALIFA_events++){
						
						if(v_CALIFA_evts_lab_energy[CALIFA_events][0] < 10) { //gammas
						//create lorentz vector in the laboratory system
						TLorentzVector temp_lorentz_lab;
						TLorentzVector temp_lorentz_cms;
						temp_lorentz_lab[0] = (v_CALIFA_evts_lab_energy[CALIFA_events][0])*sin(v_CALIFA_evts_lab_energy[CALIFA_events][1])*cos(v_CALIFA_evts_lab_energy[CALIFA_events][2]); //coordinate transform
						temp_lorentz_lab[1] = (v_CALIFA_evts_lab_energy[CALIFA_events][0])*sin(v_CALIFA_evts_lab_energy[CALIFA_events][1])*sin(v_CALIFA_evts_lab_energy[CALIFA_events][2]); //coordinate transform
						temp_lorentz_lab[2] = (v_CALIFA_evts_lab_energy[CALIFA_events][0])*cos(v_CALIFA_evts_lab_energy[CALIFA_events][1]); //coordinate transform
						temp_lorentz_lab[3] = v_CALIFA_evts_lab_energy[CALIFA_events][0];
						v_lorentz_CALIFA_evts_lab.push_back(temp_lorentz_lab);
						v_lorentz_GAMMAS_lab.push_back(temp_lorentz_lab);
						temp_lorentz_lab.Boost(-v_beta_12c);
						v_lorentz_CALIFA_evts_cms.push_back(temp_lorentz_lab);
						v_lorentz_GAMMAS_cms.push_back(temp_lorentz_lab);
						temp_lorentz_lab = TLorentzVector();
						
						}
						else {//protons
						
						TLorentzVector temp_lorentz_lab;
						temp_lorentz_lab[0] = sqrt(pow(((v_CALIFA_evts_lab_energy[CALIFA_events][0])+938.272),2)-pow(938.272,2))*sin(v_CALIFA_evts_lab_energy[CALIFA_events][1])*cos(v_CALIFA_evts_lab_energy[CALIFA_events][2]); //coordinate transform
						temp_lorentz_lab[1] = sqrt(pow(((v_CALIFA_evts_lab_energy[CALIFA_events][0])+938.272),2)-pow(938.272,2))*sin(v_CALIFA_evts_lab_energy[CALIFA_events][1])*sin(v_CALIFA_evts_lab_energy[CALIFA_events][2]); //coordinate transform
						temp_lorentz_lab[2] = sqrt(pow(((v_CALIFA_evts_lab_energy[CALIFA_events][0])+938.272),2)-pow(938.272,2))*cos(v_CALIFA_evts_lab_energy[CALIFA_events][1]); //coordinate transform
						temp_lorentz_lab[3] = v_CALIFA_evts_lab_energy[CALIFA_events][0]+938.272;
						v_lorentz_CALIFA_evts_lab.push_back(temp_lorentz_lab);
						v_lorentz_PROTONS_lab.push_back(temp_lorentz_lab);
						temp_lorentz_lab.Boost(-v_beta_12c);
						v_lorentz_CALIFA_evts_cms.push_back(temp_lorentz_lab);
						v_lorentz_PROTONS_cms.push_back(temp_lorentz_lab);
						temp_lorentz_lab = TLorentzVector();	
						}
						}
						//now sort the vectors according to the Energy.......
						//3d vectors
						sort(v_CALIFA_evts_cms_energy.begin(),v_CALIFA_evts_cms_energy.end(),sortcol_pseudo);
						sort(v_CALIFA_evts_lab_energy.begin(),v_CALIFA_evts_lab_energy.end(),sortcol_pseudo);
						//lorentz vectors
						sort(v_lorentz_CALIFA_evts_lab.begin(),v_lorentz_CALIFA_evts_lab.end(),sortcol_lorentz);
						sort(v_lorentz_CALIFA_evts_cms.begin(),v_lorentz_CALIFA_evts_cms.end(),sortcol_lorentz);
						sort(v_lorentz_GAMMAS_lab.begin(),v_lorentz_GAMMAS_lab.end(),sortcol_lorentz);
						sort(v_lorentz_GAMMAS_cms.begin(),v_lorentz_GAMMAS_cms.end(),sortcol_lorentz);
						sort(v_lorentz_PROTONS_lab.begin(),v_lorentz_PROTONS_lab.end(),sortcol_lorentz);
						sort(v_lorentz_PROTONS_cms.begin(),v_lorentz_PROTONS_cms.end(),sortcol_lorentz);

						//make restriction of having two protons (= more than 30 MeV in lab frame...)
						if (v_CALIFA_evts_lab_energy.size() >= 2 && v_CALIFA_evts_lab_energy[0][0] > 30 && v_CALIFA_evts_lab_energy[1][0] > 30){
								Double_t biggest_proton_e = v_lorentz_PROTONS_cms[0][3];
								Double_t big_proton_e = v_lorentz_PROTONS_cms[1][3];
								TLorentzVector biggest_proton_lab = v_lorentz_PROTONS_cms[0];
								TLorentzVector big_proton_lab = v_lorentz_PROTONS_cms[1];
								biggest_proton_lab.Boost(v_beta_12c);
								big_proton_lab.Boost(v_beta_12c);	
								h2_E_sum_vs_theta_sum->Fill(biggest_proton_e+big_proton_e-2*938,(biggest_proton_lab.Theta()+big_proton_lab.Theta())*(180./PI));
								//Definition of four-vector p_missing: in cms frame
								TLorentzVector p_missing = v_lorentz_PROTONS_cms[0] +v_lorentz_PROTONS_cms[1]-v_target_cms;
								//Definition of missing Energy = m_proton - e_miss (where e_miss is equal to the energy component of p_missing)
								Double_t missing_energy = 938.272 - p_missing.E();
								//Reconstructing missing mass as for src 10B (where the reconstructed mass was equal to the neutron mass)
								//here I expect a peak at around 0...
								TLorentzVector p_12C_cms(0,0,0,11177.92337806);
								TLorentzVector tl_reco = p_12C_cms + v_target_cms -  v_lorentz_PROTONS_cms[0] - v_lorentz_PROTONS_cms[1] - v_11B_cms;		
								//construct momentum of p_missing:
								Double_t test_mom = sqrt(p_missing.Px()*p_missing.Px()+p_missing.Py()*p_missing.Py()+p_missing.Pz()*p_missing.Pz());
								Double_t test_mom_11b = sqrt(v_11B_cms.Px()*v_11B_cms.Px() + v_11B_cms.Py()*v_11B_cms.Py() + v_11B_cms.Pz()*v_11B_cms.Pz());
								Double_t M_missing = tl_reco.Mag();	
								h1_missing_mass_11B->Fill(M_missing);
								h1_missing_energy_11B_p2p->Fill(missing_energy);
								h2_e_miss_vs_theta_sum_p2p_11B->Fill(missing_energy,((v_CALIFA_evts_cms_energy[0][1]+v_CALIFA_evts_cms_energy[1][1])/PI)*180.);
								//Filling histo: missing Energy(cms) vs theta1+theta2
								//Fillin histo: missing momentum (cms) where missing_mom is equal to the magnitude of the 3-vector part of p_missing
								h1_p_missing_11B->Fill(sqrt(pow(p_missing[0],2)+pow(p_missing[1],2)+pow(p_missing[2],2))); //mom_miss = sqrt(p_x^2+p_y^2+p_z^2)
								//Filling missing momentum histo  for x and z component:
								h1_p_missing_11B_z->Fill(p_missing[2]);
								h1_p_missing_11B_x->Fill(p_missing[0]);	
								//As we do just have the x,z info for the 11B fragment (MWPCs) we also need to project the p1 and p2 protons to the x,z plane for making \
								//some angular distribution analysis.... :
								TLorentzVector proton_1_x_z = v_lorentz_PROTONS_cms[0];
								//now setting the y contribution to 0:
								proton_1_x_z.SetPy(0.);
								//same for proton 2:
								TLorentzVector proton_2_x_z = v_lorentz_PROTONS_cms[1];
								proton_2_x_z.SetPy(0.);
								//p_missing in the CMS x,z plane...:
								TLorentzVector p_missing_x_z = proton_1_x_z + proton_2_x_z - v_target_cms;
								
								//Angle between p_missing_x_z (3-vec) and 11B (3-vec) in the x-z plane
								Double_t angle_pi_11B = p_missing_x_z.Angle(v_11B_cms.Vect());
								//filling histogram showing angle between p_i and 11B in CMS
								h1_angle_pi_11B->Fill(angle_pi_11B*(180/PI));
								h1_angle_pi_11B_cos->Fill(cos(angle_pi_11B));

								//making cut on the QE (p,2p) reaction: theta1 + theta2 < 90° and phi1+phi2 = 180° +- 40°
								if (v_CALIFA_evts_cms_energy[0][1]+v_CALIFA_evts_cms_energy[1][1] < PI/2. && abs(v_CALIFA_evts_cms_energy[0][2]-v_CALIFA_evts_cms_energy[1][2]) < 3.83972 \
									&& abs(v_CALIFA_evts_cms_energy[0][2]-v_CALIFA_evts_cms_energy[1][2]) > 2.44346){
									//Short Analysis for Time Stitching etc...
									Int_t entries_wr = WRMasterData->GetEntries();
									Int_t entries_wr_califa = WRCalifaData->GetEntries();
									cout << ">>>>>>>>> This is eventnr.:\t" << evtnr << endl;
									if (entries_wr_califa == 2 && entries_wr){
										wrcalifadata = new R3BWRCalifaData*[entries_wr_califa];
										wrmasterdata = new R3BWRMasterData*[entries_wr];
										wrcalifadata[0] = (R3BWRCalifaData*)WRCalifaData->At(0);
    									wrcalifadata[1] = (R3BWRCalifaData*)WRCalifaData->At(1);
										wrmasterdata[0] = (R3BWRMasterData*)WRMasterData->At(0);
										//get time difference between frst WRCALIFA and WRMaster
										h1_wr_first_wr_master->Fill((wrcalifadata[0]->GetTimeStamp())-(wrmasterdata[0]->GetTimeStamp()));
										//get time difference between second WRCALIFA and WRMaster
										h1_wr_second_wr_master->Fill((wrcalifadata[1]->GetTimeStamp())-(wrmasterdata[0]->GetTimeStamp()));

										}
									else
										cout << "NOT ENOUGH WR entries for eventnr.:\t" << evtnr << "   " << "something went wrong..." << endl;
									



									//--------------------------------------

									h2_missing_energy_vs_t1_t2_11B_p2p->Fill(missing_energy,(v_CALIFA_evts_cms_energy[0][1]+v_CALIFA_evts_cms_energy[1][1])*(180/PI));
									h1_missing_energy_11B_p2p_cut->Fill(missing_energy);	
									h1_p_missing_11B_z_cut->Fill(p_missing[2]);
									h1_p_missing_11B_y_cut->Fill(p_missing[1]);
									h1_p_missing_11B_x_cut->Fill(p_missing[0]);
									//same but with the 11B momentum in the 12C frame:
									h1_mom11B_z_cut->Fill(v_11B_cms.Pz());
									h1_mom11B_x_cut->Fill(v_11B_cms.Px());
									h1_mom11B_y_cut->Fill(v_11B_cms.Py());
									h1_mom11B_cut->Fill(sqrt(pow(v_11B_cms.Pz(),2)+pow(v_11B_cms.Px(),2)+pow(v_11B_cms.Py(),2)));
									h1_p_missing_11B_cut->Fill(sqrt(pow(p_missing[0],2)+pow(p_missing[1],2)+pow(p_missing[2],2)));
									h2_p_missing_vs_p_11B_cut->Fill(sqrt(pow(p_missing[0],2)+pow(p_missing[1],2)+pow(p_missing[2],2)),sqrt(pow(v_11B_cms[0],2)+pow(v_11B_cms[1],2)+pow(v_11B_cms[2],2)));
									h1_missing_mass_11B_p2p_cut->Fill(p_missing.M());
									//analyze if both protons hit IPhos or Barrel region...
									Double_t theta_0 = v_CALIFA_evts_lab_energy[0][1];
									Double_t theta_1 = v_CALIFA_evts_lab_energy[1][1];
									if (theta_0 < 0.75049157836 && theta_1 < 0.75049157836){
										cout << (*angle_comb)[0] << endl;
										cout << angle_comb[0] << endl;
										h1_reco_mass_11B_p2p_cut_two_iPhos->Fill(p_missing.M());
										
										h2_reco_mass_angle_comb_cut->Fill(p_missing.M(),0);
										}
									else if ((theta_0 < 0.75049157836 && theta_1 > 0.75049157836) || (theta_0 > 0.75049157836 && theta_1 < 0.75049157836)){
										cout << (*angle_comb)[1] << endl;
										cout << angle_comb[1] << endl;
										h1_reco_mass_11B_p2p_cut_iPhos_barrel->Fill(p_missing.M());
										h2_reco_mass_angle_comb_cut->Fill(p_missing.M(),1);
										if (theta_0 > theta_1){
											h2_reco_mass_11B_iphos_angle->Fill(p_missing.M(),(theta_1*180.)/PI);
											h2_reco_mass_11B_barrel_angle->Fill(p_missing.M(),(theta_0*180.)/PI);
											}
										else {
											h2_reco_mass_11B_iphos_angle->Fill(p_missing.M(),(theta_0*180.)/PI);
                                            h2_reco_mass_11B_barrel_angle->Fill(p_missing.M(),(theta_1*180.)/PI);
		
											}
										}
									else if (theta_0 > 0.75049157836 && theta_1 > 0.75049157836){
										cout << (*angle_comb)[2] << endl;
										cout << angle_comb[2] << endl;
										h1_reco_mass_11B_p2p_cut_two_barrel->Fill(p_missing.M());
										h2_reco_mass_angle_comb_cut->Fill(p_missing.M(),2);
										}
									//----------------------------------------------------
									h1_excit_energy_11B_cut->Fill((p_12C_cms + v_target_cms - v_lorentz_PROTONS_cms[0]- v_lorentz_PROTONS_cms[1]).M() - 10255.09731131);
									//here I compare the excitation energy and missing energy formula....
									h2_missing_e_2meth_11B_p2p_cut->Fill(missing_energy,(p_12C_cms + v_target_cms - v_lorentz_PROTONS_cms[0]- v_lorentz_PROTONS_cms[1]).M() - 10255.09731131);
					
									int rand_num;
                                	rand_num = rand() % 2 ;

                                	//computing P_y as Valerii does in his thesis (page 66)
                                	Double_t P_y_knocked_out = sqrt(pow(v_lorentz_PROTONS_lab[rand_num][0],2)+pow(v_lorentz_PROTONS_lab[rand_num][1],2)+pow(v_lorentz_PROTONS_lab[rand_num][2],2))*sin(v_lorentz_PROTONS_lab[rand_num].Theta())*sin(v_lorentz_PROTONS_lab[rand_num].Phi()-v_lorentz_PROTONS_lab[1-rand_num].Phi());
									Double_t P_fragment = v_11B_lab.Py(); 
									h2_corr_y_proton_fragm->Fill(P_fragment,P_y_knocked_out);
									h2_corr_x_proton_fragm->Fill(v_11B_lab.Px(),P_y_knocked_out);
									//not fully convinced from the above formula... what I get as formula P_y = Q_k*sin(theta_k)*sin(phi_k) - Q_i*sin(theta_i)*sin(phi_i) 
									//where k is the knocked out proton in the 12C and i is the target proton (at rest in laboratory frame)
									Double_t P_y_knocked_out_tj = sqrt(pow(v_lorentz_PROTONS_lab[rand_num][0],2)+pow(v_lorentz_PROTONS_lab[rand_num][1],2)+pow(v_lorentz_PROTONS_lab[rand_num][2],2))*sin(v_lorentz_PROTONS_lab[rand_num].Theta())*sin(v_lorentz_PROTONS_lab[rand_num].Phi()) -sqrt(pow(v_lorentz_PROTONS_lab[1-rand_num][0],2)+pow(v_lorentz_PROTONS_lab[1-rand_num][1],2)+pow(v_lorentz_PROTONS_lab[1-rand_num][2],2))*sin(v_lorentz_PROTONS_lab[1-rand_num].Theta())*sin(v_lorentz_PROTONS_lab[1-rand_num].Phi());
									h2_corr_y_proton_fragm_tj->Fill(P_fragment,P_y_knocked_out_tj);	
								//by setting the y-value as above for the particles to 0 (see: proton_2_x_z, p_missing_x_z,...) 
								//we get an anguar distortion. Therefore I want to use here p_missing instead and together with 
								//the angular cut compute cos(p_i and p_11B) again:
								Double_t angle_pi_11B_3d = p_missing.Angle(v_11B_cms.Vect());
								h1_angle_pi_11B_cos_cut->Fill(cos(angle_pi_11B_3d));	
								//in the paper they calculate the angle in the x-y plane, so I will do...
								TVector3 p3_missing_x_y(p_missing.Px(),p_missing.Py(),0);
								TVector3 p3_11b_x_y(v_11B_cms.Px(),v_11B_cms.Py(),0);
								h1_angle_pi_11B_cos_cut_plane->Fill(cos(p3_missing_x_y.Angle(p3_11b_x_y)));
								if (cos(p3_missing_x_y.Angle(p3_11b_x_y)) < 0){
									h2_corr_y_proton_fragm_cos_less->Fill(v_11B_lab.Py(),P_y_knocked_out);
									}
								if (cos(p3_missing_x_y.Angle(p3_11b_x_y)) > 0){
									h2_corr_y_proton_fragm_cos_more->Fill(v_11B_lab.Py(),P_y_knocked_out);
									}
								//make plot h1_angle_pi_11B_cos_cut_plane mroe restictive: 77 < theta < 85 and E_miss > -60:
								if ((v_CALIFA_evts_cms_energy[0][1]+v_CALIFA_evts_cms_energy[1][1]) > 1.343903524 && (v_CALIFA_evts_cms_energy[0][1]+v_CALIFA_evts_cms_energy[1][1]) < 1.4835298642 && missing_energy > -60.){
									h1_angle_pi_11B_cos_hcut_plane->Fill(cos(p3_missing_x_y.Angle(p3_11b_x_y)));	
									} 

								//another test:
								//check arzimuth angle of p1+p2, if positivem,give b11 an angle in y-direction of -0.03 rad, otherwise vice versa...
								Double_t p1_plus_p2_phi = (v_lorentz_PROTONS_cms[0] +v_lorentz_PROTONS_cms[1]).Phi();
									if (p1_plus_p2_phi > 0){
										TLorentzVector fake_y(0,tot_mom11B*sin(-0.03),0,0);
										TLorentzVector v_11B_cms_fake_y = (v_11B_lab + fake_y);	
										v_11B_cms_fake_y.Boost(-v_beta_12c);
										Double_t fake_angle = p_missing.Angle(v_11B_cms_fake_y.Vect());
										h1_angle_pi_11B_cos_cut_fake_y->Fill(cos(fake_angle));			
										h2_p_missing_vs_p_11B_cut_fake_y->Fill(sqrt(pow(p_missing[0],2)+pow(p_missing[1],2)+pow(p_missing[2],2)),sqrt(pow(v_11B_cms_fake_y[0],2)+pow(v_11B_cms_fake_y[1],2)+pow(v_11B_cms_fake_y[2],2)));
										}
									else if (p1_plus_p2_phi < 0){
										TLorentzVector fake_y(0,tot_mom11B*sin(0.03),0,0);
										TLorentzVector v_11B_cms_fake_y = (v_11B_lab + fake_y);
										v_11B_cms_fake_y.Boost(-v_beta_12c);
										Double_t fake_angle = p_missing.Angle(v_11B_cms_fake_y.Vect());
										h1_angle_pi_11B_cos_cut_fake_y->Fill(cos(fake_angle));	
										}
								
								
									}
								
								//-------------------------------------------------------------------------------		
								
							if (v_lorentz_GAMMAS_cms.size()){
								h1_gamma_energyE_max_val_11B->Fill(v_lorentz_GAMMAS_cms[0][3]);
								h1_gamma_energyE_max_val_11B_lab->Fill(v_lorentz_GAMMAS_lab[0][3]);
								//checking multiplicity...
								h2_gamma_energy_11B_vs_mult->Fill(v_lorentz_GAMMAS_cms[0][3],v_lorentz_GAMMAS_cms.size());
								//check angle distribution (in lab frame) vs largest energy in cm frame
								TLorentzVector highest_gamma_temp = v_lorentz_GAMMAS_cms[0];
								highest_gamma_temp.Boost(-v_beta_12c);
								h2_gamma_energy_11_vs_angle->Fill(v_lorentz_GAMMAS_cms[0][3],highest_gamma_temp.Theta()*(180./PI));
								//summing up the total energy deposited in CALIFA. As angle take the one with highest energy...
								Double_t sum_energy_gamma = highest_gamma_temp[3];
								for (int g = 1; g < v_lorentz_GAMMAS_cms.size();++g){
									TLorentzVector temp_gamma_v_lorentz = v_lorentz_GAMMAS_cms[g];
									v_lorentz_GAMMAS_cms[g].Boost(-v_beta_12c);
									sum_energy_gamma += v_lorentz_GAMMAS_cms[g][3];
								}
								Double_t doppl_corr_summ_gamma = sum_energy_gamma*(1/sqrt(1-beta_given*beta_given))*(1-beta_given*cos(highest_gamma_temp.Theta()));
								h1_gamma_energyE_sum_11B->Fill(doppl_corr_summ_gamma);
									
								//making cut on the QE (p,2p) reaction: theta1 + theta2 < 90° and phi1+phi2 = 180° +- 40°
								if (v_CALIFA_evts_cms_energy[0][1]+v_CALIFA_evts_cms_energy[1][1] < PI/2. && abs(v_CALIFA_evts_cms_energy[0][2]-v_CALIFA_evts_cms_energy[1][2]) < 3.83972  && abs(v_CALIFA_evts_cms_energy[0][2]-v_CALIFA_evts_cms_energy[1][2]) > 2.44346){
									h1_gamma_energyE_max_val_11B_cut->Fill(v_lorentz_GAMMAS_cms[0][3]);	
									h2_gamma_energy_11B_vs_mult_cut->Fill(v_lorentz_GAMMAS_cms[0][3],v_lorentz_GAMMAS_cms.size());
									h1_gamma_energyE_sum_11B_cut->Fill(doppl_corr_summ_gamma);
									//boosting back to 11B frame instead as incorrectly to the 12C frame...
									TLorentzVector tl_highest_g_11B = v_lorentz_GAMMAS_lab[0];
									tl_highest_g_11B.Boost(-sin(psi_in)*beta_precise,0,-cos(psi_in)*beta_precise);
									h1_gamma_energyE_max_val_11Bframe_cut->Fill(tl_highest_g_11B[3]);
									//cut on multiplicity < 3
									if (v_lorentz_GAMMAS_cms.size() < 3){
										h1_gamma_energyE_max_val_11B_m_cut->Fill(v_lorentz_GAMMAS_cms[0][3]);
									}
								}
								//--------------------------------------------------------------------------------------

								}	
								
							if (v_lorentz_PROTONS_lab.size() > 2){
								
								cout << "more than 3 protons, really? " << endl;
								}
							}
						
					}
                        }

//############################################################################################################
////------------END OF 11B ANALYSIS-------------------------------------------------------------------------
////##########################################################################################################

//############################################################################################################
//--------------BEGIN OF 10B ANALYSIS-------------------------------------------------------------------------
//############################################################################################################

                        if (a_q < s[0] && a_q > (s[0]-0.15) && entries_califa >= 2){  //10B
                        //psi_in = mean_psi_out_10b-slope_10b*xMW3 - angle_offs_10b;
						psi_in = u[0]-u[1]*xMW3 -u[2];
                        m = (cos(psi_out)-cos(psi_in))/(sin(psi_out)+sin(psi_in));
                        rho =  (Leff/(2*sin((psi_in+psi_out)/2)))*sqrt(pow(((m+tan(alpha_G))/(1-m*tan(alpha_G))),2)+1);
                        w = 2*abs(asin(BD/(2*rho)));
                        Double_t pathlength_precise = abs((zB-zT)/(cos(psi_in)))+rho*w+abs((zToFW-z_D)/cos(psi_out));
                        Double_t beta_precise = ((pathlength_precise/time_target_tof)*pow(10,6))/light_c;
						Double_t gamma_precise = 1/(sqrt(1-beta_precise*beta_precise));
                        Double_t a_q_precise = (pow(10,-3)*((((mag_field*rho)/(beta_precise))*1.602176634*pow(10,-19))/(1.66053906660*pow(10,-27)*light_c)))/gamma_given;

						//compute position direction of fragment after the target using the target pos. and the MW1/2 (x-pos) and MWPC3 (y-pos.)
						////--> dangerous!! set fragment y to -y, TODO: change it back after usage
						fragment_y = -(yMW3+yMW3_shift)/(pathlength_precise-760.4); //substract distance MW3 to TOFW
						fragment_x =  v_diff_mw2_mw1.X()/v_diff_mw2_mw1.Mag();
						fragment_z = v_diff_mw2_mw1.Z()/v_diff_mw2_mw1.Mag();

						//beginning after isotope has been selected:
						//TLorentzVector of 10B in the lab and cms 
						Double_t tot_mom10B = 9326.98688128*(beta_precise/(sqrt(1-(beta_precise*beta_precise)))); //now with correct mass defect
						//TODO: check this formula, it doesn't seem to be correct, as I use the lorentz factor twice...
						//Double_t energy_10B = sqrt(pow(tot_mom10B,2)+pow(10*938.272,2));
						Double_t energy_10B = (1./(sqrt(1-(beta_precise*beta_precise))))*9326.98688128; //now with correct mass defect
						TLorentzVector v_10B_lab(tot_mom10B*fragment_x,tot_mom10B*fragment_y,tot_mom10B*fragment_z,energy_10B);
						TLorentzVector v_10B_cms = v_10B_lab;
						v_10B_cms.Boost(-v_beta_12c);

						//TLorentzVector of the target in lab and cms frame
						TLorentzVector v_target_lab(0,0,0,938.272);
						TLorentzVector v_target_cms = v_target_lab;
						v_target_cms.Boost(-v_beta_12c);
						//fill entries of CALIFA in vectors, no matter if they are protons or gammas..
						//variales needed:
						vector<vector<double> > v_CALIFA_evts_lab_energy; //here I insert vectors v = (E_lab,theta_lab,phi_lab)
						vector<vector<double> > v_CALIFA_evts_cms_energy; //here I inserv vectors v = (E_cms,theta_lab,phi_lab)
						//for latest update I change to: v_CALIFA_evts_cms_energy(E_cms,theta_lab,phi_lab,Nf,Ns)
						vector<TLorentzVector> v_lorentz_CALIFA_evts_lab; //here I insert Lorentz vector v = (p_x_lab,p_y_lab,p_z_lab,E_lab)
						vector<TLorentzVector> v_lorentz_CALIFA_evts_cms; //here I insert Lorentz vector v = (p_x_cms,p_y_cms,p_z_cms,E_cms)
						vector<TLorentzVector> v_lorentz_PROTONS_lab; //here I inserrt lorentz vector v = (p_x_lab,p_y_lab,p_z_lab,E_lab) where E_lab > 30MeV
						vector<TLorentzVector> v_lorentz_PROTONS_cms; //here I inserrt lorentz vector v = (p_x_cms,p_y_cms,p_z_cms,E_cms) 
						vector<TLorentzVector> v_lorentz_GAMMAS_lab; //here I inserrt lorentz vector v = (p_x_lab,p_y_lab,p_z_lab,E_lab) where E_lab < 30MeV
						vector<TLorentzVector> v_lorentz_GAMMAS_cms; //here I inserrt lorentz vector v = (p_x_lab,p_y_lab,p_z_lab,E_cms) 
						TVector3 v_boost(0,0,beta_precise*light_c);
						
						for (Int_t j = 0.;j<entries_califa;j++){
							vector<double> v_temp(3);
							vector<double> v_temp_cms(5);
							califahitdata[j] = (R3BCalifaHitData*)CalifaHitData->At(j);
							v_temp[0] = (califahitdata[j]->GetEnergy())/1000.;
							//randomize angles to remove artifacts...	
							//if (califahitdata[j]->GetTheta() > 0.45){
                            //    Int_t rand_angl = rand() % 4500 + 1;
                            //    v_temp[1] = califahitdata[j]->GetTheta() - 0.0225 + 0.00001*rand_angl;
                            //    }
                            //else {
                            //    v_temp[1] = califahitdata[j]->GetTheta();
                            //    }	

							v_temp[1] = califahitdata[j]->GetTheta();
							v_temp[2] = califahitdata[j]->GetPhi();
							v_temp_cms[0] = v_temp[0]*(1/(sqrt(1-beta_precise*beta_precise)))*(1-beta_precise*cos(v_temp[1])); //Doppler correction
							v_temp_cms[1] = v_temp[1];
						    v_temp_cms[2] = v_temp[2];
							v_temp_cms[3] = (califahitdata[j]->GetNf());
							v_temp_cms[4] = (califahitdata[j]->GetNs());
							v_CALIFA_evts_lab_energy.push_back(v_temp);
							v_CALIFA_evts_cms_energy.push_back(v_temp_cms);
							v_temp.clear();
							v_temp_cms.clear();
							
						}
						
						
						for (int CALIFA_events = 0; CALIFA_events < entries_califa; CALIFA_events++){
						
						if(v_CALIFA_evts_lab_energy[CALIFA_events][0] < 30) { //gammas
						//create lorentz vector in the laboratory system
						TLorentzVector temp_lorentz_lab;
						TLorentzVector temp_lorentz_cms;
						temp_lorentz_lab[0] = (v_CALIFA_evts_lab_energy[CALIFA_events][0])*sin(v_CALIFA_evts_lab_energy[CALIFA_events][1])*cos(v_CALIFA_evts_lab_energy[CALIFA_events][2]); //coordinate transform
						temp_lorentz_lab[1] = (v_CALIFA_evts_lab_energy[CALIFA_events][0])*sin(v_CALIFA_evts_lab_energy[CALIFA_events][1])*sin(v_CALIFA_evts_lab_energy[CALIFA_events][2]); //coordinate transform
						temp_lorentz_lab[2] = (v_CALIFA_evts_lab_energy[CALIFA_events][0])*cos(v_CALIFA_evts_lab_energy[CALIFA_events][1]); //coordinate transform
						temp_lorentz_lab[3] = v_CALIFA_evts_lab_energy[CALIFA_events][0];
						v_lorentz_CALIFA_evts_lab.push_back(temp_lorentz_lab);
						v_lorentz_GAMMAS_lab.push_back(temp_lorentz_lab);
						//temp_lorentz_lab.Boost(v_boost);
						temp_lorentz_lab.Boost(-v_beta_12c);
						v_lorentz_CALIFA_evts_cms.push_back(temp_lorentz_lab);
						v_lorentz_GAMMAS_cms.push_back(temp_lorentz_lab);
						temp_lorentz_lab = TLorentzVector();
						
						}
						else {//protons
						
						TLorentzVector temp_lorentz_lab;
						temp_lorentz_lab[0] = sqrt(pow(((v_CALIFA_evts_lab_energy[CALIFA_events][0])+938.272),2)-pow(938.272,2))*sin(v_CALIFA_evts_lab_energy[CALIFA_events][1])*cos(v_CALIFA_evts_lab_energy[CALIFA_events][2]); //coordinate transform
						temp_lorentz_lab[1] = sqrt(pow(((v_CALIFA_evts_lab_energy[CALIFA_events][0])+938.272),2)-pow(938.272,2))*sin(v_CALIFA_evts_lab_energy[CALIFA_events][1])*sin(v_CALIFA_evts_lab_energy[CALIFA_events][2]); //coordinate transform
						temp_lorentz_lab[2] = sqrt(pow(((v_CALIFA_evts_lab_energy[CALIFA_events][0])+938.272),2)-pow(938.272,2))*cos(v_CALIFA_evts_lab_energy[CALIFA_events][1]); //coordinate transform
						temp_lorentz_lab[3] = v_CALIFA_evts_lab_energy[CALIFA_events][0]+938.272;
						v_lorentz_CALIFA_evts_lab.push_back(temp_lorentz_lab);
						v_lorentz_PROTONS_lab.push_back(temp_lorentz_lab);
						temp_lorentz_lab.Boost(-v_beta_12c);
						v_lorentz_CALIFA_evts_cms.push_back(temp_lorentz_lab);
						v_lorentz_PROTONS_cms.push_back(temp_lorentz_lab);
						temp_lorentz_lab = TLorentzVector();	
						}
						}
						//now sort the vectors according to the Energy.......
						//3d vectors
						sort(v_CALIFA_evts_cms_energy.begin(),v_CALIFA_evts_cms_energy.end(),sortcol_pseudo);
						sort(v_CALIFA_evts_lab_energy.begin(),v_CALIFA_evts_lab_energy.end(),sortcol_pseudo);
						//lorentz vectors
						sort(v_lorentz_CALIFA_evts_lab.begin(),v_lorentz_CALIFA_evts_lab.end(),sortcol_lorentz);
						sort(v_lorentz_CALIFA_evts_cms.begin(),v_lorentz_CALIFA_evts_cms.end(),sortcol_lorentz);
						sort(v_lorentz_GAMMAS_lab.begin(),v_lorentz_GAMMAS_lab.end(),sortcol_lorentz);
						sort(v_lorentz_GAMMAS_cms.begin(),v_lorentz_GAMMAS_cms.end(),sortcol_lorentz);
						sort(v_lorentz_PROTONS_lab.begin(),v_lorentz_PROTONS_lab.end(),sortcol_lorentz);
						sort(v_lorentz_PROTONS_cms.begin(),v_lorentz_PROTONS_cms.end(),sortcol_lorentz);
						
						//check the trigger pattern. For src we need entry in neuland..
						int tpat = DataCA->GetTpat();
						if (tpat == 4){
							//cout << "this is eventnr. of tpat4:\t" << evtnr << endl;
							//cout << ">>> this means, neutron has been detected in NEULAND <<<<<<<<" << endl;
						}
						//make restriction of having two protons (= more than 30 MeV in lab frame...) and no tpat ==4 (spill on + sofstart+ califa + neuland) and psi difference greater than  100 degrees
						if (v_CALIFA_evts_lab_energy.size() >= 2 && v_CALIFA_evts_lab_energy[0][0] > 30 && v_CALIFA_evts_lab_energy[1][0] > 30 && abs(v_CALIFA_evts_lab_energy[0][2] - v_CALIFA_evts_lab_energy[1][2]) > 1.745){
								h1_mom_abs_10B->Fill(sqrt(v_10B_cms.Px()*v_10B_cms.Px()+v_10B_cms.Py()*v_10B_cms.Py()+v_10B_cms.Pz()*v_10B_cms.Pz()));
								h1_mom_x_10B->Fill(v_10B_cms.Px());
								h1_mom_y_10B->Fill(v_10B_cms.Py());
								h1_mom_z_10B->Fill(v_10B_cms.Pz());
								//Definition of 12C four-vector in its cms frame:
								TLorentzVector p_12C_cms(0,0,0,11177.92337806);
								//Definition of four-vector p_missing: in cms frame
								TLorentzVector p_missing = p_12C_cms + v_target_cms -  v_lorentz_PROTONS_cms[0] - v_lorentz_PROTONS_cms[1] - v_10B_cms;
								//Definition of M_missing is the reaction missing mass (see Nature paper of Valerii)
								Double_t M_missing = p_missing.Mag();
								h1_missing_mass_10B->Fill(M_missing);	
								//Lorentz vector corresponding to the initial projectile proton in the 12C frame
								TLorentzVector p_proton_proj = v_lorentz_PROTONS_cms[0] +v_lorentz_PROTONS_cms[1]-v_target_cms;
								Double_t missing_energy = 938.272 - p_proton_proj.E();
								h2_e_miss_vs_theta_sum_p2p_10B->Fill(missing_energy,((v_CALIFA_evts_cms_energy[0][1]+v_CALIFA_evts_cms_energy[1][1])/PI)*180.);
								//Filling histo: missing momentum (cms) where missing_mom is equal to the magnitude of the 3-vector part of p_proton_proj. It's actually the momentum of the projectile proton before scattering...
								Double_t p_proton_proj_abs_val = sqrt(pow(p_proton_proj[0],2)+pow(p_proton_proj[1],2)+pow(p_proton_proj[2],2));
								h1_p_missing_10B->Fill(p_proton_proj_abs_val); //mom_miss = sqrt(p_x^2+p_y^2+p_z^2)
								
								
								//Filling histo missing_mom (CMS) of the projectile proton before QFS vs the theta angles of the two outgoing protons
								h2_missing_energy_vs_t1_t2_10B->Fill(sqrt(p_proton_proj*p_proton_proj),(v_CALIFA_evts_cms_energy[0][1]+v_CALIFA_evts_cms_energy[1][1])*(180/PI));
								
								//now we want to check the anglular distributions between the proton and neutron in the SRC\
								//therefore we need the CMS momentum of the interacting proton (= 3-vec of p_i, see calculation here...) and the CMS momentum of the \
								//do the cosine of the proton neutron pair and 10B in the x-y plane:
								TVector3 p3_missing_x_y(p_missing.Px(),p_missing.Py(),0.); //neutron vector
								TVector3 p3_p_i_x_y(p_proton_proj.Px(),p_proton_proj.Py(),0.);//initial proton vector 
								TVector3 p3_10B_x_y(v_10B_cms.Px(),v_10B_cms.Py(),0.);
								h1_angle_pi_p_n_xy_plane->Fill(cos(p3_missing_x_y.Angle(p3_p_i_x_y)));//histo p_i vs p_n in xy plane
								//histo of 10B vs src pair in xy plane
								h1_angle_src_vs_10B_cos_xy_plane->Fill(cos(((p3_p_i_x_y-p3_missing_x_y)*0.5).Angle(p3_10B_x_y)));
							    //histo of 10B vs src_pair in 3d
							    h1_angle_src_vs_10B_cos_3d->Fill(cos(((p_missing.Vect()-p_proton_proj.Vect())*0.5).Angle(v_10B_cms.Vect())));
								//use same restriction as in paper + phi_diff > 100 degree 
								if((p_proton_proj.Vect()).Mag() > 350. && ((v_CALIFA_evts_cms_energy[0][1]+v_CALIFA_evts_cms_energy[1][1])/PI)*180. > 63 && missing_energy > -110 && missing_energy < 240){
									h1_angle_src_vs_10B_cos_3d_paper_cut->Fill(cos(((p_missing.Vect()-p_proton_proj.Vect())*0.5).Angle(v_10B_cms.Vect())));	
								}

								//neutron from the SRC (= 3-vector of p_missing); put again everything into x-z plane!
								TLorentzVector p_i_x_z = v_lorentz_PROTONS_cms[0] + v_lorentz_PROTONS_cms[1] - v_target_cms;
								TLorentzVector p_missing_x_z = p_missing;
								p_missing_x_z.SetPy(0.);
								//do the x-z projection
								p_i_x_z.SetPy(0.);
								//also x-z projection for the 10B fragment
								v_10B_cms.SetPy(0.);
								//Angle between p_i (projectile proton(SRC)) and the according neutron 
								Double_t angle_pi_pn_src = p_i_x_z.Angle(p_missing_x_z.Vect());
								//Angle between src pair and 10B fragment:-> now with updated version p_rel
								Double_t angle_src_pair_10B_frag = ((p_i_x_z-p_missing_x_z)*0.5).Angle(v_10B_cms.Vect());
								//Fill histogram with angular distribution between p_i and p_n in x-z plane:
								h1_angle_pi_pn_10B->Fill(angle_pi_pn_src*(180/PI));
								h1_angle_pi_pn_10B_cos->Fill(cos(angle_pi_pn_src));
								h1_angle_src_vs_10B_cos->Fill(cos(angle_src_pair_10B_frag));
								
								TVector3 sum_mom_inner_part_cms = p_proton_proj.Vect() + p_missing.Vect() + v_10B_cms.Vect(); //sum of the momentum vectors
								//of the 12C inner patricles (should be around 0 from momentum conservation)
								h2_sum_mom_vs_cos_p_i_p_n_xy->Fill(cos(p3_missing_x_y.Angle(p3_p_i_x_y)),sum_mom_inner_part_cms.Mag());
								//Histogram theta1 vs theta2 randomizing which is particle 1 and 2
							    int rand_num;
								rand_num = rand() % 2 ;
								
								//looking at the angular distributions of the two protons for different p_missing regions
								if (p_proton_proj_abs_val < 250){
									h2_theta1_theta2_10b_lower250->Fill(((v_CALIFA_evts_cms_energy[rand_num][1])/PI)*180,((v_CALIFA_evts_cms_energy[1-rand_num][1])/PI)*180);
									h2_n_f_n_s_10b_lower250->Fill((v_CALIFA_evts_cms_energy[rand_num][3])/1000.,(v_CALIFA_evts_cms_energy[rand_num][4]/1000.));
									h2_n_f_n_s_10b_lower250->Fill((v_CALIFA_evts_cms_energy[1-rand_num][3])/1000.,(v_CALIFA_evts_cms_energy[1-rand_num][4]/1000.));
									if (tpat == 4){
										h2_theta1_theta2_10b_lower250_neu->Fill(((v_CALIFA_evts_cms_energy[rand_num][1])/PI)*180,((v_CALIFA_evts_cms_energy[1-rand_num][1])/PI)*180);
										}
									}
								if (p_proton_proj_abs_val > 250 && p_proton_proj_abs_val < 650){
									h2_theta1_theta2_10b_lower650->Fill(((v_CALIFA_evts_cms_energy[rand_num][1])/PI)*180,((v_CALIFA_evts_cms_energy[1-rand_num][1])/PI)*180);
									h2_n_f_n_s_10b_lower650->Fill((v_CALIFA_evts_cms_energy[rand_num][3])/1000.,(v_CALIFA_evts_cms_energy[rand_num][4]/1000.));
									h2_n_f_n_s_10b_lower650->Fill((v_CALIFA_evts_cms_energy[1-rand_num][3])/1000.,(v_CALIFA_evts_cms_energy[1-rand_num][4]/1000.));
									if (tpat == 4){
										h2_theta1_theta2_10b_lower650_neu->Fill(((v_CALIFA_evts_cms_energy[rand_num][1])/PI)*180,((v_CALIFA_evts_cms_energy[1-rand_num][1])/PI)*180);
										}

									}
								if (p_proton_proj_abs_val > 650){
									h2_theta1_theta2_10b_higher650->Fill(((v_CALIFA_evts_cms_energy[rand_num][1])/PI)*180,((v_CALIFA_evts_cms_energy[1-rand_num][1])/PI)*180);	
									h1_angle_pi_p_n_xy_plane_mom_cut->Fill(cos(p3_missing_x_y.Angle(p3_p_i_x_y)));
									h1_angle_src_vs_10B_cos_xy_plane_mom_cut->Fill(cos(((p3_p_i_x_y-p3_missing_x_y)*0.5).Angle(p3_10B_x_y)));
									h2_n_f_n_s_10b_higher650->Fill((v_CALIFA_evts_cms_energy[rand_num][3])/1000.,(v_CALIFA_evts_cms_energy[rand_num][4]/1000.));
									h2_n_f_n_s_10b_higher650->Fill((v_CALIFA_evts_cms_energy[1-rand_num][3])/1000.,(v_CALIFA_evts_cms_energy[1-rand_num][4]/1000.));
									if (tpat == 4){
										h2_theta1_theta2_10b_higher650_neu->Fill(((v_CALIFA_evts_cms_energy[rand_num][1])/PI)*180,((v_CALIFA_evts_cms_energy[1-rand_num][1])/PI)*180);
										}
									}
								
								h2_theta1_theta2_10b->Fill(((v_CALIFA_evts_cms_energy[rand_num][1])/PI)*180,((v_CALIFA_evts_cms_energy[1-rand_num][1])/PI)*180);
								h2_n_f_n_s_10b->Fill((v_CALIFA_evts_cms_energy[rand_num][3])/1000.,(v_CALIFA_evts_cms_energy[rand_num][4]/1000.));
								h2_n_f_n_s_10b->Fill((v_CALIFA_evts_cms_energy[1-rand_num][3])/1000.,(v_CALIFA_evts_cms_energy[1-rand_num][4]/1000.));
								if (tpat == 4){
									h2_theta1_theta2_10b_neu->Fill(((v_CALIFA_evts_cms_energy[rand_num][1])/PI)*180,((v_CALIFA_evts_cms_energy[1-rand_num][1])/PI)*180);
									}
								//making mass cut on the reconstucted neutron mass M_missing 
								if (M_missing > 850 && M_missing < 1100){
									h2_theta1_theta2_10b_m_cut->Fill(((v_CALIFA_evts_cms_energy[rand_num][1])/PI)*180,((v_CALIFA_evts_cms_energy[1-rand_num][1])/PI)*180);
									h1_p_missing_10B_m_cut->Fill(sqrt(pow(p_proton_proj[0],2)+pow(p_proton_proj[1],2)+pow(p_proton_proj[2],2)));
									//here making also cut on the reconstructed mass of the initial proton:
									if(p_proton_proj.Mag() > 844 && p_proton_proj.Mag() < 960){
										h2_theta1_theta2_10b_double_cut->Fill(((v_CALIFA_evts_cms_energy[rand_num][1])/PI)*180,((v_CALIFA_evts_cms_energy[1-rand_num][1])/PI)*180);
										} 
									}
								
								//here I calculate the mass the mass of the initial proton: p_missing = p1 + p2 -p_tg -> sqrt(pmissing^2) \
								//using also the y component
								Double_t M_mass_rec = p_proton_proj.Mag();
								h1_missing_mass_10B_rec->Fill(M_mass_rec);
								//use M_mass_rec as cut variable to check the momentum distribution...for eg 10B:
								if (M_mass_rec > 830 && M_mass_rec < 970){
									h1_p_missing_10B_p_i_cut->Fill(sqrt(pow(p_proton_proj[0],2)+pow(p_proton_proj[1],2)+pow(p_proton_proj[2],2)));
									}
								//ok, now vice versa, I use the momentum of the initial proton as cut variable. To see scr-pairs I cut on mom_p_i > 900
								if ((sqrt(pow(p_proton_proj[0],2)+pow(p_proton_proj[1],2)+pow(p_proton_proj[2],2))) > 900){
									h2_theta1_theta2_10b_miss_mom_cut->Fill(((v_CALIFA_evts_cms_energy[rand_num][1])/PI)*180,((v_CALIFA_evts_cms_energy[1-rand_num][1])/PI)*180);	
									h1_missing_mass_10B_rec_miss_cut->Fill(M_mass_rec);
									h1_missing_mass_10B_neutron_miss_cut->Fill(M_missing);
									}

							//gamma analysis for 12C(p,2p)10B reaction
							if (v_lorentz_GAMMAS_cms.size()){
                                h1_gamma_energyE_max_val_10B->Fill(v_lorentz_GAMMAS_cms[0][3]);
                                           }
				
								
							if (v_lorentz_PROTONS_lab.size() > 2){
								
								cout << "more than 3 protons, really? 10B " << endl;
								}
							}
						

                        }
//############################################################################################################
//--------------END OF 10B ANALYSIS---------------------------------------------------------------------------
//############################################################################################################

                        
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
                        h1_tof_detnr_12c->Fill(i+1);
                        h1_mw3_fy_12c->Fill(raw_mwpc3);
                        h2_tof_path_12c->Fill(pathlength_precise,time_target_tof);
                        h2_radius_path_12c->Fill(rho,time_target_tof);
                        h1_tofwall_posy_medium[i]->Fill(tof_diff_up_down);
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
//l->Add(h1_gamma_energyE_max_val_11B);
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
l->Add(h1_cos_low_gamma);
l->Add(h2_theta_sum_vs_phi_diff_11b);
l->Add(h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom);
l->Add(h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom200);
l->Add(h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom150);
l->Add(h1_m_mom_250_gamma_spec);
l->Add(h1_m_mom_200_gamma_spec);
l->Add(h1_m_mom_150_gamma_spec);
l->Add(h1_cluster_gamma_11b);
l->Add(h1_cluster_proton_11b);
l->Add(h2_radius_vs_beta_12c);
l->Add(h2_radius_vs_tof_12c);
l->Add(h1_abs_phi_diff_11b);
l->Add(h2_long_mom_p1p2_long_mom11b_large_gamma);
//for (Int_t i = 0; i < NbDets; i++){
//l->Add(h1_tofwall_posy_fast[i]);
//l->Add(h1_tofwall_posy_medium[i]);
//l->Add(h1_tofwall_posy_slow[i]);
//}
l->Add(h1_missing_energy_11B_p2p);
l->Add(h1_missing_energy_11B_p2p_cut);
l->Add(h1_missing_mass_11B_p2p_cut);
l->Add(h2_missing_energy_vs_t1_t2_11B_p2p);
l->Add(h1_p_missing_11B);
l->Add(h1_p_missing_11B_cut);
l->Add(h1_p_missing_11B_z);
l->Add(h1_p_missing_11B_z_cut);
l->Add(h1_p_missing_11B_x_cut);
l->Add(h1_p_missing_11B_y_cut);
l->Add(h1_p_missing_11B_x);
l->Add(h1_missing_mass_10B);
l->Add(h1_missing_mass_10B_neutron_miss_cut);
l->Add(h1_missing_mass_10B_rec);
l->Add(h1_missing_mass_10B_rec_miss_cut);
l->Add(h1_missing_mass_11B);
l->Add(h1_p_missing_10B);
l->Add(h2_missing_energy_vs_t1_t2_10B);
l->Add(h1_angle_pi_11B);
l->Add(h1_angle_pi_11B_cos);
l->Add(h1_angle_pi_11B_cos_cut);
l->Add(h1_angle_pi_11B_cos_cut_fake_y);
l->Add(h1_angle_pi_11B_cos_cut_plane);
l->Add(h1_angle_pi_pn_10B);
l->Add(h1_angle_pi_pn_10B_cos);
l->Add(h1_angle_src_vs_10B_cos);
l->Add(h1_gamma_energyE_max_val_11B);
l->Add(h1_gamma_energyE_max_val_11B_cut);
l->Add(h1_gamma_energyE_max_val_10B);
l->Add(h1_gamma_energyE_max_val_11B_lab);
l->Add(h2_gamma_energy_11B_vs_mult);
l->Add(h2_gamma_energy_11B_vs_mult_cut);
l->Add(h1_gamma_energyE_max_val_11B_m_cut);
l->Add(h1_gamma_energyE_sum_11B);
l->Add(h1_gamma_energyE_sum_11B_cut);
l->Add(h2_gamma_energy_11_vs_angle);
l->Add(h2_theta1_theta2_10b);
l->Add(h2_theta1_theta2_10b_m_cut);
l->Add(h2_theta1_theta2_10b_double_cut);
l->Add(h2_theta1_theta2_10b_miss_mom_cut);
l->Add(h2_p_missing_vs_p_11B_cut);
l->Add(h2_p_missing_vs_p_11B_cut_fake_y);
l->Add(h1_excit_energy_11B_cut);
l->Add(h1_p_missing_10B_m_cut);
l->Add(h1_p_missing_10B_p_i_cut);
l->Add(h1_gamma_energyE_max_val_11Bframe_cut);
l->Add(h2_corr_y_proton_fragm);
l->Add(h2_corr_x_proton_fragm);
l->Add(h2_corr_y_proton_fragm_tj);
l->Add(h2_missing_e_2meth_11B_p2p_cut);
l->Add(h1_angle_pi_11B_cos_hcut_plane);
l->Add(h2_corr_y_proton_fragm_cos_less);
l->Add(h2_corr_y_proton_fragm_cos_more);
l->Add(h2_corr_y_mws);
l->Add(h2_corr_y_mws_minusmw3);
l->Add(h2_E_sum_vs_theta_sum);
l->Add(h1_mom11B_z_cut);
l->Add(h1_mom11B_x_cut);
l->Add(h1_mom11B_y_cut);
l->Add(h1_mom11B_cut);
l->Add(h2_theta1_theta2_10b_lower250);
l->Add(h2_theta1_theta2_10b_lower650);
l->Add(h2_theta1_theta2_10b_higher650);
l->Add(h1_angle_pi_p_n_xy_plane);
l->Add(h1_angle_pi_p_n_xy_plane_mom_cut);
l->Add(h1_angle_src_vs_10B_cos_xy_plane);
l->Add(h1_angle_src_vs_10B_cos_xy_plane_mom_cut);
l->Add(h1_wr_first_wr_master);
l->Add(h1_wr_second_wr_master);
l->Add(h2_theta1_theta2_10b_neu);
l->Add(h2_theta1_theta2_10b_lower250_neu);
l->Add(h2_theta1_theta2_10b_lower650_neu);
l->Add(h2_theta1_theta2_10b_higher650_neu);
l->Add(h2_n_f_n_s_10b);
l->Add(h2_n_f_n_s_10b_lower250);
l->Add(h2_n_f_n_s_10b_lower650);
l->Add(h2_n_f_n_s_10b_higher650);
//h2_reco_mass_angle_comb_cut->LabelsDeflate("Y");
//h2_reco_mass_angle_comb_cut->LabelsOption("v");
l->Add(h2_reco_mass_angle_comb_cut);
l->Add(h1_reco_mass_11B_p2p_cut_two_iPhos);
l->Add(h1_reco_mass_11B_p2p_cut_iPhos_barrel);
l->Add(h1_reco_mass_11B_p2p_cut_two_barrel);
l->Add(h2_reco_mass_11B_iphos_angle);
l->Add(h2_reco_mass_11B_barrel_angle);
l->Add(h1_mom_abs_10B);
l->Add(h1_mom_x_10B);
l->Add(h1_mom_y_10B);
l->Add(h1_mom_z_10B);
l->Add(h2_sum_mom_vs_cos_p_i_p_n_xy);
l->Add(h1_angle_src_vs_10B_cos_3d);
l->Add(h1_angle_src_vs_10B_cos_3d_paper_cut);
l->Add(h2_e_miss_vs_theta_sum_p2p_11B);
l->Add(h2_e_miss_vs_theta_sum_p2p_10B);

l->Write("histlist", TObject::kSingleKey);

delete SofMwpc3HitData;
delete SofSciTcalData;
delete SofToFWTcalData;
delete SofMwpc0HitData;
delete SofMwpc1HitData;
delete SofMwpc2HitData;
delete SofToFWMappedData;
delete SofTwimHitData;
delete CalifaHitData;
cout << "end of process, successful" <<endl;
}



