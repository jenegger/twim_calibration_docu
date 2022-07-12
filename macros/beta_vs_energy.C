//insert all headers
#include "TTree.h"
#include "TFile.h"
#include <math.h>
#include <string>
#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TClonesArray.h"
#include "TChain.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include "/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/macros/mw_position.C"
#include "/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/macros/path_tof_calc.C"
#include "/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/macros/after_glad_pos.C"
//Assignment of the TWIM sections:
////0: Messel Down
////1: Messel Up
////2: Wixhausen Up
////3: Wixhausen Down
using namespace std;
//use geometry from S444 in 2020:TODO is that somewhat right??
//const double twim_anodes_width = 550./18.;
//const double pos_mw1_z = 0.; //chosen by myself
//const double pos_mw2_z = 660.; //in mm
//const double pos_twim_entrance = 81 + 1.5*twim_anodes_width ;
//const double twim_anode_frish_grid_dist = 110. ;//mm


//new geometry from Audrey for S455 U238
const double twim_anodes_width = 400./16; //25mm
const double pos_mw1_z = 0; //chosen by myself it's actually 300 mm from target
const double pos_mw2_z = 935 - 300; //mm
const double pos_twim_entrance = (610-200+12.5) -300;
const double twim_anode_frish_grid_dist = 110. ;//mm
vector<vector<double> > dummy_vec {
                {-1000,-1000,-1000},
                {-1000,-1000,-1000}
                };
 

void beta_vs_energy(const string& input_str){

Long64_t entries_twim = 0;
Long64_t entries_tof = 0;

string fname = string("/scratch5/ge37liw/unpacked_newroot_04_2022/s455_03_273_") + input_str + string("_unpacked.root");
string cal_fname = string("/scratch5/ge37liw/calibrated_rand_angles/s455_03_273_") + input_str + string("_calibrated.root");
TFile* file_tsplines_one = TFile::Open("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/macros/spline_position_messel_and_wixhausen.root");
TList* spline_list_one = (TList*)file_tsplines_one->Get("list_tspline");
TSpline3* spline_messel_sec0 = (TSpline3*)spline_list_one->FindObject("Position on 11th anode in music vs Calibrated Energy messel with tcut,section0_pfx");
TSpline3* spline_messel_sec1 = (TSpline3*)spline_list_one->FindObject("Position on 11th anode in music vs Calibrated Energy messel with tcut,section1_pfx");
TSpline3* spline_wix_sec2 = (TSpline3*)spline_list_one->FindObject("Position on 11th anode in music vs Calibrated Energy Wixhausen with tcut,section2_pfx");
TSpline3* spline_wix_sec3 = (TSpline3*)spline_list_one->FindObject("Position on 11th anode in music vs Calibrated Energy Wixhausen with tcut,section3_pfx");
const double mean_pos_sec0 = 226766;
const double mean_pos_sec1 = 239296;
const double mean_pos_sec2 = 274225;
const double mean_pos_sec3 = 262821;

const char* char_fname= fname.c_str();
const char* char_cal_fname = cal_fname.c_str();
TChain* chain = new TChain("evt");
chain->Reset();
chain->Add(char_fname);
TChain*  chain2 = new TChain("evt");
chain2->Reset();
chain2->Add(char_cal_fname);
chain->AddFriend(chain2);
Long64_t nevents = chain->GetEntries();

char hist_name[500];

//------ HISTOS ----------------------------------

TH2D* h2_beta_vs_energy;
sprintf(hist_name,"Beta vs Calibrated Energy");
h2_beta_vs_energy = new TH2D(hist_name,hist_name,500,0.73,0.83,950,0,450000);
h2_beta_vs_energy->GetXaxis()->SetTitle("Beta");
h2_beta_vs_energy->GetYaxis()->SetTitle("Calibrated Energy [a.u.]");
h2_beta_vs_energy->GetXaxis()->CenterTitle(true);
h2_beta_vs_energy->GetYaxis()->CenterTitle(true);
h2_beta_vs_energy->GetYaxis()->SetLabelSize(0.045);
h2_beta_vs_energy->GetYaxis()->SetTitleSize(0.045);

//do the above analysis section-wise
//section0
TH2D* h2_beta_vs_energy_sec0;
sprintf(hist_name,"Beta vs Calibrated Energy,section0");
h2_beta_vs_energy_sec0 = new TH2D(hist_name,hist_name,500,0.73,0.83,950,0,450000);
h2_beta_vs_energy_sec0->GetXaxis()->SetTitle("Beta");
h2_beta_vs_energy_sec0->GetYaxis()->SetTitle("Calibrated Energy [a.u.]");
h2_beta_vs_energy_sec0->GetXaxis()->CenterTitle(true);
h2_beta_vs_energy_sec0->GetYaxis()->CenterTitle(true);
h2_beta_vs_energy_sec0->GetYaxis()->SetLabelSize(0.045);
h2_beta_vs_energy_sec0->GetYaxis()->SetTitleSize(0.045);

//section1
TH2D* h2_beta_vs_energy_sec1;
sprintf(hist_name,"Beta vs Calibrated Energy,section1");
h2_beta_vs_energy_sec1 = new TH2D(hist_name,hist_name,500,0.73,0.83,950,0,450000);
h2_beta_vs_energy_sec1->GetXaxis()->SetTitle("Beta");
h2_beta_vs_energy_sec1->GetYaxis()->SetTitle("Calibrated Energy [a.u.]");
h2_beta_vs_energy_sec1->GetXaxis()->CenterTitle(true);
h2_beta_vs_energy_sec1->GetYaxis()->CenterTitle(true);
h2_beta_vs_energy_sec1->GetYaxis()->SetLabelSize(0.045);
h2_beta_vs_energy_sec1->GetYaxis()->SetTitleSize(0.045);

//section2
TH2D* h2_beta_vs_energy_sec2;
sprintf(hist_name,"Beta vs Calibrated Energy,section2");
h2_beta_vs_energy_sec2 = new TH2D(hist_name,hist_name,500,0.73,0.83,950,0,450000);
h2_beta_vs_energy_sec2->GetXaxis()->SetTitle("Beta");
h2_beta_vs_energy_sec2->GetYaxis()->SetTitle("Calibrated Energy [a.u.]");
h2_beta_vs_energy_sec2->GetXaxis()->CenterTitle(true);
h2_beta_vs_energy_sec2->GetYaxis()->CenterTitle(true);
h2_beta_vs_energy_sec2->GetYaxis()->SetLabelSize(0.045);
h2_beta_vs_energy_sec2->GetYaxis()->SetTitleSize(0.045);

//section3
TH2D* h2_beta_vs_energy_sec3;
sprintf(hist_name,"Beta vs Calibrated Energy,section3");
h2_beta_vs_energy_sec3 = new TH2D(hist_name,hist_name,500,0.73,0.83,950,0,450000);
h2_beta_vs_energy_sec3->GetXaxis()->SetTitle("Beta");
h2_beta_vs_energy_sec3->GetYaxis()->SetTitle("Calibrated Energy [a.u.]");
h2_beta_vs_energy_sec3->GetXaxis()->CenterTitle(true);
h2_beta_vs_energy_sec3->GetYaxis()->CenterTitle(true);
h2_beta_vs_energy_sec3->GetYaxis()->SetLabelSize(0.045);
h2_beta_vs_energy_sec3->GetYaxis()->SetTitleSize(0.045);
//end section wise analysis of beta correction


TH2D* h2_beta_vs_energy_corr;
sprintf(hist_name,"Beta vs Calibrated Energy beta corrected");
h2_beta_vs_energy_corr = new TH2D(hist_name,hist_name,500,0.73,0.83,950,0,450000);
h2_beta_vs_energy_corr->GetXaxis()->SetTitle("Beta");
h2_beta_vs_energy_corr->GetYaxis()->SetTitle("Calibrated Energy [a.u.]");
h2_beta_vs_energy_corr->GetXaxis()->CenterTitle(true);
h2_beta_vs_energy_corr->GetYaxis()->CenterTitle(true);
h2_beta_vs_energy_corr->GetYaxis()->SetLabelSize(0.045);
h2_beta_vs_energy_corr->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_beta_vs_energy_corr_sec0;
sprintf(hist_name,"Beta vs Calibrated Energy beta corrected,section0");
h2_beta_vs_energy_corr_sec0 = new TH2D(hist_name,hist_name,500,0.73,0.83,950,0,450000);
h2_beta_vs_energy_corr_sec0->GetXaxis()->SetTitle("Beta");
h2_beta_vs_energy_corr_sec0->GetYaxis()->SetTitle("Calibrated Energy [a.u.]");
h2_beta_vs_energy_corr_sec0->GetXaxis()->CenterTitle(true);
h2_beta_vs_energy_corr_sec0->GetYaxis()->CenterTitle(true);
h2_beta_vs_energy_corr_sec0->GetYaxis()->SetLabelSize(0.045);
h2_beta_vs_energy_corr_sec0->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_beta_vs_energy_corr_sec1;
sprintf(hist_name,"Beta vs Calibrated Energy beta corrected,section1");
h2_beta_vs_energy_corr_sec1 = new TH2D(hist_name,hist_name,500,0.73,0.83,950,0,450000);
h2_beta_vs_energy_corr_sec1->GetXaxis()->SetTitle("Beta");
h2_beta_vs_energy_corr_sec1->GetYaxis()->SetTitle("Calibrated Energy [a.u.]");
h2_beta_vs_energy_corr_sec1->GetXaxis()->CenterTitle(true);
h2_beta_vs_energy_corr_sec1->GetYaxis()->CenterTitle(true);
h2_beta_vs_energy_corr_sec1->GetYaxis()->SetLabelSize(0.045);
h2_beta_vs_energy_corr_sec1->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_beta_vs_energy_corr_sec2;
sprintf(hist_name,"Beta vs Calibrated Energy beta corrected,section2");
h2_beta_vs_energy_corr_sec2 = new TH2D(hist_name,hist_name,500,0.73,0.83,950,0,450000);
h2_beta_vs_energy_corr_sec2->GetXaxis()->SetTitle("Beta");
h2_beta_vs_energy_corr_sec2->GetYaxis()->SetTitle("Calibrated Energy [a.u.]");
h2_beta_vs_energy_corr_sec2->GetXaxis()->CenterTitle(true);
h2_beta_vs_energy_corr_sec2->GetYaxis()->CenterTitle(true);
h2_beta_vs_energy_corr_sec2->GetYaxis()->SetLabelSize(0.045);
h2_beta_vs_energy_corr_sec2->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_beta_vs_energy_corr_sec3;
sprintf(hist_name,"Beta vs Calibrated Energy beta corrected,section3");
h2_beta_vs_energy_corr_sec3 = new TH2D(hist_name,hist_name,500,0.73,0.83,950,0,450000);
h2_beta_vs_energy_corr_sec3->GetXaxis()->SetTitle("Beta");
h2_beta_vs_energy_corr_sec3->GetYaxis()->SetTitle("Calibrated Energy [a.u.]");
h2_beta_vs_energy_corr_sec3->GetXaxis()->CenterTitle(true);
h2_beta_vs_energy_corr_sec3->GetYaxis()->CenterTitle(true);
h2_beta_vs_energy_corr_sec3->GetYaxis()->SetLabelSize(0.045);
h2_beta_vs_energy_corr_sec3->GetYaxis()->SetTitleSize(0.045);
//TF1* beta_func = new TF1("beta_func","198059*(TMath::Power(x,-5./3.)) + 32590.6");
//const static double mean_ene = 324231.;
const static double mean_ene = 250000.;
TF1* beta_func_sec0 = new TF1("beta_func_sec0","111559*(TMath::Power(x,-5./3.)) + 35942.6");
TF1* beta_func_sec1 = new TF1("beta_func_sec1","119028*(TMath::Power(x,-5./3.)) + 43631");
TF1* beta_func_sec2 = new TF1("beta_func_sec2","107890*(TMath::Power(x,-5./3.)) + 47938.4");
TF1* beta_func_sec3 = new TF1("beta_func_sec3","110713*(TMath::Power(x,-5./3.)) + 44166.2");

TH2D* h2_x_music_vs_energy_messel;
sprintf(hist_name,"Position on 11th anode in music vs Calibrated Energy messel");
h2_x_music_vs_energy_messel = new TH2D(hist_name,hist_name,600,0.,120,950,0,450000);
h2_x_music_vs_energy_messel->GetXaxis()->SetTitle("position Messel side");
h2_x_music_vs_energy_messel->GetYaxis()->SetTitle("Calibrated Energy [a.u.]");
h2_x_music_vs_energy_messel->GetXaxis()->CenterTitle(true);
h2_x_music_vs_energy_messel->GetYaxis()->CenterTitle(true);
h2_x_music_vs_energy_messel->GetYaxis()->SetLabelSize(0.045);
h2_x_music_vs_energy_messel->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_x_music_vs_energy_messel_sec0;
sprintf(hist_name,"Position on 11th anode in music vs Calibrated Energy messel,section0");
h2_x_music_vs_energy_messel_sec0 = new TH2D(hist_name,hist_name,600,0.,120,950,0,450000);
h2_x_music_vs_energy_messel_sec0->GetXaxis()->SetTitle("position Messel side");
h2_x_music_vs_energy_messel_sec0->GetYaxis()->SetTitle("Calibrated Energy [a.u.]");
h2_x_music_vs_energy_messel_sec0->GetXaxis()->CenterTitle(true);
h2_x_music_vs_energy_messel_sec0->GetYaxis()->CenterTitle(true);
h2_x_music_vs_energy_messel_sec0->GetYaxis()->SetLabelSize(0.045);
h2_x_music_vs_energy_messel_sec0->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_x_music_vs_energy_messel_sec1;
sprintf(hist_name,"Position on 11th anode in music vs Calibrated Energy messel,section1");
h2_x_music_vs_energy_messel_sec1 = new TH2D(hist_name,hist_name,600,0.,120,950,0,450000);
h2_x_music_vs_energy_messel_sec1->GetXaxis()->SetTitle("position Messel side");
h2_x_music_vs_energy_messel_sec1->GetYaxis()->SetTitle("Calibrated Energy [a.u.]");
h2_x_music_vs_energy_messel_sec1->GetXaxis()->CenterTitle(true);
h2_x_music_vs_energy_messel_sec1->GetYaxis()->CenterTitle(true);
h2_x_music_vs_energy_messel_sec1->GetYaxis()->SetLabelSize(0.045);
h2_x_music_vs_energy_messel_sec1->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_x_music_vs_energy_messel_corr;
sprintf(hist_name,"Position on 11th anode in music vs Calibrated Energy messel corrected");
h2_x_music_vs_energy_messel_corr = new TH2D(hist_name,hist_name,600,0.,120,950,0,450000);
h2_x_music_vs_energy_messel_corr->GetXaxis()->SetTitle("position Messel side");
h2_x_music_vs_energy_messel_corr->GetYaxis()->SetTitle("Calibrated Energy [a.u.]");
h2_x_music_vs_energy_messel_corr->GetXaxis()->CenterTitle(true);
h2_x_music_vs_energy_messel_corr->GetYaxis()->CenterTitle(true);
h2_x_music_vs_energy_messel_corr->GetYaxis()->SetLabelSize(0.045);
h2_x_music_vs_energy_messel_corr->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_charge_music_corr_messel_only_beta;
sprintf(hist_name,"Energy deposited TWIM, Messel side(section 0 & 1) only beta correction");
h1_charge_music_corr_messel_only_beta = new TH1D(hist_name,hist_name,6000,0,600000);
h1_charge_music_corr_messel_only_beta->GetXaxis()->SetTitle("Cal. Energy TWIM Messel");
h1_charge_music_corr_messel_only_beta->GetYaxis()->SetTitle("Counts");
h1_charge_music_corr_messel_only_beta->GetXaxis()->CenterTitle(true);
h1_charge_music_corr_messel_only_beta->GetYaxis()->CenterTitle(true);
h1_charge_music_corr_messel_only_beta->GetYaxis()->SetLabelSize(0.045);
h1_charge_music_corr_messel_only_beta->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_charge_music_corr_messel;
sprintf(hist_name,"Energy deposited TWIM, Messel side(section 0 & 1)");
h1_charge_music_corr_messel = new TH1D(hist_name,hist_name,6000,0,600000);
h1_charge_music_corr_messel->GetXaxis()->SetTitle("Cal. Energy TWIM Messel");
h1_charge_music_corr_messel->GetYaxis()->SetTitle("Counts");
h1_charge_music_corr_messel->GetXaxis()->CenterTitle(true);
h1_charge_music_corr_messel->GetYaxis()->CenterTitle(true);
h1_charge_music_corr_messel->GetYaxis()->SetLabelSize(0.045);
h1_charge_music_corr_messel->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_charge_music_corr_sec0;
sprintf(hist_name,"Energy deposited TWIM section0, beta and pos. corrected");
h1_charge_music_corr_sec0 = new TH1D(hist_name,hist_name,6000,0,600000);
h1_charge_music_corr_sec0->GetXaxis()->SetTitle("Cal. Energy section0");
h1_charge_music_corr_sec0->GetYaxis()->SetTitle("Counts");
h1_charge_music_corr_sec0->GetXaxis()->CenterTitle(true);
h1_charge_music_corr_sec0->GetYaxis()->CenterTitle(true);
h1_charge_music_corr_sec0->GetYaxis()->SetLabelSize(0.045);
h1_charge_music_corr_sec0->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_charge_music_corr_sec1;
sprintf(hist_name,"Energy deposited TWIM section1, beta and pos. corrected");
h1_charge_music_corr_sec1 = new TH1D(hist_name,hist_name,6000,0,600000);
h1_charge_music_corr_sec1->GetXaxis()->SetTitle("Cal. Energy section1");
h1_charge_music_corr_sec1->GetYaxis()->SetTitle("Counts");
h1_charge_music_corr_sec1->GetXaxis()->CenterTitle(true);
h1_charge_music_corr_sec1->GetYaxis()->CenterTitle(true);
h1_charge_music_corr_sec1->GetYaxis()->SetLabelSize(0.045);
h1_charge_music_corr_sec1->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_charge_music_corr_sec2;
sprintf(hist_name,"Energy deposited TWIM section2, beta and pos. corrected");
h1_charge_music_corr_sec2 = new TH1D(hist_name,hist_name,6000,0,600000);
h1_charge_music_corr_sec2->GetXaxis()->SetTitle("Cal. Energy section2");
h1_charge_music_corr_sec2->GetYaxis()->SetTitle("Counts");
h1_charge_music_corr_sec2->GetXaxis()->CenterTitle(true);
h1_charge_music_corr_sec2->GetYaxis()->CenterTitle(true);
h1_charge_music_corr_sec2->GetYaxis()->SetLabelSize(0.045);
h1_charge_music_corr_sec2->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_charge_music_corr_sec3;
sprintf(hist_name,"Energy deposited TWIM section3, beta and pos. corrected");
h1_charge_music_corr_sec3 = new TH1D(hist_name,hist_name,6000,0,600000);
h1_charge_music_corr_sec3->GetXaxis()->SetTitle("Cal. Energy section3");
h1_charge_music_corr_sec3->GetYaxis()->SetTitle("Counts");
h1_charge_music_corr_sec3->GetXaxis()->CenterTitle(true);
h1_charge_music_corr_sec3->GetYaxis()->CenterTitle(true);
h1_charge_music_corr_sec3->GetYaxis()->SetLabelSize(0.045);
h1_charge_music_corr_sec3->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_charge_music_corr_messel_full_cal;
sprintf(hist_name,"Energy deposited TWIM, Messel side(section 0 & 1) now fully calibrated");
h1_charge_music_corr_messel_full_cal = new TH1D(hist_name,hist_name,6000,0,600000);
h1_charge_music_corr_messel_full_cal->GetXaxis()->SetTitle("Cal. Energy TWIM Messel");
h1_charge_music_corr_messel_full_cal->GetYaxis()->SetTitle("Counts");
h1_charge_music_corr_messel_full_cal->GetXaxis()->CenterTitle(true);
h1_charge_music_corr_messel_full_cal->GetYaxis()->CenterTitle(true);
h1_charge_music_corr_messel_full_cal->GetYaxis()->SetLabelSize(0.045);
h1_charge_music_corr_messel_full_cal->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_charge_music_corr_messel_full_cal_pol;
sprintf(hist_name,"Energy deposited TWIM, Messel side(section 0 & 1) now fully calibrated with pol");
h1_charge_music_corr_messel_full_cal_pol = new TH1D(hist_name,hist_name,6000,0,600000);
h1_charge_music_corr_messel_full_cal_pol->GetXaxis()->SetTitle("Cal. Energy TWIM Messel");
h1_charge_music_corr_messel_full_cal_pol->GetYaxis()->SetTitle("Counts");
h1_charge_music_corr_messel_full_cal_pol->GetXaxis()->CenterTitle(true);
h1_charge_music_corr_messel_full_cal_pol->GetYaxis()->CenterTitle(true);
h1_charge_music_corr_messel_full_cal_pol->GetYaxis()->SetLabelSize(0.045);
h1_charge_music_corr_messel_full_cal_pol->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_charge_music_corr_wix_only_beta;
sprintf(hist_name,"Energy deposited TWIM, Wixhausen side(section 2 & 3) only beta correction");
h1_charge_music_corr_wix_only_beta = new TH1D(hist_name,hist_name,6000,0,600000);
h1_charge_music_corr_wix_only_beta->GetXaxis()->SetTitle("Cal. Energy TWIM Wixhausen");
h1_charge_music_corr_wix_only_beta->GetYaxis()->SetTitle("Counts");
h1_charge_music_corr_wix_only_beta->GetXaxis()->CenterTitle(true);
h1_charge_music_corr_wix_only_beta->GetYaxis()->CenterTitle(true);
h1_charge_music_corr_wix_only_beta->GetYaxis()->SetLabelSize(0.045);
h1_charge_music_corr_wix_only_beta->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_charge_music_corr_wix;
sprintf(hist_name,"Energy deposited TWIM, Wixhausen side(section 2 & 3)");
h1_charge_music_corr_wix = new TH1D(hist_name,hist_name,6000,0,600000);
h1_charge_music_corr_wix->GetXaxis()->SetTitle("Cal. Energy TWIM Wixhausen");
h1_charge_music_corr_wix->GetYaxis()->SetTitle("Counts");
h1_charge_music_corr_wix->GetXaxis()->CenterTitle(true);
h1_charge_music_corr_wix->GetYaxis()->CenterTitle(true);
h1_charge_music_corr_wix->GetYaxis()->SetLabelSize(0.045);
h1_charge_music_corr_wix->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_charge_music_corr_wix_full_cal;
sprintf(hist_name,"Energy deposited TWIM, Wixhausen side(section 2 & 3) now fully calibrated");
h1_charge_music_corr_wix_full_cal = new TH1D(hist_name,hist_name,6000,0,600000);
h1_charge_music_corr_wix_full_cal->GetXaxis()->SetTitle("Cal. Energy TWIM Wixhausen");
h1_charge_music_corr_wix_full_cal->GetYaxis()->SetTitle("Counts");
h1_charge_music_corr_wix_full_cal->GetXaxis()->CenterTitle(true);
h1_charge_music_corr_wix_full_cal->GetYaxis()->CenterTitle(true);
h1_charge_music_corr_wix_full_cal->GetYaxis()->SetLabelSize(0.045);
h1_charge_music_corr_wix_full_cal->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_charge_music_corr_wix_full_cal_pol;
sprintf(hist_name,"Energy deposited TWIM, Wixhausen side(section 2 & 3) now fully calibrated with pol");
h1_charge_music_corr_wix_full_cal_pol = new TH1D(hist_name,hist_name,6000,0,600000);
h1_charge_music_corr_wix_full_cal_pol->GetXaxis()->SetTitle("Cal. Energy TWIM Wixhausen");
h1_charge_music_corr_wix_full_cal_pol->GetYaxis()->SetTitle("Counts");
h1_charge_music_corr_wix_full_cal_pol->GetXaxis()->CenterTitle(true);
h1_charge_music_corr_wix_full_cal_pol->GetYaxis()->CenterTitle(true);
h1_charge_music_corr_wix_full_cal_pol->GetYaxis()->SetLabelSize(0.045);
h1_charge_music_corr_wix_full_cal_pol->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_charge_music_corr_sum;
sprintf(hist_name,"Energy deposited in TWIM, summed up, beta and pos. calibrated");
h1_charge_music_corr_sum = new TH1D(hist_name,hist_name,6000,0,600000);
h1_charge_music_corr_sum->GetXaxis()->SetTitle("Cal. Energy TWIM");
h1_charge_music_corr_sum->GetYaxis()->SetTitle("Counts");
h1_charge_music_corr_sum->GetXaxis()->CenterTitle(true);
h1_charge_music_corr_sum->GetYaxis()->CenterTitle(true);
h1_charge_music_corr_sum->GetYaxis()->SetLabelSize(0.045);
h1_charge_music_corr_sum->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_charge_music_corr_sum_only_beta;
sprintf(hist_name,"Energy deposited in TWIM, summed up, only beta corrected");
h1_charge_music_corr_sum_only_beta = new TH1D(hist_name,hist_name,6000,0,600000);
h1_charge_music_corr_sum_only_beta->GetXaxis()->SetTitle("Cal. Energy TWIM");
h1_charge_music_corr_sum_only_beta->GetYaxis()->SetTitle("Counts");
h1_charge_music_corr_sum_only_beta->GetXaxis()->CenterTitle(true);
h1_charge_music_corr_sum_only_beta->GetYaxis()->CenterTitle(true);
h1_charge_music_corr_sum_only_beta->GetYaxis()->SetLabelSize(0.045);
h1_charge_music_corr_sum_only_beta->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_charge_music_corr_sum_only_messel;
sprintf(hist_name,"Energy deposited in TWIM, summed up, beta and pos. calibrated only messel side events");
h1_charge_music_corr_sum_only_messel = new TH1D(hist_name,hist_name,6000,0,600000);
h1_charge_music_corr_sum_only_messel->GetXaxis()->SetTitle("Cal. Energy TWIM");
h1_charge_music_corr_sum_only_messel->GetYaxis()->SetTitle("Counts");
h1_charge_music_corr_sum_only_messel->GetXaxis()->CenterTitle(true);
h1_charge_music_corr_sum_only_messel->GetYaxis()->CenterTitle(true);
h1_charge_music_corr_sum_only_messel->GetYaxis()->SetLabelSize(0.045);
h1_charge_music_corr_sum_only_messel->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_charge_music_corr_sum_only_wix;
sprintf(hist_name,"Energy deposited in TWIM, summed up, beta and pos. calibrated only wixhausen side events");
h1_charge_music_corr_sum_only_wix = new TH1D(hist_name,hist_name,6000,0,600000);
h1_charge_music_corr_sum_only_wix->GetXaxis()->SetTitle("Cal. Energy TWIM");
h1_charge_music_corr_sum_only_wix->GetYaxis()->SetTitle("Counts");
h1_charge_music_corr_sum_only_wix->GetXaxis()->CenterTitle(true);
h1_charge_music_corr_sum_only_wix->GetYaxis()->CenterTitle(true);
h1_charge_music_corr_sum_only_wix->GetYaxis()->SetLabelSize(0.045);
h1_charge_music_corr_sum_only_wix->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_charge_music_corr_sum_full_cal;
sprintf(hist_name,"Energy deposited in TWIM, summed up, beta and pos. calibrated, fully calibrated wixh, messel");
h1_charge_music_corr_sum_full_cal = new TH1D(hist_name,hist_name,6000,0,600000);
h1_charge_music_corr_sum_full_cal->GetXaxis()->SetTitle("Cal. Energy TWIM");
h1_charge_music_corr_sum_full_cal->GetYaxis()->SetTitle("Counts");
h1_charge_music_corr_sum_full_cal->GetXaxis()->CenterTitle(true);
h1_charge_music_corr_sum_full_cal->GetYaxis()->CenterTitle(true);
h1_charge_music_corr_sum_full_cal->GetYaxis()->SetLabelSize(0.045);
h1_charge_music_corr_sum_full_cal->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_charge_music_corr_sum_full_cal_pol;
sprintf(hist_name,"Energy deposited in TWIM, summed up, beta and pos. calibrated, fully calibrated wixh, messel with pol");
h1_charge_music_corr_sum_full_cal_pol = new TH1D(hist_name,hist_name,6000,0,600000);
h1_charge_music_corr_sum_full_cal_pol->GetXaxis()->SetTitle("Cal. Energy TWIM");
h1_charge_music_corr_sum_full_cal_pol->GetYaxis()->SetTitle("Counts");
h1_charge_music_corr_sum_full_cal_pol->GetXaxis()->CenterTitle(true);
h1_charge_music_corr_sum_full_cal_pol->GetYaxis()->CenterTitle(true);
h1_charge_music_corr_sum_full_cal_pol->GetYaxis()->SetLabelSize(0.045);
h1_charge_music_corr_sum_full_cal_pol->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_x_music_vs_energy_wix;
sprintf(hist_name,"Position on 11th anode in music vs Calibrated Energy Wixhausen");
h2_x_music_vs_energy_wix = new TH2D(hist_name,hist_name,600,-120,0.,950,0,450000);
h2_x_music_vs_energy_wix->GetXaxis()->SetTitle("position Messel side");
h2_x_music_vs_energy_wix->GetYaxis()->SetTitle("Calibrated Energy [a.u.]");
h2_x_music_vs_energy_wix->GetXaxis()->CenterTitle(true);
h2_x_music_vs_energy_wix->GetYaxis()->CenterTitle(true);
h2_x_music_vs_energy_wix->GetYaxis()->SetLabelSize(0.045);
h2_x_music_vs_energy_wix->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_x_music_vs_energy_wix_sec2;
sprintf(hist_name,"Position on 11th anode in music vs Calibrated Energy Wixhausen,section2");
h2_x_music_vs_energy_wix_sec2 = new TH2D(hist_name,hist_name,600,-120,0.,950,0,450000);
h2_x_music_vs_energy_wix_sec2->GetXaxis()->SetTitle("position Messel side");
h2_x_music_vs_energy_wix_sec2->GetYaxis()->SetTitle("Calibrated Energy [a.u.]");
h2_x_music_vs_energy_wix_sec2->GetXaxis()->CenterTitle(true);
h2_x_music_vs_energy_wix_sec2->GetYaxis()->CenterTitle(true);
h2_x_music_vs_energy_wix_sec2->GetYaxis()->SetLabelSize(0.045);
h2_x_music_vs_energy_wix_sec2->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_x_music_vs_energy_wix_sec3;
sprintf(hist_name,"Position on 11th anode in music vs Calibrated Energy Wixhausen,section3");
h2_x_music_vs_energy_wix_sec3 = new TH2D(hist_name,hist_name,600,-120,0.,950,0,450000);
h2_x_music_vs_energy_wix_sec3->GetXaxis()->SetTitle("position Messel side");
h2_x_music_vs_energy_wix_sec3->GetYaxis()->SetTitle("Calibrated Energy [a.u.]");
h2_x_music_vs_energy_wix_sec3->GetXaxis()->CenterTitle(true);
h2_x_music_vs_energy_wix_sec3->GetYaxis()->CenterTitle(true);
h2_x_music_vs_energy_wix_sec3->GetYaxis()->SetLabelSize(0.045);
h2_x_music_vs_energy_wix_sec3->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_x_music_vs_energy_wix_corr;
sprintf(hist_name,"Position on 11th anode in music vs Calibrated Energy Wixhausen corrected");
h2_x_music_vs_energy_wix_corr = new TH2D(hist_name,hist_name,600,-120,0.,950,0,450000);
h2_x_music_vs_energy_wix_corr->GetXaxis()->SetTitle("position Messel side");
h2_x_music_vs_energy_wix_corr->GetYaxis()->SetTitle("Calibrated Energy [a.u.]");
h2_x_music_vs_energy_wix_corr->GetXaxis()->CenterTitle(true);
h2_x_music_vs_energy_wix_corr->GetYaxis()->CenterTitle(true);
h2_x_music_vs_energy_wix_corr->GetYaxis()->SetLabelSize(0.045);
h2_x_music_vs_energy_wix_corr->GetYaxis()->SetTitleSize(0.045);

//messel side comparison fully corrected charges!!
TH2D* h2_x_energy_vs_energy_messel_corr;
sprintf(hist_name,"Charge Section 0 vs Charge Section 1 in Messel");
h2_x_energy_vs_energy_messel_corr = new TH2D(hist_name,hist_name,1000,-10,40,1000,-10,40);
h2_x_energy_vs_energy_messel_corr->GetXaxis()->SetTitle("Charge Section 0");
h2_x_energy_vs_energy_messel_corr->GetYaxis()->SetTitle("Charge Section 1");
h2_x_energy_vs_energy_messel_corr->GetXaxis()->CenterTitle(true);
h2_x_energy_vs_energy_messel_corr->GetYaxis()->CenterTitle(true);
h2_x_energy_vs_energy_messel_corr->GetYaxis()->SetLabelSize(0.045);
h2_x_energy_vs_energy_messel_corr->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_charge_sum_vs_diff_messel_corr;
sprintf(hist_name,"Sum Charge of Section 0&1 vs difference in Messel (section1-section0)");
h2_charge_sum_vs_diff_messel_corr = new TH2D(hist_name,hist_name,600,-15,15,800,10,50);
h2_charge_sum_vs_diff_messel_corr->GetXaxis()->SetTitle("Charge Difference");
h2_charge_sum_vs_diff_messel_corr->GetYaxis()->SetTitle("Charge Sum");
h2_charge_sum_vs_diff_messel_corr->GetXaxis()->CenterTitle(true);
h2_charge_sum_vs_diff_messel_corr->GetYaxis()->CenterTitle(true);
h2_charge_sum_vs_diff_messel_corr->GetYaxis()->SetLabelSize(0.045);
h2_charge_sum_vs_diff_messel_corr->GetYaxis()->SetTitleSize(0.045);

//wixhausen side comparison fully corrected charges!!
TH2D* h2_x_energy_vs_energy_wix_corr;
sprintf(hist_name,"Charge Section 2 vs Charge Section 3 in Wixhausen");
h2_x_energy_vs_energy_wix_corr = new TH2D(hist_name,hist_name,1000,-10,40,1000,-10,40);
h2_x_energy_vs_energy_wix_corr->GetXaxis()->SetTitle("Charge Section 2");
h2_x_energy_vs_energy_wix_corr->GetYaxis()->SetTitle("Charge Section 3");
h2_x_energy_vs_energy_wix_corr->GetXaxis()->CenterTitle(true);
h2_x_energy_vs_energy_wix_corr->GetYaxis()->CenterTitle(true);
h2_x_energy_vs_energy_wix_corr->GetYaxis()->SetLabelSize(0.045);
h2_x_energy_vs_energy_wix_corr->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_charge_sum_vs_diff_wix_corr;
sprintf(hist_name,"Sum Charge of Section 2&3 vs difference in Wixhausen(section3-section2)");
h2_charge_sum_vs_diff_wix_corr = new TH2D(hist_name,hist_name,600,-15,15,800,10,50);
h2_charge_sum_vs_diff_wix_corr->GetXaxis()->SetTitle("Charge Difference");
h2_charge_sum_vs_diff_wix_corr->GetYaxis()->SetTitle("Charge Sum");
h2_charge_sum_vs_diff_wix_corr->GetXaxis()->CenterTitle(true);
h2_charge_sum_vs_diff_wix_corr->GetYaxis()->CenterTitle(true);
h2_charge_sum_vs_diff_wix_corr->GetYaxis()->SetLabelSize(0.045);
h2_charge_sum_vs_diff_wix_corr->GetYaxis()->SetTitleSize(0.045);


//now fitting with Z = [0] +[1]*sqrt[E] +[2]*E
//wixhausen
TH1D* h1_one_charge_wixh_sqrt;
sprintf(hist_name,"Energy deposited in TWIM on wixhausen, beta and pos. calibrated, sqrt-fit");
h1_one_charge_wixh_sqrt = new TH1D(hist_name,hist_name,2000,0,200);
h1_one_charge_wixh_sqrt->GetXaxis()->SetTitle("Cal. Energy TWIM");
h1_one_charge_wixh_sqrt->GetYaxis()->SetTitle("Counts");
h1_one_charge_wixh_sqrt->GetXaxis()->CenterTitle(true);
h1_one_charge_wixh_sqrt->GetYaxis()->CenterTitle(true);
h1_one_charge_wixh_sqrt->GetYaxis()->SetLabelSize(0.045);
h1_one_charge_wixh_sqrt->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_sum_charge_wixh_sqrt;
sprintf(hist_name,"Energy deposited in TWIM on wixhausen SUM (for both fragments on wix), beta and pos. calibrated, sqrt-fit");
h1_sum_charge_wixh_sqrt = new TH1D(hist_name,hist_name,2000,0,200);
h1_sum_charge_wixh_sqrt->GetXaxis()->SetTitle("Cal. Energy TWIM");
h1_sum_charge_wixh_sqrt->GetYaxis()->SetTitle("Counts");
h1_sum_charge_wixh_sqrt->GetXaxis()->CenterTitle(true);
h1_sum_charge_wixh_sqrt->GetYaxis()->CenterTitle(true);
h1_sum_charge_wixh_sqrt->GetYaxis()->SetLabelSize(0.045);
h1_sum_charge_wixh_sqrt->GetYaxis()->SetTitleSize(0.045);

//now sum up without correction
TH1D* h1_sum_charge_wixh;
sprintf(hist_name,"Energy deposited in TWIM on wixhausen SUM (for both fragments on wix), beta and pos. calibrated, no fit");
h1_sum_charge_wixh = new TH1D(hist_name,hist_name,6000,0,600000);
h1_sum_charge_wixh->GetXaxis()->SetTitle("Cal. Energy TWIM");
h1_sum_charge_wixh->GetYaxis()->SetTitle("Counts");
h1_sum_charge_wixh->GetXaxis()->CenterTitle(true);
h1_sum_charge_wixh->GetYaxis()->CenterTitle(true);
h1_sum_charge_wixh->GetYaxis()->SetLabelSize(0.045);
h1_sum_charge_wixh->GetYaxis()->SetTitleSize(0.045);


//messel
TH1D* h1_one_charge_messel_sqrt;
sprintf(hist_name,"Energy deposited in TWIM on messel, beta and pos. calibrated, sqrt-fit");
h1_one_charge_messel_sqrt = new TH1D(hist_name,hist_name,2000,0,200);
h1_one_charge_messel_sqrt->GetXaxis()->SetTitle("Cal. Energy TWIM");
h1_one_charge_messel_sqrt->GetYaxis()->SetTitle("Counts");
h1_one_charge_messel_sqrt->GetXaxis()->CenterTitle(true);
h1_one_charge_messel_sqrt->GetYaxis()->CenterTitle(true);
h1_one_charge_messel_sqrt->GetYaxis()->SetLabelSize(0.045);
h1_one_charge_messel_sqrt->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_sum_charge_messel_sqrt;
sprintf(hist_name,"Energy deposited in TWIM on messel SUM (for both fragments on wix), beta and pos. calibrated, sqrt-fit");
h1_sum_charge_messel_sqrt = new TH1D(hist_name,hist_name,2000,0,200);
h1_sum_charge_messel_sqrt->GetXaxis()->SetTitle("Cal. Energy TWIM");
h1_sum_charge_messel_sqrt->GetYaxis()->SetTitle("Counts");
h1_sum_charge_messel_sqrt->GetXaxis()->CenterTitle(true);
h1_sum_charge_messel_sqrt->GetYaxis()->CenterTitle(true);
h1_sum_charge_messel_sqrt->GetYaxis()->SetLabelSize(0.045);
h1_sum_charge_messel_sqrt->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_sum_charge_messel;
sprintf(hist_name,"Energy deposited in TWIM on messel SUM (for both fragments on wix), beta and pos. calibrated, no fit");
h1_sum_charge_messel = new TH1D(hist_name,hist_name,6000,0,600000);
h1_sum_charge_messel->GetXaxis()->SetTitle("Cal. Energy TWIM");
h1_sum_charge_messel->GetYaxis()->SetTitle("Counts");
h1_sum_charge_messel->GetXaxis()->CenterTitle(true);
h1_sum_charge_messel->GetYaxis()->CenterTitle(true);
h1_sum_charge_messel->GetYaxis()->SetLabelSize(0.045);
h1_sum_charge_messel->GetYaxis()->SetTitleSize(0.045);


//now check that the sections are precalibrated... just check.... before beta and pos correction

TH1D* h1_energy_section0;
sprintf(hist_name,"Energy deposited in TWIM SECTION 0 before pos,beta cal");
h1_energy_section0 = new TH1D(hist_name,hist_name,6000,0,600000);
h1_energy_section0->GetXaxis()->SetTitle("Energy TWIM SECTION 0");
h1_energy_section0->GetYaxis()->SetTitle("Counts");
h1_energy_section0->GetXaxis()->CenterTitle(true);
h1_energy_section0->GetYaxis()->CenterTitle(true);
h1_energy_section0->GetYaxis()->SetLabelSize(0.045);
h1_energy_section0->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_energy_section1;
sprintf(hist_name,"Energy deposited in TWIM SECTION 1 before pos,beta cal");
h1_energy_section1 = new TH1D(hist_name,hist_name,6000,0,600000);
h1_energy_section1->GetXaxis()->SetTitle("Energy TWIM SECTION 1");
h1_energy_section1->GetYaxis()->SetTitle("Counts");
h1_energy_section1->GetXaxis()->CenterTitle(true);
h1_energy_section1->GetYaxis()->CenterTitle(true);
h1_energy_section1->GetYaxis()->SetLabelSize(0.045);
h1_energy_section1->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_energy_section2;
sprintf(hist_name,"Energy deposited in TWIM SECTION 2 before pos,beta cal");
h1_energy_section2 = new TH1D(hist_name,hist_name,6000,0,600000);
h1_energy_section2->GetXaxis()->SetTitle("Energy TWIM SECTION 2");
h1_energy_section2->GetYaxis()->SetTitle("Counts");
h1_energy_section2->GetXaxis()->CenterTitle(true);
h1_energy_section2->GetYaxis()->CenterTitle(true);
h1_energy_section2->GetYaxis()->SetLabelSize(0.045);
h1_energy_section2->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_energy_section3;
sprintf(hist_name,"Energy deposited in TWIM SECTION 3 before pos,beta cal");
h1_energy_section3 = new TH1D(hist_name,hist_name,6000,0,600000);
h1_energy_section3->GetXaxis()->SetTitle("Energy TWIM SECTION 3");
h1_energy_section3->GetYaxis()->SetTitle("Counts");
h1_energy_section3->GetXaxis()->CenterTitle(true);
h1_energy_section3->GetYaxis()->CenterTitle(true);
h1_energy_section3->GetYaxis()->SetLabelSize(0.045);
h1_energy_section3->GetYaxis()->SetTitleSize(0.045);

//deltax vs xal for section 0
TH2D* h2_deltax_xal_sec0[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name,"Section 0, Anode %i -> DeltaX vs Xal",i);
h2_deltax_xal_sec0[i] = new TH2D(hist_name,hist_name,1100,0,110,1000,-5,5);
h2_deltax_xal_sec0[i]->GetXaxis()->SetTitle("Xal [mm]");
h2_deltax_xal_sec0[i]->GetYaxis()->SetTitle("DeltaX (Xfit -Xal) [mm]");
h2_deltax_xal_sec0[i]->GetXaxis()->CenterTitle(true);
h2_deltax_xal_sec0[i]->GetYaxis()->CenterTitle(true);
h2_deltax_xal_sec0[i]->GetYaxis()->SetLabelSize(0.045);
h2_deltax_xal_sec0[i]->GetYaxis()->SetTitleSize(0.045);
}

//deltax vs xal for section 1
TH2D* h2_deltax_xal_sec1[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name,"Section 1, Anode %i -> DeltaX vs Xal",i);
h2_deltax_xal_sec1[i] = new TH2D(hist_name,hist_name,1100,0,110,1000,-5,5);
h2_deltax_xal_sec1[i]->GetXaxis()->SetTitle("Xal [mm]");
h2_deltax_xal_sec1[i]->GetYaxis()->SetTitle("DeltaX (Xfit -Xal) [mm]");
h2_deltax_xal_sec1[i]->GetXaxis()->CenterTitle(true);
h2_deltax_xal_sec1[i]->GetYaxis()->CenterTitle(true);
h2_deltax_xal_sec1[i]->GetYaxis()->SetLabelSize(0.045);
h2_deltax_xal_sec1[i]->GetYaxis()->SetTitleSize(0.045);
}

//deltax vs xal for section 2
TH2D* h2_deltax_xal_sec2[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name,"Section 2, Anode %i -> DeltaX vs Xal",i);
h2_deltax_xal_sec2[i] = new TH2D(hist_name,hist_name,1100,-110,0,1000,-5,5);
h2_deltax_xal_sec2[i]->GetXaxis()->SetTitle("Xal [mm]");
h2_deltax_xal_sec2[i]->GetYaxis()->SetTitle("DeltaX (Xfit -Xal) [mm]");
h2_deltax_xal_sec2[i]->GetXaxis()->CenterTitle(true);
h2_deltax_xal_sec2[i]->GetYaxis()->CenterTitle(true);
h2_deltax_xal_sec2[i]->GetYaxis()->SetLabelSize(0.045);
h2_deltax_xal_sec2[i]->GetYaxis()->SetTitleSize(0.045);
}

//deltax vs xal for section 3
TH2D* h2_deltax_xal_sec3[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name,"Section 3, Anode %i -> DeltaX vs Xal",i);
h2_deltax_xal_sec3[i] = new TH2D(hist_name,hist_name,1100,-110,0,1000,-5,5);
h2_deltax_xal_sec3[i]->GetXaxis()->SetTitle("Xal [mm]");
h2_deltax_xal_sec3[i]->GetYaxis()->SetTitle("DeltaX (Xfit -Xal) [mm]");
h2_deltax_xal_sec3[i]->GetXaxis()->CenterTitle(true);
h2_deltax_xal_sec3[i]->GetYaxis()->CenterTitle(true);
h2_deltax_xal_sec3[i]->GetYaxis()->SetLabelSize(0.045);
h2_deltax_xal_sec3[i]->GetYaxis()->SetTitleSize(0.045);
}

TH2D* h2_tofns_vs_ypos[29];
for (Int_t i = 0; i < 29; i++){
sprintf(hist_name,"TOFW pad %i tof_ns vs extrapolated y-pos.",i);
h2_tofns_vs_ypos[i] = new TH2D(hist_name,hist_name,500,-50,50,600,-300,300);
h2_tofns_vs_ypos[i]->GetXaxis()->SetTitle("tof ns");
h2_tofns_vs_ypos[i]->GetYaxis()->SetTitle("Extrap. y-pos [mm]");
h2_tofns_vs_ypos[i]->GetXaxis()->CenterTitle(true);
h2_tofns_vs_ypos[i]->GetYaxis()->CenterTitle(true);
h2_tofns_vs_ypos[i]->GetYaxis()->SetLabelSize(0.045);
h2_tofns_vs_ypos[i]->GetYaxis()->SetTitleSize(0.045);
}
TH2D* h2_tofns_vs_ypos2X2Y_LD_RU[29];
for (Int_t i = 0; i < 29; i++){
sprintf(hist_name,"TOFW pad %i tof_ns vs extrapolated y-pos. for subcase 2X2Y-LD-RU",i);
h2_tofns_vs_ypos2X2Y_LD_RU[i] = new TH2D(hist_name,hist_name,500,-50,50,600,-300,300);
h2_tofns_vs_ypos2X2Y_LD_RU[i]->GetXaxis()->SetTitle("tof ns");
h2_tofns_vs_ypos2X2Y_LD_RU[i]->GetYaxis()->SetTitle("Extrap. y-pos [mm]");
h2_tofns_vs_ypos2X2Y_LD_RU[i]->GetXaxis()->CenterTitle(true);
h2_tofns_vs_ypos2X2Y_LD_RU[i]->GetYaxis()->CenterTitle(true);
h2_tofns_vs_ypos2X2Y_LD_RU[i]->GetYaxis()->SetLabelSize(0.045);
h2_tofns_vs_ypos2X2Y_LD_RU[i]->GetYaxis()->SetTitleSize(0.045);
}
TH2D* h2_tofns_vs_ypos2X2Y_LU_RD[29];
for (Int_t i = 0; i < 29; i++){
sprintf(hist_name,"TOFW pad %i tof_ns vs extrapolated y-pos. for subcase 2X2Y-LU-RD",i);
h2_tofns_vs_ypos2X2Y_LU_RD[i] = new TH2D(hist_name,hist_name,500,-50,50,600,-300,300);
h2_tofns_vs_ypos2X2Y_LU_RD[i]->GetXaxis()->SetTitle("tof ns");
h2_tofns_vs_ypos2X2Y_LU_RD[i]->GetYaxis()->SetTitle("Extrap. y-pos [mm]");
h2_tofns_vs_ypos2X2Y_LU_RD[i]->GetXaxis()->CenterTitle(true);
h2_tofns_vs_ypos2X2Y_LU_RD[i]->GetYaxis()->CenterTitle(true);
h2_tofns_vs_ypos2X2Y_LU_RD[i]->GetYaxis()->SetLabelSize(0.045);
h2_tofns_vs_ypos2X2Y_LU_RD[i]->GetYaxis()->SetTitleSize(0.045);
}
TH2D* h2_tofns_vs_ypos2X2Y_L[29];
for (Int_t i = 0; i < 29; i++){
sprintf(hist_name,"TOFW pad %i tof_ns vs extrapolated y-pos. for subcase 2X2Y-L",i);
h2_tofns_vs_ypos2X2Y_L[i] = new TH2D(hist_name,hist_name,500,-50,50,600,-300,300);
h2_tofns_vs_ypos2X2Y_L[i]->GetXaxis()->SetTitle("tof ns");
h2_tofns_vs_ypos2X2Y_L[i]->GetYaxis()->SetTitle("Extrap. y-pos [mm]");
h2_tofns_vs_ypos2X2Y_L[i]->GetXaxis()->CenterTitle(true);
h2_tofns_vs_ypos2X2Y_L[i]->GetYaxis()->CenterTitle(true);
h2_tofns_vs_ypos2X2Y_L[i]->GetYaxis()->SetLabelSize(0.045);
h2_tofns_vs_ypos2X2Y_L[i]->GetYaxis()->SetTitleSize(0.045);
}
TH2D* h2_tofns_vs_ypos2X2Y_R[29];
for (Int_t i = 0; i < 29; i++){
sprintf(hist_name,"TOFW pad %i tof_ns vs extrapolated y-pos. for subcase 2X2Y-R",i);
h2_tofns_vs_ypos2X2Y_R[i] = new TH2D(hist_name,hist_name,500,-50,50,600,-300,300);
h2_tofns_vs_ypos2X2Y_R[i]->GetXaxis()->SetTitle("tof ns");
h2_tofns_vs_ypos2X2Y_R[i]->GetYaxis()->SetTitle("Extrap. y-pos [mm]");
h2_tofns_vs_ypos2X2Y_R[i]->GetXaxis()->CenterTitle(true);
h2_tofns_vs_ypos2X2Y_R[i]->GetYaxis()->CenterTitle(true);
h2_tofns_vs_ypos2X2Y_R[i]->GetYaxis()->SetLabelSize(0.045);
h2_tofns_vs_ypos2X2Y_R[i]->GetYaxis()->SetTitleSize(0.045);
}
//------------------------------------------------

//get the TSpines
TFile* file_tsplines = TFile::Open("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/data/output_retrieve_tspline_dt/retrieve_tspline_result.root");
TList* spline_list = (TList*)file_tsplines->Get("list_tspline");
TSpline3* spline_sec0[16];
TSpline3* spline_sec1[16];
TSpline3* spline_sec2[16];
TSpline3* spline_sec3[16];
char spline_names[500];
for (Int_t i = 0; i < 16; i++){
	sprintf(spline_names,"Section 0, Anode %i -> DeltaX vs Xal_pfx",i);
	spline_sec0[i] = (TSpline3*)spline_list->FindObject(spline_names);
	}
for (Int_t i = 0; i < 16; i++){
	sprintf(spline_names,"Section 1, Anode %i -> DeltaX vs Xal_pfx",i);
	spline_sec1[i] = (TSpline3*)spline_list->FindObject(spline_names);
	}
for (Int_t i = 0; i < 16; i++){
	sprintf(spline_names,"Section 2, Anode %i -> DeltaX vs Xal_pfx",i);
	spline_sec2[i] = (TSpline3*)spline_list->FindObject(spline_names);
	}
for (Int_t i = 0; i < 16; i++){
	sprintf(spline_names,"Section 3, Anode %i -> DeltaX vs Xal_pfx",i);
	spline_sec3[i] = (TSpline3*)spline_list->FindObject(spline_names);
	}


R3BEventHeader* DataCA = new R3BEventHeader();
TBranch* branchData = chain->GetBranch("EventHeader.");
branchData->SetAddress(&DataCA);

TClonesArray* SofTwimMappedData = new TClonesArray("R3BTwimMappedData",1);
R3BTwimHitData** softwimmappeddata;
TBranch *branchSofTwimMappedData = chain->GetBranch("TwimMappedData");
branchSofTwimMappedData->SetAddress(&SofTwimMappedData);

TClonesArray* SofMwpc1CalData = new TClonesArray("R3BMwpcCalData",5);
R3BMwpcCalData** sofmwpc1caldata;
TBranch *branchSofMwpc1CalData = chain->GetBranch("Mwpc1CalData");
branchSofMwpc1CalData->SetAddress(&SofMwpc1CalData);

TClonesArray* SofMwpc2CalData = new TClonesArray("R3BMwpcCalData",5);
R3BMwpcCalData** sofmwpc2caldata;
TBranch *branchSofMwpc2CalData = chain->GetBranch("Mwpc2CalData");
branchSofMwpc2CalData->SetAddress(&SofMwpc2CalData);

TClonesArray* SofMwpc3CalData = new TClonesArray("R3BMwpcCalData",5);
R3BMwpcCalData** sofmwpc3caldata;
TBranch *branchSofMwpc3CalData = chain->GetBranch("Mwpc3CalData");
branchSofMwpc3CalData->SetAddress(&SofMwpc3CalData);

TClonesArray* SofTofWSingleTcalData = new TClonesArray("R3BSofTofWSingleTcalData");
R3BSofTofWSingleTcalData** softofwsingletcaldata;
TBranch* branchSofTofWSingleTcalData = chain->GetBranch("SofTofWSingleTcalData");
branchSofTofWSingleTcalData->SetAddress(&SofTofWSingleTcalData);

TClonesArray* SofSciTcalData = new TClonesArray("R3BSofSciTcalData",2);
R3BSofSciTcalData** sofscitcaldata;
TBranch *branchSofSciTcalData = chain->GetBranch("SofSciTcalData");
branchSofSciTcalData->SetAddress(&SofSciTcalData);

//parameters for the raw anode energy
fstream fin;
fin.open("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/parameters/parameters_twim_anodes.csv", ios::in); 
string line, word;
vector<vector<vector<double> > > v_para_twim(4,vector<vector<double> >(16, vector<double>(2)));
getline(fin, line);
while(fin.peek()!=EOF) {
	getline(fin, line);
	stringstream s(line);
	vector<string> temp_vec;
	while (getline(s, word, ',')) {     temp_vec.push_back(word);
	}
	v_para_twim[stoi(temp_vec[0])][stoi(temp_vec[1])][0] = stod(temp_vec[2]);
	v_para_twim[stoi(temp_vec[0])][stoi(temp_vec[1])][1] = stod(temp_vec[3]);
	temp_vec.clear();
}
//parameters_for the sum_section energies
fstream fin2;
fin2.open("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/parameters/sum_anodes_parameters.csv", ios::in); 
string line2, word2;
vector<vector<double> > v_para_twim_sum;
getline(fin2, line2);
while(fin2.peek()!=EOF) {
	getline(fin2, line2);
	stringstream s(line2);
	vector<double> temp_vec;
	while (getline(s, word2, ',')) {     temp_vec.push_back(stod(word2));
	}
	v_para_twim_sum.push_back(temp_vec);
	
	temp_vec.clear();
}

//parameters for drift time raw
fstream fin3;
fin3.open("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/parameters/xal_vs_dt_params.csv", ios::in); 
string line3, word3;
vector<vector<vector<double> > > v_dt_twim(4,vector<vector<double> >(16, vector<double>(2)));
getline(fin3, line3);
while(fin3.peek()!=EOF) {
	getline(fin3, line3);
	stringstream s(line3);
	vector<string> temp_vec;
	while (getline(s, word3, ',')) {     temp_vec.push_back(word3);
	}
	v_dt_twim[stoi(temp_vec[0])][stoi(temp_vec[1])][0] = stod(temp_vec[2]);
	v_dt_twim[stoi(temp_vec[0])][stoi(temp_vec[1])][1] = stod(temp_vec[3]);
	temp_vec.clear();
}


//parameters_for TOFW calibration y position
fstream fin4;
fin4.open("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/macros/tof_params.csv", ios::in); 
string line4, word4;
vector<vector<double> > v_para_tof_y;
getline(fin4, line4);
while(fin4.peek()!=EOF) {
	getline(fin4, line4);
	stringstream s(line4);
	vector<double> temp_vec;
	while (getline(s, word4, ',')) {     temp_vec.push_back(stod(word4));
	}
	v_para_tof_y.push_back(temp_vec);
	
	temp_vec.clear();
}

//parameters for calibration of the position dependence of the energy side for the wixhausen and messel side
fstream fin5;
fin5.open("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/parameters/energy_position_calib_params.csv",ios::in);
string line5, word5;
vector<vector<double> > v_para_pos_energy;
getline(fin5,line5);
while(fin5.peek() != EOF){
	getline(fin5,line5);
	stringstream s(line5);
	vector<double> temp_vec;
	while (getline(s,word5,',')) { temp_vec.push_back(stod(word5));
	}
	v_para_pos_energy.push_back(temp_vec);
	temp_vec.clear();
}

const double slope_wix = v_para_pos_energy[0][0];
const double offset_wix = v_para_pos_energy[0][1];
const double slope_messel = v_para_pos_energy[1][0];
const double offset_messel = v_para_pos_energy[1][1];
const double wix_c0 = v_para_pos_energy[0][0];
const double wix_c1 = v_para_pos_energy[0][0];
const double wix_c2 = v_para_pos_energy[0][0];
const double wix_c3 = v_para_pos_energy[0][0];
const double wix_c4 = v_para_pos_energy[0][0];
const double wix_c5 = v_para_pos_energy[0][0];
const double wix_c6 = v_para_pos_energy[0][0];
const double wix_par0 = v_para_pos_energy[0][0];
const double wix_par1 = v_para_pos_energy[0][1];
const double wix_par2 = v_para_pos_energy[0][2];
const double messel_par0 = v_para_pos_energy[1][0];
const double messel_par1 = v_para_pos_energy[1][1];
const double messel_par2 = v_para_pos_energy[1][2];

const double messel_c0 = v_para_pos_energy[1][0];
const double messel_c1 = v_para_pos_energy[1][1];
const double messel_c2 = v_para_pos_energy[1][2];
const double messel_c3 = v_para_pos_energy[1][3];
const double messel_c4 = v_para_pos_energy[1][4];
const double messel_c5 = v_para_pos_energy[1][5];
const double messel_c6 = v_para_pos_energy[1][6];

const double sec0_c0 = v_para_pos_energy[0][0];
const double sec0_c1 = v_para_pos_energy[0][1];
const double sec0_c2 = v_para_pos_energy[0][2];

const double sec1_c0 = v_para_pos_energy[1][0];
const double sec1_c1 = v_para_pos_energy[1][1];
const double sec1_c2 = v_para_pos_energy[1][2];

const double sec2_c0 = v_para_pos_energy[2][0];
const double sec2_c1 = v_para_pos_energy[2][1];
const double sec2_c2 = v_para_pos_energy[2][2];

const double sec3_c0 = v_para_pos_energy[3][0];
const double sec3_c1 = v_para_pos_energy[3][1];
const double sec3_c2 = v_para_pos_energy[3][2];

for(Long64_t i=0;i < nevents;i++){
	Long64_t evtnr = i;
	if (i%100000==0)
		cout<<"Processing event for charge analysis "<<i<<endl;
	chain->GetEvent(i);
	entries_twim = SofTwimMappedData->GetEntries();
	entries_tof = SofTofWSingleTcalData->GetEntries();
	
	if (entries_twim > 1 && entries_tof == 2){
	//cout << "this is event with more than one TWIM entry:\t" << evtnr << endl;
		//cout << "corresponding entries:\t" << entries_twim << endl;
		Int_t nr_start;
    	nr_start = SofSciTcalData->GetEntries();
		R3BTwimMappedData** softwimmappeddata  = new R3BTwimMappedData*[entries_twim];
//		R3BMwpcCalData** sofmwpc1caldata = new R3BMwpcCalData*[entries_mw1];
//		R3BMwpcCalData** sofmwpc2caldata = new R3BMwpcCalData*[entries_mw2];
		R3BSofSciTcalData** sofscitcaldata = new R3BSofSciTcalData*[nr_start];

		//create 4 arrays with 0 as entry (= 4 sections)
		double sec_0_arr[16] = { 0.};
		double sec_1_arr[16] = { 0.};
	    double sec_2_arr[16] = { 0.};
	    double sec_3_arr[16] = { 0.};	
		//those are the arrays for the drift time
		double drift_sec_0_arr[16] = {0.};
		double drift_sec_1_arr[16] = {0.};
		double drift_sec_2_arr[16] = {0.};
		double drift_sec_3_arr[16] = {0.};
		
		double e_sum_sec_0 = 0.;
		double e_sum_sec_1 = 0.;
		double e_sum_sec_2 = 0.;
		double e_sum_sec_3 = 0.;
		
		bool single_filled_anodes = true;
		bool tref_anodes = true;
		UShort_t Tref_section0 = 0;
		UShort_t Tref_section1 = 0;
		UShort_t Tref_section2 = 0;
		UShort_t Tref_section3 = 0;
		
		//TJ just to check multiplicity of channel 16
		Int_t mult_channel_tref = 0;

		//fill the 4 arrays
		for (Int_t j = 0; j < entries_twim; j++){
			softwimmappeddata[j] = (R3BTwimMappedData*)SofTwimMappedData->At(j);
			Int_t twim_section = (softwimmappeddata[j]->GetSecID())-1;
			Int_t twim_anode = (softwimmappeddata[j]->GetAnodeID())-1;
			if (twim_section == 0 && twim_anode < 16){
				if (sec_0_arr[twim_anode] != 0){
					single_filled_anodes = false;
					continue;
					}
				sec_0_arr[twim_anode] = softwimmappeddata[j]->GetEnergy();
				drift_sec_0_arr[twim_anode] = softwimmappeddata[j]->GetTime();
				if (softwimmappeddata[j]->GetEnergy() != 0 && twim_anode != 13){  // section0, anode 13 is bad, leave out ...
					e_sum_sec_0 += (softwimmappeddata[j]->GetEnergy())*(1/v_para_twim[0][twim_anode][0]) - (v_para_twim[0][twim_anode][1]/v_para_twim[0][twim_anode][0]);
					}
				}
			else if (twim_section == 0 && twim_anode == 16){
					if (Tref_section0 != 0){
						tref_anodes = false;
						continue;
						}
					Tref_section0 = softwimmappeddata[j]->GetTime();
					mult_channel_tref++;
					}
			else if (twim_section == 1 && twim_anode < 16){
				if (sec_1_arr[twim_anode] != 0){
					single_filled_anodes = false;
					continue;
					}
				sec_1_arr[twim_anode] = softwimmappeddata[j]->GetEnergy();
				drift_sec_1_arr[twim_anode] = softwimmappeddata[j]->GetTime();
				if (softwimmappeddata[j]->GetEnergy() != 0){
					e_sum_sec_1 += (softwimmappeddata[j]->GetEnergy())*(1/v_para_twim[1][twim_anode][0]) - (v_para_twim[1][twim_anode][1]/v_para_twim[1][twim_anode][0]);
					}
				}
			else if (twim_section == 1 && twim_anode == 16){
					if (Tref_section1 != 0){
						tref_anodes = false;
						continue;
						}
				Tref_section1 = softwimmappeddata[j]->GetTime();
				
				}
			else if (twim_section == 2 && twim_anode < 16){
				if (sec_2_arr[twim_anode] != 0){
					single_filled_anodes = false;
					continue;
					}
				sec_2_arr[twim_anode] = softwimmappeddata[j]->GetEnergy();
				drift_sec_2_arr[twim_anode] = softwimmappeddata[j]->GetTime();
				if (softwimmappeddata[j]->GetEnergy() != 0){
					e_sum_sec_2 += (softwimmappeddata[j]->GetEnergy())*(1/v_para_twim[2][twim_anode][0]) - (v_para_twim[2][twim_anode][1]/v_para_twim[2][twim_anode][0]);
					}
				}
			else if (twim_section == 2 && twim_anode == 16){
					if (Tref_section2 != 0){
						tref_anodes = false;
						continue;
						}
				Tref_section2 = softwimmappeddata[j]->GetTime();
				}
			else if (twim_section == 3 && twim_anode < 16){
				if (sec_3_arr[twim_anode] != 0){
					single_filled_anodes = false;
					continue;
					}
				sec_3_arr[twim_anode] = softwimmappeddata[j]->GetEnergy();
				drift_sec_3_arr[twim_anode] = softwimmappeddata[j]->GetTime();
				if (softwimmappeddata[j]->GetEnergy() != 0){
					e_sum_sec_3 += (softwimmappeddata[j]->GetEnergy())*(1/v_para_twim[3][twim_anode][0]) - (v_para_twim[3][twim_anode][1]/v_para_twim[3][twim_anode][0]);
					}
				}		
			else if (twim_section == 3 && twim_anode == 16){
					if (Tref_section3 != 0){
						tref_anodes = false;
						continue;
						}
				Tref_section3 = softwimmappeddata[j]->GetTime();
				}
	
			}
		//check arrays if all entries are not equal 0.
		Double_t min_sec0 = sec_0_arr[0];
        Double_t min_sec1 = sec_1_arr[0];
        Double_t min_sec2 = sec_2_arr[0];
        Double_t min_sec3 = sec_3_arr[0];

        for (Int_t j = 1; j < 16; j++){
            if (sec_0_arr[j] < min_sec0){
                min_sec0 = sec_0_arr[j];
                }
            if(sec_1_arr[j] < min_sec1){
                min_sec1 = sec_1_arr[j];
                }
            if(sec_2_arr[j] < min_sec2){
                min_sec2 = sec_2_arr[j];
                }
            if(sec_3_arr[j] < min_sec3){
                min_sec3 = sec_3_arr[j];
                }
			}

		//correct energy section wise
		if (min_sec0 != 0){
			e_sum_sec_0 = e_sum_sec_0*v_para_twim_sum[0][0] +v_para_twim_sum[0][1]*16;
			}
		if (min_sec1 != 0){
			e_sum_sec_1 = e_sum_sec_1*v_para_twim_sum[1][0] +v_para_twim_sum[1][1]*16;
			}
		if (min_sec2 != 0){
            e_sum_sec_2 = e_sum_sec_2*v_para_twim_sum[2][0] +v_para_twim_sum[2][1]*16;
			}
		if (min_sec3 != 0){
            e_sum_sec_3 = e_sum_sec_3*v_para_twim_sum[3][0] +v_para_twim_sum[3][1]*16;
			}

		//FF1 is more to the left than FF2-> FF1.x > FF2.x
		//for calibration I only use events with Mw2_Xsum_nHits >1 && Mw2_Y_nHits >1
	Double_t FF1_slope;
	Double_t FF2_slope;
	Double_t FF1_offset;
	Double_t FF2_offset;
//	if (peaks_X_mw1_sum.size() > 1 && peaks_X_mw2_sum.size() > 1 && peaks_Y_mw1.size() > 1 && peaks_Y_mw2.size() > 1 && single_filled_anodes && tref_anodes && nr_start > 1){

		//---------------
		//subcase 2X2Y-D: one left down && one right down	
//		if (peaks_X_mw2_down.size() > 1 && peaks_Y_mw2[0] < -10 && peaks_Y_mw2[1] < -10 && e_sum_sec_0 > 0 && e_sum_sec_3 > 0 && peaks_Y_mw1[0] < -10 && peaks_Y_mw1[1] < -10 && e_sum_sec_1 == 0 && e_sum_sec_2 == 0   ){
//			//cout << "subcase 2X2Y-D: one left down && one right down" << endl;
//			double arr_xal_sec_0[16] = {-1000.};
//			double arr_xal_sec_3[16] = {-1000.};
//			for (Int_t j = 0; j < 16; j++){
//				Double_t anode_z_pos = pos_twim_entrance + j*twim_anodes_width;
//				arr_xal_sec_0[j] = v_dt_twim[0][j][0]*(drift_sec_0_arr[j]-Tref_section0)+v_dt_twim[0][j][1];
//				arr_xal_sec_3[j] = v_dt_twim[3][j][0]*(drift_sec_3_arr[j]-Tref_section3) +v_dt_twim[3][j][1];
//				}
//			Int_t n_fitpoints = 10;
//			Double_t x[10], y_sec0[10], y_sec3[10];
//			for ( Int_t j = 3; j < 13; j++){
//				x[j-3] = pos_twim_entrance + j*twim_anodes_width;
//				y_sec0[j-3] = arr_xal_sec_0[j];
//				y_sec3[j-3] = arr_xal_sec_3[j];
//				}
//			
//			TGraph* gr_sec0 = new TGraph(n_fitpoints,x,y_sec0);
//			TGraph* gr_sec3 = new TGraph(n_fitpoints,x,y_sec3);
//			TF1* fit_sec0 = new TF1("fit_sec0","[0]*x+ [1]");
//			TF1* fit_sec3 = new TF1("fit_sec3","[0]*x+ [1]");
//			gr_sec0->Fit(fit_sec0,"Q");
//			gr_sec3->Fit(fit_sec3,"Q");
//			// Get the fit parameters for section 0
//			Double_t slope_sec0 = fit_sec0->GetParameter(0);
//			Double_t offset_sec0 = fit_sec0->GetParameter(1);
//
//			// Get the fit parameters for section 3
//			Double_t slope_sec3 = fit_sec3->GetParameter(0);
//			Double_t offset_sec3 = fit_sec3->GetParameter(1);
//
//			//Get the deltax values
//			 double arr_deltax_sec_0[16] = {-1000.};
//			 double arr_deltax_sec_3[16] = {-1000.};
//			for ( Int_t j = 0; j < 16; j++){
//				arr_deltax_sec_0[j] = (slope_sec0*(pos_twim_entrance + j*twim_anodes_width) + offset_sec0) - arr_xal_sec_0[j];
//				arr_deltax_sec_3[j] = (slope_sec3*(pos_twim_entrance + j*twim_anodes_width) + offset_sec3) - arr_xal_sec_3[j];
//				h2_deltax_xal_sec0[j]->Fill(arr_xal_sec_0[j],arr_deltax_sec_0[j]);
//				h2_deltax_xal_sec2[j]->Fill(arr_xal_sec_3[j],arr_deltax_sec_3[j]);
//				}
//			delete  gr_sec0;
//			delete  gr_sec3;
//			delete fit_sec0;
//			delete fit_sec3;	
//			}
		//subcase 2X2Y-U: one left up && one right up
//		if (peaks_X_mw2_up.size() > 1 && peaks_Y_mw2[0] > 10 && peaks_Y_mw2[1] > 10 && e_sum_sec_1 > 0 && e_sum_sec_2 > 0 && peaks_Y_mw1[0] > 10 && peaks_Y_mw1[1] > 10 && e_sum_sec_0 == 0 && e_sum_sec_3 == 0){
//		//	cout << "subcase 2X2Y-U: one left up && one right up" << endl;
//			double arr_xal_sec_1[16] = {-1000.};
//			double arr_xal_sec_2[16] = {-1000.};
//			for (Int_t j = 0; j < 16; j++){
//				Double_t anode_z_pos = pos_twim_entrance + j*twim_anodes_width;
//				arr_xal_sec_1[j] = v_dt_twim[1][j][0]*(drift_sec_1_arr[j]-Tref_section1) + v_dt_twim[1][j][1];	
//				arr_xal_sec_2[j] = v_dt_twim[2][j][0]*(drift_sec_2_arr[j]-Tref_section2) +v_dt_twim[2][j][1];	
//				}
//			//just use the first 10 anode dt values to make the fit
//			Int_t n_fitpoints = 10;
//			Double_t x[10], y_sec1[10], y_sec2[10];
//			for ( Int_t j = 3; j < 13; j++){
//				x[j-3] = pos_twim_entrance + j*twim_anodes_width;
//				y_sec1[j-3] = arr_xal_sec_1[j];
//				y_sec2[j-3] = arr_xal_sec_2[j];
//				}
//			TGraph* gr_sec1 = new TGraph(n_fitpoints,x,y_sec1);
//			TGraph* gr_sec2 = new TGraph(n_fitpoints,x,y_sec2);
//			TF1* fit_sec1 = new TF1("fit_sec1","[0]*x+ [1]");
//			TF1* fit_sec2 = new TF1("fit_sec2","[0]*x+ [1]");
//			gr_sec1->Fit(fit_sec1,"Q");
//			gr_sec2->Fit(fit_sec2,"Q");
//			// Get the fit parameters for section 1
//			Double_t slope_sec1 = fit_sec1->GetParameter(0);
//			Double_t offset_sec1 = fit_sec1->GetParameter(1);
//
//			// Get the fit parameters for section 2
//			Double_t slope_sec2 = fit_sec2->GetParameter(0);
//			Double_t offset_sec2 = fit_sec2->GetParameter(1);
//
//			//Get the deltax values
//			 double arr_deltax_sec_1[16] = {-1000.};
//			 double arr_deltax_sec_2[16] = {-1000.};
//			for ( Int_t j = 0; j < 16; j++){
//				arr_deltax_sec_1[j] = (slope_sec1*(pos_twim_entrance + j*twim_anodes_width) + offset_sec1) - arr_xal_sec_1[j];
//				arr_deltax_sec_2[j] = (slope_sec2*(pos_twim_entrance + j*twim_anodes_width) + offset_sec2) - arr_xal_sec_2[j];
//				h2_deltax_xal_sec1[j]->Fill(arr_xal_sec_1[j],arr_deltax_sec_1[j]);
//				h2_deltax_xal_sec2[j]->Fill(arr_xal_sec_2[j],arr_deltax_sec_2[j]);
//				}
//			delete  gr_sec1;
//			delete  gr_sec2;
//			delete fit_sec1;
//			delete fit_sec2;	
//			}
		//subcase with one  FF down and one  FF up 
//		if (peaks_X_mw2_down.size() > 0 && peaks_X_mw2_up.size() > 0  && peaks_X_mw1_down.size() > 0 && peaks_X_mw1_up.size() > 0 && ((peaks_Y_mw2[0] <= 0 && peaks_Y_mw2[1] > 0) || (peaks_Y_mw2[0] > 0 && peaks_Y_mw2[1] <= 0))  && ((peaks_Y_mw1[0] <= 0 && peaks_Y_mw1[1] > 0) || (peaks_Y_mw1[0] > 0 && peaks_Y_mw1[1] <= 0))){
			//initialize start detector info
			vector<double> v_start_times;
			for (Int_t i = 0; i < 2; i++){
				sofscitcaldata[i] = (R3BSofSciTcalData*)SofSciTcalData->At(i);
				v_start_times.push_back(sofscitcaldata[i]->GetRawTimeNs());
				}
			Double_t start_time = 0.5*(v_start_times[0]-v_start_times[1]); //TODO: substract time rom start to target
			
			//subsubcase 2X2Y-LD-RU: one left down && one right up
			if (e_sum_sec_0 > 0 && e_sum_sec_2 > 0  && min_sec0 != 0 && min_sec2 != 0&& e_sum_sec_1 == 0 && e_sum_sec_3 == 0 && (abs(sec_0_arr[15]- sec_0_arr[0])) < 3000 && (abs(sec_2_arr[15]- sec_2_arr[0])) < 3000){
			//cout << "subsubcase 2X2Y-LD-RU: one left down && one right up" << endl;
			double arr_xal_sec_0[16] = {-1000.};
			double arr_xal_sec_2[16] = {-1000.};
			for (Int_t j = 0; j < 16; j++){
				Double_t anode_z_pos = pos_twim_entrance + j*twim_anodes_width;
				arr_xal_sec_0[j] = v_dt_twim[0][j][0]*(drift_sec_0_arr[j]-Tref_section0)+v_dt_twim[0][j][1];	
				arr_xal_sec_2[j] = v_dt_twim[2][j][0]*(drift_sec_2_arr[j]-Tref_section2)+v_dt_twim[2][j][1];	
				}
			//just use the first 10 anode dt values to make the fit
			Int_t n_fitpoints = 10;
			Double_t x[10], y_sec0[10], y_sec2[10];
			for ( Int_t j = 3; j < 13; j++){
				x[j-3] = pos_twim_entrance + j*twim_anodes_width;
				y_sec0[j-3] = arr_xal_sec_0[j] + (spline_sec0[j]->Eval(arr_xal_sec_0[j]));
				y_sec2[j-3] = arr_xal_sec_2[j] + (spline_sec2[j]->Eval(arr_xal_sec_2[j]));
				}
			TGraph* gr_sec0 = new TGraph(n_fitpoints,x,y_sec0);
			TGraph* gr_sec2 = new TGraph(n_fitpoints,x,y_sec2);
			TF1* fit_sec0 = new TF1("fit_sec0","[0]*x+ [1]");
			TF1* fit_sec2 = new TF1("fit_sec2","[0]*x+ [1]");
			gr_sec0->Fit(fit_sec0,"Q");
			gr_sec2->Fit(fit_sec2,"Q");
			// Get the fit parameters for section 0
			Double_t slope_sec0 = fit_sec0->GetParameter(0);
			Double_t offset_sec0 = fit_sec0->GetParameter(1);

			// Get the fit parameters for section 2
			Double_t slope_sec2 = fit_sec2->GetParameter(0);
			Double_t offset_sec2 = fit_sec2->GetParameter(1);
			//Calculate psi_in (from slope)
			Double_t psi_in_sec0 = atan(slope_sec0);
			Double_t psi_in_sec2 = atan(slope_sec2);

			Double_t mw1_x_sec0 = slope_sec0*4.5 + offset_sec0;
			Double_t mw1_x_sec2 = slope_sec2*4.5 + offset_sec2;
			//cout << "eventnr:\t" << evtnr << endl;
			vector<vector<double> > tof_and_y = tof_and_y_calibrated(SofMwpc3CalData,SofTofWSingleTcalData); 
			//cout << "the two times:\t" << tof_and_y[0][0] << "   " << tof_and_y[1][0] << endl;
			if (tof_and_y != dummy_vec){
			vector<double> pathlength_sec0 = path_tof_calc(psi_in_sec0,mw1_x_sec0,tof_and_y[0][2],tof_and_y[0][0],start_time);
			vector<double> pathlength_sec2 = path_tof_calc(psi_in_sec2,mw1_x_sec2,tof_and_y[1][2],tof_and_y[1][0],start_time);
			Double_t full_path_sec0 = sqrt(pow(pathlength_sec0[0],2)+pow(0.1*(tof_and_y[0][1]),2));
			Double_t full_path_sec2 = sqrt(pow(pathlength_sec2[0],2)+pow(0.1*(tof_and_y[1][1]),2));
			Double_t time_of_flight_sec0 = pathlength_sec0[1];
			Double_t time_of_flight_sec2 = pathlength_sec2[1];
			Double_t beta0 = ((full_path_sec0/time_of_flight_sec0)*pow(10,7))/299792458;
			Double_t beta2 = ((full_path_sec2/time_of_flight_sec2)*pow(10,7))/299792458;
			h2_beta_vs_energy->Fill(beta0,e_sum_sec_0);
			h2_beta_vs_energy_sec0->Fill(beta0,e_sum_sec_0);
			h2_beta_vs_energy_corr->Fill(beta0,e_sum_sec_0*mean_ene/(beta_func_sec0->Eval(beta0)));
			h2_beta_vs_energy_corr_sec0->Fill(beta0,e_sum_sec_0*mean_ene/(beta_func_sec0->Eval(beta0)));
			h2_beta_vs_energy->Fill(beta2,e_sum_sec_2);
			h2_beta_vs_energy_sec2->Fill(beta2,e_sum_sec_2);
			h2_beta_vs_energy_corr->Fill(beta2,e_sum_sec_2*mean_ene/(beta_func_sec2->Eval(beta2)));
			h2_beta_vs_energy_corr_sec2->Fill(beta2,e_sum_sec_2*mean_ene/(beta_func_sec2->Eval(beta2)));
			h2_x_music_vs_energy_messel->Fill(y_sec0[7],e_sum_sec_0*mean_ene/(beta_func_sec0->Eval(beta0)));
			h2_x_music_vs_energy_messel_sec0->Fill(y_sec0[7],e_sum_sec_0*mean_ene/(beta_func_sec0->Eval(beta0)));
			h2_x_music_vs_energy_wix->Fill(y_sec2[7],e_sum_sec_2*mean_ene/(beta_func_sec2->Eval(beta2)));
			h2_x_music_vs_energy_wix_sec2->Fill(y_sec2[7],e_sum_sec_2*mean_ene/(beta_func_sec2->Eval(beta2)));
			//now comes the final energy correction with positioning tspline
			Double_t beta_corr_energy_sec0 = e_sum_sec_0*mean_ene/(beta_func_sec0->Eval(beta0));
			Double_t beta_corr_energy_sec2 = e_sum_sec_2*mean_ene/(beta_func_sec2->Eval(beta2));
			Double_t pos_corr_energy_sec0 = beta_corr_energy_sec0*mean_pos_sec0/(spline_messel_sec0->Eval(y_sec0[7]));
			Double_t pos_corr_energy_sec2 = beta_corr_energy_sec2*mean_pos_sec2/(spline_wix_sec2->Eval(y_sec2[7]));
			h2_x_music_vs_energy_wix_corr->Fill(y_sec2[7],pos_corr_energy_sec2);
			h2_x_music_vs_energy_messel_corr->Fill(y_sec0[7],pos_corr_energy_sec0);
			h1_charge_music_corr_messel->Fill(pos_corr_energy_sec0);
			h1_charge_music_corr_wix->Fill(pos_corr_energy_sec2);
			//best possible calibration....
			h1_charge_music_corr_sec0->Fill(pos_corr_energy_sec0);
			h1_charge_music_corr_sec2->Fill(pos_corr_energy_sec2);
			h1_charge_music_corr_sum->Fill(pos_corr_energy_sec0+pos_corr_energy_sec2);
			h1_charge_music_corr_messel_full_cal->Fill(pos_corr_energy_sec0*slope_messel + offset_messel);
			h1_charge_music_corr_wix_full_cal->Fill(pos_corr_energy_sec2*slope_wix + offset_wix);
			h1_charge_music_corr_sum_full_cal->Fill(pos_corr_energy_sec0*slope_messel + offset_messel+pos_corr_energy_sec2*slope_wix + offset_wix);
			h1_charge_music_corr_wix_full_cal_pol->Fill(messel_c0 +messel_c1*pos_corr_energy_sec2+messel_c2*pow(pos_corr_energy_sec2,2)+messel_c3*pow(pos_corr_energy_sec2,3)+messel_c4*pow(pos_corr_energy_sec2,4)+messel_c5*pow(pos_corr_energy_sec2,5)+messel_c6*pow(pos_corr_energy_sec2,6));
			h1_charge_music_corr_messel_full_cal_pol->Fill(pos_corr_energy_sec0);
			h1_charge_music_corr_sum_full_cal_pol->Fill(messel_c0 +messel_c1*pos_corr_energy_sec2+messel_c2*pow(pos_corr_energy_sec2,2)+messel_c3*pow(pos_corr_energy_sec2,3)+messel_c4*pow(pos_corr_energy_sec2,4)+messel_c5*pow(pos_corr_energy_sec2,5)+messel_c6*pow(pos_corr_energy_sec2,6)+pos_corr_energy_sec0);
			//plots only with beta correction:
			h1_charge_music_corr_messel_only_beta->Fill(beta_corr_energy_sec0);
			h1_charge_music_corr_wix_only_beta->Fill(beta_corr_energy_sec2);
			h1_charge_music_corr_sum_only_beta->Fill(beta_corr_energy_sec0+beta_corr_energy_sec2);
			//plot with sqrt fit
			h1_one_charge_wixh_sqrt->Fill(wix_par0+ wix_par1*sqrt(pos_corr_energy_sec2)+wix_par2*pos_corr_energy_sec2);
			h1_one_charge_messel_sqrt->Fill(messel_par0+ messel_par1*sqrt(pos_corr_energy_sec0)+ messel_par2*pos_corr_energy_sec0);
 			//check precalibrated state
 			h1_energy_section0->Fill(e_sum_sec_0);
 			h1_energy_section2->Fill(e_sum_sec_2);
			
			//if (pos_corr_energy_sec0 < 10 || pos_corr_energy_sec2 < 10){
			//	cout << "this is eventnr:\t" << evtnr << endl; 
			//	cout << "e_sum_sec_0:\t" << e_sum_sec_0 << endl;
			//	cout << "e_sum_sec_2:\t" << e_sum_sec_2 << endl;
			//	cout << "beta_corr_energy_sec0:\t" << beta_corr_energy_sec0 << endl;
			//	cout << "beta_corr_energy_sec2:\t" << beta_corr_energy_sec2 << endl;
			//	cout << "pos_corr_energy_sec0:\t" << pos_corr_energy_sec0 << endl;
			//	cout << "pos_corr_energy_sec2:\t" << pos_corr_energy_sec2 << endl;
			//	cout << "spine eval wix:\t" << spline_wix->Eval(y_sec2[7]) << endl;
			//	cout << "spline eval messel:\t" << spline_messel->Eval(y_sec0[7]) << endl;
			//	}
			
				}
			delete  gr_sec0;
			delete  gr_sec2;
			delete fit_sec0;
			delete fit_sec2;	
			}
			//subsubcase 2X2Y-LU-RD: one left up && one right down
			if (e_sum_sec_1 > 0 && e_sum_sec_3 > 0  && min_sec1 != 0 && min_sec3 != 0 && e_sum_sec_0 == 0 && e_sum_sec_2 == 0 && (abs(sec_1_arr[15]- sec_1_arr[0])) < 3000 && (abs(sec_3_arr[15]- sec_3_arr[0])) < 3000){
			//	cout << "subsubcase 2X2Y-LU-RD: one left up && one right down" << endl;
			double arr_xal_sec_1[16] = {-1000.};
			double arr_xal_sec_3[16] = {-1000.};
			for (Int_t j = 0; j < 16; j++){
				Double_t anode_z_pos = pos_twim_entrance + j*twim_anodes_width;
				arr_xal_sec_1[j] = v_dt_twim[1][j][0]*(drift_sec_1_arr[j]-Tref_section1) + v_dt_twim[1][j][1];	
				arr_xal_sec_3[j] = v_dt_twim[3][j][0]*(drift_sec_3_arr[j]-Tref_section3)+v_dt_twim[3][j][1];	
				}
			//just use the first 10 anode dt values to make the fit
			Int_t n_fitpoints = 10;
			Double_t x[10], y_sec1[10], y_sec3[10];
			for ( Int_t j = 3; j < 13; j++){
				x[j-3] = pos_twim_entrance + j*twim_anodes_width;
				y_sec1[j-3] = arr_xal_sec_1[j] + (spline_sec1[j]->Eval(arr_xal_sec_1[j]));
				y_sec3[j-3] = arr_xal_sec_3[j] + (spline_sec3[j]->Eval(arr_xal_sec_3[j]));
				}
			TGraph* gr_sec1 = new TGraph(n_fitpoints,x,y_sec1);
			TGraph* gr_sec3 = new TGraph(n_fitpoints,x,y_sec3);
			TF1* fit_sec1 = new TF1("fit_sec1","[0]*x+ [1]");
			TF1* fit_sec3 = new TF1("fit_sec3","[0]*x+ [1]");
			gr_sec1->Fit(fit_sec1,"Q");
			gr_sec3->Fit(fit_sec3,"Q");
			// Get the fit parameters for section 1
			Double_t slope_sec1 = fit_sec1->GetParameter(0);
			Double_t offset_sec1 = fit_sec1->GetParameter(1);

			// Get the fit parameters for section 2
			Double_t slope_sec3 = fit_sec3->GetParameter(0);
			Double_t offset_sec3 = fit_sec3->GetParameter(1);
			//Calculate psi_in (from slope)
			Double_t psi_in_sec1 = atan(slope_sec1);
			Double_t psi_in_sec3 = atan(slope_sec3);

			Double_t mw1_x_sec1 = slope_sec1*4.5 + offset_sec1;
			Double_t mw1_x_sec3 = slope_sec3*4.5 + offset_sec3;
			//cout << "eventnr:\t" << evtnr << endl;
			vector<vector<double> > tof_and_y = tof_and_y_calibrated(SofMwpc3CalData,SofTofWSingleTcalData);
            if (tof_and_y != dummy_vec){
            vector<double> pathlength_sec1 = path_tof_calc(psi_in_sec1,mw1_x_sec1,tof_and_y[1][2],tof_and_y[1][0],start_time);
            vector<double> pathlength_sec3 = path_tof_calc(psi_in_sec3,mw1_x_sec3,tof_and_y[0][2],tof_and_y[0][0],start_time);
            Double_t full_path_sec1 = sqrt(pow(pathlength_sec1[0],2)+pow(0.1*(tof_and_y[1][1]),2));
            Double_t full_path_sec3 = sqrt(pow(pathlength_sec3[0],2)+pow(0.1*(tof_and_y[0][1]),2));
            Double_t time_of_flight_sec1 = pathlength_sec1[1];
            Double_t time_of_flight_sec3 = pathlength_sec3[1];
			//cout << "time of flight\t" << time_of_flight_sec1<< endl;
			//cout << "pathlength\t" << full_path_sec1 << endl;
			Double_t beta3 = ((full_path_sec3/time_of_flight_sec3)*pow(10,7))/299792458;
			Double_t beta1 = ((full_path_sec1/time_of_flight_sec1)*pow(10,7))/299792458;	
			h2_beta_vs_energy->Fill(beta3,e_sum_sec_3);
			h2_beta_vs_energy_sec3->Fill(beta3,e_sum_sec_3);
			h2_beta_vs_energy_corr->Fill(beta3,e_sum_sec_3*mean_ene/(beta_func_sec3->Eval(beta3)));
			h2_beta_vs_energy_corr_sec3->Fill(beta3,e_sum_sec_3*mean_ene/(beta_func_sec3->Eval(beta3)));
			h2_beta_vs_energy->Fill(beta1,e_sum_sec_1);
			h2_beta_vs_energy_sec1->Fill(beta1,e_sum_sec_1);
			h2_beta_vs_energy_corr->Fill(beta1,e_sum_sec_1*mean_ene/(beta_func_sec1->Eval(beta1)));
			h2_beta_vs_energy_corr_sec1->Fill(beta1,e_sum_sec_1*mean_ene/(beta_func_sec1->Eval(beta1)));
			h2_x_music_vs_energy_messel->Fill(y_sec1[7],e_sum_sec_1*mean_ene/(beta_func_sec1->Eval(beta1)));
			h2_x_music_vs_energy_messel_sec1->Fill(y_sec1[7],e_sum_sec_1*mean_ene/(beta_func_sec1->Eval(beta1)));
			h2_x_music_vs_energy_wix->Fill(y_sec3[7],e_sum_sec_3*mean_ene/(beta_func_sec3->Eval(beta3)));
			h2_x_music_vs_energy_wix_sec3->Fill(y_sec3[7],e_sum_sec_3*mean_ene/(beta_func_sec3->Eval(beta3)));
			//now comes the final energy correction with positioning tspline
			Double_t beta_corr_energy_sec3 = e_sum_sec_3*mean_ene/(beta_func_sec3->Eval(beta3));
			Double_t beta_corr_energy_sec1 = e_sum_sec_1*mean_ene/(beta_func_sec1->Eval(beta1));
			Double_t pos_corr_energy_sec3 = beta_corr_energy_sec3*mean_pos_sec3/(spline_wix_sec3->Eval(y_sec3[7]));
			Double_t pos_corr_energy_sec1 = beta_corr_energy_sec1*mean_pos_sec1/(spline_messel_sec1->Eval(y_sec1[7]));
			h2_x_music_vs_energy_messel_corr->Fill(y_sec1[7],pos_corr_energy_sec1);
			h2_x_music_vs_energy_wix_corr->Fill(y_sec3[7],pos_corr_energy_sec3);
			h1_charge_music_corr_messel->Fill(pos_corr_energy_sec1);
			h1_charge_music_corr_wix->Fill(pos_corr_energy_sec3);
			//best possible calibration ....
			h1_charge_music_corr_sec1->Fill(pos_corr_energy_sec1);
			h1_charge_music_corr_sec3->Fill(pos_corr_energy_sec3);
			h1_charge_music_corr_sum->Fill(pos_corr_energy_sec1+pos_corr_energy_sec3);
			h1_charge_music_corr_messel_full_cal->Fill(pos_corr_energy_sec1*slope_messel + offset_messel);
			h1_charge_music_corr_wix_full_cal->Fill(pos_corr_energy_sec3*slope_wix + offset_wix);
			h1_charge_music_corr_sum_full_cal->Fill(pos_corr_energy_sec1*slope_messel + offset_messel+pos_corr_energy_sec3*slope_wix + offset_wix);
			h1_charge_music_corr_wix_full_cal_pol->Fill(messel_c0 +messel_c1*pos_corr_energy_sec3+messel_c2*pow(pos_corr_energy_sec3,2)+messel_c3*pow(pos_corr_energy_sec3,3)+messel_c4*pow(pos_corr_energy_sec3,4)+messel_c5*pow(pos_corr_energy_sec3,5)+messel_c6*pow(pos_corr_energy_sec3,6));
			h1_charge_music_corr_messel_full_cal_pol->Fill(pos_corr_energy_sec1);
			h1_charge_music_corr_sum_full_cal_pol->Fill(messel_c0 +messel_c1*pos_corr_energy_sec3+messel_c2*pow(pos_corr_energy_sec3,2)+messel_c3*pow(pos_corr_energy_sec3,3)+messel_c4*pow(pos_corr_energy_sec3,4)+messel_c5*pow(pos_corr_energy_sec3,5)+messel_c6*pow(pos_corr_energy_sec3,6)+pos_corr_energy_sec1);
			//plots only with beta correction:
			h1_charge_music_corr_messel_only_beta->Fill(beta_corr_energy_sec1);
			h1_charge_music_corr_wix_only_beta->Fill(beta_corr_energy_sec3);
			h1_charge_music_corr_sum_only_beta->Fill(beta_corr_energy_sec1+beta_corr_energy_sec3);	
			//plot with sqrt fit
			h1_one_charge_wixh_sqrt->Fill(wix_par0+ wix_par1*sqrt(pos_corr_energy_sec3)+wix_par2*pos_corr_energy_sec3);
			h1_one_charge_messel_sqrt->Fill(messel_par0+ messel_par1*sqrt(pos_corr_energy_sec1)+ messel_par2*pos_corr_energy_sec1);
 			//check precalibrated state
 			h1_energy_section1->Fill(e_sum_sec_1);
 			h1_energy_section3->Fill(e_sum_sec_3);
			//if (pos_corr_energy_sec1 < 10 || pos_corr_energy_sec3 < 10){
			//	cout << "this is eventnr:\t" << evtnr << endl; 
			//	cout << "e_sum_sec_1:\t" << e_sum_sec_1 << endl;
			//	cout << "e_sum_sec_3:\t" << e_sum_sec_3 << endl;
			//	cout << "beta_corr_energy_sec1:\t" << beta_corr_energy_sec1 << endl;
			//	cout << "beta_corr_energy_sec3:\t" << beta_corr_energy_sec3 << endl;
			//	cout << "pos_corr_energy_sec1:\t" << pos_corr_energy_sec1 << endl;
			//	cout << "pos_corr_energy_sec3:\t" << pos_corr_energy_sec3 << endl;
			//	cout << "spine eval wix:\t" << spline_wix->Eval(y_sec3[7]) << endl;
			//	cout << "spline eval messel:\t" << spline_messel->Eval(y_sec1[7]) << endl;
			//	}
			
			
                }
			delete  gr_sec1;
			delete  gr_sec3;
			delete fit_sec1;
			delete fit_sec3;	
			}
			//subsubcase 2X2Y-L: one left down && one left up
			if (e_sum_sec_0 > 0 && e_sum_sec_1 > 0 && min_sec0 != 0 && min_sec1 != 0  && e_sum_sec_2 == 0 && e_sum_sec_3 == 0 && (abs(sec_0_arr[15]- sec_0_arr[0])) < 3000 && (abs(sec_1_arr[15]- sec_1_arr[0])) < 3000){
			//		cout << "subsubcase 2X2Y-L: one left down && one left up" << endl;
			double arr_xal_sec_0[16] = {-1000.};
			double arr_xal_sec_1[16] = {-1000.};
			for (Int_t j = 0; j < 16; j++){
				Double_t anode_z_pos = pos_twim_entrance + j*twim_anodes_width;
				arr_xal_sec_0[j] = v_dt_twim[0][j][0]*(drift_sec_0_arr[j]-Tref_section0)+v_dt_twim[0][j][1];	
				arr_xal_sec_1[j] = v_dt_twim[1][j][0]*(drift_sec_1_arr[j]-Tref_section1)+v_dt_twim[1][j][1];	
				}
			//just use the first 10 anode dt values to make the fit
			Int_t n_fitpoints = 10;
			Double_t x[10], y_sec0[10], y_sec1[10];
			for ( Int_t j = 3; j < 13; j++){
				x[j-3] = pos_twim_entrance + j*twim_anodes_width;
				y_sec0[j-3] = arr_xal_sec_0[j] + (spline_sec0[j]->Eval(arr_xal_sec_0[j]));
				y_sec1[j-3] = arr_xal_sec_1[j] + (spline_sec1[j]->Eval(arr_xal_sec_1[j]));
				}
			TGraph* gr_sec0 = new TGraph(n_fitpoints,x,y_sec0);
			TGraph* gr_sec1 = new TGraph(n_fitpoints,x,y_sec1);
			TF1* fit_sec0 = new TF1("fit_sec0","[0]*x+ [1]");
			TF1* fit_sec1 = new TF1("fit_sec1","[0]*x+ [1]");
			gr_sec0->Fit(fit_sec0,"Q");
			gr_sec1->Fit(fit_sec1,"Q");
			// Get the fit parameters for section 0
			Double_t slope_sec0 = fit_sec0->GetParameter(0);
			Double_t offset_sec0 = fit_sec0->GetParameter(1);

			// Get the fit parameters for section 1
			Double_t slope_sec1 = fit_sec1->GetParameter(0);
			Double_t offset_sec1 = fit_sec1->GetParameter(1);
			
			//Calculate psi_in (from slope)
			Double_t psi_in_sec0 = atan(slope_sec0);
			Double_t psi_in_sec1 = atan(slope_sec1);

			Double_t mw1_x_sec0 = slope_sec0*4.5 + offset_sec0;
			Double_t mw1_x_sec1 = slope_sec1*4.5 + offset_sec1;
			//cout << "eventnr:\t" << evtnr << endl;
			vector<vector<double> > tof_and_y = tof_and_y_calibrated(SofMwpc3CalData,SofTofWSingleTcalData);
            if (tof_and_y != dummy_vec){
            vector<double> pathlength_sec0 = path_tof_calc(psi_in_sec0,mw1_x_sec0,tof_and_y[0][2],tof_and_y[0][0],start_time);
            vector<double> pathlength_sec1 = path_tof_calc(psi_in_sec1,mw1_x_sec1,tof_and_y[1][2],tof_and_y[1][0],start_time);
            Double_t full_path_sec0 = sqrt(pow(pathlength_sec0[0],2)+pow(0.1*(tof_and_y[0][1]),2));
            Double_t full_path_sec1 = sqrt(pow(pathlength_sec1[0],2)+pow(0.1*(tof_and_y[1][1]),2));
            Double_t time_of_flight_sec0 = pathlength_sec0[1];
            Double_t time_of_flight_sec1 = pathlength_sec1[1];
	    		Double_t beta0 = ((full_path_sec0/time_of_flight_sec0)*pow(10,7))/299792458;
			Double_t beta1 = ((full_path_sec1/time_of_flight_sec1)*pow(10,7))/299792458;
			h2_beta_vs_energy->Fill(beta0,e_sum_sec_0);
			h2_beta_vs_energy_sec0->Fill(beta0,e_sum_sec_0);
			h2_beta_vs_energy_corr->Fill(beta0,e_sum_sec_0*mean_ene/(beta_func_sec0->Eval(beta0)));
			h2_beta_vs_energy_corr_sec0->Fill(beta0,e_sum_sec_0*mean_ene/(beta_func_sec0->Eval(beta0)));
			h2_beta_vs_energy->Fill(beta1,e_sum_sec_1);
			h2_beta_vs_energy_sec1->Fill(beta1,e_sum_sec_1);
			h2_beta_vs_energy_corr->Fill(beta1,e_sum_sec_1*mean_ene/(beta_func_sec1->Eval(beta1)));
			h2_beta_vs_energy_corr_sec1->Fill(beta1,e_sum_sec_1*mean_ene/(beta_func_sec1->Eval(beta1)));
			h2_x_music_vs_energy_messel->Fill(y_sec0[7],e_sum_sec_0*mean_ene/(beta_func_sec0->Eval(beta0)));
			h2_x_music_vs_energy_messel->Fill(y_sec1[7],e_sum_sec_1*mean_ene/(beta_func_sec1->Eval(beta1)));
			h2_x_music_vs_energy_messel_sec0->Fill(y_sec0[7],e_sum_sec_0*mean_ene/(beta_func_sec0->Eval(beta0)));
			h2_x_music_vs_energy_messel_sec1->Fill(y_sec1[7],e_sum_sec_1*mean_ene/(beta_func_sec1->Eval(beta1)));
			//now comes the final energy correction with positioning tspline
			Double_t beta_corr_energy_sec0 = e_sum_sec_0*mean_ene/(beta_func_sec0->Eval(beta0));
			Double_t beta_corr_energy_sec1 = e_sum_sec_1*mean_ene/(beta_func_sec1->Eval(beta1));
			Double_t pos_corr_energy_sec0 = beta_corr_energy_sec0*mean_pos_sec0/(spline_messel_sec0->Eval(y_sec0[7]));
			Double_t pos_corr_energy_sec1 = beta_corr_energy_sec1*mean_pos_sec1/(spline_messel_sec1->Eval(y_sec1[7]));
			Double_t charge_sec0 = sec0_c0 + sec0_c1*sqrt(pos_corr_energy_sec0) +sec0_c2*pos_corr_energy_sec0; 
			Double_t charge_sec1 = sec1_c0 + sec1_c1*sqrt(pos_corr_energy_sec1) +sec1_c2*pos_corr_energy_sec1;
			h2_x_energy_vs_energy_messel_corr->Fill(charge_sec0,charge_sec1);
			h2_charge_sum_vs_diff_messel_corr->Fill(charge_sec1-charge_sec0,charge_sec1+charge_sec0);
			h2_x_music_vs_energy_messel_corr->Fill(y_sec0[7],pos_corr_energy_sec0);
			h2_x_music_vs_energy_messel_corr->Fill(y_sec1[7],pos_corr_energy_sec1);
			h1_charge_music_corr_messel->Fill(pos_corr_energy_sec0);
			h1_charge_music_corr_messel->Fill(pos_corr_energy_sec1);
			//best possible calibration....
			h1_charge_music_corr_sec0->Fill(pos_corr_energy_sec0);
			h1_charge_music_corr_sec1->Fill(pos_corr_energy_sec1);
			h1_charge_music_corr_sum->Fill(pos_corr_energy_sec0+pos_corr_energy_sec1);
			h1_charge_music_corr_messel_full_cal->Fill(pos_corr_energy_sec0*slope_messel + offset_messel);
			h1_charge_music_corr_messel_full_cal->Fill(pos_corr_energy_sec1*slope_messel + offset_messel);
			h1_charge_music_corr_sum_full_cal->Fill(pos_corr_energy_sec0*slope_messel + offset_messel+pos_corr_energy_sec1*slope_messel + offset_messel);
			h1_charge_music_corr_messel_full_cal_pol->Fill(pos_corr_energy_sec0);
			h1_charge_music_corr_messel_full_cal_pol->Fill(pos_corr_energy_sec1);
			h1_charge_music_corr_sum_full_cal_pol->Fill(pos_corr_energy_sec0+pos_corr_energy_sec1);
			h1_charge_music_corr_sum_only_messel->Fill(pos_corr_energy_sec0+pos_corr_energy_sec1);
			//plots only with beta correction:
			h1_charge_music_corr_messel_only_beta->Fill(beta_corr_energy_sec0);
			h1_charge_music_corr_messel_only_beta->Fill(beta_corr_energy_sec1);
			h1_charge_music_corr_sum_only_beta->Fill(beta_corr_energy_sec0+beta_corr_energy_sec1);
			//plot with sqrt fit
			h1_one_charge_messel_sqrt->Fill(messel_par0+ messel_par1*sqrt(pos_corr_energy_sec0)+ messel_par2*pos_corr_energy_sec0);
			h1_one_charge_messel_sqrt->Fill(messel_par0+ messel_par1*sqrt(pos_corr_energy_sec1)+ messel_par2*pos_corr_energy_sec1);
			h1_sum_charge_messel_sqrt->Fill(messel_par0+ messel_par1*sqrt(pos_corr_energy_sec0)+ messel_par2*pos_corr_energy_sec0+messel_par0+ messel_par1*sqrt(pos_corr_energy_sec1)+ messel_par2*pos_corr_energy_sec1);
			h1_sum_charge_messel->Fill(pos_corr_energy_sec0+pos_corr_energy_sec1);
 			//check precalibrated state
 			h1_energy_section1->Fill(e_sum_sec_1);
 			h1_energy_section0->Fill(e_sum_sec_0);
			//if (pos_corr_energy_sec0 < 10 || pos_corr_energy_sec1 < 10){
			//	cout << "this is eventnr:\t" << evtnr << endl; 
			//	cout << "e_sum_sec_0:\t" << e_sum_sec_0 << endl;
			//	cout << "e_sum_sec_1:\t" << e_sum_sec_1 << endl;
			//	cout << "beta_corr_energy_sec0:\t" << beta_corr_energy_sec0 << endl;
			//	cout << "beta_corr_energy_sec1:\t" << beta_corr_energy_sec1 << endl;
			//	cout << "pos_corr_energy_sec0:\t" << pos_corr_energy_sec0 << endl;
			//	cout << "pos_corr_energy_sec1:\t" << pos_corr_energy_sec1 << endl;
			//	cout << "spine eval messel:\t" <<  spline_messel->Eval(y_sec0[7]) << endl;
			//	cout << "spline eval messel:\t" << spline_messel->Eval(y_sec1[7]) << endl;
			//	}

                }
			delete  gr_sec0;
			delete  gr_sec1;
			delete fit_sec0;
			delete fit_sec1;	
				}
			//subsubcase 2X2Y-R: one right down && one right up
			if (e_sum_sec_2 > 0 && e_sum_sec_3 > 0 && min_sec2 != 0 && min_sec3 != 0 && e_sum_sec_0 == 0 && e_sum_sec_1 == 0 && (abs(sec_2_arr[15]- sec_2_arr[0])) < 3000 && (abs(sec_3_arr[15]- sec_3_arr[0])) < 3000){
			//	cout << "subsubcase 2X2Y-R: one right down && one right up" << endl;

			double arr_xal_sec_2[16] = {-1000.};
			double arr_xal_sec_3[16] = {-1000.};
			for (Int_t j = 0; j < 16; j++){
				Double_t anode_z_pos = pos_twim_entrance + j*twim_anodes_width;
				arr_xal_sec_2[j] = v_dt_twim[2][j][0]*(drift_sec_2_arr[j]-Tref_section2) + v_dt_twim[2][j][1];	
				arr_xal_sec_3[j] = v_dt_twim[3][j][0]*(drift_sec_3_arr[j]-Tref_section3) +v_dt_twim[3][j][1];	
				}


			//just use the first 10 anode dt values to make the fit
			Int_t n_fitpoints = 10;
			Double_t x[10], y_sec2[10], y_sec3[10];
			for ( Int_t j = 3; j < 13; j++){
				x[j-3] = pos_twim_entrance + j*twim_anodes_width;
				y_sec2[j-3] = arr_xal_sec_2[j] + (spline_sec2[j]->Eval(arr_xal_sec_2[j]));
				y_sec3[j-3] = arr_xal_sec_3[j] + (spline_sec3[j]->Eval(arr_xal_sec_3[j]));
				}
			TGraph* gr_sec2 = new TGraph(n_fitpoints,x,y_sec2);
			TGraph* gr_sec3 = new TGraph(n_fitpoints,x,y_sec3);
			TF1* fit_sec2 = new TF1("fit_sec2","[0]*x+ [1]");
			TF1* fit_sec3 = new TF1("fit_sec3","[0]*x+ [1]");
			gr_sec2->Fit(fit_sec2,"Q");
			gr_sec3->Fit(fit_sec3,"Q");
			// Get the fit parameters for section 2
			Double_t slope_sec2 = fit_sec2->GetParameter(0);
			Double_t offset_sec2 = fit_sec2->GetParameter(1);

			// Get the fit parameters for section 3
			Double_t slope_sec3 = fit_sec3->GetParameter(0);
			Double_t offset_sec3 = fit_sec3->GetParameter(1);
			
			//Calculate psi_in (from slope)
			Double_t psi_in_sec2 = atan(slope_sec2);
			Double_t psi_in_sec3 = atan(slope_sec3);

			Double_t mw1_x_sec2 = slope_sec2*4.5 + offset_sec2;
			Double_t mw1_x_sec3 = slope_sec3*4.5 + offset_sec3;
			//cout << "eventnr:\t" << evtnr << endl;
			if (abs(psi_in_sec2) > 1){
				cout << "eventnr:\t" << evtnr << endl;
			for (int ll = 3; ll < 13;ll++){
				cout << "anode in mm  " << arr_xal_sec_2[ll] << endl;
				}
			}
			vector<vector<double> > tof_and_y = tof_and_y_calibrated(SofMwpc3CalData,SofTofWSingleTcalData);
            if (tof_and_y != dummy_vec){
            vector<double> pathlength_sec2 = path_tof_calc(psi_in_sec2,mw1_x_sec2,tof_and_y[1][2],tof_and_y[1][0],start_time);
            vector<double> pathlength_sec3 = path_tof_calc(psi_in_sec3,mw1_x_sec3,tof_and_y[0][2],tof_and_y[0][0],start_time);
            Double_t full_path_sec2 = sqrt(pow(pathlength_sec2[0],2)+pow(0.1*(tof_and_y[1][1]),2));
            Double_t full_path_sec3 = sqrt(pow(pathlength_sec3[0],2)+pow(0.1*(tof_and_y[0][1]),2));
            Double_t time_of_flight_sec2 = pathlength_sec2[1];
            Double_t time_of_flight_sec3 = pathlength_sec3[1];
	    		Double_t beta2 = ((full_path_sec2/time_of_flight_sec2)*pow(10,7))/299792458;
			Double_t beta3 = ((full_path_sec3/time_of_flight_sec3)*pow(10,7))/299792458;
			h2_beta_vs_energy->Fill(beta2,e_sum_sec_2);
			h2_beta_vs_energy_sec2->Fill(beta2,e_sum_sec_2);
			h2_beta_vs_energy_corr->Fill(beta2,e_sum_sec_2*mean_ene/(beta_func_sec2->Eval(beta2)));
			h2_beta_vs_energy_corr_sec2->Fill(beta2,e_sum_sec_2*mean_ene/(beta_func_sec2->Eval(beta2)));
			h2_beta_vs_energy->Fill(beta3,e_sum_sec_3);
			h2_beta_vs_energy_sec3->Fill(beta3,e_sum_sec_3);
			h2_beta_vs_energy_corr->Fill(beta3,e_sum_sec_3*mean_ene/(beta_func_sec3->Eval(beta3)));
			h2_beta_vs_energy_corr_sec3->Fill(beta3,e_sum_sec_3*mean_ene/(beta_func_sec3->Eval(beta3)));
			h2_x_music_vs_energy_wix->Fill(y_sec2[7],e_sum_sec_2*mean_ene/(beta_func_sec2->Eval(beta2)));
			h2_x_music_vs_energy_wix->Fill(y_sec3[7],e_sum_sec_3*mean_ene/(beta_func_sec3->Eval(beta3)));
			h2_x_music_vs_energy_wix_sec2->Fill(y_sec2[7],e_sum_sec_2*mean_ene/(beta_func_sec2->Eval(beta2)));
			h2_x_music_vs_energy_wix_sec3->Fill(y_sec3[7],e_sum_sec_3*mean_ene/(beta_func_sec3->Eval(beta3)));
			//now comes the final energy correction with positioning tspline
			Double_t beta_corr_energy_sec3 = e_sum_sec_3*mean_ene/(beta_func_sec3->Eval(beta3));
			Double_t beta_corr_energy_sec2 = e_sum_sec_2*mean_ene/(beta_func_sec2->Eval(beta2));
			Double_t pos_corr_energy_sec3 = beta_corr_energy_sec3*mean_pos_sec3/(spline_wix_sec3->Eval(y_sec3[7]));
			Double_t pos_corr_energy_sec2 = beta_corr_energy_sec2*mean_pos_sec2/(spline_wix_sec2->Eval(y_sec2[7]));
			Double_t charge_sec2 = sec2_c0 + sec2_c1*sqrt(pos_corr_energy_sec2) +sec2_c2*pos_corr_energy_sec2;
			Double_t charge_sec3 = sec3_c0 + sec3_c1*sqrt(pos_corr_energy_sec3) +sec3_c2*pos_corr_energy_sec3;
			h2_x_energy_vs_energy_wix_corr->Fill(charge_sec2,charge_sec3);
			h2_charge_sum_vs_diff_wix_corr->Fill(charge_sec3-charge_sec2,charge_sec2+charge_sec3);
			

			h2_x_music_vs_energy_wix_corr->Fill(y_sec2[7],pos_corr_energy_sec2);
			h2_x_music_vs_energy_wix_corr->Fill(y_sec3[7],pos_corr_energy_sec3);
			h1_charge_music_corr_wix->Fill(pos_corr_energy_sec2);
			h1_charge_music_corr_wix->Fill(pos_corr_energy_sec3);
			//best possible calibration ....
			h1_charge_music_corr_sec2->Fill(pos_corr_energy_sec2);
			h1_charge_music_corr_sec3->Fill(pos_corr_energy_sec3);

			h1_charge_music_corr_sum->Fill(pos_corr_energy_sec3+pos_corr_energy_sec2);
			h1_charge_music_corr_wix_full_cal->Fill(pos_corr_energy_sec2*slope_wix + offset_wix);
			h1_charge_music_corr_wix_full_cal->Fill(pos_corr_energy_sec3*slope_wix + offset_wix);
			h1_charge_music_corr_sum_full_cal->Fill(pos_corr_energy_sec2*slope_wix + offset_wix+pos_corr_energy_sec3*slope_wix + offset_wix);
			h1_charge_music_corr_wix_full_cal_pol->Fill(messel_c0 +messel_c1*pos_corr_energy_sec2+messel_c2*pow(pos_corr_energy_sec2,2)+messel_c3*pow(pos_corr_energy_sec2,3)+messel_c4*pow(pos_corr_energy_sec2,4)+messel_c5*pow(pos_corr_energy_sec2,5)+messel_c6*pow(pos_corr_energy_sec2,6));
			h1_charge_music_corr_wix_full_cal_pol->Fill(messel_c0 +messel_c1*pos_corr_energy_sec3+messel_c2*pow(pos_corr_energy_sec3,2)+messel_c3*pow(pos_corr_energy_sec3,3)+messel_c4*pow(pos_corr_energy_sec3,4)+messel_c5*pow(pos_corr_energy_sec3,5)+messel_c6*pow(pos_corr_energy_sec3,6));
			h1_charge_music_corr_sum_full_cal_pol->Fill(messel_c0 +messel_c1*pos_corr_energy_sec2+messel_c2*pow(pos_corr_energy_sec2,2)+messel_c3*pow(pos_corr_energy_sec2,3)+messel_c4*pow(pos_corr_energy_sec2,4)+messel_c5*pow(pos_corr_energy_sec2,5)+messel_c6*pow(pos_corr_energy_sec2,6)+messel_c0 +messel_c1*pos_corr_energy_sec3+messel_c2*pow(pos_corr_energy_sec3,2)+messel_c3*pow(pos_corr_energy_sec3,3)+messel_c4*pow(pos_corr_energy_sec3,4)+messel_c5*pow(pos_corr_energy_sec3,5)+messel_c6*pow(pos_corr_energy_sec3,6));
			h1_charge_music_corr_sum_only_wix->Fill(pos_corr_energy_sec2+pos_corr_energy_sec3);
			//plots only with beta correction:
			h1_charge_music_corr_wix_only_beta->Fill(beta_corr_energy_sec2);
			h1_charge_music_corr_wix_only_beta->Fill(beta_corr_energy_sec3);
			h1_charge_music_corr_sum_only_beta->Fill(beta_corr_energy_sec2+beta_corr_energy_sec3);
			//plot with sqrt fit
			h1_one_charge_wixh_sqrt->Fill(wix_par0+ wix_par1*sqrt(pos_corr_energy_sec3)+wix_par2*pos_corr_energy_sec3);
			h1_one_charge_wixh_sqrt->Fill(wix_par0+ wix_par1*sqrt(pos_corr_energy_sec2)+wix_par2*pos_corr_energy_sec2);
			h1_sum_charge_wixh_sqrt->Fill(wix_par0+ wix_par1*sqrt(pos_corr_energy_sec3)+wix_par2*pos_corr_energy_sec3+wix_par0+ wix_par1*sqrt(pos_corr_energy_sec2)+wix_par2*pos_corr_energy_sec2);
			h1_sum_charge_wixh->Fill(pos_corr_energy_sec3+pos_corr_energy_sec2);
 			//check precalibrated state
 			h1_energy_section2->Fill(e_sum_sec_2);
 			h1_energy_section3->Fill(e_sum_sec_3);
			
			//if (pos_corr_energy_sec2 < 10 || pos_corr_energy_sec3 < 10){
			//	cout << "this is eventnr:\t" << evtnr << endl; 
			//	cout << "e_sum_sec_2:\t" << e_sum_sec_2 << endl;
			//	cout << "e_sum_sec_3:\t" << e_sum_sec_3 << endl;
			//	cout << "beta_corr_energy_sec2:\t" << beta_corr_energy_sec2 << endl;
			//	cout << "beta_corr_energy_sec3:\t" << beta_corr_energy_sec3 << endl;
			//	cout << "pos_corr_energy_sec2:\t" << pos_corr_energy_sec2 << endl;
			//	cout << "pos_corr_energy_sec3:\t" << pos_corr_energy_sec3 << endl;
			//	cout << "spine eval wix:\t" << spline_wix->Eval(y_sec2[7]) << endl;
			//	cout << "spline eval wix:\t" << spline_messel->Eval(y_sec3[7]) << endl;
			//	}

			//cout << "pahtlength\t" << full_path_sec2 << endl;
			//cout << "time of flight\t" << time_of_flight_sec2 << endl;
			//cout << "beta" << ((full_path_sec2/time_of_flight_sec2)*pow(10,7))/299792458 << endl;
			//cout << "corresponding energy:\t" << e_sum_sec_3 << endl;		
                }

			delete  gr_sec2;
			delete  gr_sec3;
			delete fit_sec2;
			delete fit_sec3;	
				}

			//}
		//}

			
		delete [] softwimmappeddata;
//		delete [] sofmwpc1caldata;
//		delete [] sofmwpc2caldata; 
		delete [] sofscitcaldata;
		}
	}
char f_out_name[500];
sprintf(f_out_name,"/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/data/beta_corr/beta_vs_energy_subrun_%s.root",input_str.c_str());
TFile * f = new TFile(f_out_name,"RECREATE");
TList *l = new TList();
//for (Int_t i = 0; i < 16; i++){
//	l->Add(h2_deltax_xal_sec0[i]);
//	}
//for (Int_t i = 0; i < 16; i++){
//	l->Add(h2_deltax_xal_sec1[i]);
//	}
//for (Int_t i = 0; i < 16; i++){
//	l->Add(h2_deltax_xal_sec2[i]);
//	}
//for (Int_t i = 0; i < 16; i++){
//	l->Add(h2_deltax_xal_sec3[i]);
//	}
//for (Int_t i = 0; i < 29; i++){
//	l->Add(h2_tofns_vs_ypos[i]);
//	}
//for (Int_t i = 0; i < 29; i++){
//	l->Add(h2_tofns_vs_ypos2X2Y_LD_RU[i]);
//	}
//for (Int_t i = 0; i < 29; i++){
//	l->Add(h2_tofns_vs_ypos2X2Y_LU_RD[i]);
//	}
//for (Int_t i = 0; i < 29; i++){
//	l->Add(h2_tofns_vs_ypos2X2Y_L[i]);
//	}
//for (Int_t i = 0; i < 29; i++){
//	l->Add(h2_tofns_vs_ypos2X2Y_R[i]);
//	}
l->Add(h2_beta_vs_energy);
l->Add(h2_beta_vs_energy_corr_sec0);
l->Add(h2_beta_vs_energy_corr_sec1);
l->Add(h2_beta_vs_energy_corr_sec2);
l->Add(h2_beta_vs_energy_corr_sec3);
l->Add(h2_beta_vs_energy_corr);
l->Add(h2_x_music_vs_energy_messel);
l->Add(h2_x_music_vs_energy_wix);
l->Add(h2_x_music_vs_energy_wix_corr);
l->Add(h2_x_music_vs_energy_messel_corr);
l->Add(h1_charge_music_corr_wix);
l->Add(h1_charge_music_corr_messel);	
l->Add(h1_charge_music_corr_sum);
l->Add(h1_charge_music_corr_wix_full_cal);
l->Add(h1_charge_music_corr_messel_full_cal);
l->Add(h1_charge_music_corr_sum_full_cal);
l->Add(h1_charge_music_corr_wix_full_cal_pol);
l->Add(h1_charge_music_corr_messel_full_cal_pol);
l->Add(h1_charge_music_corr_sum_full_cal_pol);
l->Add(h1_charge_music_corr_sum_only_wix);
l->Add(h1_charge_music_corr_sum_only_messel);
l->Add(h1_charge_music_corr_messel_only_beta);
l->Add(h1_charge_music_corr_wix_only_beta);
l->Add(h1_charge_music_corr_sum_only_beta);
l->Add(h1_one_charge_wixh_sqrt);
l->Add(h1_sum_charge_wixh_sqrt);
l->Add(h1_one_charge_messel_sqrt);
l->Add(h1_sum_charge_messel_sqrt);
l->Add(h1_sum_charge_wixh);
l->Add(h1_sum_charge_messel);
l->Add(h1_energy_section0);
l->Add(h1_energy_section1);
l->Add(h1_energy_section2);
l->Add(h1_energy_section3);
l->Add(h2_beta_vs_energy_sec0);
l->Add(h2_beta_vs_energy_sec1);
l->Add(h2_beta_vs_energy_sec2);
l->Add(h2_beta_vs_energy_sec3);
l->Add(h2_x_music_vs_energy_messel_sec0);
l->Add(h2_x_music_vs_energy_messel_sec1);
l->Add(h2_x_music_vs_energy_wix_sec2);
l->Add(h2_x_music_vs_energy_wix_sec3);
l->Add(h1_charge_music_corr_sec0);
l->Add(h1_charge_music_corr_sec1);
l->Add(h1_charge_music_corr_sec2);
l->Add(h1_charge_music_corr_sec3);
l->Add(h2_x_energy_vs_energy_wix_corr);
l->Add(h2_charge_sum_vs_diff_wix_corr);
l->Add(h2_x_energy_vs_energy_messel_corr);
l->Add(h2_charge_sum_vs_diff_messel_corr);
l->Write("histlist", TObject::kSingleKey);
}
