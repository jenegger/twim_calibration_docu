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
#include "/scratch8/ge37liw/workingspace/exp_s455/my_macros/t_pattern.C"
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
                {-1000,-1000,-1000,-1000},
                {-1000,-1000,-1000,-1000}
                };
 

void u238_primary_analysis(const string& input_str){

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

//keep all mean values same
const double mean_pos_sec0 = 226766;
const double mean_pos_sec1 = 226766;
const double mean_pos_sec2 = 226766;
const double mean_pos_sec3 = 226766;
vector<double> dummy_vec{-1000,-1000,-1000,-1000};

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

const static double mean_ene = 250000.;

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
const double mag_field = 2470*0.0006527728074785267;
const double light_c = 299792458.;
const double gamma_given = 1.68385; //equals to 637 AMeV

TF1* beta_func_sec0 = new TF1("beta_func_sec0","122298*(TMath::Power(x,-5./3.)) + 9355.66");
TF1* beta_func_sec1 = new TF1("beta_func_sec1","122298*(TMath::Power(x,-5./3.)) + 9355.66");
TF1* beta_func_sec2 = new TF1("beta_func_sec2","99421.5*(TMath::Power(x,-5./3.)) + 32183.6");
TF1* beta_func_sec3 = new TF1("beta_func_sec3","110607*(TMath::Power(x,-5./3.)) + 16275.1");


//SOME HISTOS
TH2D* h2_charge_vs_pos_only_sec0;
sprintf(hist_name,"Charge vs position in TWIM, when primary 238U only through sec0");
h2_charge_vs_pos_only_sec0 = new TH2D(hist_name,hist_name,500,60,110,440,-110,110);
h2_charge_vs_pos_only_sec0->GetXaxis()->SetTitle("Charge");
h2_charge_vs_pos_only_sec0->GetYaxis()->SetTitle("Position TWIM [mm]");
h2_charge_vs_pos_only_sec0->GetXaxis()->CenterTitle(true);
h2_charge_vs_pos_only_sec0->GetYaxis()->CenterTitle(true);
h2_charge_vs_pos_only_sec0->GetYaxis()->SetLabelSize(0.045);
h2_charge_vs_pos_only_sec0->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_charge_vs_pos_only_sec1;
sprintf(hist_name,"Charge vs position in TWIM, when primary 238U only through sec1");
h2_charge_vs_pos_only_sec1 = new TH2D(hist_name,hist_name,500,60,110,440,-110,110);
h2_charge_vs_pos_only_sec1->GetXaxis()->SetTitle("Charge");
h2_charge_vs_pos_only_sec1->GetYaxis()->SetTitle("Position TWIM [mm]");
h2_charge_vs_pos_only_sec1->GetXaxis()->CenterTitle(true);
h2_charge_vs_pos_only_sec1->GetYaxis()->CenterTitle(true);
h2_charge_vs_pos_only_sec1->GetYaxis()->SetLabelSize(0.045);
h2_charge_vs_pos_only_sec1->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_charge_vs_pos_sec0_sec1;
sprintf(hist_name,"Charge vs position in TWIM, when primary 238U through sec0&1");
h2_charge_vs_pos_sec0_sec1 = new TH2D(hist_name,hist_name,500,60,110,440,-110,110);
h2_charge_vs_pos_sec0_sec1->GetXaxis()->SetTitle("Charge");
h2_charge_vs_pos_sec0_sec1->GetYaxis()->SetTitle("Position TWIM [mm]");
h2_charge_vs_pos_sec0_sec1->GetXaxis()->CenterTitle(true);
h2_charge_vs_pos_sec0_sec1->GetYaxis()->CenterTitle(true);
h2_charge_vs_pos_sec0_sec1->GetYaxis()->SetLabelSize(0.045);
h2_charge_vs_pos_sec0_sec1->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_charge_vs_pos_only_sec2;
sprintf(hist_name,"Charge vs position in TWIM, when primary 238U only through sec2");
h2_charge_vs_pos_only_sec2 = new TH2D(hist_name,hist_name,500,60,110,440,-110,110);
h2_charge_vs_pos_only_sec2->GetXaxis()->SetTitle("Charge");
h2_charge_vs_pos_only_sec2->GetYaxis()->SetTitle("Position TWIM [mm]");
h2_charge_vs_pos_only_sec2->GetXaxis()->CenterTitle(true);
h2_charge_vs_pos_only_sec2->GetYaxis()->CenterTitle(true);
h2_charge_vs_pos_only_sec2->GetYaxis()->SetLabelSize(0.045);
h2_charge_vs_pos_only_sec2->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_charge_vs_pos_only_sec3;
sprintf(hist_name,"Charge vs position in TWIM, when primary 238U only through sec3");
h2_charge_vs_pos_only_sec3 = new TH2D(hist_name,hist_name,500,60,110,440,-110,110);
h2_charge_vs_pos_only_sec3->GetXaxis()->SetTitle("Charge");
h2_charge_vs_pos_only_sec3->GetYaxis()->SetTitle("Position TWIM [mm]");
h2_charge_vs_pos_only_sec3->GetXaxis()->CenterTitle(true);
h2_charge_vs_pos_only_sec3->GetYaxis()->CenterTitle(true);
h2_charge_vs_pos_only_sec3->GetYaxis()->SetLabelSize(0.045);
h2_charge_vs_pos_only_sec3->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_charge_vs_pos_sec2_sec3;
sprintf(hist_name,"Charge vs position in TWIM, when primary 238U through sec2&3");
h2_charge_vs_pos_sec2_sec3 = new TH2D(hist_name,hist_name,500,60,110,440,-110,110);
h2_charge_vs_pos_sec2_sec3->GetXaxis()->SetTitle("Charge");
h2_charge_vs_pos_sec2_sec3->GetYaxis()->SetTitle("Position TWIM [mm]");
h2_charge_vs_pos_sec2_sec3->GetXaxis()->CenterTitle(true);
h2_charge_vs_pos_sec2_sec3->GetYaxis()->CenterTitle(true);
h2_charge_vs_pos_sec2_sec3->GetYaxis()->SetLabelSize(0.045);
h2_charge_vs_pos_sec2_sec3->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_charge_vs_pos;
sprintf(hist_name,"Charge vs position in TWIM, with 238U primary");
h2_charge_vs_pos = new TH2D(hist_name,hist_name,500,60,110,440,-110,110);
h2_charge_vs_pos->GetXaxis()->SetTitle("Charge");
h2_charge_vs_pos->GetYaxis()->SetTitle("Position TWIM [mm]");
h2_charge_vs_pos->GetXaxis()->CenterTitle(true);
h2_charge_vs_pos->GetYaxis()->CenterTitle(true);
h2_charge_vs_pos->GetYaxis()->SetLabelSize(0.045);
h2_charge_vs_pos->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_mw3_x_pos;
sprintf(hist_name,"x position of primary U238 on MW3 after glad");
h1_mw3_x_pos = new TH1D(hist_name,hist_name,800,-400,400);
h1_mw3_x_pos->GetXaxis()->SetTitle("x_pos_mw3 [mm]");
h1_mw3_x_pos->GetYaxis()->SetTitle("counts");
h1_mw3_x_pos->GetXaxis()->CenterTitle(true);
h1_mw3_x_pos->GetYaxis()->CenterTitle(true);
h1_mw3_x_pos->GetYaxis()->SetLabelSize(0.045);
h1_mw3_x_pos->GetYaxis()->SetTitleSize(0.045);

//END OF HISTOS
for(Long64_t i=0;i < nevents;i++){
	Long64_t evtnr = i;
	if (i%100000==0)
		cout<<"Processing event for charge analysis "<<i<<endl;
	chain->GetEvent(i);
	int trigger_pattern = t_pattern(DataCA);
	//cout << "trigger pattern:\t" << trigger_pattern <<  endl;
	//if (trigger_pattern == 1)  cout << "event with trigger pattern 1:\t" << evtnr << endl;
	entries_twim = SofTwimMappedData->GetEntries();
	entries_tof = SofTofWSingleTcalData->GetEntries();
	
	if (entries_twim > 1 && entries_tof == 1){
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
		double e_sum_no_cal_sec0 = e_sum_sec_0;	
		double e_sum_no_cal_sec1 = e_sum_sec_1;	
		double e_sum_no_cal_sec2 = e_sum_sec_2;	
		double e_sum_no_cal_sec3 = e_sum_sec_3;	
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

			vector<double> v_start_times;
			for (Int_t i = 0; i < 2; i++){
				sofscitcaldata[i] = (R3BSofSciTcalData*)SofSciTcalData->At(i);
				v_start_times.push_back(sofscitcaldata[i]->GetRawTimeNs());
				}
			Double_t start_time = 0.5*(v_start_times[0]-v_start_times[1]); //TODO: substract time rom start to target
			//238U primary beam analysis
			//subcase primary beam goes through sec0 or sec1
			if((e_sum_sec_0 > 0 || e_sum_sec_1 > 0) && e_sum_no_cal_sec2 == 0  && e_sum_no_cal_sec3 == 0){
				double arr_xal_comb[16] = {-1000.};
				double xal_sec0 = -1000;
				double xal_sec1 = -1000;
				
				for (Int_t j = 0; j < 16; j++){
					if (Tref_section0 != 0 && drift_sec_0_arr[j] != 0){
						xal_sec0 = v_dt_twim[0][j][0]*(drift_sec_0_arr[j]-Tref_section0)+v_dt_twim[0][j][1];
						}
					if (Tref_section1 != 0 && drift_sec_1_arr[j] != 0){
						xal_sec1 = v_dt_twim[1][j][0]*(drift_sec_1_arr[j]-Tref_section1)+v_dt_twim[1][j][1];
						}
					if (xal_sec0 == -1000. && xal_sec1 != -1000.){
						arr_xal_comb[j] = xal_sec1 + (spline_sec1[j]->Eval(xal_sec1));
						}	
					if (xal_sec1 == -1000. && xal_sec0 != -1000){
						arr_xal_comb[j] = xal_sec0 + (spline_sec0[j]->Eval(xal_sec0));
						}
					if (xal_sec0 == -1000. && xal_sec1 == -1000.){
						arr_xal_comb[j] = -1000.;
						}
					if (xal_sec0 != -1000 && xal_sec1 != -1000.){
						if (sec_0_arr[j] > sec_1_arr[j]){
							arr_xal_comb[j] = xal_sec0 + (spline_sec0[j]->Eval(xal_sec0));
							}
						if (sec_0_arr[j] < sec_1_arr[j]){
							arr_xal_comb[j] = xal_sec1 + (spline_sec1[j]->Eval(xal_sec1));
							}
						}
					xal_sec0 = -1000;
					xal_sec1 = -1000;
					}	
				Int_t n_fitpoints = 10;
				vector<double> vec_x;
				vector<double> vec_y;
				for ( Int_t j = 3; j < 13; j++){
					if (arr_xal_comb[j] != -1000.){
						vec_x.push_back(pos_twim_entrance + j*twim_anodes_width);
						vec_y.push_back(110. - arr_xal_comb[j]);
						}
						
					}	
				Int_t size_arr = vec_x.size();
				Double_t x[size_arr],y[size_arr];
				for (Int_t j = 0; j < vec_x.size(); j++){
					x[j] = vec_x[j];
					y[j] = vec_y[j];		
					}
				TGraph* gr_u238 = new TGraph(size_arr,x,y);
				TF1* fit_position = new TF1("fit_position","[0]*x+ [1]");
				gr_u238->Fit(fit_position,"Q");
				Double_t slope_u238 = fit_position->GetParameter(0);
				Double_t offset_u238 = fit_position->GetParameter(1);
				Double_t psi_in = atan(slope_u238);
				Double_t mw1_x = slope_u238*4.5 + offset_u238;
				delete gr_u238;
				delete fit_position;
				vector<double> tof_and_y = single_tof_and_y_calibrated(SofMwpc3CalData,SofTofWSingleTcalData);
				if (tof_and_y != dummy_vec){
					h1_mw3_x_pos->Fill(tof_and_y[2]);
					vector<double> pathlength = path_tof_calc(psi_in,mw1_x,tof_and_y[2],tof_and_y[0],start_time);
					Double_t paddle = tof_and_y[3];
					Double_t y_pos = 0.1*(tof_and_y[1]);
					Double_t full_path = sqrt(pow(pathlength[0],2)+pow(0.1*(tof_and_y[1]),2));
					Double_t time_of_flight = pathlength[1];
					Double_t rho = pathlength[2];
					Double_t beta = ((full_path/time_of_flight)*pow(10,7))/299792458;
					Double_t gamma = 1/(sqrt(1-beta*beta));
					Double_t a_q = (pow(10,-3)*((((mag_field*rho)/(beta))*1.602176634*pow(10,-19))/(1.66053906660*pow(10,-27)*light_c)))/gamma;
					Double_t beta_corr_energy_sec0 = e_sum_sec_0*mean_ene/(beta_func_sec0->Eval(beta));
					Double_t beta_corr_energy_sec1 = e_sum_sec_1*mean_ene/(beta_func_sec1->Eval(beta));
					Double_t pos_corr_energy_sec0 = beta_corr_energy_sec0*mean_pos_sec0/(spline_messel_sec0->Eval(y[7]));		
					//cout << "beta corrected energy sec0:\t" << beta_corr_energy_sec0 << endl;
					//cout << "position for sec0:\t" << y[7] << endl;
					//cout << "position corrected energy sec0:\t" << pos_corr_energy_sec0 << endl;
					Double_t pos_corr_energy_sec1 = beta_corr_energy_sec1*mean_pos_sec1/(spline_messel_sec1->Eval(y[7]));		
					Double_t charge_sec0 = sec0_c0 + sec0_c1*sqrt(pos_corr_energy_sec0) +sec0_c2*pos_corr_energy_sec0;
					Double_t charge_sec1 = sec1_c0 + sec1_c1*sqrt(pos_corr_energy_sec1) +sec1_c2*pos_corr_energy_sec1;
					Double_t final_charge0 = charge_sec0 + 30.;
					//cout << "charge 0:\t" << final_charge0 << endl;
					Double_t final_charge1 = charge_sec1 + 29.;
					//cout << "charge 1:\t" << final_charge1 << endl;
					Double_t sum_charge;
					//cout << "e_sum_sec_0:\t" << e_sum_sec_0 << endl;
					//cout << "e_sum_sec_1:\t" << e_sum_sec_1 << endl;
					//cout << "charge0:\t" << final_charge0 << endl;
					//cout << "charge1:\t" << final_charge1 << endl;
					if (e_sum_sec_0 == 0){
						sum_charge = final_charge1;
						h2_charge_vs_pos_only_sec1->Fill(sum_charge,y[7]);
						h2_charge_vs_pos->Fill(sum_charge,y[7]);
						}
					if (e_sum_sec_1 == 0){
						sum_charge = final_charge0;
						h2_charge_vs_pos_only_sec0->Fill(sum_charge,y[7]);
						h2_charge_vs_pos->Fill(sum_charge,y[7]);
						}
					else{
						sum_charge = sqrt(pow(final_charge0,2)+pow(final_charge1,2));
						h2_charge_vs_pos->Fill(sum_charge,y[7]);
						h2_charge_vs_pos_sec0_sec1->Fill(sum_charge,y[7]);
						}
						
					//if (sum_charge > 80) cout << "sum charge:\t" << sum_charge << endl;


					}
				}

			//subcase primary beam goes trhough sec2 or sec3
			if((e_sum_sec_2 > 0 || e_sum_sec_3 > 0) && e_sum_no_cal_sec0 == 0  && e_sum_no_cal_sec1 == 0){
				double arr_xal_comb[16] = {-1000.};
				double xal_sec2 = -1000;
				double xal_sec3 = -1000;
				
				for (Int_t j = 0; j < 16; j++){
					if (Tref_section2 != 0 && drift_sec_2_arr[j] != 0){
						xal_sec2 = v_dt_twim[2][j][0]*(drift_sec_2_arr[j]-Tref_section2)+v_dt_twim[2][j][1];
						}
					if (Tref_section3 != 0 && drift_sec_3_arr[j] != 0){
						xal_sec3 = v_dt_twim[3][j][0]*(drift_sec_3_arr[j]-Tref_section3)+v_dt_twim[3][j][1];
						}
					if (xal_sec2 == -1000. && xal_sec3 != -1000.){
						arr_xal_comb[j] = xal_sec3 + (spline_sec3[j]->Eval(xal_sec3));
						}	
					if (xal_sec3 == -1000. && xal_sec2 != -1000){
						arr_xal_comb[j] = xal_sec2 + (spline_sec2[j]->Eval(xal_sec2));
						}
					if (xal_sec2 == -1000. && xal_sec3 == -1000.){
						arr_xal_comb[j] = -1000.;
						}
					if (xal_sec2 != -1000 && xal_sec3 != -1000.){
						if (sec_2_arr[j] > sec_3_arr[j]){
							arr_xal_comb[j] = xal_sec2 + (spline_sec2[j]->Eval(xal_sec2));
							}
						if (sec_2_arr[j] < sec_3_arr[j]){
							arr_xal_comb[j] = xal_sec3 + (spline_sec3[j]->Eval(xal_sec3));
							}
						}
					xal_sec2 = -1000;
					xal_sec3 = -1000;
					}	
				Int_t n_fitpoints = 10;
				vector<double> vec_x;
				vector<double> vec_y;
				for ( Int_t j = 3; j < 13; j++){
					if (arr_xal_comb[j] != -1000.){
						vec_x.push_back(pos_twim_entrance + j*twim_anodes_width);
						vec_y.push_back(-110. - arr_xal_comb[j]);
						}
						
					}	
				Int_t size_arr = vec_x.size();
				Double_t x[size_arr],y[size_arr];
				for (Int_t j = 0; j < vec_x.size(); j++){
					x[j] = vec_x[j];
					y[j] = vec_y[j];		
					}
				TGraph* gr_u238 = new TGraph(size_arr,x,y);
				TF1* fit_position = new TF1("fit_position","[0]*x+ [1]");
				gr_u238->Fit(fit_position,"Q");
				Double_t slope_u238 = fit_position->GetParameter(0);
				Double_t offset_u238 = fit_position->GetParameter(1);
				Double_t psi_in = atan(slope_u238);
				Double_t mw1_x = slope_u238*4.5 + offset_u238;
				delete gr_u238;
				delete fit_position;
				vector<double> tof_and_y = single_tof_and_y_calibrated(SofMwpc3CalData,SofTofWSingleTcalData);
				if (tof_and_y != dummy_vec){
					h1_mw3_x_pos->Fill(tof_and_y[2]);
					vector<double> pathlength = path_tof_calc(psi_in,mw1_x,tof_and_y[2],tof_and_y[0],start_time);
					Double_t paddle = tof_and_y[3];
					Double_t y_pos = 0.1*(tof_and_y[1]);
					Double_t full_path = sqrt(pow(pathlength[0],2)+pow(0.1*(tof_and_y[1]),2));
					Double_t time_of_flight = pathlength[1];
					Double_t rho = pathlength[2];
					Double_t beta = ((full_path/time_of_flight)*pow(10,7))/299792458;
					Double_t gamma = 1/(sqrt(1-beta*beta));
					Double_t a_q = (pow(10,-3)*((((mag_field*rho)/(beta))*1.602176634*pow(10,-19))/(1.66053906660*pow(10,-27)*light_c)))/gamma;
					Double_t beta_corr_energy_sec2 = e_sum_sec_2*mean_ene/(beta_func_sec2->Eval(beta));
					Double_t beta_corr_energy_sec3 = e_sum_sec_3*mean_ene/(beta_func_sec3->Eval(beta));
					Double_t pos_corr_energy_sec2 = beta_corr_energy_sec2*mean_pos_sec2/(spline_wix_sec2->Eval(y[7]));		
					Double_t pos_corr_energy_sec3 = beta_corr_energy_sec3*mean_pos_sec3/(spline_wix_sec3->Eval(y[7]));		
					Double_t charge_sec2 = sec2_c0 + sec2_c1*sqrt(pos_corr_energy_sec2) +sec2_c2*pos_corr_energy_sec2;
					Double_t charge_sec3 = sec3_c0 + sec3_c1*sqrt(pos_corr_energy_sec3) +sec3_c2*pos_corr_energy_sec3;
					Double_t final_charge2 = charge_sec2 + 28.;
					Double_t final_charge3 = charge_sec3 + 28.;
					Double_t sum_charge;
					//cout << "e_sum_sec_2:\t" << e_sum_sec_2 << endl;
					//cout << "e_sum_sec_3:\t" << e_sum_sec_3 << endl;
					//cout << "charge2:\t" << final_charge2 << endl;
					//cout << "charge3:\t" << final_charge3 << endl;
					if (e_sum_sec_2 == 0){
						sum_charge = final_charge3;
						h2_charge_vs_pos_only_sec3->Fill(sum_charge,y[7]);
						h2_charge_vs_pos->Fill(sum_charge,y[7]);
						}
					if (e_sum_sec_3 == 0){
						sum_charge = final_charge2;
						h2_charge_vs_pos_only_sec3->Fill(sum_charge,y[7]);
						h2_charge_vs_pos->Fill(sum_charge,y[7]);
						}
					else{
						sum_charge = sqrt(pow(final_charge2,2)+pow(final_charge3,2));
						h2_charge_vs_pos->Fill(sum_charge,y[7]);
						h2_charge_vs_pos_sec2_sec3->Fill(sum_charge,y[7]);
						}
					}

				}


			
		delete [] softwimmappeddata;
		delete [] sofscitcaldata;
		}
	}
char f_out_name[500];
sprintf(f_out_name,"/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/data/beta_corr/primary_beam_subrun_%s.root",input_str.c_str());
TFile * f = new TFile(f_out_name,"RECREATE");
TList *l = new TList();
l->Add(h2_charge_vs_pos_only_sec0);
l->Add(h2_charge_vs_pos_only_sec1);
l->Add(h2_charge_vs_pos_only_sec2);
l->Add(h2_charge_vs_pos_only_sec3);
l->Add(h2_charge_vs_pos);
l->Add(h2_charge_vs_pos_sec0_sec1);
l->Add(h2_charge_vs_pos_sec2_sec3);
l->Add(h1_mw3_x_pos);
l->Write("histlist", TObject::kSingleKey);
}
