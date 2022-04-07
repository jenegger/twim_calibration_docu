//insert headers

#include "TTree.h"
#include "TFile.h"
#include <math.h>
#include <string>
#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TClonesArray.h"
#include "TChain.h"
#include "/scratch8/ge37liw/workingspace/exp_s455/my_macros/t_pattern.C"
//Assignment of the TWIM sections:
////0: Messel Down
////1: Messel Up
////2: Wixhausen Up
////3: Wixhausen Down
using namespace std;
//use geometry from S444 in 2020:TODO is that somewhat right??
const double twim_anodes_width = 400./16.;
const double pos_mw1_z = 0.; //chosen by myself
const double pos_mw2_z = 935 - 300; //in mm
const double pos_twim_entrance = (610-200+12.5) -300;
const double twim_anode_frish_grid_dist = 110. ;//mm

void delta_x_drift( const string& input_str){


Long64_t entries_twim = 0;
Long64_t entries_mw1 = 0;
Long64_t entries_mw2 = 0;

string fname = string("/scratch8/ge37liw/workingspace/exp_s455/data/unpacked/s455_03_273_") + input_str + string("_unpacked.root");
string cal_fname = string("/scratch8/ge37liw/workingspace/exp_s455/data/calibrated/with_mw1/s455_03_273_") + input_str + string("_calibrated.root");
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

//----- HISTOS -------------------------
//deltax vs xal for section 0
TH2D* h2_deltax_xal_sec0[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name,"Section 0, Anode %i -> DeltaX vs Xal",i);
h2_deltax_xal_sec0[i] = new TH2D(hist_name,hist_name,1250,-50,200,500,-5,5);
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
h2_deltax_xal_sec1[i] = new TH2D(hist_name,hist_name,1250,-50,200,500,-5,5);
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
h2_deltax_xal_sec2[i] = new TH2D(hist_name,hist_name,1250,-200,50,500,-5,5);
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
h2_deltax_xal_sec3[i] = new TH2D(hist_name,hist_name,1250,-200,50,500,-5,5);
h2_deltax_xal_sec3[i]->GetXaxis()->SetTitle("Xal [mm]");
h2_deltax_xal_sec3[i]->GetYaxis()->SetTitle("DeltaX (Xfit -Xal) [mm]");
h2_deltax_xal_sec3[i]->GetXaxis()->CenterTitle(true);
h2_deltax_xal_sec3[i]->GetYaxis()->CenterTitle(true);
h2_deltax_xal_sec3[i]->GetYaxis()->SetLabelSize(0.045);
h2_deltax_xal_sec3[i]->GetYaxis()->SetTitleSize(0.045);
}
//--------------------------------------


R3BEventHeader* DataCA = new R3BEventHeader();
TBranch* branchData = chain->GetBranch("EventHeader.");
branchData->SetAddress(&DataCA);

TClonesArray* SofTwimMappedData = new TClonesArray("R3BSofTwimMappedData",2);
R3BSofTwimHitData** softwimmappeddata;
TBranch *branchSofTwimMappedData = chain->GetBranch("TwimMappedData");
branchSofTwimMappedData->SetAddress(&SofTwimMappedData);


//parameters for the raw anode energy
fstream fin;
fin.open("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration/parameters_twim_anodes.csv", ios::in); 
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
//parameters for drift time raw
fstream fin2;
fin2.open("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration/xal_vs_dt_params.csv", ios::in); 
string line2, word2;
vector<vector<vector<double> > > v_dt_twim(4,vector<vector<double> >(16, vector<double>(2)));
getline(fin2, line2);
while(fin2.peek()!=EOF) {
	getline(fin2, line2);
	stringstream s(line2);
	vector<string> temp_vec;
	while (getline(s, word2, ',')) {     temp_vec.push_back(word2);
	}
	v_dt_twim[stoi(temp_vec[0])][stoi(temp_vec[1])][0] = stod(temp_vec[2]);
	v_dt_twim[stoi(temp_vec[0])][stoi(temp_vec[1])][1] = stod(temp_vec[3]);
	temp_vec.clear();
}
//

for(Long64_t i=0;i < nevents;i++){
	Long64_t evtnr = i;
	if (i%100000==0)
		cout<<"Processing event for charge analysis "<<i<<endl;
	chain->GetEvent(i);
	entries_twim = SofTwimMappedData->GetEntries();
//	entries_mw1 = SofMwpc1HitData->GetEntries();
//	entries_mw2 = SofMwpc2HitData->GetEntries();
	int trigger_pattern = t_pattern(DataCA);
//	if (entries_twim > 1 && (trigger_pattern == 2 || trigger_pattern == 4 || trigger_pattern == 8 || trigger_pattern ==10)){
if (entries_twim > 1 && trigger_pattern > 1){
		R3BSofTwimMappedData** softwimmappeddata  = new R3BSofTwimMappedData*[entries_twim];
		

		
			
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
	bool tref_anodes = true;
	bool single_filled_anodes = true;
    UShort_t Tref_section0 = 0;
    UShort_t Tref_section1 = 0;
    UShort_t Tref_section2 = 0;
    UShort_t Tref_section3 = 0;
	
		//fill the 4 arrays
	for (Int_t j = 0; j < entries_twim; j++){
		softwimmappeddata[j] = (R3BSofTwimMappedData*)SofTwimMappedData->At(j);
		Int_t twim_section = softwimmappeddata[j]->GetSecID();
		Int_t twim_anode = softwimmappeddata[j]->GetAnodeID();

		if (twim_section == 0 && twim_anode < 16){
			sec_0_arr[twim_anode] = softwimmappeddata[j]->GetEnergy();
			if (drift_sec_0_arr[twim_anode] != 0){
				single_filled_anodes = false;
				continue;
				}
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
			//Double_t temp_tref_sec0 = softwimmappeddata[j]->GetTime();
			//if (temp_tref_sec0 > 0){
			//	v_tref0.push_back(temp_tref_sec0);
			//}
			}
		else if (twim_section == 1 && twim_anode < 16){
			sec_1_arr[twim_anode] = softwimmappeddata[j]->GetEnergy();
			if (drift_sec_1_arr[twim_anode] != 0){
				single_filled_anodes = false;
				continue;
				}
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
			//Double_t temp_tref_sec1 = softwimmappeddata[j]->GetTime();
			//if (temp_tref_sec1 > 0){
			//	v_tref1.push_back(temp_tref_sec1);
			//}
			}
		else if (twim_section == 2 && twim_anode < 16){
			sec_2_arr[twim_anode] = softwimmappeddata[j]->GetEnergy();
			if (drift_sec_2_arr[twim_anode] != 0){
				single_filled_anodes = false;
				continue;
				}
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
			//Double_t temp_tref_sec2 = softwimmappeddata[j]->GetTime();
			//if (temp_tref_sec2 > 0){
			//	v_tref2.push_back(temp_tref_sec2);
			//}
			}
		else if (twim_section == 3 && twim_anode < 16){
			sec_3_arr[twim_anode] = softwimmappeddata[j]->GetEnergy();
			if (drift_sec_3_arr[twim_anode] != 0){
				single_filled_anodes = false;
				continue;
				}
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
			//Double_t temp_tref_sec3 = softwimmappeddata[j]->GetTime();
			//if (temp_tref_sec3 > 0){
			//	v_tref3.push_back(temp_tref_sec3);
			//}
			}

		}//end of filling arrays

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
		if (tref_anodes && single_filled_anodes){	
		//LD-RU (section 0 and 2)
		if (min_sec0 > 0 && min_sec2 > 0 && (e_sum_sec_0 > 0) && (e_sum_sec_2 > 0) && (e_sum_sec_1 == 0) && (e_sum_sec_3 == 0)){
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
			for ( Int_t j = 1; j < 11; j++){
				x[j-1] = pos_twim_entrance + j*twim_anodes_width;
				y_sec0[j-1] = arr_xal_sec_0[j];
				y_sec2[j-1] = arr_xal_sec_2[j];
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

			//Get the deltax values
			 double arr_deltax_sec_0[16] = {-1000.};
			 double arr_deltax_sec_2[16] = {-1000.};
			for ( Int_t j = 0; j < 16; j++){
				arr_deltax_sec_0[j] = (slope_sec0*(pos_twim_entrance + j*twim_anodes_width) + offset_sec0) - arr_xal_sec_0[j];
				arr_deltax_sec_2[j] = (slope_sec2*(pos_twim_entrance + j*twim_anodes_width) + offset_sec2) - arr_xal_sec_2[j];
				h2_deltax_xal_sec0[j]->Fill(arr_xal_sec_0[j],arr_deltax_sec_0[j]);
				h2_deltax_xal_sec2[j]->Fill(arr_xal_sec_2[j],arr_deltax_sec_2[j]);
				}
			delete  gr_sec0;
			delete  gr_sec2;
			delete fit_sec0;
			delete fit_sec2;	
			}

		//LD-RD (section 0 and 3)
		if (min_sec0 > 0 && min_sec3 > 0 && (e_sum_sec_0 > 0) && (e_sum_sec_3 > 0) && (e_sum_sec_1 == 0) && (e_sum_sec_2 == 0)){
			double arr_xal_sec_0[16] = {-1000.};
			double arr_xal_sec_3[16] = {-1000.};
			for (Int_t j = 0; j < 16; j++){
				Double_t anode_z_pos = pos_twim_entrance + j*twim_anodes_width;
				arr_xal_sec_0[j] = v_dt_twim[0][j][0]*(drift_sec_0_arr[j]-Tref_section0) + v_dt_twim[0][j][1];	
				arr_xal_sec_3[j] = v_dt_twim[3][j][0]*(drift_sec_3_arr[j]-Tref_section3) +v_dt_twim[3][j][1];	
				}


			//just use the first 10 anode dt values to make the fit
			Int_t n_fitpoints = 10;
			Double_t x[10], y_sec0[10], y_sec3[10];
			for ( Int_t j = 1; j < 11; j++){
				x[j-1] = pos_twim_entrance + j*twim_anodes_width;
				y_sec0[j-1] = arr_xal_sec_0[j];
				y_sec3[j-1] = arr_xal_sec_3[j];
				}
			TGraph* gr_sec0 = new TGraph(n_fitpoints,x,y_sec0);
			TGraph* gr_sec3 = new TGraph(n_fitpoints,x,y_sec3);
			TF1* fit_sec0 = new TF1("fit_sec0","[0]*x+ [1]");
			TF1* fit_sec3 = new TF1("fit_sec3","[0]*x+ [1]");
			gr_sec0->Fit(fit_sec0,"Q");
			gr_sec3->Fit(fit_sec3,"Q");
			// Get the fit parameters for section 0
			Double_t slope_sec0 = fit_sec0->GetParameter(0);
			Double_t offset_sec0 = fit_sec0->GetParameter(1);

			// Get the fit parameters for section 3
			Double_t slope_sec3 = fit_sec3->GetParameter(0);
			Double_t offset_sec3 = fit_sec3->GetParameter(1);

			//Get the deltax values
			 double arr_deltax_sec_0[16] = {-1000.};
			 double arr_deltax_sec_3[16] = {-1000.};
			for ( Int_t j = 0; j < 16; j++){
				arr_deltax_sec_0[j] = (slope_sec0*(pos_twim_entrance + j*twim_anodes_width) + offset_sec0) - arr_xal_sec_0[j];
				arr_deltax_sec_3[j] = (slope_sec3*(pos_twim_entrance + j*twim_anodes_width) + offset_sec3) - arr_xal_sec_3[j];
				h2_deltax_xal_sec0[j]->Fill(arr_xal_sec_0[j],arr_deltax_sec_0[j]);
				h2_deltax_xal_sec3[j]->Fill(arr_xal_sec_3[j],arr_deltax_sec_3[j]);
				}
			delete  gr_sec0;
			delete  gr_sec3;
			delete fit_sec0;
			delete fit_sec3;	

			}

		//LU-RU (section 1 and 2)
		if (min_sec1 > 0 && min_sec2 > 0 && (e_sum_sec_1 > 0) && (e_sum_sec_2 > 0) && (e_sum_sec_0 == 0) && (e_sum_sec_3 == 0)){
			double arr_xal_sec_1[16] = {-1000.};
			double arr_xal_sec_2[16] = {-1000.};
			for (Int_t j = 0; j < 16; j++){
				Double_t anode_z_pos = pos_twim_entrance + j*twim_anodes_width;
				arr_xal_sec_1[j] = v_dt_twim[1][j][0]*(drift_sec_1_arr[j]-Tref_section1) + v_dt_twim[1][j][1];	
				arr_xal_sec_2[j] = v_dt_twim[2][j][0]*(drift_sec_2_arr[j]-Tref_section2) +v_dt_twim[2][j][1];	
				}
			//test tj
			for (Int_t s = 0; s < 16; s++){
				if (arr_xal_sec_1[s] < 0){
					cout << "Important event:\t" << evtnr << endl;
					cout << "Corresponding drift time:\t" << drift_sec_1_arr[s]-Tref_section1 << endl;
					}
				}
			//just use the first 10 anode dt values to make the fit
			Int_t n_fitpoints = 10;
			Double_t x[10], y_sec1[10], y_sec2[10];
			for ( Int_t j = 1; j < 11; j++){
				x[j-1] = pos_twim_entrance + j*twim_anodes_width;
				y_sec1[j-1] = arr_xal_sec_1[j];
				y_sec2[j-1] = arr_xal_sec_2[j];
				}
			TGraph* gr_sec1 = new TGraph(n_fitpoints,x,y_sec1);
			TGraph* gr_sec2 = new TGraph(n_fitpoints,x,y_sec2);
			TF1* fit_sec1 = new TF1("fit_sec1","[0]*x+ [1]");
			TF1* fit_sec2 = new TF1("fit_sec2","[0]*x+ [1]");
			gr_sec1->Fit(fit_sec1,"Q");
			gr_sec2->Fit(fit_sec2,"Q");
			// Get the fit parameters for section 1
			Double_t slope_sec1 = fit_sec1->GetParameter(0);
			Double_t offset_sec1 = fit_sec1->GetParameter(1);

			// Get the fit parameters for section 2
			Double_t slope_sec2 = fit_sec2->GetParameter(0);
			Double_t offset_sec2 = fit_sec2->GetParameter(1);

			//Get the deltax values
			 double arr_deltax_sec_1[16] = {-1000.};
			 double arr_deltax_sec_2[16] = {-1000.};
			for ( Int_t j = 0; j < 16; j++){
				arr_deltax_sec_1[j] = (slope_sec1*(pos_twim_entrance + j*twim_anodes_width) + offset_sec1) - arr_xal_sec_1[j];
				arr_deltax_sec_2[j] = (slope_sec2*(pos_twim_entrance + j*twim_anodes_width) + offset_sec2) - arr_xal_sec_2[j];
				h2_deltax_xal_sec1[j]->Fill(arr_xal_sec_1[j],arr_deltax_sec_1[j]);
				h2_deltax_xal_sec2[j]->Fill(arr_xal_sec_2[j],arr_deltax_sec_2[j]);
				}
			delete  gr_sec1;
			delete  gr_sec2;
			delete fit_sec1;
			delete fit_sec2;	
			}

		//LU-RD (section 1 and 3)
		if (min_sec1 > 0 && min_sec3 > 0 && (e_sum_sec_1 > 0) && (e_sum_sec_3 > 0) && (e_sum_sec_0 == 0) && (e_sum_sec_2 == 0)){
			double arr_xal_sec_1[16] = {-1000.};
			double arr_xal_sec_3[16] = {-1000.};
			for (Int_t j = 0; j < 16; j++){
				Double_t anode_z_pos = pos_twim_entrance + j*twim_anodes_width;
				arr_xal_sec_1[j] = v_dt_twim[1][j][0]*(drift_sec_1_arr[j]-Tref_section1) + v_dt_twim[1][j][1];	
				arr_xal_sec_3[j] = v_dt_twim[3][j][0]*(drift_sec_3_arr[j]-Tref_section3)+v_dt_twim[3][j][1];	
				}
			//test tj
			for (Int_t s = 0; s < 16; s++){
				if (arr_xal_sec_1[s] < 0){
					cout << "Important event:\t" << evtnr << endl;
					cout << "Corresponding drift time:\t" << drift_sec_1_arr[s]-Tref_section1 << endl;
					}
				}
			//just use the first 10 anode dt values to make the fit
			Int_t n_fitpoints = 10;
			Double_t x[10], y_sec1[10], y_sec3[10];
			for ( Int_t j = 1; j < 11; j++){
				x[j-1] = pos_twim_entrance + j*twim_anodes_width;
				y_sec1[j-1] = arr_xal_sec_1[j];
				y_sec3[j-1] = arr_xal_sec_3[j];
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

			//Get the deltax values
			 double arr_deltax_sec_1[16] = {-1000.};
			 double arr_deltax_sec_3[16] = {-1000.};
			for ( Int_t j = 0; j < 16; j++){
				arr_deltax_sec_1[j] = (slope_sec1*(pos_twim_entrance + j*twim_anodes_width) + offset_sec1) - arr_xal_sec_1[j];
				arr_deltax_sec_3[j] = (slope_sec3*(pos_twim_entrance + j*twim_anodes_width) + offset_sec3) - arr_xal_sec_3[j];
				h2_deltax_xal_sec1[j]->Fill(arr_xal_sec_1[j],arr_deltax_sec_1[j]);
				h2_deltax_xal_sec3[j]->Fill(arr_xal_sec_3[j],arr_deltax_sec_3[j]);
				}
			delete  gr_sec1;
			delete  gr_sec3;
			delete fit_sec1;
			delete fit_sec3;	

			}	
			}
		delete [] softwimmappeddata;	
	}
}
char f_out_name[500];
sprintf(f_out_name,"deltax_xal_drift_plots_subrun_%s.root",input_str.c_str());
TFile * f = new TFile(f_out_name,"RECREATE");
TList *l = new TList();
for (Int_t i = 0; i < 16; i++){
	l->Add(h2_deltax_xal_sec0[i]);
	}
for (Int_t i = 0; i < 16; i++){
	l->Add(h2_deltax_xal_sec1[i]);
	}
for (Int_t i = 0; i < 16; i++){
	l->Add(h2_deltax_xal_sec2[i]);
	}
for (Int_t i = 0; i < 16; i++){
	l->Add(h2_deltax_xal_sec3[i]);
	}
l->Write("histlist", TObject::kSingleKey);
}

