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
 

void tof_calibration(const string& input_str){

Long64_t entries_twim = 0;
Long64_t entries_mw1 = 0;
Long64_t entries_mw2 = 0;

string fname = string("/scratch5/ge37liw/unpacked_newroot_04_2022/s455_03_273_") + input_str + string("_unpacked.root");
string cal_fname = string("/scratch5/ge37liw/calibrated_rand_angles/s455_03_273_") + input_str + string("_calibrated.root");
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

for(Long64_t i=0;i < nevents;i++){
	Long64_t evtnr = i;
	if (i%100000==0)
		cout<<"Processing event for charge analysis "<<i<<endl;
	chain->GetEvent(i);
	entries_twim = SofTwimMappedData->GetEntries();
	entries_mw1 = SofMwpc1CalData->GetEntries();
	entries_mw2 = SofMwpc2CalData->GetEntries();
	if (entries_twim > 1 && entries_mw1 > 1 && entries_mw2 > 1){
		////cout << "this is event with more than one TWIM entry:\t" << evtnr << endl;
		//cout << "corresponding entries:\t" << entries_twim << endl;
		Int_t nr_start;
    	nr_start = SofSciTcalData->GetEntries();
		Int_t nr_tof;
		nr_tof = SofTofWSingleTcalData->GetEntries();
		
		R3BTwimMappedData** softwimmappeddata  = new R3BTwimMappedData*[entries_twim];
		R3BMwpcCalData** sofmwpc1caldata = new R3BMwpcCalData*[entries_mw1];
		R3BMwpcCalData** sofmwpc2caldata = new R3BMwpcCalData*[entries_mw2];
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

		//new section to explicitly calculate positions on MW1 and MW2
		vector<double> vec_X_mw1_up(64,0); //-> 64 vertical pads
		vector<double> vec_X_mw1_down(64,0);
		vector<double> vec_X_mw1_sum(64,0);
		vector<double> vec_X_mw2_up(64,0);
		vector<double> vec_X_mw2_down(64,0);
		vector<double> vec_X_mw2_sum(64,0);
		vector<double> vec_Y_mw1(40,0); //-> 40 horizontal pads
		vector<double> vec_Y_mw2(40,0); //-> 40 horizontal pads
		
		for (Int_t i = 0; i < entries_mw1; i++){
			sofmwpc1caldata[i] = (R3BMwpcCalData*)SofMwpc1CalData->At(i);
			Int_t planeId = sofmwpc1caldata[i]->GetPlane();
			Int_t padId = sofmwpc1caldata[i]->GetPad()-1;
			Double_t charge = sofmwpc1caldata[i]->GetQ();
			if (planeId == 2){  //-> plane down
				vec_X_mw1_down[padId] = charge;
				vec_X_mw1_sum[padId] += charge;
				}
			if  (planeId == 1){ //->plane up
				vec_X_mw1_up[padId] = charge;
				vec_X_mw1_sum[padId] += charge;
				}
			if (planeId == 3){  //plane Y
				vec_Y_mw1[padId] = charge;
				}
			}
		for (Int_t i = 0; i < entries_mw2; i++){
			sofmwpc2caldata[i] = (R3BMwpcCalData*)SofMwpc2CalData->At(i);
			Int_t planeId = sofmwpc2caldata[i]->GetPlane();
			Int_t padId = sofmwpc2caldata[i]->GetPad()-1;
			Double_t charge = sofmwpc2caldata[i]->GetQ();
			if (planeId == 2){  //-> plane down
				vec_X_mw2_down[padId] = charge;
				vec_X_mw2_sum[padId] += charge;
				}
			if  (planeId == 1){ //->plane up
				vec_X_mw2_up[padId] = charge;
				vec_X_mw2_sum[padId] += charge;
				}
			if (planeId == 3){  //plane Y
				vec_Y_mw2[padId] = charge;
				}

			}
		//get the peaks for mw1
		vector<double> peaks_X_mw1_up;
		vector<double> peaks_X_mw1_down;
		vector<double> peaks_X_mw1_sum;

		vector<double> peaks_X_mw2_up;
		vector<double> peaks_X_mw2_down;
		vector<double> peaks_X_mw2_sum;
		
		vector<double> peaks_Y_mw1;
		vector<double> peaks_Y_mw2;
		
		for (Int_t i = 0; i < 2;i++){
			Int_t max_element_index_X_mw1_up = max_element(vec_X_mw1_up.begin(),vec_X_mw1_up.end())-vec_X_mw1_up.begin();
			Double_t max_element_X_mw1_up = *max_element(vec_X_mw1_up.begin(),vec_X_mw1_up.end());
			Int_t max_element_index_X_mw1_down = max_element(vec_X_mw1_down.begin(),vec_X_mw1_down.end())-vec_X_mw1_down.begin();
			Double_t max_element_X_mw1_down = *max_element(vec_X_mw1_down.begin(),vec_X_mw1_down.end());
			Int_t max_element_index_X_mw1_sum = max_element(vec_X_mw1_sum.begin(),vec_X_mw1_sum.end())-vec_X_mw1_sum.begin();
			Double_t max_element_X_mw1_sum = *max_element(vec_X_mw1_sum.begin(),vec_X_mw1_sum.end());
			Int_t max_element_index_X_mw2_up = max_element(vec_X_mw2_up.begin(),vec_X_mw2_up.end())-vec_X_mw2_up.begin();
			Double_t max_element_X_mw2_up = *max_element(vec_X_mw2_up.begin(),vec_X_mw2_up.end());
			Int_t max_element_index_X_mw2_down = max_element(vec_X_mw2_down.begin(),vec_X_mw2_down.end())-vec_X_mw2_down.begin();
			Double_t max_element_X_mw2_down = *max_element(vec_X_mw2_down.begin(),vec_X_mw2_down.end());
			Int_t max_element_index_X_mw2_sum =  max_element(vec_X_mw2_sum.begin(),vec_X_mw2_sum.end())-vec_X_mw2_sum.begin();
			Double_t max_element_X_mw2_sum = *max_element(vec_X_mw2_sum.begin(),vec_X_mw2_sum.end());
			Int_t max_element_index_Y_mw1 = max_element(vec_Y_mw1.begin(),vec_Y_mw1.end())-vec_Y_mw1.begin();
			Double_t max_element_Y_mw1 = *max_element(vec_Y_mw1.begin(),vec_Y_mw1.end());
			Int_t max_element_index_Y_mw2 = max_element(vec_Y_mw2.begin(),vec_Y_mw2.end())-vec_Y_mw2.begin();
			Double_t max_element_Y_mw2 = *max_element(vec_Y_mw2.begin(),vec_Y_mw2.end());
			//mw1_x_up
			if (max_element_X_mw1_up > 200){
				Double_t q_left;
				Double_t q_right;
				Double_t q_max = vec_X_mw1_up[max_element_index_X_mw1_up];
				if (max_element_index_X_mw1_up == 0) q_left=1;
				else if (max_element_index_X_mw1_up == 63) q_right=1;
				else {
				q_left = vec_X_mw1_up[max_element_index_X_mw1_up-1];
				q_right = vec_X_mw1_up[max_element_index_X_mw1_up +1];
				if (q_left == 0) q_left=1;
				if(q_right==0) q_right=1;
				}
				peaks_X_mw1_up.push_back(get_mw12_pos_X(q_max,max_element_index_X_mw1_up,q_left,q_right));
				if (max_element_index_X_mw1_up == 0){
					vec_X_mw1_up[max_element_index_X_mw1_up] = 0;
					vec_X_mw1_up[max_element_index_X_mw1_up+1] = 0;
					}
				else if (max_element_index_X_mw1_up == 63){
					vec_X_mw1_up[max_element_index_X_mw1_up] = 0;
					vec_X_mw1_up[max_element_index_X_mw1_up-1] = 0;	
					}
				else {
					vec_X_mw1_up[max_element_index_X_mw1_up] = 0;
					vec_X_mw1_up[max_element_index_X_mw1_up-1] = 0;
					vec_X_mw1_up[max_element_index_X_mw1_up+1] = 0;
					}
				}

			//mw1_x_down
			if (max_element_X_mw1_down > 200){
				Double_t q_left;
				Double_t q_right;
				Double_t q_max = vec_X_mw1_down[max_element_index_X_mw1_down];
				if (max_element_index_X_mw1_down == 0) q_left=1;
				else if (max_element_index_X_mw1_down == 63) q_right=1;
				else {
				q_left = vec_X_mw1_down[max_element_index_X_mw1_down-1];
				q_right = vec_X_mw1_down[max_element_index_X_mw1_down +1];
				if (q_left == 0) q_left=1;
				if(q_right==0) q_right=1;
				}
				peaks_X_mw1_down.push_back(get_mw12_pos_X(q_max,max_element_index_X_mw1_down,q_left,q_right));
				if (max_element_index_X_mw1_down == 0){
					vec_X_mw1_down[max_element_index_X_mw1_down] = 0;
					vec_X_mw1_down[max_element_index_X_mw1_down+1] = 0;
					}
				else if (max_element_index_X_mw1_down == 63){
					vec_X_mw1_down[max_element_index_X_mw1_down] = 0;
					vec_X_mw1_down[max_element_index_X_mw1_down-1] = 0;	
					}
				else {
					vec_X_mw1_down[max_element_index_X_mw1_down] = 0;
					vec_X_mw1_down[max_element_index_X_mw1_down-1] = 0;
					vec_X_mw1_down[max_element_index_X_mw1_down+1] = 0;
					}
				}

			//mw1_x_sum
			if (max_element_X_mw1_sum > 200){
				Double_t q_left;
				Double_t q_right;
				Double_t q_max = vec_X_mw1_sum[max_element_index_X_mw1_sum];
				if (max_element_index_X_mw1_sum == 0) q_left=1;
				else if (max_element_index_X_mw1_sum == 63) q_right=1;
				else {
				q_left = vec_X_mw1_sum[max_element_index_X_mw1_sum-1];
				q_right = vec_X_mw1_sum[max_element_index_X_mw1_sum +1];
				if (q_left == 0) q_left=1;
				if(q_right==0) q_right=1;
				}
				peaks_X_mw1_sum.push_back(get_mw12_pos_X(q_max,max_element_index_X_mw1_sum,q_left,q_right));
				if (max_element_index_X_mw1_sum == 0){
					vec_X_mw1_sum[max_element_index_X_mw1_sum] = 0;
					vec_X_mw1_sum[max_element_index_X_mw1_sum+1] = 0;
					}
				else if (max_element_index_X_mw1_sum == 63){
					vec_X_mw1_sum[max_element_index_X_mw1_sum] = 0;
					vec_X_mw1_sum[max_element_index_X_mw1_sum-1] = 0;	
					}
				else {
					vec_X_mw1_sum[max_element_index_X_mw1_sum] = 0;
					vec_X_mw1_sum[max_element_index_X_mw1_sum-1] = 0;
					vec_X_mw1_sum[max_element_index_X_mw1_sum+1] = 0;
					}
				}

			//mw2_x_up
			if (max_element_X_mw2_up > 200){
				Double_t q_left;
				Double_t q_right;
				Double_t q_max = vec_X_mw2_up[max_element_index_X_mw2_up];
				if (max_element_index_X_mw2_up == 0) q_left=1;
				else if (max_element_index_X_mw2_up == 63) q_right=1;
				else {
				q_left = vec_X_mw2_up[max_element_index_X_mw2_up-1];
				q_right = vec_X_mw2_up[max_element_index_X_mw2_up +1];
				if (q_left == 0) q_left=1;
				if(q_right==0) q_right=1;
				}
				peaks_X_mw2_up.push_back(get_mw12_pos_X(q_max,max_element_index_X_mw2_up,q_left,q_right));
				if (max_element_index_X_mw2_up == 0){
					vec_X_mw2_up[max_element_index_X_mw2_up] = 0;
					vec_X_mw2_up[max_element_index_X_mw2_up+1] = 0;
					}
				else if (max_element_index_X_mw2_up == 63){
					vec_X_mw2_up[max_element_index_X_mw2_up] = 0;
					vec_X_mw2_up[max_element_index_X_mw2_up-1] = 0;	
					}
				else {
					vec_X_mw2_up[max_element_index_X_mw2_up] = 0;
					vec_X_mw2_up[max_element_index_X_mw2_up-1] = 0;
					vec_X_mw2_up[max_element_index_X_mw2_up+1] = 0;
					}
				}

			//mw2_x_down
			if (max_element_X_mw2_down > 200){
				Double_t q_left;
				Double_t q_right;
				Double_t q_max = vec_X_mw2_down[max_element_index_X_mw2_down];
				if (max_element_index_X_mw2_down == 0) q_left=1;
				else if (max_element_index_X_mw2_down == 63) q_right=1;
				else {
				q_left = vec_X_mw2_down[max_element_index_X_mw2_down-1];
				q_right = vec_X_mw2_down[max_element_index_X_mw2_down +1];
				if (q_left == 0) q_left=1;
				if(q_right==0) q_right=1;
				}
				peaks_X_mw2_down.push_back(get_mw12_pos_X(q_max,max_element_index_X_mw2_down,q_left,q_right));
				if (max_element_index_X_mw2_down == 0){
					vec_X_mw2_down[max_element_index_X_mw2_down] = 0;
					vec_X_mw2_down[max_element_index_X_mw2_down+1] = 0;
					}
				else if (max_element_index_X_mw2_down == 63){
					vec_X_mw2_down[max_element_index_X_mw2_down] = 0;
					vec_X_mw2_down[max_element_index_X_mw2_down-1] = 0;	
					}
				else {
					vec_X_mw2_down[max_element_index_X_mw2_down] = 0;
					vec_X_mw2_down[max_element_index_X_mw2_down-1] = 0;
					vec_X_mw2_down[max_element_index_X_mw2_down+1] = 0;
					}
				}

			//mw2_x_sum
			if (max_element_X_mw2_sum > 200){
				Double_t q_left;
				Double_t q_right;
				Double_t q_max = vec_X_mw2_sum[max_element_index_X_mw2_sum];
				if (max_element_index_X_mw2_sum == 0) q_left=1;
				else if (max_element_index_X_mw2_sum == 63) q_right=1;
				else {
				q_left = vec_X_mw2_sum[max_element_index_X_mw2_sum-1];
				q_right = vec_X_mw2_sum[max_element_index_X_mw2_sum +1];
				if (q_left == 0) q_left=1;
				if(q_right==0) q_right=1;
				}
				peaks_X_mw2_sum.push_back(get_mw12_pos_X(q_max,max_element_index_X_mw2_sum,q_left,q_right));
				if (max_element_index_X_mw2_sum == 0){
					vec_X_mw2_sum[max_element_index_X_mw2_sum] = 0;
					vec_X_mw2_sum[max_element_index_X_mw2_sum+1] = 0;
					}
				else if (max_element_index_X_mw2_sum == 63){
					vec_X_mw2_sum[max_element_index_X_mw2_sum] = 0;
					vec_X_mw2_sum[max_element_index_X_mw2_sum-1] = 0;	
					}
				else {
					vec_X_mw2_sum[max_element_index_X_mw2_sum] = 0;
					vec_X_mw2_sum[max_element_index_X_mw2_sum-1] = 0;
					vec_X_mw2_sum[max_element_index_X_mw2_sum+1] = 0;
					}
				}


			//mw2_Y
			if (max_element_Y_mw2 > 200){
				Double_t q_left;
				Double_t q_right;
				Double_t q_max = vec_Y_mw2[max_element_index_Y_mw2];
				if (max_element_index_Y_mw2 == 0) q_left=1;
				else if (max_element_index_Y_mw2 == 63) q_right=1;
				else {
				q_left = vec_Y_mw2[max_element_index_Y_mw2-1];
				q_right = vec_Y_mw2[max_element_index_Y_mw2 +1];
				if (q_left == 0) q_left=1;
				if(q_right==0) q_right=1;
				}
				peaks_Y_mw2.push_back(get_mw12_pos_Y(q_max,max_element_index_Y_mw2,q_left,q_right));
				if (max_element_index_Y_mw2 == 0){
					vec_Y_mw2[max_element_index_Y_mw2] = 0;
					vec_Y_mw2[max_element_index_Y_mw2+1] = 0;
					}
				else if (max_element_index_Y_mw2 == 63){
					vec_Y_mw2[max_element_index_Y_mw2] = 0;
					vec_Y_mw2[max_element_index_Y_mw2-1] = 0;	
					}
				else {
					vec_Y_mw2[max_element_index_Y_mw2] = 0;
					vec_Y_mw2[max_element_index_Y_mw2-1] = 0;
					vec_Y_mw2[max_element_index_Y_mw2+1] = 0;
					}
				}


			//mw1_Y
			if (max_element_Y_mw1 > 200){
				Double_t q_left;
				Double_t q_right;
				Double_t q_max = vec_Y_mw1[max_element_index_Y_mw1];
				if (max_element_index_Y_mw1 == 0) q_left=1;
				else if (max_element_index_Y_mw1 == 63) q_right=1;
				else {
				q_left = vec_Y_mw1[max_element_index_Y_mw1-1];
				q_right = vec_Y_mw1[max_element_index_Y_mw1 +1];
				if (q_left == 0) q_left=1;
				if(q_right==0) q_right=1;
				}
				peaks_Y_mw1.push_back(get_mw12_pos_Y(q_max,max_element_index_Y_mw1,q_left,q_right));
				if (max_element_index_Y_mw1 == 0){
					vec_Y_mw1[max_element_index_Y_mw1] = 0;
					vec_Y_mw1[max_element_index_Y_mw1+1] = 0;
					}
				else if (max_element_index_Y_mw1 == 63){
					vec_Y_mw1[max_element_index_Y_mw1] = 0;
					vec_Y_mw1[max_element_index_Y_mw1-1] = 0;	
					}
				else {
					vec_Y_mw1[max_element_index_Y_mw1] = 0;
					vec_Y_mw1[max_element_index_Y_mw1-1] = 0;
					vec_Y_mw1[max_element_index_Y_mw1+1] = 0;
					}
				}

			
			}
		
		//this is the new part 01.02.2022
		//FF1 is more to the left than FF2-> FF1.x > FF2.x
		//for calibration I only use events with Mw2_Xsum_nHits >1 && Mw2_Y_nHits >1
	Double_t FF1_slope;
	Double_t FF2_slope;
	Double_t FF1_offset;
	Double_t FF2_offset;
	if (peaks_X_mw1_sum.size() > 1 && peaks_X_mw2_sum.size() > 1 && peaks_Y_mw1.size() > 1 && peaks_Y_mw2.size() > 1 && single_filled_anodes && tref_anodes && nr_start > 1 && nr_tof ==2){

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
		if (peaks_X_mw2_down.size() > 0 && peaks_X_mw2_up.size() > 0  && peaks_X_mw1_down.size() > 0 && peaks_X_mw1_up.size() > 0 && ((peaks_Y_mw2[0] <= 0 && peaks_Y_mw2[1] > 0) || (peaks_Y_mw2[0] > 0 && peaks_Y_mw2[1] <= 0))  && ((peaks_Y_mw1[0] <= 0 && peaks_Y_mw1[1] > 0) || (peaks_Y_mw1[0] > 0 && peaks_Y_mw1[1] <= 0))){
			//initialize start detector info
			vector<double> v_start_times;
			for (Int_t i = 0; i < 2; i++){
				sofscitcaldata[i] = (R3BSofSciTcalData*)SofSciTcalData->At(i);
				v_start_times.push_back(sofscitcaldata[i]->GetRawTimeNs());
				}
			Double_t start_time = 0.5*(v_start_times[0]+v_start_times[1]); //TODO: substract time rom start to target
			
			//subsubcase 2X2Y-LD-RU: one left down && one right up
			if (peaks_X_mw1_down[0] > 0 && peaks_X_mw2_down[0] > 0 && peaks_X_mw1_up[0] < 0 && peaks_X_mw2_up[0] < 0 && e_sum_sec_0 > 0 && e_sum_sec_2 > 0 && e_sum_sec_1 == 0 && e_sum_sec_3 == 0 && (abs(sec_0_arr[15]- sec_0_arr[0])) < 3000 && (abs(sec_2_arr[15]- sec_2_arr[0])) < 3000){
				//cout << "subsubcase 2X2Y-LD-RU: one left down && one right up" << endl;
				cout <<"before"<< endl;
			double arr_xal_sec_0[16] = {-1000.};
			double arr_xal_sec_2[16] = {-1000.};
			for (Int_t j = 0; j < 16; j++){
				Double_t anode_z_pos = pos_twim_entrance + j*twim_anodes_width;
				arr_xal_sec_0[j] = v_dt_twim[0][j][0]*(drift_sec_0_arr[j]-Tref_section0)+v_dt_twim[0][j][1];	
				arr_xal_sec_2[j] = v_dt_twim[2][j][0]*(drift_sec_2_arr[j]-Tref_section2)+v_dt_twim[2][j][1];	
				}
			cout << "in the middle" << endl;
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
			cout << "after fitting" << endl;
			//Calculate psi_in (from slope)
			Double_t psi_in_sec0 = atan(slope_sec0);
			Double_t psi_in_sec2 = atan(slope_sec2);
			cout << "after atan" << endl;

			Double_t mw1_x_sec0 = slope_sec0*4.5 + offset_sec0;
			Double_t mw1_x_sec2 = slope_sec2*4.5 + offset_sec2;
			cout << "Entries mw3" << SofMwpc3CalData->GetEntriesFast() << endl;
			cout << "Entries tofsingletcal" << SofTofWSingleTcalData->GetEntriesFast() << endl;
			vector<vector<double> > tof_mw_values = after_glad_pos(SofMwpc3CalData,SofTofWSingleTcalData);
			cout << "mw3 positions" << tof_mw_values[0][0] << "   " << tof_mw_values[1][0] << endl;
			vector<double> pathlength_sec0 = path_tof_calc(psi_in_sec0,mw1_x_sec0,tof_mw_values[0][0],tof_mw_values[0][2],start_time);
			vector<double> pathlength_sec2 = path_tof_calc(psi_in_sec2,mw1_x_sec2,tof_mw_values[1][0],tof_mw_values[1][2],start_time);
			cout << "after pathlength computation" << endl;
			//now interpolate y-position on TOF from MW1_y and MW2_y
			sort(peaks_Y_mw2.begin(),peaks_Y_mw2.end());
			sort(peaks_Y_mw1.begin(),peaks_Y_mw1.end());
			Double_t slope_Y_sec0 = (peaks_Y_mw2[0]-peaks_Y_mw1[0])/(93.5-4.5);
			Double_t slope_Y_sec2 = (peaks_Y_mw2[1]-peaks_Y_mw1[1])/(93.5-4.5);
			Double_t offset_Y_sec0 = peaks_Y_mw2[0] - 93.5*slope_Y_sec0;
			Double_t offset_Y_sec2 = peaks_Y_mw2[1] - 93.5*slope_Y_sec2;
			Double_t y_extrap_sec0 = (pathlength_sec0[0]-73.75)*slope_Y_sec0+offset_Y_sec0; //offset of 73.75 because target position is at -73.75
			Double_t y_extrap_sec2 = (pathlength_sec2[0]-73.75)*slope_Y_sec2+offset_Y_sec2; 
			//now the final tof calibration
			//for section 0: tof_mw_values[0][1],y_extrap_sec0
			//for section 2: tof_mw_values[1][1],y_extrap_sec2
			cout << "after" << endl;

			delete  gr_sec0;
			delete  gr_sec2;
			delete fit_sec0;
			delete fit_sec2;	
				}
			//subsubcase 2X2Y-LU-RD: one left up && one right down
			if (peaks_X_mw1_down[0] < 0 && peaks_X_mw2_down[0] < 0 && peaks_X_mw1_up[0] > 0 && peaks_X_mw2_up[0] > 0 && e_sum_sec_1 > 0 && e_sum_sec_3 > 0 && e_sum_sec_0 == 0 && e_sum_sec_2 == 0){
				//cout << "subsubcase 2X2Y-LU-RD: one left up && one right down" << endl;
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
				y_sec1[j-3] = arr_xal_sec_1[j];
				y_sec3[j-3] = arr_xal_sec_3[j];
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
			//subsubcase 2X2Y-L: one left down && one left up
			if (peaks_X_mw1_down[0] > 0 && peaks_X_mw2_down[0] > 0 && peaks_X_mw1_up[0] > 0 && peaks_X_mw2_up[0] > 0 && e_sum_sec_0 > 0 && e_sum_sec_1 > 0 && e_sum_sec_2 == 0 && e_sum_sec_3 == 0 ){
				//cout << "subsubcase 2X2Y-L: one left down && one left up" << endl;
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
				y_sec0[j-3] = arr_xal_sec_0[j];
				y_sec1[j-3] = arr_xal_sec_1[j];
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

			//Get the deltax values
			 double arr_deltax_sec_0[16] = {-1000.};
			 double arr_deltax_sec_1[16] = {-1000.};
			for ( Int_t j = 0; j < 16; j++){
				arr_deltax_sec_0[j] = (slope_sec0*(pos_twim_entrance + j*twim_anodes_width) + offset_sec0) - arr_xal_sec_0[j];
				arr_deltax_sec_1[j] = (slope_sec1*(pos_twim_entrance + j*twim_anodes_width) + offset_sec1) - arr_xal_sec_1[j];
				h2_deltax_xal_sec0[j]->Fill(arr_xal_sec_0[j],arr_deltax_sec_0[j]);
				h2_deltax_xal_sec1[j]->Fill(arr_xal_sec_1[j],arr_deltax_sec_1[j]);
				}
			delete  gr_sec0;
			delete  gr_sec1;
			delete fit_sec0;
			delete fit_sec1;	
				}
			//subsubcase 2X2Y-R: one right down && one right up
			if (peaks_X_mw1_down[0] < 0 && peaks_X_mw2_down[0] < 0 && peaks_X_mw1_up[0] < 0 && peaks_X_mw2_up[0] < 0 && e_sum_sec_2 > 0 && e_sum_sec_3 > 0 && e_sum_sec_0 == 0 && e_sum_sec_1 == 0){
				//cout << "subsubcase 2X2Y-R: one right down && one right up" << endl;

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
				y_sec2[j-3] = arr_xal_sec_2[j];
				y_sec3[j-3] = arr_xal_sec_3[j];
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

			//Get the deltax values
			 double arr_deltax_sec_2[16] = {-1000.};
			 double arr_deltax_sec_3[16] = {-1000.};
			for ( Int_t j = 0; j < 16; j++){
				arr_deltax_sec_2[j] = (slope_sec2*(pos_twim_entrance + j*twim_anodes_width) + offset_sec2) - arr_xal_sec_2[j];
				arr_deltax_sec_3[j] = (slope_sec3*(pos_twim_entrance + j*twim_anodes_width) + offset_sec3) - arr_xal_sec_3[j];
				h2_deltax_xal_sec2[j]->Fill(arr_xal_sec_2[j],arr_deltax_sec_2[j]);
				h2_deltax_xal_sec3[j]->Fill(arr_xal_sec_3[j],arr_deltax_sec_3[j]);
				}
			delete  gr_sec2;
			delete  gr_sec3;
			delete fit_sec2;
			delete fit_sec3;	
				}

			}
		}

			
		delete [] softwimmappeddata;
		delete [] sofmwpc1caldata;
		delete [] sofmwpc2caldata; 
		delete [] sofscitcaldata;
		}
	}
char f_out_name[500];
sprintf(f_out_name,"/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/data/output_delta_x_final/foo_subrun_%s.root",input_str.c_str());
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
