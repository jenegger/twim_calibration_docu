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
//Assignment of the TWIM sections:
////0: Messel Down
////1: Messel Up
////2: Wixhausen Up
////3: Wixhausen Down
using namespace std;
void twim_sum_energy(const string& input_str){

Long64_t entries_twim = 0;
string fname = string("/scratch8/ge37liw/workingspace/exp_s455/data/unpacked/s455_03_273_") + input_str + string("_unpacked.root");
string cal_fname = string("/scratch8/ge37liw/workingspace/exp_s455/data/calibrated/s455_03_273_") + input_str + string("_calibrated.root");
const char* char_fname= fname.c_str();
const char* char_cal_fname = cal_fname.c_str();
TChain* chain = new TChain("evt");
chain->Reset();
chain->Add(char_fname);
TChain*  chain2 = new TChain("evt");
chain2->Reset();
chain2->Add(char_cal_fname);
chain->AddFriend(chain2);
chain->Print("all");
Long64_t nevents = chain->GetEntries();
cout << "total number of entries:\t" << nevents << endl;
cout << "total number of entries from calibrated:\t"<< chain2->GetEntries() << endl; 

char hist_name[500];
//HISTOS-------------------------
TH2D* h2_e_sum_tof_sec0;
sprintf(hist_name,"Fission Fragment  Messel side (section0) 238U");
h2_e_sum_tof_sec0 = new TH2D(hist_name,hist_name,1000,35.2,38,10000,0,50000);
h2_e_sum_tof_sec0->GetXaxis()->SetTitle("Time Start to TOFW [ns]");
h2_e_sum_tof_sec0->GetYaxis()->SetTitle("Energy");
h2_e_sum_tof_sec0->GetXaxis()->CenterTitle(true);
h2_e_sum_tof_sec0->GetYaxis()->CenterTitle(true);
h2_e_sum_tof_sec0->GetXaxis()->SetLabelSize(0.045);
h2_e_sum_tof_sec0->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_e_sum_tof_sec0_corr;
sprintf(hist_name,"Corrected Fission Fragment  Messel side (section0) 238U");
h2_e_sum_tof_sec0_corr = new TH2D(hist_name,hist_name,1000,35.2,38,10000,0,50000);
h2_e_sum_tof_sec0_corr->GetXaxis()->SetTitle("Time Start to TOFW [ns]");
h2_e_sum_tof_sec0_corr->GetYaxis()->SetTitle("Energy");
h2_e_sum_tof_sec0_corr->GetXaxis()->CenterTitle(true);
h2_e_sum_tof_sec0_corr->GetYaxis()->CenterTitle(true);
h2_e_sum_tof_sec0_corr->GetXaxis()->SetLabelSize(0.045);
h2_e_sum_tof_sec0_corr->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_e_sum_tof_sec1;
sprintf(hist_name,"Fission Fragment  Messel side (section1) 238U");
h2_e_sum_tof_sec1 = new TH2D(hist_name,hist_name,1000,35.2,38,10000,0,50000);
h2_e_sum_tof_sec1->GetXaxis()->SetTitle("Time Start to TOFW [ns]");
h2_e_sum_tof_sec1->GetYaxis()->SetTitle("Energy");
h2_e_sum_tof_sec1->GetXaxis()->CenterTitle(true);
h2_e_sum_tof_sec1->GetYaxis()->CenterTitle(true);
h2_e_sum_tof_sec1->GetXaxis()->SetLabelSize(0.045);
h2_e_sum_tof_sec1->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_e_sum_tof_sec2;
sprintf(hist_name,"Fission Fragment  Wixhausen side (section2) 238U");
h2_e_sum_tof_sec2 = new TH2D(hist_name,hist_name,1000,35.2,38,10000,0,50000);
h2_e_sum_tof_sec2->GetXaxis()->SetTitle("Time Start to TOFW [ns]");
h2_e_sum_tof_sec2->GetYaxis()->SetTitle("Energy");
h2_e_sum_tof_sec2->GetXaxis()->CenterTitle(true);
h2_e_sum_tof_sec2->GetYaxis()->CenterTitle(true);
h2_e_sum_tof_sec2->GetXaxis()->SetLabelSize(0.045);
h2_e_sum_tof_sec2->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_e_sum_tof_sec3;
sprintf(hist_name,"Fission Fragment  Wixhausen side (section3) 238U");
h2_e_sum_tof_sec3 = new TH2D(hist_name,hist_name,1000,35.2,38,10000,0,50000);
h2_e_sum_tof_sec3->GetXaxis()->SetTitle("Time Start to TOFW [ns]");
h2_e_sum_tof_sec3->GetYaxis()->SetTitle("Energy");
h2_e_sum_tof_sec3->GetXaxis()->CenterTitle(true);
h2_e_sum_tof_sec3->GetYaxis()->CenterTitle(true);
h2_e_sum_tof_sec3->GetXaxis()->SetLabelSize(0.045);
h2_e_sum_tof_sec3->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_e_sum_tof_sec2_corr;
sprintf(hist_name,"Corrected Fission Fragment  Wixhausen side (section2) 238U");
h2_e_sum_tof_sec2_corr = new TH2D(hist_name,hist_name,1000,35.2,38,10000,0,50000);
h2_e_sum_tof_sec2_corr->GetXaxis()->SetTitle("Time Start to TOFW [ns]");
h2_e_sum_tof_sec2_corr->GetYaxis()->SetTitle("Energy");
h2_e_sum_tof_sec2_corr->GetXaxis()->CenterTitle(true);
h2_e_sum_tof_sec2_corr->GetYaxis()->CenterTitle(true);
h2_e_sum_tof_sec2_corr->GetXaxis()->SetLabelSize(0.045);
h2_e_sum_tof_sec2_corr->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_e_sum_sec0;
sprintf(hist_name,"TWIM Energy Section 0 238U");
h1_e_sum_sec0 = new TH1D(hist_name,hist_name,10000,0,50000);
h1_e_sum_sec0->GetXaxis()->SetTitle("TWIM Energy");
h1_e_sum_sec0->GetYaxis()->SetTitle("Counts");
h1_e_sum_sec0->GetXaxis()->CenterTitle(true);
h1_e_sum_sec0->GetYaxis()->CenterTitle(true);
h1_e_sum_sec0->GetYaxis()->SetLabelSize(0.045);
h1_e_sum_sec0->GetYaxis()->SetTitleSize(0.045);


TH1D* h1_e_sum_sec1;
sprintf(hist_name,"TWIM Energy Section 1 238U");
h1_e_sum_sec1 = new TH1D(hist_name,hist_name,10000,0,50000);
h1_e_sum_sec1->GetXaxis()->SetTitle("TWIM Energy");
h1_e_sum_sec1->GetYaxis()->SetTitle("Counts");
h1_e_sum_sec1->GetXaxis()->CenterTitle(true);
h1_e_sum_sec1->GetYaxis()->CenterTitle(true);
h1_e_sum_sec1->GetYaxis()->SetLabelSize(0.045);
h1_e_sum_sec1->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_e_sum_sec2;
sprintf(hist_name,"TWIM Energy Section 2 238U");
h1_e_sum_sec2 = new TH1D(hist_name,hist_name,10000,0,50000);
h1_e_sum_sec2->GetXaxis()->SetTitle("TWIM Energy");
h1_e_sum_sec2->GetYaxis()->SetTitle("Counts");
h1_e_sum_sec2->GetXaxis()->CenterTitle(true);
h1_e_sum_sec2->GetYaxis()->CenterTitle(true);
h1_e_sum_sec2->GetYaxis()->SetLabelSize(0.045);
h1_e_sum_sec2->GetYaxis()->SetTitleSize(0.045);


TH1D* h1_e_sum_sec3;
sprintf(hist_name,"TWIM Energy Section 3 238U");
h1_e_sum_sec3 = new TH1D(hist_name,hist_name,10000,0,50000);
h1_e_sum_sec3->GetXaxis()->SetTitle("TWIM Energy");
h1_e_sum_sec3->GetYaxis()->SetTitle("Counts");
h1_e_sum_sec3->GetXaxis()->CenterTitle(true);
h1_e_sum_sec3->GetYaxis()->CenterTitle(true);
h1_e_sum_sec3->GetYaxis()->SetLabelSize(0.045);
h1_e_sum_sec3->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_e_sum_sec0_range;
sprintf(hist_name,"TWIM Energy Section 0 238U in range 36.6-36.9");
h1_e_sum_sec0_range = new TH1D(hist_name,hist_name,10000,0,50000);
h1_e_sum_sec0_range->GetXaxis()->SetTitle("TWIM Energy");
h1_e_sum_sec0_range->GetYaxis()->SetTitle("Counts");
h1_e_sum_sec0_range->GetXaxis()->CenterTitle(true);
h1_e_sum_sec0_range->GetYaxis()->CenterTitle(true);
h1_e_sum_sec0_range->GetYaxis()->SetLabelSize(0.045);
h1_e_sum_sec0_range->GetYaxis()->SetTitleSize(0.045);


TH1D* h1_e_sum_sec1_range;
sprintf(hist_name,"TWIM Energy Section 1 238U in range 36.6-36.9");
h1_e_sum_sec1_range = new TH1D(hist_name,hist_name,10000,0,50000);
h1_e_sum_sec1_range->GetXaxis()->SetTitle("TWIM Energy");
h1_e_sum_sec1_range->GetYaxis()->SetTitle("Counts");
h1_e_sum_sec1_range->GetXaxis()->CenterTitle(true);
h1_e_sum_sec1_range->GetYaxis()->CenterTitle(true);
h1_e_sum_sec1_range->GetYaxis()->SetLabelSize(0.045);
h1_e_sum_sec1_range->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_e_sum_sec2_range;
sprintf(hist_name,"TWIM Energy Section 2 238U in range 36.6-36.9");
h1_e_sum_sec2_range = new TH1D(hist_name,hist_name,10000,0,50000);
h1_e_sum_sec2_range->GetXaxis()->SetTitle("TWIM Energy");
h1_e_sum_sec2_range->GetYaxis()->SetTitle("Counts");
h1_e_sum_sec2_range->GetXaxis()->CenterTitle(true);
h1_e_sum_sec2_range->GetYaxis()->CenterTitle(true);
h1_e_sum_sec2_range->GetYaxis()->SetLabelSize(0.045);
h1_e_sum_sec2_range->GetYaxis()->SetTitleSize(0.045);


TH1D* h1_e_sum_sec3_range;
sprintf(hist_name,"TWIM Energy Section 3 238U in range 36.6-36.9");
h1_e_sum_sec3_range = new TH1D(hist_name,hist_name,10000,0,50000);
h1_e_sum_sec3_range->GetXaxis()->SetTitle("TWIM Energy");
h1_e_sum_sec3_range->GetYaxis()->SetTitle("Counts");
h1_e_sum_sec3_range->GetXaxis()->CenterTitle(true);
h1_e_sum_sec3_range->GetYaxis()->CenterTitle(true);
h1_e_sum_sec3_range->GetYaxis()->SetLabelSize(0.045);
h1_e_sum_sec3_range->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_y_diff_tof;
sprintf(hist_name,"y_tof section2 (RU) -y_tofw section 0 (LD)");
h1_y_diff_tof = new TH1D(hist_name,hist_name,500,-500,500);
h1_y_diff_tof->GetXaxis()->SetTitle("TWIM Energy");
h1_y_diff_tof->GetYaxis()->SetTitle("Counts");
h1_y_diff_tof->GetXaxis()->CenterTitle(true);
h1_y_diff_tof->GetYaxis()->CenterTitle(true);
h1_y_diff_tof->GetYaxis()->SetLabelSize(0.045);
h1_y_diff_tof->GetYaxis()->SetTitleSize(0.045);
//-------------------------------


R3BEventHeader* DataCA = new R3BEventHeader();
TBranch* branchData = chain->GetBranch("EventHeader.");
branchData->SetAddress(&DataCA);

TClonesArray* SofTwimMappedData = new TClonesArray("R3BSofTwimMappedData",2);
R3BSofTwimHitData** softwimmappeddata;
TBranch *branchSofTwimMappedData = chain->GetBranch("TwimMappedData");
branchSofTwimMappedData->SetAddress(&SofTwimMappedData);

TClonesArray* SofSciTcalData = new TClonesArray("R3BSofSciTcalData",2);
R3BSofSciTcalData** sofscitcaldata;
TBranch *branchSofSciTcalData = chain->GetBranch("SofSciTcalData");
branchSofSciTcalData->SetAddress(&SofSciTcalData);

TClonesArray* SofToFWHitData = new TClonesArray("R3BSofTofWHitData",2);
R3BSofTofWHitData** softofwhitdata;
TBranch *branchSofToFWHitData = chain->GetBranch("TofWHitData");
branchSofToFWHitData->SetAddress(&SofToFWHitData);


fstream fin;
fin.open("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration/parameters_twim_anodes.csv", ios::in); //FIXME: use right file...
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


for(Long64_t i=0;i< nevents;i++){
	Long64_t evtnr = i;
	if (i%100000==0)
		cout<<"Processing event for charge analysis "<<i<<endl;
	chain->GetEvent(i);
	entries_twim = SofTwimMappedData->GetEntries();
	if (entries_twim > 31){
		//cout << "this is event with more than one TWIM entry:\t" << evtnr << endl;
		//cout << "corresponding entries:\t" << entries_twim << endl;
		
		R3BSofTwimMappedData** softwimmappeddata  = new R3BSofTwimMappedData*[entries_twim];
	
		//create 4 arrays with 0 as entry (= 4 sections)
		double sec_0_arr[16] = { 0.};
		double sec_1_arr[16] = { 0.};
	    double sec_2_arr[16] = { 0.};
	    double sec_3_arr[16] = { 0.};	
		double e_sum_sec_0 = 0.;
		double e_sum_sec_1 = 0.;
		double e_sum_sec_2 = 0.;
		double e_sum_sec_3 = 0.;

		
		//fill the 4 arrays
		for (Int_t j = 0; j < entries_twim; j++){
			softwimmappeddata[j] = (R3BSofTwimMappedData*)SofTwimMappedData->At(j);
			Int_t twim_section = softwimmappeddata[j]->GetSecID();
			Int_t twim_anode = softwimmappeddata[j]->GetAnodeID();
	
			if (twim_section == 0 && twim_anode < 16){
				sec_0_arr[twim_anode] = softwimmappeddata[j]->GetEnergy();
				if (softwimmappeddata[j]->GetEnergy() != 0 && twim_anode != 13){  // section0, anode 13 is bad, leave out ...
					e_sum_sec_0 += (softwimmappeddata[j]->GetEnergy())*(1/v_para_twim[0][twim_anode][0]) - (v_para_twim[0][twim_anode][1]/v_para_twim[0][twim_anode][0]);
					}
				}
			else if (twim_section == 1 && twim_anode < 16){
				sec_1_arr[twim_anode] = softwimmappeddata[j]->GetEnergy();
				if (softwimmappeddata[j]->GetEnergy() != 0){
					e_sum_sec_1 += (softwimmappeddata[j]->GetEnergy())*(1/v_para_twim[1][twim_anode][0]) - (v_para_twim[1][twim_anode][1]/v_para_twim[1][twim_anode][0]);
					}
				}
			else if (twim_section == 2 && twim_anode < 16){
				sec_2_arr[twim_anode] = softwimmappeddata[j]->GetEnergy();
				if (softwimmappeddata[j]->GetEnergy() != 0){
					e_sum_sec_2 += (softwimmappeddata[j]->GetEnergy())*(1/v_para_twim[2][twim_anode][0]) - (v_para_twim[2][twim_anode][1]/v_para_twim[2][twim_anode][0]);
					}
				}
			else if (twim_section == 3 && twim_anode < 16){
				sec_3_arr[twim_anode] = softwimmappeddata[j]->GetEnergy();
				if (softwimmappeddata[j]->GetEnergy() != 0){
					e_sum_sec_3 += (softwimmappeddata[j]->GetEnergy())*(1/v_para_twim[3][twim_anode][0]) - (v_para_twim[3][twim_anode][1]/v_para_twim[3][twim_anode][0]);
					}
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
		
		//LD-RU (section 0 and 2)
		if ((min_sec0 > 1) && (min_sec2 > 1) && (e_sum_sec_1 == 0) && (e_sum_sec_3 == 0)){
			Int_t nr_start = SofSciTcalData->GetEntries();
			Int_t nr_tof = SofToFWHitData->GetEntries();
			if (nr_tof == 2  && nr_start > 1){
				R3BSofTofWHitData** softofwhitdata = new R3BSofTofWHitData*[nr_tof];
				R3BSofSciTcalData** sofscitcaldata = new R3BSofSciTcalData*[nr_start];
				vector<double> v_start_times;
				vector<double> v_tof_x;
				vector<double> v_tof_y;
				vector<double> v_tof_times;
			  for (Int_t i = 0; i < 2; i++){
					sofscitcaldata[i] = (R3BSofSciTcalData*)SofSciTcalData->At(i);
			  		softofwhitdata[i] = (R3BSofTofWHitData*)SofToFWHitData->At(i);
					v_start_times.push_back(sofscitcaldata[i]->GetRawTimeNs());
					v_tof_x.push_back(softofwhitdata[i]->GetX());
					v_tof_y.push_back(softofwhitdata[i]->GetY());
					v_tof_times.push_back(softofwhitdata[i]->GetTof());	
			  	}
				Double_t start_time = 0.5*(v_start_times[0]+v_start_times[1]);
				Double_t left_time; //Messel
				Double_t right_time; //Wixhausen
				Double_t left_y;
				Double_t right_y;
				cout << v_tof_times[0] << "   " << v_tof_times[1] << endl;
				if (v_tof_x[0] < v_tof_x[1]){
					right_time = v_tof_times[0];
					left_y = v_tof_y[0];
					left_time = v_tof_times[1];
					right_y = v_tof_y[1];
					}
				if(v_tof_x[0] > v_tof_x[1]){
					right_time = v_tof_times[1];
					left_y = v_tof_y[1];
					left_time = v_tof_times[0];
					right_y = v_tof_y[0];
					}
				if (left_time > 35.2 && left_time < 38 && right_time > 35.2 && right_time < 38){
					//now fill histos with sum energy vs time and the 1D histo for the energy					
					h2_e_sum_tof_sec0->Fill(left_time,e_sum_sec_0/16.);
					h2_e_sum_tof_sec2->Fill(right_time,e_sum_sec_2/16.);
					h1_y_diff_tof->Fill(right_y-left_y);
					}
				if (left_time > 36.9 && left_time < 37){
					h1_e_sum_sec0->Fill(e_sum_sec_0/16.);
					}
				if (right_time > 36.9 && right_time < 37){
					h1_e_sum_sec2->Fill(e_sum_sec_2/16.);
					}
				if (left_time > 36.6 && left_time < 36.9){
					h1_e_sum_sec0_range->Fill(e_sum_sec_0/16.);
					}
				if (right_time > 36.6 && right_time < 36.9){
					h1_e_sum_sec2_range->Fill(e_sum_sec_2/16.);
					}
				//corrected version as I feel unsecure if it is correctly reconstructed...
				Double_t corr_left_time;
				Double_t corr_right_time;
				Double_t corr_left_y;
				Double_t corr_right_y;
				if (v_tof_y[0] < v_tof_y[1]){
					corr_left_time = v_tof_times[0];
					corr_right_time = v_tof_times[1];
					}
				if (v_tof_y[0] > v_tof_y[1]){
					corr_left_time = v_tof_times[1];
					corr_right_time = v_tof_times[0];
					}
				if (corr_left_time > 35.2 && corr_left_time < 38 && corr_right_time > 35.2 && corr_right_time < 38){
					h2_e_sum_tof_sec0_corr->Fill(left_time,e_sum_sec_0/16.);
					h2_e_sum_tof_sec2_corr->Fill(right_time,e_sum_sec_2/16.);
					}
				}

			}
		//LD-RD (section 0 and 3)
		if ((min_sec0 > 1) && (min_sec3 > 1) && (e_sum_sec_1 == 0) && (e_sum_sec_2 == 0)){
			Int_t nr_start = SofSciTcalData->GetEntries();
			Int_t nr_tof = SofToFWHitData->GetEntries();
			if (nr_tof == 2  && nr_start > 1){
				R3BSofTofWHitData** softofwhitdata = new R3BSofTofWHitData*[nr_tof];
				R3BSofSciTcalData** sofscitcaldata = new R3BSofSciTcalData*[nr_start];
				vector<double> v_start_times;
				vector<double> v_tof_x;
				vector<double> v_tof_times;
			  for (Int_t i = 0; i < 2; i++){
					sofscitcaldata[i] = (R3BSofSciTcalData*)SofSciTcalData->At(i);
			  		softofwhitdata[i] = (R3BSofTofWHitData*)SofToFWHitData->At(i);
					v_start_times.push_back(sofscitcaldata[i]->GetRawTimeNs());
					v_tof_x.push_back(softofwhitdata[i]->GetX());
					v_tof_times.push_back(softofwhitdata[i]->GetTof());	
			  	}
				Double_t start_time = 0.5*(v_start_times[0]+v_start_times[1]);
				Double_t left_time; //Messel
				Double_t right_time; //Wixhausen
				//cout << v_tof_times[0] << "   " << v_tof_times[1] << endl;
				//cout << e_sum_sec_0 << "   " << e_sum_sec_3 << endl;
				if (v_tof_x[0] < v_tof_x[1]){
					right_time = v_tof_times[0];
					left_time = v_tof_times[1];
					}
				if(v_tof_x[0] > v_tof_x[1]){
					right_time = v_tof_times[1];
					left_time = v_tof_times[0];
					}
				if (left_time > 35.2 && left_time < 38 && right_time > 35.2 && right_time < 38){
					//now fill histos with sum energy vs time and the 1D histo for the energy					
					h2_e_sum_tof_sec0->Fill(left_time,e_sum_sec_0/16.);
					h2_e_sum_tof_sec3->Fill(right_time,e_sum_sec_3/16.);
					}
				if (left_time > 36.9 && left_time < 37){
					h1_e_sum_sec0->Fill(e_sum_sec_0/16.);		
					}
				if (right_time > 36.9 && right_time < 37){
					h1_e_sum_sec3->Fill(e_sum_sec_3/16.);
					}
				if (left_time > 36.6 && left_time < 36.9){
					h1_e_sum_sec0_range->Fill(e_sum_sec_0/16.);
					}
				if (right_time > 36.6 && right_time < 36.9){
					h1_e_sum_sec3_range->Fill(e_sum_sec_3/16.);
					}
				}
			}
		//LU-RU (section 1 and 2)
		if ((min_sec1 > 1) && (min_sec2 > 1) && (e_sum_sec_0 == 0) && (e_sum_sec_3 == 0)){
			Int_t nr_start = SofSciTcalData->GetEntries();
			Int_t nr_tof = SofToFWHitData->GetEntries();
			if (nr_tof == 2  && nr_start > 1){
				R3BSofTofWHitData** softofwhitdata = new R3BSofTofWHitData*[nr_tof];
				R3BSofSciTcalData** sofscitcaldata = new R3BSofSciTcalData*[nr_start];
				vector<double> v_start_times;
				vector<double> v_tof_x;
				vector<double> v_tof_times;
			  for (Int_t i = 0; i < 2; i++){
					sofscitcaldata[i] = (R3BSofSciTcalData*)SofSciTcalData->At(i);
			  		softofwhitdata[i] = (R3BSofTofWHitData*)SofToFWHitData->At(i);
					v_start_times.push_back(sofscitcaldata[i]->GetRawTimeNs());
					v_tof_x.push_back(softofwhitdata[i]->GetX());
					v_tof_times.push_back(softofwhitdata[i]->GetTof());	
			  	}
				Double_t start_time = 0.5*(v_start_times[0]+v_start_times[1]);
				Double_t left_time; //Messel
				Double_t right_time; //Wixhausen
				if (v_tof_x[0] < v_tof_x[1]){
					right_time = v_tof_times[0];
					left_time = v_tof_times[1];
					}
				if(v_tof_x[0] > v_tof_x[1]){
					right_time = v_tof_times[1];
					left_time = v_tof_times[0];
					}
				if (left_time > 35.2 && left_time < 38 && right_time > 35.2 && right_time < 38){
					//now fill histos with sum energy vs time and the 1D histo for the energy					
					h2_e_sum_tof_sec1->Fill(left_time,e_sum_sec_1/16.);	
					h2_e_sum_tof_sec2->Fill(right_time,e_sum_sec_2/16.);
					}
				if (left_time > 36.9 && left_time < 37){
					h1_e_sum_sec1->Fill(e_sum_sec_1/16.);
					}
				if (right_time > 36.9 && right_time < 37){
					h1_e_sum_sec2->Fill(e_sum_sec_2/16.);
					}
				if (left_time > 36.6 && left_time < 36.9){
					h1_e_sum_sec1_range->Fill(e_sum_sec_1/16.);
					}
				if (right_time > 36.6 && right_time < 36.9){
					h1_e_sum_sec2_range->Fill(e_sum_sec_2/16.);
					}
				}
			}
		//LU-RD (section 1 and 3)
		if ((min_sec1 > 1) && (min_sec3 > 1) && (e_sum_sec_0 == 0) && (e_sum_sec_2 == 0)){
			Int_t nr_start = SofSciTcalData->GetEntries();
			Int_t nr_tof = SofToFWHitData->GetEntries();
			if (nr_tof == 2  && nr_start > 1){
				R3BSofTofWHitData** softofwhitdata = new R3BSofTofWHitData*[nr_tof];
				R3BSofSciTcalData** sofscitcaldata = new R3BSofSciTcalData*[nr_start];
				vector<double> v_start_times;
				vector<double> v_tof_x;
				vector<double> v_tof_times;
			  for (Int_t i = 0; i < 2; i++){
					sofscitcaldata[i] = (R3BSofSciTcalData*)SofSciTcalData->At(i);
			  		softofwhitdata[i] = (R3BSofTofWHitData*)SofToFWHitData->At(i);
					v_start_times.push_back(sofscitcaldata[i]->GetRawTimeNs());
					v_tof_x.push_back(softofwhitdata[i]->GetX());
					v_tof_times.push_back(softofwhitdata[i]->GetTof());	
			  	}
				Double_t start_time = 0.5*(v_start_times[0]+v_start_times[1]);
				Double_t left_time; //Messel
				Double_t right_time; //Wixhausen
				if (v_tof_x[0] < v_tof_x[1]){
					right_time = v_tof_times[0];
					left_time = v_tof_times[1];
					}
				if(v_tof_x[0] > v_tof_x[1]){
					right_time = v_tof_times[1];
					left_time = v_tof_times[0];
					}
				if (left_time > 35.2 && left_time < 38 && right_time > 35.2 && right_time < 38){
					//now fill histos with sum energy vs time and the 1D histo for the energy					
					h2_e_sum_tof_sec1->Fill(left_time,e_sum_sec_1/16.);
					h2_e_sum_tof_sec3->Fill(right_time,e_sum_sec_3/16.);
					}
				if (left_time > 36.9 && left_time < 37){
					h1_e_sum_sec1->Fill(e_sum_sec_1/16.);
					}
				if (right_time > 36.9 && right_time < 37){
					h1_e_sum_sec3->Fill(e_sum_sec_3/16.);
					}
				if (left_time > 36.6 && left_time < 36.9){
					h1_e_sum_sec1_range->Fill(e_sum_sec_1/16.);
					}
				if (right_time > 36.6 && right_time < 36.9){
					h1_e_sum_sec3_range->Fill(e_sum_sec_3/16.);
					}
				}
			}
		}

//Select events with two TWIM hits: LD-RU, LD-RD,LU-RU,LU-RD

//Then check that all 16 channels of the selected sections have signal

//Then limit to specific ToF range

	}
char f_out_name[500];
sprintf(f_out_name,"test_esum_plots_subrun_%s.root",input_str.c_str());
TFile * f = new TFile(f_out_name,"RECREATE");
TList *l = new TList();
l->Add(h2_e_sum_tof_sec0);
l->Add(h2_e_sum_tof_sec1);
l->Add(h2_e_sum_tof_sec2);
l->Add(h2_e_sum_tof_sec3);
l->Add(h1_e_sum_sec0);
l->Add(h1_e_sum_sec1);
l->Add(h1_e_sum_sec2);
l->Add(h1_e_sum_sec3);
l->Add(h1_e_sum_sec0_range);
l->Add(h1_e_sum_sec1_range);
l->Add(h1_e_sum_sec2_range);
l->Add(h1_e_sum_sec3_range);
l->Add(h1_y_diff_tof);
l->Add(h2_e_sum_tof_sec0_corr);
l->Add(h2_e_sum_tof_sec2_corr);
l->Write("histlist", TObject::kSingleKey);
}

