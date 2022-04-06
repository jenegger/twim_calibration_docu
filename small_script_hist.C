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
#include <stdlib.h>

using namespace std;
void small_script_hist(const string& input_str){

char hist_name[500];
//histograms for section 0
TH2D* h2_twim_sec0_energy_anode[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name, "Section0: Energy AnodeID %i vs Energy ref. Anode 5",i);
h2_twim_sec0_energy_anode[i] = new TH2D(hist_name,hist_name,1000,0,50000,1000,0,50000);
h2_twim_sec0_energy_anode[i]->GetXaxis()->SetTitle("Energy AnodeID 5");
h2_twim_sec0_energy_anode[i]->GetYaxis()->SetTitle("Energy AnodeID %i");
h2_twim_sec0_energy_anode[i]->GetXaxis()->CenterTitle(true);
h2_twim_sec0_energy_anode[i]->GetYaxis()->CenterTitle(true);
h2_twim_sec0_energy_anode[i]->GetYaxis()->SetLabelSize(0.045);
h2_twim_sec0_energy_anode[i]->GetYaxis()->SetTitleSize(0.045);
}
//histograms for section 1
TH2D* h2_twim_sec1_energy_anode[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name, "Section1: Energy AnodeID %i vs Energy ref. Anode 5",i);
h2_twim_sec1_energy_anode[i] = new TH2D(hist_name,hist_name,1000,0,50000,1000,0,50000);
h2_twim_sec1_energy_anode[i]->GetXaxis()->SetTitle("Energy AnodeID 5");
h2_twim_sec1_energy_anode[i]->GetYaxis()->SetTitle("Energy AnodeID %i");
h2_twim_sec1_energy_anode[i]->GetXaxis()->CenterTitle(true);
h2_twim_sec1_energy_anode[i]->GetYaxis()->CenterTitle(true);
h2_twim_sec1_energy_anode[i]->GetYaxis()->SetLabelSize(0.045);
h2_twim_sec1_energy_anode[i]->GetYaxis()->SetTitleSize(0.045);
}
//histograms for section 2
TH2D* h2_twim_sec2_energy_anode[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name, "Section2: Energy AnodeID %i vs Energy ref. Anode 5",i);
h2_twim_sec2_energy_anode[i] = new TH2D(hist_name,hist_name,1000,0,50000,1000,0,50000);
h2_twim_sec2_energy_anode[i]->GetXaxis()->SetTitle("Energy AnodeID 5");
h2_twim_sec2_energy_anode[i]->GetYaxis()->SetTitle("Energy AnodeID %i");
h2_twim_sec2_energy_anode[i]->GetXaxis()->CenterTitle(true);
h2_twim_sec2_energy_anode[i]->GetYaxis()->CenterTitle(true);
h2_twim_sec2_energy_anode[i]->GetYaxis()->SetLabelSize(0.045);
h2_twim_sec2_energy_anode[i]->GetYaxis()->SetTitleSize(0.045);
}
//histograms for section 3
TH2D* h2_twim_sec3_energy_anode[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name, "Section3: Energy AnodeID %i vs Energy ref. Anode 5",i);
h2_twim_sec3_energy_anode[i] = new TH2D(hist_name,hist_name,1000,0,50000,1000,0,50000);
h2_twim_sec3_energy_anode[i]->GetXaxis()->SetTitle("Energy AnodeID 5");
h2_twim_sec3_energy_anode[i]->GetYaxis()->SetTitle("Energy AnodeID %i");
h2_twim_sec3_energy_anode[i]->GetXaxis()->CenterTitle(true);
h2_twim_sec3_energy_anode[i]->GetYaxis()->CenterTitle(true);
h2_twim_sec3_energy_anode[i]->GetYaxis()->SetLabelSize(0.045);
h2_twim_sec3_energy_anode[i]->GetYaxis()->SetTitleSize(0.045);
}

Long64_t entries_califa = 0;
string fname = string("/scratch8/ge37liw/workingspace/exp_s455/data/unpacked/s455_03_273_") + input_str + string("_unpacked.root");
const char* char_fname= fname.c_str();
TChain* chain = new TChain("evt");
chain->Reset();
chain->Add(char_fname);
Long64_t nevents = chain->GetEntries();
cout << "total number of entries:\t" << nevents << endl;

R3BEventHeader* DataCA = new R3BEventHeader();
TBranch* branchData = chain->GetBranch("EventHeader.");
branchData->SetAddress(&DataCA);

TClonesArray* SofTwimMappedData = new TClonesArray("R3BSofTwimMappedData",2);
R3BSofTwimHitData** softwimmappeddata;
TBranch *branchSofTwimMappedData = chain->GetBranch("TwimMappedData");
branchSofTwimMappedData->SetAddress(&SofTwimMappedData);

for(Long64_t i = 0; i < nevents; i++){
    Long64_t evtnr = i;
    if (i%100000==0)
        cout<<"Processing event for charge analysis "<<i<<endl;
    chain->GetEvent(i);
	
	Long_t entries_twim = SofTwimMappedData->GetEntriesFast();
//	if (entries_twim > 0){
//		cout << "entries of twim:\t " << entries_twim << "  for eventnr:\t" << evtnr << endl;
//	}
	R3BSofTwimMappedData** softwimmappeddata  = new R3BSofTwimMappedData*[entries_twim];

	//create 4 arrays with 0 as entry (= 4 sections)
	double sec_0_arr[16] = { 0.};
	double sec_1_arr[16] = { 0.};
    double sec_2_arr[16] = { 0.};
    double sec_3_arr[16] = { 0.};	
	
	//fill the 4 arrays
	for (Int_t j = 0; j < entries_twim; j++){
		softwimmappeddata[j] = (R3BSofTwimMappedData*)SofTwimMappedData->At(j);
		Int_t twim_section = softwimmappeddata[j]->GetSecID();
		Int_t twim_anode = softwimmappeddata[j]->GetAnodeID();

		if (twim_section == 0 && twim_anode < 16){
			sec_0_arr[twim_anode] = softwimmappeddata[j]->GetEnergy();
			}
		else if (twim_section == 1 && twim_anode < 16){
			sec_1_arr[twim_anode] = softwimmappeddata[j]->GetEnergy();
			}
		else if (twim_section == 2 && twim_anode < 16){
			sec_2_arr[twim_anode] = softwimmappeddata[j]->GetEnergy();
			}
		else if (twim_section == 3 && twim_anode < 16){
			sec_3_arr[twim_anode] = softwimmappeddata[j]->GetEnergy();
			}		
//		else {
//			cout << "WARNING, no TWIM Section selected!" <<  "for eventnr:\t" << evtnr << endl;
//			cout << "anodeID:\t" << twim_anode << "  and energy" << softwimmappeddata[j]->GetEnergy() << endl;
//			cout << "SECTION is:\t" << twim_section << "END----------" << endl;
//			}



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
		
	if (min_sec0 > 1){
		for (Int_t j = 0; j < 16;j++){
			h2_twim_sec0_energy_anode[j]->Fill(sec_0_arr[5],sec_0_arr[j]);				
			}
		}		

	if (min_sec1 > 1){
		for (Int_t j = 0; j < 16;j++){
			h2_twim_sec1_energy_anode[j]->Fill(sec_1_arr[5],sec_1_arr[j]);				
			}
		}		

	if (min_sec2 > 1){
		for (Int_t j = 0; j < 16;j++){
			h2_twim_sec2_energy_anode[j]->Fill(sec_2_arr[5],sec_2_arr[j]);				
			}
		}		

	if (min_sec3 > 1){
		for (Int_t j = 0; j < 16;j++){
			h2_twim_sec3_energy_anode[j]->Fill(sec_3_arr[5],sec_3_arr[j]);				
			}
		}		
}

char f_out_name[500];
sprintf(f_out_name,"section_plots_subrun_%s.root",input_str.c_str());
TFile * f = new TFile(f_out_name,"RECREATE");
//fill histos in TList
TList *l = new TList();
gStyle->SetOptFit(1);
//section 0
for(Int_t i = 0; i < 16; i++){
	l->Add(h2_twim_sec0_energy_anode[i]);
}
//section 1
for(Int_t i = 0; i < 16; i++){
    l->Add(h2_twim_sec1_energy_anode[i]);
}
//section 2
for(Int_t i = 0; i < 16; i++){
    l->Add(h2_twim_sec2_energy_anode[i]);
}
//section 3
for(Int_t i = 0; i < 16; i++){
    l->Add(h2_twim_sec3_energy_anode[i]);
}
l->Write("histlist", TObject::kSingleKey);
}
