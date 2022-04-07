#include <iostream>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include <math.h>
#include "TH1F.h"
#include "TH2F.h"
using namespace std;
char hist_name[500];
char fname[500];
char f_out_name[500];

void calibrate_mw(char const count_i[50]){
Int_t entries_mw0 = 0.;
Int_t entries_mw1 = 0.;
Int_t entries_mw2 = 0.;
Int_t entries_mw3 = 0.;

TH1F* h1_mw0_x;

sprintf(hist_name, "MWPC0 f.X");
h1_mw0_x = new TH1F(hist_name,hist_name,1200,-300,300);
h1_mw0_x->GetXaxis()->SetTitle("MWPC0.fX [mm]");
h1_mw0_x->GetYaxis()->SetTitle("Counts");
h1_mw0_x->GetXaxis()->CenterTitle(true);
h1_mw0_x->GetYaxis()->CenterTitle(true);
h1_mw0_x->GetYaxis()->SetLabelSize(0.045);
h1_mw0_x->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_mw0_y;

sprintf(hist_name, "MWPC0 f.Y");
h1_mw0_y = new TH1F(hist_name,hist_name,1200,-300,300);
h1_mw0_y->GetXaxis()->SetTitle("MWPC0.fY [mm]");
h1_mw0_y->GetYaxis()->SetTitle("Counts");
h1_mw0_y->GetXaxis()->CenterTitle(true);
h1_mw0_y->GetYaxis()->CenterTitle(true);
h1_mw0_y->GetYaxis()->SetLabelSize(0.045);
h1_mw0_y->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_mw1_x;

sprintf(hist_name, "MWPC1 f.X");
h1_mw1_x = new TH1F(hist_name,hist_name,1200,-300,300);
h1_mw1_x->GetXaxis()->SetTitle("MWPC1.fX [mm]");
h1_mw1_x->GetYaxis()->SetTitle("Counts");
h1_mw1_x->GetXaxis()->CenterTitle(true);
h1_mw1_x->GetYaxis()->CenterTitle(true);
h1_mw1_x->GetYaxis()->SetLabelSize(0.045);
h1_mw1_x->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_mw1_y;

sprintf(hist_name, "MWPC1 f.Y");
h1_mw1_y = new TH1F(hist_name,hist_name,1200,-300,300);
h1_mw1_y->GetXaxis()->SetTitle("MWPC1.fY [mm]");
h1_mw1_y->GetYaxis()->SetTitle("Counts");
h1_mw1_y->GetXaxis()->CenterTitle(true);
h1_mw1_y->GetYaxis()->CenterTitle(true);
h1_mw1_y->GetYaxis()->SetLabelSize(0.045);
h1_mw1_y->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_mw2_x;

sprintf(hist_name, "MWPC2 f.X");
h1_mw2_x = new TH1F(hist_name,hist_name,1200,-300,300);
h1_mw2_x->GetXaxis()->SetTitle("MWPC2.fX [mm]");
h1_mw2_x->GetYaxis()->SetTitle("Counts");
h1_mw2_x->GetXaxis()->CenterTitle(true);
h1_mw2_x->GetYaxis()->CenterTitle(true);
h1_mw2_x->GetYaxis()->SetLabelSize(0.045);
h1_mw2_x->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_mw2_y;

sprintf(hist_name, "MWPC2 f.Y");
h1_mw2_y = new TH1F(hist_name,hist_name,1200,-300,300);
h1_mw2_y->GetXaxis()->SetTitle("MWPC2.fY [mm]");
h1_mw2_y->GetYaxis()->SetTitle("Counts");
h1_mw2_y->GetXaxis()->CenterTitle(true);
h1_mw2_y->GetYaxis()->CenterTitle(true);
h1_mw2_y->GetYaxis()->SetLabelSize(0.045);
h1_mw2_y->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_mw3_x;

sprintf(hist_name, "MWPC3 f.X");
h1_mw3_x = new TH1F(hist_name,hist_name,2000,-500,500);
h1_mw3_x->GetXaxis()->SetTitle("MWPC3.fX [mm]");
h1_mw3_x->GetYaxis()->SetTitle("Counts");
h1_mw3_x->GetXaxis()->CenterTitle(true);
h1_mw3_x->GetYaxis()->CenterTitle(true);
h1_mw3_x->GetYaxis()->SetLabelSize(0.045);
h1_mw3_x->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_mw3_y;

sprintf(hist_name, "MWPC3 f.Y");
h1_mw3_y = new TH1F(hist_name,hist_name,2000,-500,500);
h1_mw3_y->GetXaxis()->SetTitle("MWPC3.fY [mm]");
h1_mw3_y->GetYaxis()->SetTitle("Counts");
h1_mw3_y->GetXaxis()->CenterTitle(true);
h1_mw3_y->GetYaxis()->CenterTitle(true);
h1_mw3_y->GetYaxis()->SetLabelSize(0.045);
h1_mw3_y->GetYaxis()->SetTitleSize(0.045);

sprintf(fname,"/scratch8/ge37liw/workingspace/data/root_files/all_ts_unpack/sweep_target/ts_cone_cluster_%s_zero_cm_y_corr.root",count_i);
sprintf(f_out_name,"/home/ge37liw/plots/histos/zero_with_miss/fast_macro/calibration_mwpcs_x_y_%s.root", count_i);
TFile * f = new TFile(f_out_name,"RECREATE");

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

Long64_t nevents = chain->GetEntries();
for(Long64_t i=0;i< nevents;i++){
    Long64_t evtnr = i;
    if (i%100000==0)
        cout<<"Processing event "<<i<<endl;
	chain->GetEvent(i);
	entries_mw0 = SofMwpc0HitData->GetEntries();
    entries_mw1 = SofMwpc1HitData->GetEntries();
    entries_mw2 = SofMwpc2HitData->GetEntries();
    entries_mw3 = SofMwpc3HitData->GetEntries();
	
	if (entries_mw0 == 1 && entries_mw1 == 1 && entries_mw2 == 1 && entries_mw3 == 1){
		sofmwpc0hitdata = new R3BSofMwpcHitData*[1];
        sofmwpc1hitdata = new R3BSofMwpcHitData*[1];
        sofmwpc2hitdata = new R3BSofMwpcHitData*[1];
        sofmwpc3hitdata = new R3BSofMwpcHitData*[1];
		sofmwpc0hitdata[0] = (R3BSofMwpcHitData*)SofMwpc0HitData->At(0);
        sofmwpc1hitdata[0] = (R3BSofMwpcHitData*)SofMwpc1HitData->At(0);
        sofmwpc2hitdata[0] = (R3BSofMwpcHitData*)SofMwpc2HitData->At(0);
        sofmwpc3hitdata[0] = (R3BSofMwpcHitData*)SofMwpc3HitData->At(0);
		Double_t xMW0 = sofmwpc0hitdata[0]->GetX();
        Double_t xMW1 = sofmwpc1hitdata[0]->GetX();
        Double_t xMW2 = sofmwpc2hitdata[0]->GetX();
        Double_t xMW3 = sofmwpc3hitdata[0]->GetX();
        Double_t yMW0 = sofmwpc0hitdata[0]->GetY();
        Double_t yMW1 = sofmwpc1hitdata[0]->GetY();
        Double_t yMW2 = sofmwpc2hitdata[0]->GetY();
        Double_t yMW3 = sofmwpc3hitdata[0]->GetY();
		
		//fill histos
		h1_mw0_x->Fill(xMW0);
		h1_mw0_y->Fill(yMW0);
		h1_mw1_x->Fill(xMW1);
		h1_mw1_y->Fill(yMW1);
		h1_mw2_x->Fill(xMW2);
		h1_mw2_y->Fill(yMW2);
		h1_mw3_x->Fill(xMW3);
		h1_mw3_y->Fill(yMW3);
		
		}

}
TList *l = new TList();
l->Add(h1_mw0_x);
l->Add(h1_mw0_y);
l->Add(h1_mw1_x);
l->Add(h1_mw1_y);
l->Add(h1_mw2_x);
l->Add(h1_mw2_y);
l->Add(h1_mw3_x);
l->Add(h1_mw3_y);
l->Write("histlist", TObject::kSingleKey);
cout << "end of process, successful" <<endl;
}
