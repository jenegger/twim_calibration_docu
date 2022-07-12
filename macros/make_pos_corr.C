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
#include <iostream>
#include <fstream>


#include <stdlib.h>
using namespace std;
void make_pos_corr(const string& input_str){
string fname = string(input_str);
char hist_name[500];
TH2D* h2_energy_vs_pos_messel_tspline_sec0;
sprintf(hist_name, "Position on 11th anode in music vs Calibrated Energy messel with tcut,section0");
h2_energy_vs_pos_messel_tspline_sec0 = new TH2D(hist_name,hist_name,600,0.,120,950,0,450000);
h2_energy_vs_pos_messel_tspline_sec0->GetXaxis()->SetTitle("pos[mm]");
h2_energy_vs_pos_messel_tspline_sec0->GetYaxis()->SetTitle("Energy");
h2_energy_vs_pos_messel_tspline_sec0->GetXaxis()->CenterTitle(true);
h2_energy_vs_pos_messel_tspline_sec0->GetYaxis()->CenterTitle(true);
h2_energy_vs_pos_messel_tspline_sec0->GetYaxis()->SetLabelSize(0.045);
h2_energy_vs_pos_messel_tspline_sec0->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_energy_vs_pos_messel_tspline_sec1;
sprintf(hist_name, "Position on 11th anode in music vs Calibrated Energy messel with tcut,section1");
h2_energy_vs_pos_messel_tspline_sec1 = new TH2D(hist_name,hist_name,600,0.,120,950,0,450000);
h2_energy_vs_pos_messel_tspline_sec1->GetXaxis()->SetTitle("pos[mm]");
h2_energy_vs_pos_messel_tspline_sec1->GetYaxis()->SetTitle("Energy");
h2_energy_vs_pos_messel_tspline_sec1->GetXaxis()->CenterTitle(true);
h2_energy_vs_pos_messel_tspline_sec1->GetYaxis()->CenterTitle(true);
h2_energy_vs_pos_messel_tspline_sec1->GetYaxis()->SetLabelSize(0.045);
h2_energy_vs_pos_messel_tspline_sec1->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_energy_vs_pos_wix_tspline_sec2;
sprintf(hist_name, "Position on 11th anode in music vs Calibrated Energy Wixhausen with tcut,section2");
h2_energy_vs_pos_wix_tspline_sec2 = new TH2D(hist_name,hist_name,600,-120,0.,950,0,450000);
h2_energy_vs_pos_wix_tspline_sec2->GetXaxis()->SetTitle("pos[mm]");
h2_energy_vs_pos_wix_tspline_sec2->GetYaxis()->SetTitle("Energy");
h2_energy_vs_pos_wix_tspline_sec2->GetXaxis()->CenterTitle(true);
h2_energy_vs_pos_wix_tspline_sec2->GetYaxis()->CenterTitle(true);
h2_energy_vs_pos_wix_tspline_sec2->GetYaxis()->SetLabelSize(0.045);
h2_energy_vs_pos_wix_tspline_sec2->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_energy_vs_pos_wix_tspline_sec3;
sprintf(hist_name, "Position on 11th anode in music vs Calibrated Energy Wixhausen with tcut,section3");
h2_energy_vs_pos_wix_tspline_sec3 = new TH2D(hist_name,hist_name,600,-120,0.,950,0,450000);
h2_energy_vs_pos_wix_tspline_sec3->GetXaxis()->SetTitle("pos[mm]");
h2_energy_vs_pos_wix_tspline_sec3->GetYaxis()->SetTitle("Energy");
h2_energy_vs_pos_wix_tspline_sec3->GetXaxis()->CenterTitle(true);
h2_energy_vs_pos_wix_tspline_sec3->GetYaxis()->CenterTitle(true);
h2_energy_vs_pos_wix_tspline_sec3->GetYaxis()->SetLabelSize(0.045);
h2_energy_vs_pos_wix_tspline_sec3->GetYaxis()->SetTitleSize(0.045);

const char* char_fname = fname.c_str(); 
TFile* file_input(TFile::Open(char_fname ,"READ"));
TFile* file_tcut_sec0(TFile::Open("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/data/beta_corr/section0_cut.root","READ"));
TFile* file_tcut_sec1(TFile::Open("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/data/beta_corr/section1_cut.root","READ"));
TFile* file_tcut_sec2(TFile::Open("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/data/beta_corr/section2_cut.root","READ"));
TFile* file_tcut_sec3(TFile::Open("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/data/beta_corr/section3_cut.root","READ"));

TList* mylist = (TList*)file_input->Get("histlist");

TCutG* mycut_sec0= (TCutG*)file_tcut_sec0->Get("section0_cut");
TCutG* mycut_sec1= (TCutG*)file_tcut_sec1->Get("section1_cut");
TCutG* mycut_sec2= (TCutG*)file_tcut_sec2->Get("section2_cut");
TCutG* mycut_sec3= (TCutG*)file_tcut_sec3->Get("section3_cut");

//mycut->SetVarX("x");
//mycut->SetVarY("y");
sprintf(hist_name,"Position on 11th anode in music vs Calibrated Energy messel,section0");
TH2D* h2_energy_vs_pos_messel_sec0 = (TH2D*)mylist->FindObject(hist_name);
sprintf(hist_name,"Position on 11th anode in music vs Calibrated Energy messel,section1");
TH2D* h2_energy_vs_pos_messel_sec1 = (TH2D*)mylist->FindObject(hist_name);
sprintf(hist_name,"Position on 11th anode in music vs Calibrated Energy Wixhausen,section2");
TH2D* h2_energy_vs_pos_wix_sec2 = (TH2D*)mylist->FindObject(hist_name);
sprintf(hist_name,"Position on 11th anode in music vs Calibrated Energy Wixhausen,section3");
TH2D* h2_energy_vs_pos_wix_sec3 = (TH2D*)mylist->FindObject(hist_name);
for(Int_t k = 1; k < 600; k++){
	for(Int_t l = 1; l < 950; l++){
		Double_t y_val_messel_sec0 = h2_energy_vs_pos_messel_sec0->GetYaxis()->GetBinCenter(l);
		Double_t x_val_messel_sec0 = h2_energy_vs_pos_messel_sec0->GetXaxis()->GetBinCenter(k);
		if (mycut_sec0->IsInside(x_val_messel_sec0,y_val_messel_sec0)){
			Double_t  bin_cont = h2_energy_vs_pos_messel_sec0->GetBinContent(k,l);
			h2_energy_vs_pos_messel_tspline_sec0->SetBinContent(k,l,bin_cont);
			}
		Double_t y_val_messel_sec1 = h2_energy_vs_pos_messel_sec1->GetYaxis()->GetBinCenter(l);
		Double_t x_val_messel_sec1 = h2_energy_vs_pos_messel_sec1->GetXaxis()->GetBinCenter(k);
		if (mycut_sec1->IsInside(x_val_messel_sec1,y_val_messel_sec1)){
			Double_t  bin_cont = h2_energy_vs_pos_messel_sec1->GetBinContent(k,l);
			h2_energy_vs_pos_messel_tspline_sec1->SetBinContent(k,l,bin_cont);
			}
		Double_t y_val_wix_sec2 = h2_energy_vs_pos_wix_sec2->GetYaxis()->GetBinCenter(l);
		Double_t x_val_wix_sec2 = h2_energy_vs_pos_wix_sec2->GetXaxis()->GetBinCenter(k);
		if(mycut_sec2->IsInside(x_val_wix_sec2,y_val_wix_sec2)){
			Double_t  bin_cont = h2_energy_vs_pos_wix_sec2->GetBinContent(k,l);
			h2_energy_vs_pos_wix_tspline_sec2->SetBinContent(k,l,bin_cont);
			}
		Double_t y_val_wix_sec3 = h2_energy_vs_pos_wix_sec3->GetYaxis()->GetBinCenter(l);
		Double_t x_val_wix_sec3 = h2_energy_vs_pos_wix_sec3->GetXaxis()->GetBinCenter(k);
		if(mycut_sec3->IsInside(x_val_wix_sec3,y_val_wix_sec3)){
			Double_t  bin_cont = h2_energy_vs_pos_wix_sec3->GetBinContent(k,l);
			h2_energy_vs_pos_wix_tspline_sec3->SetBinContent(k,l,bin_cont);
			}

		}
	}
TH1D* h1_profile_energy_vs_pos_messel_tspline_sec0 = (TH1D*)h2_energy_vs_pos_messel_tspline_sec0->ProfileX();
TH1D* h1_profile_energy_vs_pos_messel_tspline_sec1 = (TH1D*)h2_energy_vs_pos_messel_tspline_sec1->ProfileX();
TH1D* h1_profile_energy_vs_pos_wix_tspline_sec2 = (TH1D*)h2_energy_vs_pos_wix_tspline_sec2->ProfileX();
TH1D* h1_profile_energy_vs_pos_wix_tspline_sec3 = (TH1D*)h2_energy_vs_pos_wix_tspline_sec3->ProfileX();
TSpline3* spline_messel_sec0 = new TSpline3(h1_profile_energy_vs_pos_messel_tspline_sec0);
TSpline3* spline_messel_sec1 = new TSpline3(h1_profile_energy_vs_pos_messel_tspline_sec1);
TSpline3* spline_wix_sec2 = new TSpline3(h1_profile_energy_vs_pos_wix_tspline_sec2);
TSpline3* spline_wix_sec3 = new TSpline3(h1_profile_energy_vs_pos_wix_tspline_sec3);
cout << "eval of tspline messel sec3 " << spline_messel_sec0->Eval(60) << endl;



Double_t mean_messel_sec0 = h2_energy_vs_pos_messel_tspline_sec0->GetMean(2);
Double_t mean_messel_sec1 = h2_energy_vs_pos_messel_tspline_sec1->GetMean(2);
Double_t mean_wix_sec2 = h2_energy_vs_pos_wix_tspline_sec2->GetMean(2);
Double_t mean_wix_sec3 = h2_energy_vs_pos_wix_tspline_sec3->GetMean(2);
cout << "the mean value on messel side SECTION 0 is:\t" << mean_messel_sec0 << endl;
cout << "the mean value on messel side SECTION 1 is:\t" << mean_messel_sec1 << endl;
cout << "the mean value on wixhausen side SECTION 2 is:\t" << mean_wix_sec2 << endl;
cout << "the mean value on wixhausen side SECTION 3 is:\t" << mean_wix_sec3 << endl;
//h2_energy_vs_pos_wix_tspline->Draw("colz");
//spline_wix->Draw("same");
char f_out_name[500];
sprintf(f_out_name,"spline_position_messel_and_wixhausen.root");
TFile * f = new TFile(f_out_name,"RECREATE");
TList *list_spline = new TList();
list_spline->Add(spline_messel_sec0);
list_spline->Add(spline_messel_sec1);
list_spline->Add(spline_wix_sec2);
list_spline->Add(spline_wix_sec3);
list_spline->Write("list_tspline", TObject::kSingleKey);
}
