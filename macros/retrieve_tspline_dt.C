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
void retrieve_tspline_dt(const string& input_str){
string fname = string(input_str);
const char* char_fname= fname.c_str();
char hist_name[500];
vector<vector<double> > vec_para;
//Define histograms to fit:
TH2D* h2_deltax_vs_xal_sec0[16];
TH2D* h2_deltax_vs_xal_sec1[16];
TH2D* h2_deltax_vs_xal_sec2[16];
TH2D* h2_deltax_vs_xal_sec3[16];

TH1D* h1_deltax_y_proj_sec0[16];
TH1D* h1_deltax_y_proj_sec1[16];
TH1D* h1_deltax_y_proj_sec2[16];
TH1D* h1_deltax_y_proj_sec3[16];

//par file for the standard deviation of X_cal:
#include <ctime>
time_t *rawtime = new time_t;
struct tm * timeinfo;
time(rawtime);
timeinfo = localtime(rawtime);
ofstream par_file;
par_file.open("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/parameters/xcal_std_dev_params.csv");

par_file << "#FILE created by retrieve_tspline_dt.C at time:" << timeinfo->tm_hour << ":" << timeinfo->tm_min << "Date:\t" << timeinfo->tm_mday << "/" << timeinfo->tm_mon << "/" << timeinfo->tm_year << endl;
TFile* file_input(TFile::Open(char_fname ,"READ"));
file_input->ls();


//section 0

TH1D* h1_deltax_profile_sec0[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name, "Section0: TProfile for DeltaX vs Xal for AnodeID %i ",i);
h1_deltax_profile_sec0[i] = new TH1D(hist_name,hist_name,1250,-50,200);
h1_deltax_profile_sec0[i]->GetXaxis()->SetTitle("Xal [mm]");
h1_deltax_profile_sec0[i]->GetYaxis()->SetTitle("DeltaX [mm]");
h1_deltax_profile_sec0[i]->GetXaxis()->CenterTitle(true);
h1_deltax_profile_sec0[i]->GetYaxis()->CenterTitle(true);
h1_deltax_profile_sec0[i]->GetYaxis()->SetLabelSize(0.045);
h1_deltax_profile_sec0[i]->GetYaxis()->SetTitleSize(0.045);
}

TH2D* h2_delta_x_vs_xcal_sec0[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name,"Section 0, Anode %i -> Xfit - Xcal vs Xcal (using TSpline3)",i);
h2_delta_x_vs_xcal_sec0[i] = new TH2D(hist_name,hist_name,1250,-50,200,500,-5,5);
h2_delta_x_vs_xcal_sec0[i]->GetXaxis()->SetTitle("Xcal [mm]");
h2_delta_x_vs_xcal_sec0[i]->GetYaxis()->SetTitle("DeltaX (Xfit -Xal) [mm]");
h2_delta_x_vs_xcal_sec0[i]->GetXaxis()->CenterTitle(true);
h2_delta_x_vs_xcal_sec0[i]->GetYaxis()->CenterTitle(true);
h2_delta_x_vs_xcal_sec0[i]->GetYaxis()->SetLabelSize(0.045);
h2_delta_x_vs_xcal_sec0[i]->GetYaxis()->SetTitleSize(0.045);
}

//section 1

TH1D* h1_deltax_profile_sec1[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name, "Section1: TProfile for DeltaX vs Xal for AnodeID %i ",i);
h1_deltax_profile_sec1[i] = new TH1D(hist_name,hist_name,1250,-50,200);
h1_deltax_profile_sec1[i]->GetXaxis()->SetTitle("Xal [mm]");
h1_deltax_profile_sec1[i]->GetYaxis()->SetTitle("DeltaX [mm]");
h1_deltax_profile_sec1[i]->GetXaxis()->CenterTitle(true);
h1_deltax_profile_sec1[i]->GetYaxis()->CenterTitle(true);
h1_deltax_profile_sec1[i]->GetYaxis()->SetLabelSize(0.045);
h1_deltax_profile_sec1[i]->GetYaxis()->SetTitleSize(0.045);
}

TH2D* h2_delta_x_vs_xcal_sec1[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name,"Section 1, Anode %i -> Xfit - Xcal vs Xcal (using TSpline3)",i);
h2_delta_x_vs_xcal_sec1[i] = new TH2D(hist_name,hist_name,1250,-50,200,500,-5,5);
h2_delta_x_vs_xcal_sec1[i]->GetXaxis()->SetTitle("Xcal [mm]");
h2_delta_x_vs_xcal_sec1[i]->GetYaxis()->SetTitle("DeltaX (Xfit -Xal) [mm]");
h2_delta_x_vs_xcal_sec1[i]->GetXaxis()->CenterTitle(true);
h2_delta_x_vs_xcal_sec1[i]->GetYaxis()->CenterTitle(true);
h2_delta_x_vs_xcal_sec1[i]->GetYaxis()->SetLabelSize(0.045);
h2_delta_x_vs_xcal_sec1[i]->GetYaxis()->SetTitleSize(0.045);
}

//section 2

TH1D* h1_deltax_profile_sec2[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name, "Section2: TProfile for DeltaX vs Xal for AnodeID %i ",i);
h1_deltax_profile_sec2[i] = new TH1D(hist_name,hist_name,1250,-200,50);
h1_deltax_profile_sec2[i]->GetXaxis()->SetTitle("Xal [mm]");
h1_deltax_profile_sec2[i]->GetYaxis()->SetTitle("DeltaX [mm]");
h1_deltax_profile_sec2[i]->GetXaxis()->CenterTitle(true);
h1_deltax_profile_sec2[i]->GetYaxis()->CenterTitle(true);
h1_deltax_profile_sec2[i]->GetYaxis()->SetLabelSize(0.045);
h1_deltax_profile_sec2[i]->GetYaxis()->SetTitleSize(0.045);
}

TH2D* h2_delta_x_vs_xcal_sec2[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name,"Section 2, Anode %i -> Xfit - Xcal vs Xcal (using TSpline3)",i);
h2_delta_x_vs_xcal_sec2[i] = new TH2D(hist_name,hist_name,1250,-200,50,500,-5,5);
h2_delta_x_vs_xcal_sec2[i]->GetXaxis()->SetTitle("Xcal [mm]");
h2_delta_x_vs_xcal_sec2[i]->GetYaxis()->SetTitle("DeltaX (Xfit -Xal) [mm]");
h2_delta_x_vs_xcal_sec2[i]->GetXaxis()->CenterTitle(true);
h2_delta_x_vs_xcal_sec2[i]->GetYaxis()->CenterTitle(true);
h2_delta_x_vs_xcal_sec2[i]->GetYaxis()->SetLabelSize(0.045);
h2_delta_x_vs_xcal_sec2[i]->GetYaxis()->SetTitleSize(0.045);
}
//section 3

TH1D* h1_deltax_profile_sec3[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name, "Section3: TProfile for DeltaX vs Xal for AnodeID %i ",i);
h1_deltax_profile_sec3[i] = new TH1D(hist_name,hist_name,1250,-200,50);
h1_deltax_profile_sec3[i]->GetXaxis()->SetTitle("Xal [mm]");
h1_deltax_profile_sec3[i]->GetYaxis()->SetTitle("DeltaX [mm]");
h1_deltax_profile_sec3[i]->GetXaxis()->CenterTitle(true);
h1_deltax_profile_sec3[i]->GetYaxis()->CenterTitle(true);
h1_deltax_profile_sec3[i]->GetYaxis()->SetLabelSize(0.045);
h1_deltax_profile_sec3[i]->GetYaxis()->SetTitleSize(0.045);
}

TH2D* h2_delta_x_vs_xcal_sec3[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name,"Section 3, Anode %i -> Xfit - Xcal vs Xcal (using TSpline3)",i);
h2_delta_x_vs_xcal_sec3[i] = new TH2D(hist_name,hist_name,1250,-200,50,500,-5,5);
h2_delta_x_vs_xcal_sec3[i]->GetXaxis()->SetTitle("Xcal [mm]");
h2_delta_x_vs_xcal_sec3[i]->GetYaxis()->SetTitle("DeltaX (Xfit -Xal) [mm]");
h2_delta_x_vs_xcal_sec3[i]->GetXaxis()->CenterTitle(true);
h2_delta_x_vs_xcal_sec3[i]->GetYaxis()->CenterTitle(true);
h2_delta_x_vs_xcal_sec3[i]->GetYaxis()->SetLabelSize(0.045);
h2_delta_x_vs_xcal_sec3[i]->GetYaxis()->SetTitleSize(0.045);
}
char f_out_name[500];
char f_splines[500];
sprintf(f_out_name,"/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/data/output_retrieve_tspline_dt/retrieve_tspline_result.root");
TFile * f = new TFile(f_out_name,"RECREATE");
//
//TCutG *gcut1 = new TCutG("gcut1",7);

//section 0
TList* MyHistList0= (TList*)file_input->Get("histlist");
TCanvas sec0_c("sec0_c","Section0: TProfile for DeltaX vs Xal for Anodes");
sec0_c.Divide(4,4);
for (Int_t i = 0; i < 16; i++){
	sec0_c.cd(i+1);
	gPad->SetLogz(i+1);
	gStyle->SetOptStat(1);
	gStyle->SetOptFit(111);
	sprintf(hist_name,"Section 0, Anode %i -> DeltaX vs Xal",i); 
	h2_deltax_vs_xal_sec0[i] = (TH2D*)MyHistList0->FindObject(hist_name);
	h1_deltax_profile_sec0[i] = (TH1D*)h2_deltax_vs_xal_sec0[i]->ProfileX();
	TSpline3* splines_sec0 = new TSpline3(h1_deltax_profile_sec0[i],0,100);
	
	//filling histo
	for(Int_t k = 0; k < 1100; k++){ //bins in x range
		for(Int_t l = 0; l < 1000; l++){ //bins in y range
			Double_t delta_x = h2_deltax_vs_xal_sec0[i]->GetYaxis()->GetBinCenter(l);
			Double_t x_al = h2_deltax_vs_xal_sec0[i]->GetXaxis()->GetBinCenter(k);
			Double_t  bin_cont = h2_deltax_vs_xal_sec0[i]->GetBinContent(k,l);
			Double_t x_fit = delta_x + x_al;
			Double_t x_cal = x_al + (splines_sec0->Eval(x_al));
			Double_t new_y_val = x_fit - x_cal;
			TAxis *xaxis = h2_delta_x_vs_xcal_sec0[i]->GetXaxis();
			TAxis *yaxis = h2_delta_x_vs_xcal_sec0[i]->GetYaxis();
			Int_t binx = xaxis->FindBin(x_cal);
			Int_t biny = yaxis->FindBin(new_y_val);
			h2_delta_x_vs_xcal_sec0[i]->SetBinContent(binx,biny,bin_cont);
			}
		}
	h2_deltax_vs_xal_sec0[i]->Draw("colz");
	splines_sec0->Draw("same");
	sprintf(f_splines,"/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/parameters/spline_sec_0_anode_%i.xml",i);	
	splines_sec0->SaveAs(f_splines);

}
sec0_c.Modified();
sec0_c.Update();
sec0_c.Print("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/data/output_retrieve_tspline_dt/deltax_vs_al_sec0.pdf");
sec0_c.Write();

TCanvas tspline_corr_sec0_c("tspline_corr_sec0_c","Section0: Xfit -Xcal vs Xcal ");
tspline_corr_sec0_c.Divide(4,4);
for (Int_t i = 0; i < 16; i++){
	tspline_corr_sec0_c.cd(i+1);
	gPad->SetLogz(i+1);
	gStyle->SetOptStat(1);
	gStyle->SetOptFit(111);
	h2_delta_x_vs_xcal_sec0[i]->Draw("colz");
	h1_deltax_y_proj_sec0[i] = h2_delta_x_vs_xcal_sec0[i]->ProjectionY();
	h1_deltax_y_proj_sec0[i]->GetXaxis()->SetRangeUser(-0.1,0.1);
	}
tspline_corr_sec0_c.Modified();
tspline_corr_sec0_c.Update();
tspline_corr_sec0_c.Print("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/data/output_retrieve_tspline_dt/deltax_vs_xcal_sec0.pdf");
tspline_corr_sec0_c.Write();

TCanvas resolution_sec0("resolution_sec0","Section0:Projection Deltax with gauss fit");
resolution_sec0.Divide(4,4);
for (Int_t i = 0; i < 16; i++){
	resolution_sec0.cd(i+1);
	gStyle->SetOptFit(111);
	h1_deltax_y_proj_sec0[i]->Draw();
	TF1* my_func = new TF1("my_func","gaus");
	h1_deltax_y_proj_sec0[i]->Fit("my_func","Q","",-0.1,0.1);
	my_func->Draw("same");
	par_file << 0 << "," << i << ","  << my_func->GetParameter(2) << endl;
	}
resolution_sec0.Modified();
resolution_sec0.Update();
resolution_sec0.Print("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/data/output_retrieve_tspline_dt/resolution_sec0.pdf");
resolution_sec0.Write();


//section 1


TCanvas sec1_c("sec1_c","Section1: TProfile for DeltaX vs Xal for Anodes");
sec1_c.Divide(4,4);
for (Int_t i = 0; i < 16; i++){
	sec1_c.cd(i+1);
	gPad->SetLogz(i+1);
	gStyle->SetOptStat(1);
	gStyle->SetOptFit(111);
	sprintf(hist_name,"Section 1, Anode %i -> DeltaX vs Xal",i); 
	h2_deltax_vs_xal_sec1[i] = (TH2D*)MyHistList0->FindObject(hist_name);
	h1_deltax_profile_sec1[i] = (TH1D*)h2_deltax_vs_xal_sec1[i]->ProfileX();
	TSpline3* splines_sec1 = new TSpline3(h1_deltax_profile_sec1[i]);
	
	//filling histo
	for(Int_t k = 0; k < 1100; k++){ //bins in x range
		for(Int_t l = 0; l < 1000; l++){ //bins in y range
			Double_t delta_x = h2_deltax_vs_xal_sec1[i]->GetYaxis()->GetBinCenter(l);
			Double_t x_al = h2_deltax_vs_xal_sec1[i]->GetXaxis()->GetBinCenter(k);
			Double_t  bin_cont = h2_deltax_vs_xal_sec1[i]->GetBinContent(k,l);
			Double_t x_fit = delta_x + x_al;
			Double_t x_cal = x_al + (splines_sec1->Eval(x_al));
			Double_t new_y_val = x_fit - x_cal;
			TAxis *xaxis = h2_delta_x_vs_xcal_sec1[i]->GetXaxis();
			TAxis *yaxis = h2_delta_x_vs_xcal_sec1[i]->GetYaxis();
			Int_t binx = xaxis->FindBin(x_cal);
			Int_t biny = yaxis->FindBin(new_y_val);
			h2_delta_x_vs_xcal_sec1[i]->SetBinContent(binx,biny,bin_cont);
			}
		}
	h2_deltax_vs_xal_sec1[i]->Draw("colz");
	splines_sec1->Draw("same");
	sprintf(f_splines,"/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/parameters/spline_sec_1_anode_%i.root",i);	
	splines_sec1->SaveAs(f_splines);

}
sec1_c.Modified();
sec1_c.Update();
sec1_c.Print("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/data/output_retrieve_tspline_dt/deltax_vs_al_sec1.pdf");
sec1_c.Write();

TCanvas tspline_corr_sec1_c("tspline_corr_sec1_c","Section1: Xfit -Xcal vs Xcal ");
tspline_corr_sec1_c.Divide(4,4);
for (Int_t i = 0; i < 16; i++){
	tspline_corr_sec1_c.cd(i+1);
	gPad->SetLogz(i+1);
	gStyle->SetOptStat(1);
	gStyle->SetOptFit(111);
	h2_delta_x_vs_xcal_sec1[i]->Draw("colz");
	h1_deltax_y_proj_sec1[i] = h2_delta_x_vs_xcal_sec1[i]->ProjectionY();
	h1_deltax_y_proj_sec1[i]->GetXaxis()->SetRangeUser(-0.1,0.1);
	}
tspline_corr_sec1_c.Modified();
tspline_corr_sec1_c.Update();
tspline_corr_sec1_c.Print("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/data/output_retrieve_tspline_dt/deltax_vs_xcal_sec1.pdf");
tspline_corr_sec1_c.Write();

TCanvas resolution_sec1("resolution_sec1","Section1:Projection Deltax with gauss fit");
resolution_sec1.Divide(4,4);
for (Int_t i = 0; i < 16; i++){
	resolution_sec1.cd(i+1);
	gStyle->SetOptFit(111);
	h1_deltax_y_proj_sec1[i]->Draw();
	TF1* my_func = new TF1("my_func","gaus");
	h1_deltax_y_proj_sec1[i]->Fit("my_func","Q","",-0.1,0.1);
	my_func->Draw("same");
	par_file << 1 << "," << i << ","  << my_func->GetParameter(2) << endl;
	}
resolution_sec1.Modified();
resolution_sec1.Update();
resolution_sec1.Print("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/data/output_retrieve_tspline_dt/resolution_sec1.pdf");
resolution_sec1.Write();


//section 2

TCanvas sec2_c("sec2_c","Section2: TProfile for DeltaX vs Xal for Anodes");
sec2_c.Divide(4,4);
for (Int_t i = 0; i < 16; i++){
	sec2_c.cd(i+1);
	gPad->SetLogz(i+1);
	gStyle->SetOptStat(1);
	gStyle->SetOptFit(111);
	sprintf(hist_name,"Section 2, Anode %i -> DeltaX vs Xal",i); 
	h2_deltax_vs_xal_sec2[i] = (TH2D*)MyHistList0->FindObject(hist_name);
	h1_deltax_profile_sec2[i] = (TH1D*)h2_deltax_vs_xal_sec2[i]->ProfileX();
	TSpline3* splines_sec2 = new TSpline3(h1_deltax_profile_sec2[i]);
	
	//filling histo
	for(Int_t k = 0; k < 1100; k++){ //bins in x range
		for(Int_t l = 0; l < 1000; l++){ //bins in y range
			Double_t delta_x = h2_deltax_vs_xal_sec2[i]->GetYaxis()->GetBinCenter(l);
			Double_t x_al = h2_deltax_vs_xal_sec2[i]->GetXaxis()->GetBinCenter(k);
			Double_t  bin_cont = h2_deltax_vs_xal_sec2[i]->GetBinContent(k,l);
			Double_t x_fit = delta_x + x_al;
			Double_t x_cal = x_al + (splines_sec2->Eval(x_al));
			Double_t new_y_val = x_fit - x_cal;
			TAxis *xaxis = h2_delta_x_vs_xcal_sec2[i]->GetXaxis();
			TAxis *yaxis = h2_delta_x_vs_xcal_sec2[i]->GetYaxis();
			Int_t binx = xaxis->FindBin(x_cal);
			Int_t biny = yaxis->FindBin(new_y_val);
			h2_delta_x_vs_xcal_sec2[i]->SetBinContent(binx,biny,bin_cont);
			}
		}
	h2_deltax_vs_xal_sec2[i]->Draw("colz");
	splines_sec2->Draw("same");
	sprintf(f_splines,"/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/parameters/spline_sec_2_anode_%i.root",i);	
	splines_sec2->SaveAs(f_splines);

}
sec2_c.Modified();
sec2_c.Update();
sec2_c.Print("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/data/output_retrieve_tspline_dt/deltax_vs_al_sec2.pdf");
sec2_c.Write();

TCanvas tspline_corr_sec2_c("tspline_corr_sec2_c","Section2: Xfit -Xcal vs Xcal ");
tspline_corr_sec2_c.Divide(4,4);
for (Int_t i = 0; i < 16; i++){
	tspline_corr_sec2_c.cd(i+1);
	gPad->SetLogz(i+1);
	gStyle->SetOptStat(1);
	gStyle->SetOptFit(111);
	h2_delta_x_vs_xcal_sec2[i]->Draw("colz");
	h1_deltax_y_proj_sec2[i] = h2_delta_x_vs_xcal_sec2[i]->ProjectionY();
	h1_deltax_y_proj_sec2[i]->GetXaxis()->SetRangeUser(-0.1,0.1);
	}
tspline_corr_sec2_c.Modified();
tspline_corr_sec2_c.Update();
tspline_corr_sec2_c.Print("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/data/output_retrieve_tspline_dt/deltax_vs_xcal_sec2.pdf");
tspline_corr_sec2_c.Write();

TCanvas resolution_sec2("resolution_sec2","Section2:Projection Deltax with gauss fit");
resolution_sec2.Divide(4,4);
for (Int_t i = 0; i < 16; i++){
	resolution_sec2.cd(i+1);
	gStyle->SetOptFit(111);
	h1_deltax_y_proj_sec2[i]->Draw();
	TF1* my_func = new TF1("my_func","gaus");
	h1_deltax_y_proj_sec2[i]->Fit("my_func","Q","",-0.1,0.1);
	my_func->Draw("same");
	par_file << 2 << "," << i << ","  << my_func->GetParameter(2) << endl;
	}
resolution_sec2.Modified();
resolution_sec2.Update();
resolution_sec2.Print("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/data/output_retrieve_tspline_dt/resolution_sec2.pdf");
resolution_sec2.Write();


//section 3


TCanvas sec3_c("sec3_c","Section3: TProfile for DeltaX vs Xal for Anodes");
sec3_c.Divide(4,4);
for (Int_t i = 0; i < 16; i++){
	sec3_c.cd(i+1);
	gPad->SetLogz(i+1);
	gStyle->SetOptStat(1);
	gStyle->SetOptFit(111);
	sprintf(hist_name,"Section 3, Anode %i -> DeltaX vs Xal",i); 
	h2_deltax_vs_xal_sec3[i] = (TH2D*)MyHistList0->FindObject(hist_name);
	h1_deltax_profile_sec3[i] = (TH1D*)h2_deltax_vs_xal_sec3[i]->ProfileX();
	TSpline3* splines_sec3 = new TSpline3(h1_deltax_profile_sec3[i]);
	
	//filling histo
	for(Int_t k = 0; k < 1100; k++){ //bins in x range
		for(Int_t l = 0; l < 1000; l++){ //bins in y range
			Double_t delta_x = h2_deltax_vs_xal_sec3[i]->GetYaxis()->GetBinCenter(l);
			Double_t x_al = h2_deltax_vs_xal_sec3[i]->GetXaxis()->GetBinCenter(k);
			Double_t  bin_cont = h2_deltax_vs_xal_sec3[i]->GetBinContent(k,l);
			Double_t x_fit = delta_x + x_al;
			Double_t x_cal = x_al + (splines_sec3->Eval(x_al));
			Double_t new_y_val = x_fit - x_cal;
			TAxis *xaxis = h2_delta_x_vs_xcal_sec3[i]->GetXaxis();
			TAxis *yaxis = h2_delta_x_vs_xcal_sec3[i]->GetYaxis();
			Int_t binx = xaxis->FindBin(x_cal);
			Int_t biny = yaxis->FindBin(new_y_val);
			h2_delta_x_vs_xcal_sec3[i]->SetBinContent(binx,biny,bin_cont);
			}
		}
	h2_deltax_vs_xal_sec3[i]->Draw("colz");
	splines_sec3->Draw("same");
	sprintf(f_splines,"/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/parameters/spline_sec_3_anode_%i.root",i);	
	splines_sec3->SaveAs(f_splines);

}
sec3_c.Modified();
sec3_c.Update();
sec3_c.Print("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/data/output_retrieve_tspline_dt/deltax_vs_al_sec3.pdf");
sec3_c.Write();

TCanvas tspline_corr_sec3_c("tspline_corr_sec3_c","Section3: Xfit -Xcal vs Xcal ");
tspline_corr_sec3_c.Divide(4,4);
for (Int_t i = 0; i < 16; i++){
	tspline_corr_sec3_c.cd(i+1);
	gPad->SetLogz(i+1);
	gStyle->SetOptStat(1);
	gStyle->SetOptFit(111);
	h2_delta_x_vs_xcal_sec3[i]->Draw("colz");
	h1_deltax_y_proj_sec3[i] = h2_delta_x_vs_xcal_sec3[i]->ProjectionY();
	h1_deltax_y_proj_sec3[i]->GetXaxis()->SetRangeUser(-0.1,0.1);
	}
tspline_corr_sec3_c.Modified();
tspline_corr_sec3_c.Update();
tspline_corr_sec3_c.Print("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/data/output_retrieve_tspline_dt/deltax_vs_xcal_sec3.pdf");
tspline_corr_sec3_c.Write();

TCanvas resolution_sec3("resolution_sec3","Section3:Projection Deltax with gauss fit");
resolution_sec3.Divide(4,4);
for (Int_t i = 0; i < 16; i++){
	resolution_sec3.cd(i+1);
	gStyle->SetOptFit(111);
	h1_deltax_y_proj_sec3[i]->Draw();
	TF1* my_func = new TF1("my_func","gaus");
	h1_deltax_y_proj_sec3[i]->Fit("my_func","Q","",-0.1,0.1);
	my_func->Draw("same");
	par_file << 3 << "," << i << ","  << my_func->GetParameter(2) << endl;
	}
resolution_sec3.Modified();
resolution_sec3.Update();
resolution_sec3.Print("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/data/output_retrieve_tspline_dt/resolution_sec3.pdf");
resolution_sec3.Write();

TList *l = new TList();
for (Int_t i = 0; i < 16; i++){
l->Add(h2_delta_x_vs_xcal_sec0[i]);
}
for (Int_t i = 0; i < 16; i++){
l->Add(h2_delta_x_vs_xcal_sec1[i]);
}
for (Int_t i = 0; i < 16; i++){
l->Add(h2_delta_x_vs_xcal_sec2[i]);
}
for (Int_t i = 0; i < 16; i++){
l->Add(h2_delta_x_vs_xcal_sec3[i]);
}
l->Write("histlist", TObject::kSingleKey);
par_file.close();
}
