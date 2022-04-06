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
#include <iostream>
#include <fstream>


#include <stdlib.h>
using namespace std;
void retrieve_drift_fits(const string& input_str){
string fname = string(input_str);
const char* char_fname= fname.c_str();
char hist_name[500];

//--- HISTOS ----------------

//drift time vs Xal histograms for section 0
TH2D* h2_xal_vs_DTraw_sec0[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name,"Section 0, Anode %i -> Xal vs DTraw",i);
h2_xal_vs_DTraw_sec0[i] = new TH2D(hist_name,hist_name,325,0,32500,130,-120,10);
h2_xal_vs_DTraw_sec0[i]->GetXaxis()->SetTitle("Drift Time");
h2_xal_vs_DTraw_sec0[i]->GetYaxis()->SetTitle("Xal (extrapol. from MW1/2)");
h2_xal_vs_DTraw_sec0[i]->GetXaxis()->CenterTitle(true);
h2_xal_vs_DTraw_sec0[i]->GetYaxis()->CenterTitle(true);
h2_xal_vs_DTraw_sec0[i]->GetYaxis()->SetLabelSize(0.045);
h2_xal_vs_DTraw_sec0[i]->GetYaxis()->SetTitleSize(0.045);
}


TH2D* h2_xal_vs_DTraw_sec0_after_fit[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name,"Section 0, Anode %i -> Xal vs DTraw after fit",i);
h2_xal_vs_DTraw_sec0_after_fit[i] = new TH2D(hist_name,hist_name,325,0,32500,130,-10,120);
h2_xal_vs_DTraw_sec0_after_fit[i]->GetXaxis()->SetTitle("Drift Time");
h2_xal_vs_DTraw_sec0_after_fit[i]->GetYaxis()->SetTitle("Xal (extrapol. from MW1/2)");
h2_xal_vs_DTraw_sec0_after_fit[i]->GetXaxis()->CenterTitle(true);
h2_xal_vs_DTraw_sec0_after_fit[i]->GetYaxis()->CenterTitle(true);
h2_xal_vs_DTraw_sec0_after_fit[i]->GetYaxis()->SetLabelSize(0.045);
h2_xal_vs_DTraw_sec0_after_fit[i]->GetYaxis()->SetTitleSize(0.045);
}

//drift time vs Xal histograms for section 1

TH2D* h2_xal_vs_DTraw_sec1[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name,"Section 1, Anode %i -> Xal vs DTraw",i);
h2_xal_vs_DTraw_sec1[i] = new TH2D(hist_name,hist_name,325,0,32500,400,-200,200);
h2_xal_vs_DTraw_sec1[i]->GetXaxis()->SetTitle("Drift Time");
h2_xal_vs_DTraw_sec1[i]->GetYaxis()->SetTitle("Xal (extrapol. from MW1/2)");
h2_xal_vs_DTraw_sec1[i]->GetXaxis()->CenterTitle(true);
h2_xal_vs_DTraw_sec1[i]->GetYaxis()->CenterTitle(true);
h2_xal_vs_DTraw_sec1[i]->GetYaxis()->SetLabelSize(0.045);
h2_xal_vs_DTraw_sec1[i]->GetYaxis()->SetTitleSize(0.045);
}


TH2D* h2_xal_vs_DTraw_sec1_after_fit[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name,"Section 1, Anode %i -> Xal vs DTraw after fit",i);
h2_xal_vs_DTraw_sec1_after_fit[i] = new TH2D(hist_name,hist_name,325,0,32500,130,-10,120);
h2_xal_vs_DTraw_sec1_after_fit[i]->GetXaxis()->SetTitle("Drift Time");
h2_xal_vs_DTraw_sec1_after_fit[i]->GetYaxis()->SetTitle("Xal (extrapol. from MW1/2)");
h2_xal_vs_DTraw_sec1_after_fit[i]->GetXaxis()->CenterTitle(true);
h2_xal_vs_DTraw_sec1_after_fit[i]->GetYaxis()->CenterTitle(true);
h2_xal_vs_DTraw_sec1_after_fit[i]->GetYaxis()->SetLabelSize(0.045);
h2_xal_vs_DTraw_sec1_after_fit[i]->GetYaxis()->SetTitleSize(0.045);
}

//drift time vs Xal histograms for section 2

TH2D* h2_xal_vs_DTraw_sec2[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name,"Section 2, Anode %i -> Xal vs DTraw",i);
h2_xal_vs_DTraw_sec2[i] = new TH2D(hist_name,hist_name,325,0,32500,400,-200,200);
h2_xal_vs_DTraw_sec2[i]->GetXaxis()->SetTitle("Drift Time");
h2_xal_vs_DTraw_sec2[i]->GetYaxis()->SetTitle("Xal (extrapol. from MW1/2)");
h2_xal_vs_DTraw_sec2[i]->GetXaxis()->CenterTitle(true);
h2_xal_vs_DTraw_sec2[i]->GetYaxis()->CenterTitle(true);
h2_xal_vs_DTraw_sec2[i]->GetYaxis()->SetLabelSize(0.045);
h2_xal_vs_DTraw_sec2[i]->GetYaxis()->SetTitleSize(0.045);
}

TH2D* h2_xal_vs_DTraw_sec2_after_fit[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name,"Section 2, Anode %i -> Xal vs DTraw after fit",i);
h2_xal_vs_DTraw_sec2_after_fit[i] = new TH2D(hist_name,hist_name,325,0,32500,130,-120,10);
h2_xal_vs_DTraw_sec2_after_fit[i]->GetXaxis()->SetTitle("Drift Time");
h2_xal_vs_DTraw_sec2_after_fit[i]->GetYaxis()->SetTitle("Xal (extrapol. from MW1/2)");
h2_xal_vs_DTraw_sec2_after_fit[i]->GetXaxis()->CenterTitle(true);
h2_xal_vs_DTraw_sec2_after_fit[i]->GetYaxis()->CenterTitle(true);
h2_xal_vs_DTraw_sec2_after_fit[i]->GetYaxis()->SetLabelSize(0.045);
h2_xal_vs_DTraw_sec2_after_fit[i]->GetYaxis()->SetTitleSize(0.045);
}
//drift time vs Xal histograms for section 3

TH2D* h2_xal_vs_DTraw_sec3[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name,"Section 3, Anode %i -> Xal vs DTraw",i);
h2_xal_vs_DTraw_sec3[i] = new TH2D(hist_name,hist_name,325,0,32500,400,-200,200);
h2_xal_vs_DTraw_sec3[i]->GetXaxis()->SetTitle("Drift Time");
h2_xal_vs_DTraw_sec3[i]->GetYaxis()->SetTitle("Xal (extrapol. from MW1/2)");
h2_xal_vs_DTraw_sec3[i]->GetXaxis()->CenterTitle(true);
h2_xal_vs_DTraw_sec3[i]->GetYaxis()->CenterTitle(true);
h2_xal_vs_DTraw_sec3[i]->GetYaxis()->SetLabelSize(0.045);
h2_xal_vs_DTraw_sec3[i]->GetYaxis()->SetTitleSize(0.045);
}
TH2D* h2_xal_vs_DTraw_sec3_after_fit[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name,"Section 3, Anode %i -> Xal vs DTraw after fit",i);
h2_xal_vs_DTraw_sec3_after_fit[i] = new TH2D(hist_name,hist_name,325,0,32500,130,-120,10);
h2_xal_vs_DTraw_sec3_after_fit[i]->GetXaxis()->SetTitle("Drift Time");
h2_xal_vs_DTraw_sec3_after_fit[i]->GetYaxis()->SetTitle("Xal (extrapol. from MW1/2)");
h2_xal_vs_DTraw_sec3_after_fit[i]->GetXaxis()->CenterTitle(true);
h2_xal_vs_DTraw_sec3_after_fit[i]->GetYaxis()->CenterTitle(true);
h2_xal_vs_DTraw_sec3_after_fit[i]->GetYaxis()->SetLabelSize(0.045);
h2_xal_vs_DTraw_sec3_after_fit[i]->GetYaxis()->SetTitleSize(0.045);
}


//---------------------------
//file to store alignment parameters
#include <ctime>
time_t *rawtime = new time_t;
struct tm * timeinfo;
time(rawtime);
timeinfo = localtime(rawtime);
ofstream par_file;
par_file.open("xal_vs_dt_params.csv");

par_file << "#FILE created by retrieve_drift_fits.C at time:" << timeinfo->tm_hour << ":" << timeinfo->tm_min << "Date:\t" << timeinfo->tm_mday << "/" << timeinfo->tm_mon << "/" << timeinfo->tm_year << endl;
TFile* file_input(TFile::Open(char_fname ,"READ"));
file_input->ls();

char f_out_name[500];
sprintf(f_out_name,"combined_dt_fit_plots.root");
TFile * f = new TFile(f_out_name,"RECREATE");

//section 0
TList* MyHistList0= (TList*)file_input->Get("histlist");
TCanvas sec0_c("sec0_c","Section0: Xal vs DTraw");
sec0_c.Divide(4,4);
for (Int_t i = 0; i < 16; i++){
	sec0_c.cd(i+1);
	gPad->SetLogz(i+1);
	gStyle->SetOptStat(1);
	gStyle->SetOptFit(111);
	sprintf(hist_name,"Section 0, Anode %i -> Xal vs DTraw",i);
	h2_xal_vs_DTraw_sec0[i] = (TH2D*)MyHistList0->FindObject(hist_name);
	Double_t mean = h2_xal_vs_DTraw_sec0[i]->GetMean(1);
	Double_t st_dev = h2_xal_vs_DTraw_sec0[i]->GetStdDev(1);
	TF1* fit_example_sec0 = new TF1("fit_example_sec0","[0]*x+ [1]",mean-1.2*st_dev,mean+1.2*st_dev);
	h2_xal_vs_DTraw_sec0[i]->Fit(fit_example_sec0, "R");
	Double_t slope = fit_example_sec0->GetParameter(0);
	Double_t offset = fit_example_sec0->GetParameter(1);
//	h2_xal_vs_DTraw_sec0[i]->Draw("colz");
	//fit_example_sec0->Draw("same");
	//use fit parameters to make new fit
	for (Int_t k = 0; k < 1625;k++){
		for (Int_t l = 0; l < 650;l++){
			Double_t xal = h2_xal_vs_DTraw_sec0[i]->GetYaxis()->GetBinCenter(l);
			Double_t dt_raw = h2_xal_vs_DTraw_sec0[i]->GetXaxis()->GetBinCenter(k);
			Double_t  bin_cont = h2_xal_vs_DTraw_sec0[i]->GetBinContent(k,l);
			if ( xal < dt_raw*slope+offset+5 && xal > dt_raw*slope+offset-5){
				TAxis *xaxis = h2_xal_vs_DTraw_sec0_after_fit[i]->GetXaxis();
				TAxis *yaxis = h2_xal_vs_DTraw_sec0_after_fit[i]->GetYaxis();
				Int_t binx = xaxis->FindBin(dt_raw);
				Int_t biny = yaxis->FindBin(xal);
				h2_xal_vs_DTraw_sec0_after_fit[i]->SetBinContent(binx,biny,bin_cont);
				}	
			}	
		}
	delete fit_example_sec0;
	h2_xal_vs_DTraw_sec0_after_fit[i]->Draw("colz");
	TF1* final_fit_sec0 = new TF1("final_fit_sec0","[0]*x+ [1]",mean -2*st_dev,mean+2*st_dev);
	final_fit_sec0->SetLineWidth(1);
	h2_xal_vs_DTraw_sec0_after_fit[i]->Fit(final_fit_sec0);
	final_fit_sec0->Draw("same");
	par_file << 0 << "," << i << "," << slope << "," << offset << endl;
}
sec0_c.Modified();
sec0_c.Update();
sec0_c.Print("section_0_drift_fit.pdf");
sec0_c.Write();

//section 1
TCanvas sec1_c("sec1_c","Section1: Xal vs DTraw");
sec1_c.Divide(4,4);
for (Int_t i = 0; i < 16; i++){
	sec1_c.cd(i+1);
	gPad->SetLogz(i+1);
	gStyle->SetOptStat(1);
	gStyle->SetOptFit(111);
	sprintf(hist_name,"Section 1, Anode %i -> Xal vs DTraw",i);
	h2_xal_vs_DTraw_sec1[i] = (TH2D*)MyHistList0->FindObject(hist_name);
	Double_t mean = h2_xal_vs_DTraw_sec1[i]->GetMean(1);
	Double_t st_dev = h2_xal_vs_DTraw_sec1[i]->GetStdDev(1);
	TF1* fit_example_sec1 = new TF1("fit_example_sec1","[0]*x+ [1]",mean-1.2*st_dev,mean+1.2*st_dev);
	h2_xal_vs_DTraw_sec1[i]->Fit(fit_example_sec1, "R");
	Double_t slope = fit_example_sec1->GetParameter(0);
	Double_t offset = fit_example_sec1->GetParameter(1);
	//h2_xal_vs_DTraw_sec1[i]->Draw("colz");
	//fit_example_sec1->Draw("same");
	//use fit parameters to make new fit
	for (Int_t k = 0; k < 1625;k++){
		for (Int_t l = 0; l < 650;l++){
			Double_t xal = h2_xal_vs_DTraw_sec1[i]->GetYaxis()->GetBinCenter(l);
			Double_t dt_raw = h2_xal_vs_DTraw_sec1[i]->GetXaxis()->GetBinCenter(k);
			Double_t  bin_cont = h2_xal_vs_DTraw_sec1[i]->GetBinContent(k,l);
			if ( xal < dt_raw*slope+offset+5 && xal > dt_raw*slope+offset-5){
				TAxis *xaxis = h2_xal_vs_DTraw_sec1_after_fit[i]->GetXaxis();
				TAxis *yaxis = h2_xal_vs_DTraw_sec1_after_fit[i]->GetYaxis();
				Int_t binx = xaxis->FindBin(dt_raw);
				Int_t biny = yaxis->FindBin(xal);
				h2_xal_vs_DTraw_sec1_after_fit[i]->SetBinContent(binx,biny,bin_cont);
				}	
			}	
		}
	delete fit_example_sec1;
	h2_xal_vs_DTraw_sec1_after_fit[i]->Draw("colz");
	TF1* final_fit_sec1 = new TF1("final_fit_sec1","[0]*x+ [1]",mean -2*st_dev,mean+2*st_dev);
	final_fit_sec1->SetLineWidth(1);
	h2_xal_vs_DTraw_sec1_after_fit[i]->Fit(final_fit_sec1);
	final_fit_sec1->Draw("same");
	par_file << 1 << "," << i << "," << slope << "," << offset << endl;
}
sec1_c.Modified();
sec1_c.Update();
sec1_c.Print("section_1_drift_fit.pdf");
sec1_c.Write();


//section 2
TCanvas sec2_c("sec2_c","Section2: Xal vs DTraw");
sec2_c.Divide(4,4);
for (Int_t i = 0; i < 16; i++){
	sec2_c.cd(i+1);
	gPad->SetLogz(i+1);
	gStyle->SetOptStat(1);
	gStyle->SetOptFit(111);
	sprintf(hist_name,"Section 2, Anode %i -> Xal vs DTraw",i);
	h2_xal_vs_DTraw_sec2[i] = (TH2D*)MyHistList0->FindObject(hist_name);
	Double_t mean = h2_xal_vs_DTraw_sec2[i]->GetMean(1);
	Double_t st_dev = h2_xal_vs_DTraw_sec2[i]->GetStdDev(1);
	TF1* fit_example_sec2 = new TF1("fit_example_sec2","[0]*x+ [1]",mean-1.2*st_dev,mean+1.2*st_dev);
	h2_xal_vs_DTraw_sec2[i]->Fit(fit_example_sec2, "R");
	Double_t slope = fit_example_sec2->GetParameter(0);
	Double_t offset = fit_example_sec2->GetParameter(1);
	//h2_xal_vs_DTraw_sec2[i]->Draw("colz");
	//fit_example_sec2->Draw("same");
	//use fit parameters to make new fit
	for (Int_t k = 0; k < 1625;k++){
		for (Int_t l = 0; l < 650;l++){
			Double_t xal = h2_xal_vs_DTraw_sec2[i]->GetYaxis()->GetBinCenter(l);
			Double_t dt_raw = h2_xal_vs_DTraw_sec2[i]->GetXaxis()->GetBinCenter(k);
			Double_t  bin_cont = h2_xal_vs_DTraw_sec2[i]->GetBinContent(k,l);
			if ( xal < dt_raw*slope+offset+5 && xal > dt_raw*slope+offset-5){
				TAxis *xaxis = h2_xal_vs_DTraw_sec2_after_fit[i]->GetXaxis();
				TAxis *yaxis = h2_xal_vs_DTraw_sec2_after_fit[i]->GetYaxis();
				Int_t binx = xaxis->FindBin(dt_raw);
				Int_t biny = yaxis->FindBin(xal);
				h2_xal_vs_DTraw_sec2_after_fit[i]->SetBinContent(binx,biny,bin_cont);
				}	
			}	
		}
	delete fit_example_sec2;
	h2_xal_vs_DTraw_sec2_after_fit[i]->Draw("colz");
	TF1* final_fit_sec2 = new TF1("final_fit_sec2","[0]*x+ [1]",mean -2*st_dev,mean+2*st_dev);
	final_fit_sec2->SetLineWidth(1);
	h2_xal_vs_DTraw_sec2_after_fit[i]->Fit(final_fit_sec2);
	final_fit_sec2->Draw("same");
	par_file << 2 << "," << i << "," << slope << "," << offset << endl;
}
sec2_c.Modified();
sec2_c.Update();
sec2_c.Print("section_2_drift_fit.pdf");
sec2_c.Write();


//section 3
TCanvas sec3_c("sec3_c","Section3: Xal vs DTraw");
sec3_c.Divide(4,4);
for (Int_t i = 0; i < 16; i++){
	sec3_c.cd(i+1);
	gPad->SetLogz(i+1);
	gStyle->SetOptStat(1);
	gStyle->SetOptFit(111);
	sprintf(hist_name,"Section 3, Anode %i -> Xal vs DTraw",i);
	h2_xal_vs_DTraw_sec3[i] = (TH2D*)MyHistList0->FindObject(hist_name);
	Double_t mean = h2_xal_vs_DTraw_sec3[i]->GetMean(1);
	Double_t st_dev = h2_xal_vs_DTraw_sec3[i]->GetStdDev(1);
	TF1* fit_example_sec3 = new TF1("fit_example_sec3","[0]*x+ [1]",mean-1.2*st_dev,mean+1.2*st_dev);
	h2_xal_vs_DTraw_sec3[i]->Fit(fit_example_sec3, "R");
	Double_t slope = fit_example_sec3->GetParameter(0);
	Double_t offset = fit_example_sec3->GetParameter(1);
	//h2_xal_vs_DTraw_sec3[i]->Draw("colz");
	//fit_example_sec3->Draw("same");
	//use fit parameters to make new fit
	for (Int_t k = 0; k < 1625;k++){
		for (Int_t l = 0; l < 650;l++){
			Double_t xal = h2_xal_vs_DTraw_sec3[i]->GetYaxis()->GetBinCenter(l);
			Double_t dt_raw = h2_xal_vs_DTraw_sec3[i]->GetXaxis()->GetBinCenter(k);
			Double_t  bin_cont = h2_xal_vs_DTraw_sec3[i]->GetBinContent(k,l);
			if ( xal < dt_raw*slope+offset+5 && xal > dt_raw*slope+offset-5){
				TAxis *xaxis = h2_xal_vs_DTraw_sec3_after_fit[i]->GetXaxis();
				TAxis *yaxis = h2_xal_vs_DTraw_sec3_after_fit[i]->GetYaxis();
				Int_t binx = xaxis->FindBin(dt_raw);
				Int_t biny = yaxis->FindBin(xal);
				h2_xal_vs_DTraw_sec3_after_fit[i]->SetBinContent(binx,biny,bin_cont);
				}	
			}	
		}
	delete fit_example_sec3;
	h2_xal_vs_DTraw_sec3_after_fit[i]->Draw("colz");
	TF1* final_fit_sec3 = new TF1("final_fit_sec3","[0]*x+ [1]",mean -2*st_dev,mean+2*st_dev);
	final_fit_sec3->SetLineWidth(1);
	h2_xal_vs_DTraw_sec3_after_fit[i]->Fit(final_fit_sec3);
	final_fit_sec3->Draw("same");
	par_file << 3 << "," << i << "," << slope << "," << offset << endl;
}
sec3_c.Modified();
sec3_c.Update();
sec3_c.Print("section_3_drift_fit.pdf");
sec3_c.Write();

par_file.close();
}
