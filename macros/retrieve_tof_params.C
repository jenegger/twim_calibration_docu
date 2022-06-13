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
void retrieve_tof_params(const string& input_str){
string fname = string(input_str);
const char* char_fname= fname.c_str();
TFile* file_input(TFile::Open(char_fname ,"READ"));
TList* MyHistList0= (TList*)file_input->Get("histlist");
char f_out_name[500];
sprintf(f_out_name,"fit_y_tofw.root");
TFile * f = new TFile(f_out_name,"RECREATE");
char hist_name[500];
TH2D* h2_tof_vs_y_pos[29];
TH2D* h2_tof_vs_y_pos_fit[29];
//for (Int_t i = 0; i < 29; i++){
//	sprintf(hist_name,"TOFW pad %i tof_ns vs extrapolated y-pos.",i);
//	h2_tof_vs_y_pos[i] = new TH2D(hist_name,hist_name,500,-50,50,600,-300,300);
//	h2_tof_vs_y_pos[i]->GetXaxis()->SetTitle("tof ns");
//	h2_tof_vs_y_pos[i]->GetYaxis()->SetTitle("Extrap. y-pos [mm]");
//	h2_tof_vs_y_pos[i]->GetXaxis()->CenterTitle(true);
//	h2_tof_vs_y_pos[i]->GetYaxis()->CenterTitle(true);
//	h2_tof_vs_y_pos[i]->GetYaxis()->SetLabelSize(0.045);
//	h2_tof_vs_y_pos[i]->GetYaxis()->SetTitleSize(0.045);
//}
//file to store  parameters
#include <ctime>
time_t *rawtime = new time_t;
struct tm * timeinfo;
time(rawtime);
timeinfo = localtime(rawtime); 
ofstream par_file;
par_file.open("tof_params.csv");
par_file << "#FILE created by retrieve_tof_params  at time:" << timeinfo->tm_hour << ":" << timeinfo->tm_min << "Date:\t" << timeinfo->tm_mday << "/" << timeinfo->tm_mon << "/" << timeinfo->tm_year << endl;
for (Int_t i = 0; i < 29; i++){
	sprintf(hist_name,"FIT TOFW pad %i tof_ns vs extrapolated y-pos.",i);
	h2_tof_vs_y_pos_fit[i] = new TH2D(hist_name,hist_name,500,-50,50,600,-300,300);
	h2_tof_vs_y_pos_fit[i]->GetXaxis()->SetTitle("tof ns");
	h2_tof_vs_y_pos_fit[i]->GetYaxis()->SetTitle("Extrap. y-pos [mm]");
	h2_tof_vs_y_pos_fit[i]->GetXaxis()->CenterTitle(true);
	h2_tof_vs_y_pos_fit[i]->GetYaxis()->CenterTitle(true);
	h2_tof_vs_y_pos_fit[i]->GetYaxis()->SetLabelSize(0.045);
	h2_tof_vs_y_pos_fit[i]->GetYaxis()->SetTitleSize(0.045);
}
TCanvas tof_canvas("tof_canvas","TOF_pad[i] vs extrapolated y-pos");
tof_canvas.Divide(5,6);
for (Int_t i = 0; i < 29; i++){
	tof_canvas.cd(i+1);
	gPad->SetLogz(i+1);
	gStyle->SetOptStat(1);
	gStyle->SetOptFit(111);
	gStyle->SetLineWidth(1);
	sprintf(hist_name,"TOFW pad %i tof_ns vs extrapolated y-pos.",i);
	h2_tof_vs_y_pos[i] = (TH2D*)MyHistList0->FindObject(hist_name);
	TF1* fit_tof_pads = new TF1("fit_tof_pads","[0]*x+ [1]");
	h2_tof_vs_y_pos[i]->Fit(fit_tof_pads);
	Double_t slope = fit_tof_pads->GetParameter(0);
	Double_t offset = fit_tof_pads->GetParameter(1);
	for(Int_t k = 1; k < 500; k++){
		for(Int_t l = 1; l < 600; l++){
			Double_t y_val = h2_tof_vs_y_pos[i]->GetYaxis()->GetBinCenter(l);
			Double_t x_val = h2_tof_vs_y_pos[i]->GetXaxis()->GetBinCenter(k);
			if (abs((x_val*slope+offset)-y_val) < 40){
				Double_t  bin_cont = h2_tof_vs_y_pos[i]->GetBinContent(k,l);
				h2_tof_vs_y_pos_fit[i]->SetBinContent(k,l,bin_cont);
				}
			else continue;
			}
		}
	TF1* final_fit = new TF1("final_fit","[0]*x+ [1]");
	h2_tof_vs_y_pos_fit[i]->Fit(final_fit);
	h2_tof_vs_y_pos_fit[i]->Draw("colz");
	final_fit->Draw("same");
	Double_t par0 = final_fit->GetParameter(0);
	Double_t par1 = final_fit->GetParameter(1);
	par_file << i << "," << par0 << "," << par1 << endl;
	}
	tof_canvas.Modified();
	tof_canvas.Update();
	tof_canvas.Print("tof_calibration_y_pos.pdf");
	tof_canvas.Write();
	par_file.close();
}
