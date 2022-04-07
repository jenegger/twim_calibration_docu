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
void retrieve_fits_hist(const string& input_str){
string fname = string(input_str);
const char* char_fname= fname.c_str();
char hist_name[500];
vector<vector<double> > vec_para;
//Define histograms to fit:
TH2D* h2_sec0_energy_anode[16];
TH2D* h2_sec1_energy_anode[16];
TH2D* h2_sec2_energy_anode[16];
TH2D* h2_sec3_energy_anode[16];

TH2D* h2_test;
sprintf(hist_name,"histogram with cut...");
h2_test = new TH2D(hist_name,hist_name,1000,0,50000,1000,0,50000);
h2_test->GetXaxis()->SetTitle("Energy form anode 5");
h2_test->GetYaxis()->SetTitle("Energy from last anode");
h2_test->GetXaxis()->CenterTitle(true);
h2_test->GetYaxis()->CenterTitle(true);
h2_test->GetXaxis()->SetLabelSize(0.045);
h2_test->GetYaxis()->SetTitleSize(0.045);


TH2D* h2_test2;
sprintf(hist_name,"only with tcutg");
h2_test2 = new TH2D(hist_name,hist_name,1000,0,50000,1000,0,50000);
h2_test2->GetXaxis()->SetTitle("Energy form anode 5");
h2_test2->GetYaxis()->SetTitle("Energy from last anode");
h2_test2->GetXaxis()->CenterTitle(true);
h2_test2->GetYaxis()->CenterTitle(true);
h2_test2->GetXaxis()->SetLabelSize(0.045);
h2_test2->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_twim_sec0_fit[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name, "Section0: Fit Energy AnodeID %i vs Energy ref. Anode 5",i);
h2_twim_sec0_fit[i] = new TH2D(hist_name,hist_name,1000,0,50000,1000,0,50000);
h2_twim_sec0_fit[i]->GetXaxis()->SetTitle("Energy AnodeID 5");
h2_twim_sec0_fit[i]->GetYaxis()->SetTitle("Energy AnodeID %i");
h2_twim_sec0_fit[i]->GetXaxis()->CenterTitle(true);
h2_twim_sec0_fit[i]->GetYaxis()->CenterTitle(true);
h2_twim_sec0_fit[i]->GetYaxis()->SetLabelSize(0.045);
h2_twim_sec0_fit[i]->GetYaxis()->SetTitleSize(0.045);
}

TH2D* h2_twim_sec1_fit[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name, "Section1: Fit Energy AnodeID %i vs Energy ref. Anode 5",i);
h2_twim_sec1_fit[i] = new TH2D(hist_name,hist_name,1000,0,50000,1000,0,50000);
h2_twim_sec1_fit[i]->GetXaxis()->SetTitle("Energy AnodeID 5");
h2_twim_sec1_fit[i]->GetYaxis()->SetTitle("Energy AnodeID %i");
h2_twim_sec1_fit[i]->GetXaxis()->CenterTitle(true);
h2_twim_sec1_fit[i]->GetYaxis()->CenterTitle(true);
h2_twim_sec1_fit[i]->GetYaxis()->SetLabelSize(0.045);
h2_twim_sec1_fit[i]->GetYaxis()->SetTitleSize(0.045);
}

TH2D* h2_twim_sec2_fit[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name, "Section2: Fit Energy AnodeID %i vs Energy ref. Anode 5",i);
h2_twim_sec2_fit[i] = new TH2D(hist_name,hist_name,1000,0,50000,1000,0,50000);
h2_twim_sec2_fit[i]->GetXaxis()->SetTitle("Energy AnodeID 5");
h2_twim_sec2_fit[i]->GetYaxis()->SetTitle("Energy AnodeID %i");
h2_twim_sec2_fit[i]->GetXaxis()->CenterTitle(true);
h2_twim_sec2_fit[i]->GetYaxis()->CenterTitle(true);
h2_twim_sec2_fit[i]->GetYaxis()->SetLabelSize(0.045);
h2_twim_sec2_fit[i]->GetYaxis()->SetTitleSize(0.045);
}

TH2D* h2_twim_sec3_fit[16];
for (Int_t i = 0; i < 16; i++){
sprintf(hist_name, "Section3: Fit Energy AnodeID %i vs Energy ref. Anode 5",i);
h2_twim_sec3_fit[i] = new TH2D(hist_name,hist_name,1000,0,50000,1000,0,50000);
h2_twim_sec3_fit[i]->GetXaxis()->SetTitle("Energy AnodeID 5");
h2_twim_sec3_fit[i]->GetYaxis()->SetTitle("Energy AnodeID %i");
h2_twim_sec3_fit[i]->GetXaxis()->CenterTitle(true);
h2_twim_sec3_fit[i]->GetYaxis()->CenterTitle(true);
h2_twim_sec3_fit[i]->GetYaxis()->SetLabelSize(0.045);
h2_twim_sec3_fit[i]->GetYaxis()->SetTitleSize(0.045);
}

//file to store alignment parameters
#include <ctime>
time_t *rawtime = new time_t;
struct tm * timeinfo;
time(rawtime);
timeinfo = localtime(rawtime); 
ofstream par_file;
par_file.open("parameters_twim_anodes.csv");
par_file << "#FILE created by retrieve_fits_hist.C at time:" << timeinfo->tm_hour << ":" << timeinfo->tm_min << "Date:\t" << timeinfo->tm_mday << "/" << timeinfo->tm_mon << "/" << timeinfo->tm_year << endl;



TFile* file_input(TFile::Open(char_fname ,"READ"));
file_input->ls();
vector<double>v_params_temp(2);
	
char f_out_name[500];
sprintf(f_out_name,"final_fit_plots.root");
TFile * f = new TFile(f_out_name,"RECREATE");

TCutG *gcut1 = new TCutG("gcut1",7);

//section 0
TList* MyHistList0= (TList*)file_input->Get("histlist");
TCanvas sec0_c("sec0_c","Section0: Energy AnodeID[i] vs ref Energy AnodeID = 5");
sec0_c.Divide(4,4);
for (Int_t i = 0; i < 16; i++){
	sec0_c.cd(i+1);
	gPad->SetLogz(i+1);
	gStyle->SetOptStat(1);
	gStyle->SetOptFit(111);
	sprintf(hist_name,"Section0: Energy AnodeID %i vs Energy ref. Anode 5",i); 
	h2_sec0_energy_anode[i] = (TH2D*)MyHistList0->FindObject(hist_name);
	Double_t mean = h2_sec0_energy_anode[i]->GetMean(1);
	Double_t st_dev = h2_sec0_energy_anode[i]->GetStdDev(1);
	TF1* fit_example_sec0 = new TF1("fit_example_sec0","[0]*x+ [1]",mean-st_dev,mean+st_dev);
	h2_sec0_energy_anode[i]->Fit(fit_example_sec0, "R");
	Double_t slope = fit_example_sec0->GetParameter(0);
	Double_t offset = fit_example_sec0->GetParameter(1);
	gcut1->SetVarX("x");
	gcut1->SetVarY("y");
	gcut1->SetPoint(0,0.,0.);
	gcut1->SetPoint(1,((5000-offset)/slope),0.);
	gcut1->SetPoint(2,50000,slope*50000+offset-5000);
	gcut1->SetPoint(3,50000.,50000.);
	gcut1->SetPoint(4,((45000-offset)/slope),50000.);
	gcut1->SetPoint(5,0,offset+5000);
	gcut1->SetPoint(6,0.,0.);
	v_params_temp[0] = slope;
	v_params_temp[1] = offset;
	vec_para.push_back(v_params_temp);
	delete fit_example_sec0;
	//fill new histo with cutg:
	for(Int_t k = 1; k < 1000; k++){
		for(Int_t l = 1; l < 1000; l++){
			Double_t y_val = h2_sec0_energy_anode[i]->GetYaxis()->GetBinCenter(l);
			Double_t x_val = h2_sec0_energy_anode[i]->GetXaxis()->GetBinCenter(k);
			if (gcut1->IsInside(x_val,y_val)){
				Double_t  bin_cont = h2_sec0_energy_anode[i]->GetBinContent(k,l);
				h2_twim_sec0_fit[i]->SetBinContent(k,l,bin_cont);
				}
			else
				continue;
			}
		}
	Double_t mean_final = h2_twim_sec0_fit[i]->GetMean(1);
	Double_t st_dev_final =h2_twim_sec0_fit[i]->GetStdDev(1);
	TF1* fit_example_single = new TF1("fit_example_single","[0]*x+ [1]",mean_final-st_dev_final,mean_final+st_dev_final);
	h2_twim_sec0_fit[i]->Fit(fit_example_single, "R");
	h2_twim_sec0_fit[i]->Draw("colz");
	fit_example_single->Draw("same");
	//get the fit parameters
	Double_t par0 = fit_example_single->GetParameter(0);
	Double_t par1 = fit_example_single->GetParameter(1);
	//fill the parameters into file, first row = section, second row= anode, third row = par0, fourth row = par1
	par_file << 0 << "," << i << "," << par0 << "," << par1 << endl;
	
}
sec0_c.Modified();
sec0_c.Update();
sec0_c.Print("section_0_channels.pdf");
sec0_c.Write();


//section1
TCanvas sec1_c("sec1_c","Section1: Energy AnodeID[i] vs ref Energy AnodeID = 5");
sec1_c.Divide(4,4);
for (Int_t i = 0; i < 16; i++){
	sec1_c.cd(i+1);
	gPad->SetLogz(i+1);
	gStyle->SetOptStat(1);
	gStyle->SetOptFit(111);
	sprintf(hist_name,"Section1: Energy AnodeID %i vs Energy ref. Anode 5",i); 
	h2_sec1_energy_anode[i] = (TH2D*)MyHistList0->FindObject(hist_name);
	Double_t mean = h2_sec1_energy_anode[i]->GetMean(1);
	Double_t st_dev = h2_sec1_energy_anode[i]->GetStdDev(1);
	TF1* fit_example_sec1 = new TF1("fit_example_sec1","[0]*x+ [1]",mean-st_dev,mean+st_dev);
	h2_sec1_energy_anode[i]->Fit(fit_example_sec1, "R");
	Double_t slope = fit_example_sec1->GetParameter(0);
	Double_t offset = fit_example_sec1->GetParameter(1);
	gcut1->SetVarX("x");
	gcut1->SetVarY("y");
	gcut1->SetPoint(0,0.,0.);
	gcut1->SetPoint(1,((5000-offset)/slope),0.);
	gcut1->SetPoint(2,50000,slope*50000+offset-5000);
	gcut1->SetPoint(3,50000.,50000.);
	gcut1->SetPoint(4,((45000-offset)/slope),50000.);
	gcut1->SetPoint(5,0,offset+5000);
	gcut1->SetPoint(6,0.,0.);
	v_params_temp[0] = slope;
	v_params_temp[1] = offset;
	vec_para.push_back(v_params_temp);
	delete fit_example_sec1;
	//fill new histo with cutg:
	for(Int_t k = 1; k < 1000; k++){
		for(Int_t l = 1; l < 1000; l++){
			Double_t y_val = h2_sec1_energy_anode[i]->GetYaxis()->GetBinCenter(l);
			Double_t x_val = h2_sec1_energy_anode[i]->GetXaxis()->GetBinCenter(k);
			if (gcut1->IsInside(x_val,y_val)){
				Double_t  bin_cont = h2_sec1_energy_anode[i]->GetBinContent(k,l);
				h2_twim_sec1_fit[i]->SetBinContent(k,l,bin_cont);
				}
			else
				continue;
			}
		}
	Double_t mean_final = h2_twim_sec1_fit[i]->GetMean(1);
	Double_t st_dev_final =h2_twim_sec1_fit[i]->GetStdDev(1);
	TF1* fit_example_single = new TF1("fit_example_single","[0]*x+ [1]",mean_final-st_dev_final,mean_final+st_dev_final);
	h2_twim_sec1_fit[i]->Fit(fit_example_single, "R");
	h2_twim_sec1_fit[i]->Draw("colz");
	fit_example_single->Draw("same");
	//get the fit parameters
	Double_t par0 = fit_example_single->GetParameter(0);
	Double_t par1 = fit_example_single->GetParameter(1);
	//fill the parameters into file, first row = section, second row= anode, third row = par0, fourth row = par1
	par_file << 1 << "," << i << "," << par0 << "," << par1 << endl;
}
sec1_c.Modified();
sec1_c.Update();
sec1_c.Print("section_1_channels.pdf");
sec1_c.Write();

//section2

TCanvas sec2_c("sec2_c","Section2: Energy AnodeID[i] vs ref Energy AnodeID = 5");
sec2_c.Divide(4,4);
for (Int_t i = 0; i < 16; i++){
	sec2_c.cd(i+1);
	gPad->SetLogz(i+1);
	gStyle->SetOptStat(1);
	gStyle->SetOptFit(111);
	sprintf(hist_name,"Section2: Energy AnodeID %i vs Energy ref. Anode 5",i); 
	h2_sec2_energy_anode[i] = (TH2D*)MyHistList0->FindObject(hist_name);
	Double_t mean = h2_sec2_energy_anode[i]->GetMean(1);
	Double_t st_dev = h2_sec2_energy_anode[i]->GetStdDev(1);
	TF1* fit_example_sec2 = new TF1("fit_example_sec2","[0]*x+ [1]",mean-st_dev,mean+st_dev);
	h2_sec2_energy_anode[i]->Fit(fit_example_sec2, "R");
	Double_t slope = fit_example_sec2->GetParameter(0);
	Double_t offset = fit_example_sec2->GetParameter(1);
	gcut1->SetVarX("x");
	gcut1->SetVarY("y");
	gcut1->SetPoint(0,0.,0.);
	gcut1->SetPoint(1,((5000-offset)/slope),0.);
	gcut1->SetPoint(2,50000,slope*50000+offset-5000);
	gcut1->SetPoint(3,50000.,50000.);
	gcut1->SetPoint(4,((45000-offset)/slope),50000.);
	gcut1->SetPoint(5,0,offset+5000);
	gcut1->SetPoint(6,0.,0.);
	v_params_temp[0] = slope;
	v_params_temp[1] = offset;
	vec_para.push_back(v_params_temp);
	delete fit_example_sec2;
	//fill new histo with cutg:
	for(Int_t k = 1; k < 1000; k++){
		for(Int_t l = 1; l < 1000; l++){
			Double_t y_val = h2_sec2_energy_anode[i]->GetYaxis()->GetBinCenter(l);
			Double_t x_val = h2_sec2_energy_anode[i]->GetXaxis()->GetBinCenter(k);
			if (gcut1->IsInside(x_val,y_val)){
				Double_t  bin_cont = h2_sec2_energy_anode[i]->GetBinContent(k,l);
				h2_twim_sec2_fit[i]->SetBinContent(k,l,bin_cont);
				}
			else
				continue;
			}
		}
	Double_t mean_final = h2_twim_sec2_fit[i]->GetMean(1);
	Double_t st_dev_final =h2_twim_sec2_fit[i]->GetStdDev(1);
	TF1* fit_example_single = new TF1("fit_example_single","[0]*x+ [1]",mean_final-st_dev_final,mean_final+st_dev_final);
	h2_twim_sec2_fit[i]->Fit(fit_example_single, "R");
	h2_twim_sec2_fit[i]->Draw("colz");
	fit_example_single->Draw("same");
	//get the fit parameters
	Double_t par0 = fit_example_single->GetParameter(0);
	Double_t par1 = fit_example_single->GetParameter(1);
	//fill the parameters into file, first row = section, second row= anode, third row = par0, fourth row = par1
	par_file << 2 << "," << i << "," << par0 << "," << par1 << endl;
}
sec2_c.Modified();
sec2_c.Update();
sec2_c.Print("section_2_channels.pdf");
sec2_c.Write();

//section3

TCanvas sec3_c("sec3_c","Section3: Energy AnodeID[i] vs ref Energy AnodeID = 5");
sec3_c.Divide(4,4);
for (Int_t i = 0; i < 16; i++){
	sec3_c.cd(i+1);
	gPad->SetLogz(i+1);
	gStyle->SetOptStat(1);
	gStyle->SetOptFit(111);
	sprintf(hist_name,"Section3: Energy AnodeID %i vs Energy ref. Anode 5",i); 
	h2_sec3_energy_anode[i] = (TH2D*)MyHistList0->FindObject(hist_name);
	Double_t mean = h2_sec3_energy_anode[i]->GetMean(1);
	Double_t st_dev = h2_sec3_energy_anode[i]->GetStdDev(1);
	TF1* fit_example_sec3 = new TF1("fit_example_sec3","[0]*x+ [1]",mean-st_dev,mean+st_dev);
	h2_sec3_energy_anode[i]->Fit(fit_example_sec3, "R");
	Double_t slope = fit_example_sec3->GetParameter(0);
	Double_t offset = fit_example_sec3->GetParameter(1);
	gcut1->SetVarX("x");
	gcut1->SetVarY("y");
	gcut1->SetPoint(0,0.,0.);
	gcut1->SetPoint(1,((5000-offset)/slope),0.);
	gcut1->SetPoint(2,50000,slope*50000+offset-5000);
	gcut1->SetPoint(3,50000.,50000.);
	gcut1->SetPoint(4,((45000-offset)/slope),50000.);
	gcut1->SetPoint(5,0,offset+5000);
	gcut1->SetPoint(6,0.,0.);
	v_params_temp[0] = slope;
	v_params_temp[1] = offset;
	vec_para.push_back(v_params_temp);
	delete fit_example_sec3;
	//fill new histo with cutg:
	for(Int_t k = 1; k < 1000; k++){
		for(Int_t l = 1; l < 1000; l++){
			Double_t y_val = h2_sec3_energy_anode[i]->GetYaxis()->GetBinCenter(l);
			Double_t x_val = h2_sec3_energy_anode[i]->GetXaxis()->GetBinCenter(k);
			if (gcut1->IsInside(x_val,y_val)){
				Double_t  bin_cont = h2_sec3_energy_anode[i]->GetBinContent(k,l);
				h2_twim_sec3_fit[i]->SetBinContent(k,l,bin_cont);
				}
			else
				continue;
			}
		}
	Double_t mean_final = h2_twim_sec3_fit[i]->GetMean(1);
	Double_t st_dev_final =h2_twim_sec3_fit[i]->GetStdDev(1);
	TF1* fit_example_single = new TF1("fit_example_single","[0]*x+ [1]",mean_final-st_dev_final,mean_final+st_dev_final);
	h2_twim_sec3_fit[i]->Fit(fit_example_single, "R");
	h2_twim_sec3_fit[i]->Draw("colz");
	fit_example_single->Draw("same");
	//get the fit parameters
	Double_t par0 = fit_example_single->GetParameter(0);
	Double_t par1 = fit_example_single->GetParameter(1);
	//fill the parameters into file, first row = section, second row= anode, third row = par0, fourth row = par1
	par_file << 3 << "," << i << "," << par0 << "," << par1 << endl;
}
sec3_c.Modified();
sec3_c.Update();
sec3_c.Print("section_3_channels.pdf");
sec3_c.Write();

par_file.close();
}
