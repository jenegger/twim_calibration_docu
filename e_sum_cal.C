//how to find peak of 1D histo:

#include "TH1F.h"
#include "TH2F.h"
#include "TClonesArray.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include<bits/stdc++.h>
using namespace std;

void e_sum_cal(const string& input_str){
string fname = string(input_str);
const char* char_fname= fname.c_str();
char hist_name[500];
#include <ctime>
time_t *rawtime = new time_t;
struct tm * timeinfo;
time(rawtime);
timeinfo = localtime(rawtime);
ofstream par_file;
par_file.open("sum_anodes_parameters.csv");
par_file << "#FILE created by e_sum_cal.C at time:" << timeinfo->tm_hour << ":" << timeinfo->tm_min << "Date:\t" << timeinfo->tm_mday << "/" << timeinfo->tm_mon << "/" << timeinfo->tm_year << endl;

TFile* file_input(TFile::Open(char_fname ,"READ"));
TList* MyHistList0= (TList*)file_input->Get("histlist");

TH1D* h1_e_sum_0;
sprintf(hist_name,"TWIM Energy Section 0 peaks");
h1_e_sum_0 = new TH1D(hist_name,hist_name,10000,0,50000);
h1_e_sum_0->GetXaxis()->SetTitle("TWIM Energy");
h1_e_sum_0->GetYaxis()->SetTitle("Counts");
h1_e_sum_0->GetXaxis()->CenterTitle(true);
h1_e_sum_0->GetYaxis()->CenterTitle(true);
h1_e_sum_0->GetYaxis()->SetLabelSize(0.045);
h1_e_sum_0->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_e_sum_1;
sprintf(hist_name,"TWIM Energy Section 1 peaks");
h1_e_sum_1 = new TH1D(hist_name,hist_name,10000,0,50000);
h1_e_sum_1->GetXaxis()->SetTitle("TWIM Energy");
h1_e_sum_1->GetYaxis()->SetTitle("Counts");
h1_e_sum_1->GetXaxis()->CenterTitle(true);
h1_e_sum_1->GetYaxis()->CenterTitle(true);
h1_e_sum_1->GetYaxis()->SetLabelSize(0.045);
h1_e_sum_1->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_e_sum_2;
sprintf(hist_name,"TWIM Energy Section 2 peaks");
h1_e_sum_2 = new TH1D(hist_name,hist_name,10000,0,50000);
h1_e_sum_2->GetXaxis()->SetTitle("TWIM Energy");
h1_e_sum_2->GetYaxis()->SetTitle("Counts");
h1_e_sum_2->GetXaxis()->CenterTitle(true);
h1_e_sum_2->GetYaxis()->CenterTitle(true);
h1_e_sum_2->GetYaxis()->SetLabelSize(0.045);
h1_e_sum_2->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_e_sum_3;
sprintf(hist_name,"TWIM Energy Section 3 peaks");
h1_e_sum_3 = new TH1D(hist_name,hist_name,10000,0,50000);
h1_e_sum_3->GetXaxis()->SetTitle("TWIM Energy");
h1_e_sum_3->GetYaxis()->SetTitle("Counts");
h1_e_sum_3->GetXaxis()->CenterTitle(true);
h1_e_sum_3->GetYaxis()->CenterTitle(true);
h1_e_sum_3->GetYaxis()->SetLabelSize(0.045);
h1_e_sum_3->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name,"TWIM Energy Section 0 238U");
h1_e_sum_0 = (TH1D*)MyHistList0->FindObject(hist_name); 
h1_e_sum_0->Rebin();
h1_e_sum_0->Rebin();
h1_e_sum_0->Rebin();

sprintf(hist_name,"TWIM Energy Section 1 238U");
h1_e_sum_1 = (TH1D*)MyHistList0->FindObject(hist_name); 
h1_e_sum_1->Rebin();
h1_e_sum_1->Rebin();
h1_e_sum_1->Rebin();

sprintf(hist_name,"TWIM Energy Section 2 238U");
h1_e_sum_2 = (TH1D*)MyHistList0->FindObject(hist_name); 
h1_e_sum_2->Rebin();
h1_e_sum_2->Rebin();
h1_e_sum_2->Rebin();

sprintf(hist_name,"TWIM Energy Section 3 238U");
h1_e_sum_3 = (TH1D*)MyHistList0->FindObject(hist_name); 
h1_e_sum_3->Rebin();
h1_e_sum_3->Rebin();
h1_e_sum_3->Rebin();


TSpectrum *s0 = new TSpectrum();
const int Npeaks0 = s0->Search(h1_e_sum_0);
TSpectrum *s1 = new TSpectrum();
const int Npeaks1 = s1->Search(h1_e_sum_1);
TSpectrum *s2 = new TSpectrum();
const int Npeaks2 = s2->Search(h1_e_sum_2);
TSpectrum *s3 = new TSpectrum();
const int Npeaks3 = s3->Search(h1_e_sum_3);


double *pos0 = s0->GetPositionX();
double *pos1 = s1->GetPositionX();
double *pos2 = s2->GetPositionX();
double *pos3 = s3->GetPositionX();

vector<double>v_peaks_0;
vector<double>v_peaks_1;
vector<double>v_peaks_2;
vector<double>v_peaks_3;

for (int i = 0; i < Npeaks0; ++i) {
	int bin = h1_e_sum_0->FindBin(pos0[i]);
	v_peaks_0.push_back(h1_e_sum_0->GetBinLowEdge(bin));
	}

for (int i = 0; i < Npeaks1; ++i) {
	int bin = h1_e_sum_1->FindBin(pos1[i]);
	v_peaks_1.push_back(h1_e_sum_1->GetBinLowEdge(bin));
	}

for (int i = 0; i < Npeaks2; ++i) {
	int bin = h1_e_sum_2->FindBin(pos2[i]);
	v_peaks_2.push_back(h1_e_sum_2->GetBinLowEdge(bin));
	}
for (int i = 0; i < Npeaks3; ++i) {
	int bin = h1_e_sum_3->FindBin(pos3[i]);
	v_peaks_3.push_back(h1_e_sum_3->GetBinLowEdge(bin));
	}

//h1_e_sum_0->Draw();
//h1_e_sum_2->Draw();
//h1_e_sum_1->Draw();
//h1_e_sum_0->Draw();
//double *pos = s->GetPositionX();
//cout << pos[4] << endl;
//cout << "... Number of peaks: " << Npeaks << endl;
//cout << "... Peak positions: " << endl;
//for (int i = 0; i < Npeaks; ++i) {
//	int bin = h1_e_sum_0->FindBin(pos[i]);
//	cout << h1_e_sum_0->GetBinLowEdge(bin);
//	if (i != Npeaks-1)  cout << ", ";
//	else  cout << endl;
//	}
sort(v_peaks_0.begin(),v_peaks_0.end());
sort(v_peaks_1.begin(),v_peaks_1.end());
sort(v_peaks_2.begin(),v_peaks_2.end());
sort(v_peaks_3.begin(),v_peaks_3.end());

cout << "this is pos of first peak 4thchamber:\t" << v_peaks_3[0] << endl;
cout << "this is pos of first peak 3rd chamber:\t" << v_peaks_2[0] << endl;
cout << "this is pos of first peak 2nd chamber:\t" << v_peaks_1[0] << endl;
cout << "this is pos of first peak 1st chamber:\t" << v_peaks_0[0] << endl;
cout << "diff first two peaks,chamber 1:\t" << v_peaks_0[1] - v_peaks_0[0] << endl;
cout << "diff first two peaks,chamber 2:\t" << v_peaks_1[1] - v_peaks_1[0] << endl;
cout << "diff first two peaks,chamber 3:\t" << v_peaks_2[1] - v_peaks_2[0] << endl;
cout << "diff first two peaks,chamber 4:\t" << v_peaks_3[1] - v_peaks_3[0] << endl;
Double_t distance_peaks = (v_peaks_0[1] - v_peaks_0[0] + v_peaks_1[1] - v_peaks_1[0] + v_peaks_2[1] - v_peaks_2[0] + v_peaks_3[1] - v_peaks_3[0])/4.;
vector<double> v_first_entries_sec{v_peaks_0[0],v_peaks_1[0],v_peaks_2[0],v_peaks_3[0]};
Int_t max_first_peak_sec = max_element(v_first_entries_sec.begin(),v_first_entries_sec.end()) - v_first_entries_sec.begin();
cout << "section with maximum peak:\t" << max_first_peak_sec << endl;
vector<vector<double> > v_all_sec;
v_all_sec.push_back(v_peaks_0);
v_all_sec.push_back(v_peaks_1);
v_all_sec.push_back(v_peaks_2);
v_all_sec.push_back(v_peaks_3);

vector<double> reduce_size{-1,-1,-1,-1};
for (int j = 0; j < 4; j++){
	bool peak_search = true;
	int p = 0;
	while (peak_search){
		Double_t diff_peaks = abs(v_all_sec[max_first_peak_sec][0] - v_all_sec[j][p]);
		if (diff_peaks < distance_peaks){
			peak_search = false;
			reduce_size[j] = p;	
			}		
		p+=1;
		}
	}
cout << "hello in front of first resizement" <<  endl;
//now resize all the vectors
vector<int> size_data_arr{-1,-1,-1,-1};
//section 0
for (int i = 0; i < reduce_size[0]; i++){
	v_peaks_0.erase(v_peaks_0.begin());
	}
size_data_arr[0] = v_peaks_0.size();
//section 1
for (int i = 0; i < reduce_size[1]; i++){
	v_peaks_1.erase(v_peaks_1.begin());
    }
size_data_arr[1] = v_peaks_1.size();
//section 2
for (int i = 0; i < reduce_size[2]; i++){
	v_peaks_2.erase(v_peaks_2.begin());
    }
size_data_arr[2] = v_peaks_2.size();
//section 3
for (int i = 0; i < reduce_size[3]; i++){
	v_peaks_3.erase(v_peaks_3.begin());
    }
size_data_arr[3] = v_peaks_3.size();

Int_t min_data_size = *min_element(size_data_arr.begin(),size_data_arr.end());
cout << "min element size:\t" <<  min_data_size << endl;

for (int i = 0; i < 4; i++){
	cout << "size of section" << i << "is:\t" << size_data_arr[i] << endl;
}
//final resizement
//section 0
cout << "first section size" << v_peaks_0.size() << endl;
cout << "size diff for first section :\t" << v_peaks_0.size()-min_data_size << endl;
if (v_peaks_0.size()-min_data_size > 0){
for (int i = 0; i < (v_peaks_0.size()-min_data_size); i++){
	v_peaks_0.erase(v_peaks_0.end()-1);
	}
	}
//section 1
if ((v_peaks_1.size()-min_data_size) > 0){
for (int i = 0; i < (v_peaks_1.size()-min_data_size); i++){
	v_peaks_1.erase(v_peaks_1.end()-1);
	}
	}
//section 2
if ((v_peaks_2.size()-min_data_size) > 0){
for (int i = 0; i < (v_peaks_2.size()-min_data_size); i++){
	v_peaks_2.erase(v_peaks_2.end()-1);
	}
	}
//section 3
if ((v_peaks_3.size()-min_data_size) > 0){
for (int i = 0; i < (v_peaks_3.size()-min_data_size); i++){
	v_peaks_3.erase(v_peaks_3.end()-1);
	}
	}



for (int i = 0; i < min_data_size; i++){
cout << "sec 0 " <<	v_peaks_0[i] << endl;
cout << "sec 1 " << v_peaks_1[i] << endl;
cout << "sec 2 " << v_peaks_2[i] << endl;
cout << "sec 3 " << v_peaks_3[i] << endl;
}
//just for test

TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,500,300);
c1->Divide(2,2);
//fit section 0
c1->cd(1);
Double_t x[100], y[100];
for (int j = 0; j < min_data_size; j++){
	x[j] = v_peaks_0[j];
	y[j] = v_peaks_0[j];	
}
TGraph* gr0 = new TGraph(min_data_size,x,y);
TF1 *fit0 = new TF1("fit0","[0]+x*[1]");
gr0->Fit("fit0");
gr0->Draw("AC*");
Double_t slope0 =  1./(fit0->GetParameter(1));
Double_t offset0 = -(fit0->GetParameter(0))/(fit0->GetParameter(1));
par_file << slope0 << "," << offset0 << endl;

//fit section 1
c1->cd(2);
for (int j = 0; j < min_data_size; j++){
	x[j] = v_peaks_0[j];
	y[j] = v_peaks_1[j];
}
TGraph* gr1 = new TGraph(min_data_size,x,y);
TF1 *fit1 = new TF1("fit1","[0]+x*[1]");
gr1->Fit("fit1");
gr1->Draw("AC*");
Double_t slope1 =  1./(fit1->GetParameter(1));
Double_t offset1 = -(fit1->GetParameter(0))/(fit1->GetParameter(1));
par_file << slope1 << "," << offset1 << endl;

//section 2
c1->cd(3);
for (int j = 0; j < min_data_size; j++){
	x[j] = v_peaks_0[j];
	y[j] = v_peaks_2[j];
}
TGraph* gr2 = new TGraph(min_data_size,x,y);
TF1 *fit2 = new TF1("fit2","[0]+x*[1]");
gr2->Fit("fit2");
gr2->Draw("AC*");
Double_t slope2 =  1./(fit2->GetParameter(1));
Double_t offset2 = -(fit2->GetParameter(0))/(fit2->GetParameter(1));
par_file << slope2 << "," << offset2 << endl;

//section 3
c1->cd(4);
for (int j = 0; j < min_data_size; j++){
	x[j] = v_peaks_0[j];
	y[j] = v_peaks_3[j];
}
TGraph* gr3 = new TGraph(min_data_size,x,y);
TF1 *fit3 = new TF1("fit3","[0]+x*[1]");
gr3->Fit("fit3");
gr3->Draw("AC*");
Double_t slope3 =  1./(fit3->GetParameter(1));
Double_t offset3 = -(fit3->GetParameter(0))/(fit3->GetParameter(1));
par_file << slope3 << "," << offset3 << endl;

cout << "some parameters" << fit3->GetParameter(0) << "   " << fit3->GetParameter(1) << endl;
par_file.close();

}

