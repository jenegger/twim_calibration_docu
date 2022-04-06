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

using namespace std;

void resolution_plots(const string& input_str){
fstream fin;
fin.open("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration/xcal_std_dev_params.csv", ios::in);
string line, word;
vector<vector<double> > v_sigma(4,vector<double>(16));
getline(fin, line);
while(fin.peek()!=EOF) {
	getline(fin, line);
	stringstream s(line);
	vector<string> temp_vec;
	while (getline(s, word, ',')) {     temp_vec.push_back(word);
	}
	v_sigma[stoi(temp_vec[0])][stoi(temp_vec[1])] = stod(temp_vec[2]);
} 

Double_t res_0[16];
Double_t res_1[16];
Double_t res_2[16];
Double_t res_3[16];
Double_t counting[16];
Double_t const_val[16];
for (Int_t i = 0; i < 16; i++){
	res_0[i] = v_sigma[0][i];
	res_1[i] = v_sigma[1][i];
	res_2[i] = v_sigma[2][i];
	res_3[i] = v_sigma[3][i];
	counting[i] = i;
	const_val[i] = 0.03;
	}
//create graphs
Int_t n = 16;
cout << "res,sec0,anode5:\t" << res_0[4] << endl;
TGraph *gr_sec0 = new TGraph(n,counting,res_0);
gr_sec0->SetTitle("Section 0, Sigma Xfit-Xcorrected ");
TGraph *gr_sec1 = new TGraph(n,counting,res_1);
gr_sec1->SetTitle("Section 1, Sigma Xfit-Xcorrected ");
TGraph *gr_sec2 = new TGraph(n,counting,res_2);
gr_sec2->SetTitle("Section 2, Sigma Xfit-Xcorrected ");
TGraph *gr_sec3 = new TGraph(n,counting,res_3);
gr_sec3->SetTitle("Section 3, Sigma Xfit-Xcorrected ");

TGraph* gr_straight_line = new TGraph(n,counting,const_val);
gr_straight_line->SetLineColor(2);

TCanvas canvas1("canvas1","Sigma Xfit-Xcorrected ");
canvas1.Divide(2,2);
canvas1.cd(1);
gr_sec0->Draw("AP*");
gr_straight_line->Draw("C");
canvas1.cd(2);
gr_sec1->Draw("AP*");
gr_straight_line->Draw("C");
canvas1.cd(3);
gr_sec2->Draw("AP*");
gr_straight_line->Draw("C");
canvas1.cd(4);
gr_sec3->Draw("AP*");
gr_straight_line->Draw("C");
canvas1.Modified();
canvas1.Update();
canvas1.Print("resolution_xfit_x_corr.pdf");
}


