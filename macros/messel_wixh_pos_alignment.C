//includes
#include "TH1F.h"
#include "TH2F.h"
#include "TClonesArray.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include<bits/stdc++.h>
using namespace std;
void messel_wixh_pos_alignment(const string& input_str){
string fname = string(input_str);
const char* char_fname= fname.c_str();
char hist_name[500];
#include <ctime>
time_t *rawtime = new time_t;
struct tm * timeinfo;
time(rawtime);
timeinfo = localtime(rawtime);
ofstream par_file;
par_file.open("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/parameters/energy_position_calib_params.csv");
par_file << "#FILE created by messel_wixh_alignment.C at time:" << timeinfo->tm_hour << ":" << timeinfo->tm_min << "Date:\t" << timeinfo->tm_mday << "/" << timeinfo->tm_mon << "/" << timeinfo->tm_year << endl;
TFile* file_input(TFile::Open(char_fname ,"READ"));
TList* MyHistList0= (TList*)file_input->Get("histlist");
TH1D* h1_e_sum_wix_sec2;
TH1D* h1_e_sum_wix_sec3;
TH1D* h1_e_sum_messel_sec0;
TH1D* h1_e_sum_messel_sec1;


sprintf(hist_name,"Energy deposited TWIM section2, beta and pos. corrected");
h1_e_sum_wix_sec2 = (TH1D*)MyHistList0->FindObject(hist_name);
h1_e_sum_wix_sec2->Rebin();
h1_e_sum_wix_sec2->Rebin();
h1_e_sum_wix_sec2->Rebin();
sprintf(hist_name,"Energy deposited TWIM section0, beta and pos. corrected");
h1_e_sum_messel_sec0 = (TH1D*)MyHistList0->FindObject(hist_name);
h1_e_sum_messel_sec0->Rebin();
h1_e_sum_messel_sec0->Rebin();
h1_e_sum_messel_sec0->Rebin();
sprintf(hist_name,"Energy deposited TWIM section3, beta and pos. corrected");
h1_e_sum_wix_sec3 = (TH1D*)MyHistList0->FindObject(hist_name);
h1_e_sum_wix_sec3->Rebin();
h1_e_sum_wix_sec3->Rebin();
h1_e_sum_wix_sec3->Rebin();
sprintf(hist_name,"Energy deposited TWIM section1, beta and pos. corrected");
h1_e_sum_messel_sec1 = (TH1D*)MyHistList0->FindObject(hist_name);
h1_e_sum_messel_sec1->Rebin();
h1_e_sum_messel_sec1->Rebin();
h1_e_sum_messel_sec1->Rebin();

TSpectrum *s_wix_sec2 = new TSpectrum();
const int Npeaks_wix_sec2 = s_wix_sec2->Search(h1_e_sum_wix_sec2);
TSpectrum *s_wix_sec3 = new TSpectrum();
const int Npeaks_wix_sec3 = s_wix_sec3->Search(h1_e_sum_wix_sec3);
TSpectrum *s_messel_sec0 = new TSpectrum();
const int Npeaks_messel_sec0 = s_messel_sec0->Search(h1_e_sum_messel_sec0);
TSpectrum *s_messel_sec1 = new TSpectrum();
const int Npeaks_messel_sec1 = s_messel_sec1->Search(h1_e_sum_messel_sec1);

double *pos_wix_sec2 = s_wix_sec2->GetPositionX();
double *pos_wix_sec3 = s_wix_sec3->GetPositionX();
double *pos_messel_sec0 = s_messel_sec0->GetPositionX();
double *pos_messel_sec1 = s_messel_sec1->GetPositionX();

vector<double>v_peaks_wix_sec2;
vector<double>v_peaks_wix_sec3;
vector<double>v_peaks_messel_sec0;
vector<double>v_peaks_messel_sec1;

for (int i = 0; i < Npeaks_wix_sec2; ++i) {
	int bin = h1_e_sum_wix_sec2->FindBin(pos_wix_sec2[i]);
	v_peaks_wix_sec2.push_back(h1_e_sum_wix_sec2->GetBinLowEdge(bin));
	}
for (int i = 0; i < Npeaks_wix_sec3; ++i) {
	int bin = h1_e_sum_wix_sec3->FindBin(pos_wix_sec3[i]);
	v_peaks_wix_sec3.push_back(h1_e_sum_wix_sec3->GetBinLowEdge(bin));
	}
for (int i = 0; i < Npeaks_messel_sec0; ++i){
	int bin = h1_e_sum_messel_sec0->FindBin(pos_messel_sec0[i]);
	v_peaks_messel_sec0.push_back(h1_e_sum_messel_sec0->GetBinLowEdge(bin));
	}
for (int i = 0; i < Npeaks_messel_sec1; ++i){
	int bin = h1_e_sum_messel_sec1->FindBin(pos_messel_sec1[i]);
	v_peaks_messel_sec1.push_back(h1_e_sum_messel_sec1->GetBinLowEdge(bin));
	}
sort(v_peaks_wix_sec2.begin(),v_peaks_wix_sec2.end());
sort(v_peaks_wix_sec3.begin(),v_peaks_wix_sec3.end());
sort(v_peaks_messel_sec0.begin(),v_peaks_messel_sec0.end());
sort(v_peaks_messel_sec1.begin(),v_peaks_messel_sec1.end());
cout << "This is position of first peak on messel side,section0:\t" << v_peaks_messel_sec0[0] << endl;
cout << "This is position of first peak on messel side,section1:\t" << v_peaks_messel_sec1[0] << endl;
cout << "This is position of first peak on wixhausen side,section2:\t" << v_peaks_wix_sec2[0] << endl;
cout << "This is position of first peak on wixhausen side,section3:\t" << v_peaks_wix_sec3[0] << endl;

for (int i = 0; i < 5; i++){
    v_peaks_wix_sec2.erase(v_peaks_wix_sec2.end()-1);
    v_peaks_wix_sec3.erase(v_peaks_wix_sec3.end()-1);
    v_peaks_messel_sec0.erase(v_peaks_messel_sec0.end()-1);
    v_peaks_messel_sec1.erase(v_peaks_messel_sec1.end()-1);
    }

for (int i = 0; i < 5; i++){
    v_peaks_wix_sec2.erase(v_peaks_wix_sec2.end()-1);
    v_peaks_wix_sec3.erase(v_peaks_wix_sec3.end()-1);
    v_peaks_messel_sec0.erase(v_peaks_messel_sec0.end()-1);
    v_peaks_messel_sec1.erase(v_peaks_messel_sec1.end()-1);
    }


Double_t distance_peaks_wix = (v_peaks_wix_sec2[1]-v_peaks_wix_sec2[0]+v_peaks_wix_sec3[1]-v_peaks_wix_sec3[0])/2.;
vector<double> v_first_entries_wix{v_peaks_wix_sec2[0],v_peaks_wix_sec3[0]};
cout << "first entry sec2 " << v_first_entries_wix[0] << endl;
cout << "first entry sec3" << v_first_entries_wix[1] << endl;
Int_t max_first_peak_wix = max_element(v_first_entries_wix.begin(),v_first_entries_wix.end()) - v_first_entries_wix.begin();
vector<vector<double> > v_full_wix;
v_full_wix.push_back(v_peaks_wix_sec2);
v_full_wix.push_back(v_peaks_wix_sec3);
vector<double> reduce_size{-1,-1};
for (int j = 0; j < 2; j++){
	bool peak_search = true;
	int p = 0;
	while (peak_search){
		Double_t diff_peaks = abs(v_full_wix[max_first_peak_wix][0] - v_full_wix[j][p]);
		if (diff_peaks < distance_peaks_wix){
			peak_search = false;
			reduce_size[j] = p;
			}
		p+=1;
		}
	}
cout << "reduce size:\t" << reduce_size[0] << "    " << reduce_size[1] << endl;
for (int i = 0; i<reduce_size[0]; i++){
	v_peaks_wix_sec2.erase(v_peaks_wix_sec2.begin());
	}

for (int i = 0; i<reduce_size[1]; i++){
	v_peaks_wix_sec3.erase(v_peaks_wix_sec3.begin());
	}
vector<int> size_data_arr{-1,-1};
size_data_arr[0] = v_peaks_wix_sec2.size();
size_data_arr[1] = v_peaks_wix_sec3.size();
cout << "first element of sec2: \t" << v_peaks_wix_sec2[0] << endl;
cout << "first element of sec3:\t" << v_peaks_wix_sec3[0] << endl;
cout << "size sec2 " << v_peaks_wix_sec2.size() << "   and of section 3:\t" <<v_peaks_wix_sec3.size() << endl;
if (size_data_arr[0] > size_data_arr[1]){
	cout <<  "FOOO" << endl;
	for (int i = 0; i < (size_data_arr[0]-size_data_arr[1]); i++){
		v_peaks_wix_sec2.erase(v_peaks_wix_sec2.end()-1);
		}
	cout << "SECTION2" << endl;
	for (int i = 0; i < size_data_arr[1]; i++){
		cout << v_peaks_wix_sec2[i] << endl;	
		}
	cout << "SECTION3" << endl;
	for (int i = 0; i < size_data_arr[1]; i++){
		cout << v_peaks_wix_sec3[i] << endl;	
		}
	}
if (size_data_arr[0] < size_data_arr[1]){
	cout << "BAAARR" << endl;
	for (int i = 0; i < (size_data_arr[1]-size_data_arr[0]); i++){
		v_peaks_wix_sec3.erase(v_peaks_wix_sec3.end()-1);
		}
	cout << "SECTION2" << endl;
	for (int i = 0; i < size_data_arr[0]; i++){
		cout << v_peaks_wix_sec2[i] << endl;	
		}
	cout << "SECTION3" << endl;
	for (int i = 0; i < size_data_arr[0]; i++){
		cout << v_peaks_wix_sec3[i] << endl;	
		}
	}
if (size_data_arr[0] == size_data_arr[1]){
	cout << "THE LENGTH IS THE SAME; CONTINUE::::"<< endl;
	}


cout << "......section0.........." << endl;
for (int i = 0; i < v_peaks_messel_sec0.size(); i++){
cout << v_peaks_messel_sec0[i] << endl;
}
cout << "................................" << endl;
cout << "......section1.........." << endl;
for (int i = 0; i < v_peaks_messel_sec1.size(); i++){
cout << v_peaks_messel_sec1[i] << endl;
}
cout << "................................" << endl;
Double_t distance_peaks_messel = 0.6*(v_peaks_messel_sec0[1]-v_peaks_messel_sec0[0]+v_peaks_messel_sec1[1]-v_peaks_messel_sec1[0])/2.;
cout << "distance between peaks:\t" << distance_peaks_messel << endl;
vector<double> v_first_entries_messel{v_peaks_messel_sec0[0],v_peaks_messel_sec1[0]};
cout << "first entry sec2 " << v_first_entries_messel[0] << endl;
cout << "first entry sec3" << v_first_entries_messel[1] << endl;
Int_t max_first_peak_messel = max_element(v_first_entries_messel.begin(),v_first_entries_messel.end()) - v_first_entries_messel.begin();
vector<vector<double> > v_full_messel;
v_full_messel.push_back(v_peaks_messel_sec0);
v_full_messel.push_back(v_peaks_messel_sec1);
vector<double> reduce_size_messel{-1,-1};
for (int j = 0; j < 2; j++){
	bool peak_search = true;
	int p = 0;
	while (peak_search){
		Double_t diff_peaks = abs(v_full_messel[max_first_peak_messel][0] - v_full_messel[j][p]);
		if (diff_peaks < distance_peaks_messel){
			peak_search = false;
			reduce_size_messel[j] = p;
			}
		p+=1;
		}
	}
cout << "reduce size:\t" << reduce_size_messel[0] << "    " << reduce_size_messel[1] << endl;
for (int i = 0; i<reduce_size_messel[0]; i++){
	v_peaks_messel_sec0.erase(v_peaks_messel_sec0.begin());
	}

for (int i = 0; i<reduce_size_messel[1]; i++){
	v_peaks_messel_sec1.erase(v_peaks_messel_sec1.begin());
	}
vector<int> size_data_arr_messel{-1,-1};
size_data_arr_messel[0] = v_peaks_messel_sec0.size();
size_data_arr_messel[1] = v_peaks_messel_sec1.size();
cout << "first element of sec0: \t" << v_peaks_messel_sec0[0] << endl;
cout << "first element of sec1:\t" << v_peaks_messel_sec1[0] << endl;
cout << "size sec0 " << v_peaks_messel_sec0.size() << "   and of section 1:\t" <<v_peaks_messel_sec1.size() << endl;
if (size_data_arr_messel[0] > size_data_arr_messel[1]){
	cout <<  "FOOO" << endl;
	for (int i = 0; i < (size_data_arr_messel[0]-size_data_arr_messel[1]); i++){
		v_peaks_messel_sec0.erase(v_peaks_messel_sec0.end()-1);
		}
	cout << "SECTION0" << endl;
	for (int i = 0; i < size_data_arr_messel[1]; i++){
		cout << v_peaks_messel_sec0[i] << endl;	
		}
	cout << "SECTION1" << endl;
	for (int i = 0; i < size_data_arr_messel[1]; i++){
		cout << v_peaks_messel_sec1[i] << endl;	
		}
	}
if (size_data_arr_messel[0] < size_data_arr_messel[1]){
	cout << "BAAARR" << endl;
	for (int i = 0; i < (size_data_arr_messel[1]-size_data_arr_messel[0]); i++){
		v_peaks_messel_sec1.erase(v_peaks_messel_sec1.end()-1);
		}
	cout << "SECTION0" << endl;
	for (int i = 0; i < size_data_arr_messel[0]; i++){
		cout << v_peaks_messel_sec0[i] << endl;	
		}
	cout << "SECTION1" << endl;
	for (int i = 0; i < size_data_arr_messel[0]; i++){
		cout << v_peaks_messel_sec1[i] << endl;	
		}
	}
if (size_data_arr_messel[0] == size_data_arr_messel[1]){
	cout << "THE LENGTH IS THE SAME; CONTINUE::::"<< endl;
	}
//Double_t distance_peaks = (v_peaks_wix[0] - v_peaks_wix[1] + v_peaks_messel[0] + v_peaks_messel[1])/2.;
//vector<double> v_first_entries_sec{v_peaks_wix[0],v_peaks_messel[0]};
//Int_t max_first_peak_sec = max_element(v_first_entries_sec.begin(),v_first_entries_sec.end()) - v_first_entries_sec.begin();
//vector<vector<double> > v_all_sec;
//v_all_sec.push_back(v_peaks_wix);
//v_all_sec.push_back(v_peaks_messel);
//vector<double> reduce_size{-1,-1};
//for (int j = 0; j < 2; j++){
//    bool peak_search = true;
//    int p = 0;
//    while (peak_search){
//        Double_t diff_peaks = abs(v_all_sec[max_first_peak_sec][0] - v_all_sec[j][p]);
//        if (diff_peaks < distance_peaks){
//            peak_search = false;
//            reduce_size[j] = p;
//            }
//        p+=1;
//        }
//    }
//vector<int> size_data_arr{-1,-1};
////wixhausen side
//for (int i = 0; i < reduce_size[0]; i++){
//	v_peaks_wix.erase(v_peaks_wix.begin());
//	}
//
//size_data_arr[0] = v_peaks_wix.size();
//
//for (int i = 0; i < reduce_size[1]; i++){
//	v_peaks_messel.erase(v_peaks_messel.begin());
//	}
//size_data_arr[1] = v_peaks_messel.size();
//
//Int_t min_data_size = *min_element(size_data_arr.begin(),size_data_arr.end());
//if ((v_peaks_wix.size() -min_data_size) > 0){
//	for (int i = 0; i < (v_peaks_wix.size()-min_data_size); i++){
//		v_peaks_wix.erase(v_peaks_wix.end()-1);
//		}
//	}
//
//if ((v_peaks_messel.size() -min_data_size) > 0){
//	for (int i = 0; i < (v_peaks_messel.size()-min_data_size); i++){
//		v_peaks_messel.erase(v_peaks_messel.end()-1);
//		}
//	}
//for (int i = 0; i < 8; i++){
//	v_peaks_messel_sec0.erase(v_peaks_messel_sec0.begin());
//	v_peaks_messel_sec1.erase(v_peaks_messel_sec1.begin());
//	v_peaks_wix_sec2.erase(v_peaks_wix_sec2.begin());
//	v_peaks_wix_sec3.erase(v_peaks_wix_sec3.begin());
//	}
//for (int i = 0; i < 8; i++){
//	v_peaks_wix_sec2.erase(v_peaks_wix_sec2.end()-1);
//	v_peaks_wix_sec3.erase(v_peaks_wix_sec3.end()-1);
//	v_peaks_messel_sec0.erase(v_peaks_messel_sec0.end()-1);
//	v_peaks_messel_sec1.erase(v_peaks_messel_sec1.end()-1);
//	}
//min_data_size = v_peaks_wix.size();
//
//vector<int>charge_number;
//for (int i = 0; i < min_data_size; i++){
//	charge_number.push_back(i);
//	cout << "peak wixhausen:\t" << v_peaks_wix[i] << endl;
//	cout << "peak messel:\t" << v_peaks_messel[i] << endl;
//	}
//
vector<int> charge_number_sec0;
Double_t y_sec0[v_peaks_messel_sec0.size()],x_sec0[v_peaks_messel_sec0.size()];
cout << "**********Section0************" << endl;
for (int i = 0; i< v_peaks_messel_sec0.size(); i++){
	charge_number_sec0.push_back(i);	
	y_sec0[i] = i;	
	x_sec0[i] = v_peaks_messel_sec0[i];
	cout << "Peaknumber:\t" << i << "\t" << "Peakposition:\t" << v_peaks_messel_sec0[i] << endl;
	}
vector<int> charge_number_sec1;
Double_t y_sec1[v_peaks_messel_sec1.size()],x_sec1[v_peaks_messel_sec1.size()];
cout << "**********Section1************" << endl;
for (int i = 0; i< v_peaks_messel_sec1.size(); i++){
	charge_number_sec1.push_back(i);	
	y_sec1[i] = i;
	x_sec1[i] = v_peaks_messel_sec1[i];
	cout << "Peaknumber:\t" << i << "\t" << "Peakposition:\t" << v_peaks_messel_sec1[i] << endl;
	}
vector<int> charge_number_sec2;
Double_t y_sec2[v_peaks_wix_sec2.size()],x_sec2[v_peaks_wix_sec2.size()];
cout << "**********Section2************" << endl;
for (int i = 0; i< v_peaks_wix_sec2.size(); i++){
	charge_number_sec2.push_back(i);	
	y_sec2[i] = i;
	x_sec2[i] = v_peaks_wix_sec2[i];
	cout << "Peaknumber:\t" << i << "\t" << "Peakposition:\t" << v_peaks_wix_sec2[i] << endl;
	}
vector<int> charge_number_sec3;
Double_t y_sec3[v_peaks_wix_sec3.size()],x_sec3[v_peaks_wix_sec3.size()];
cout << "**********Section2************" << endl;
for (int i = 0; i< v_peaks_wix_sec3.size(); i++){
	charge_number_sec3.push_back(i);	
	y_sec3[i] = i;
	x_sec3[i] = v_peaks_wix_sec3[i];
	cout << "Peaknumber:\t" << i << "\t" << "Peakposition:\t" << v_peaks_wix_sec3[i] << endl;
	}
TCanvas *c1 = new TCanvas("c1","Canvas for Fit");
c1->Divide(2,2);
c1->cd(1);
TGraph* graph_sec0 = new TGraph(v_peaks_messel_sec0.size(),x_sec0,y_sec0);
graph_sec0->SetTitle("Alignment, Section 0");
TF1* fit_sec0 = new TF1("fit_sec0","[0]+[1]*sqrt(x)+[2]*x");
graph_sec0->Fit("fit_sec0");
graph_sec0->Draw("AP*");
par_file << fit_sec0->GetParameter(0) << "," << fit_sec0->GetParameter(1) << "," << fit_sec0->GetParameter(2) << endl;
c1->cd(2);
TGraph* graph_sec1 = new TGraph(v_peaks_messel_sec1.size(),x_sec1,y_sec1);
graph_sec1->SetTitle("Alignment, Section 1");
TF1* fit_sec1 = new TF1("fit_sec1","[0]+[1]*sqrt(x)+[2]*x");
graph_sec1->Fit("fit_sec1");
graph_sec1->Draw("AP*");
par_file << fit_sec1->GetParameter(0) << "," << fit_sec1->GetParameter(1) << "," << fit_sec1->GetParameter(2) << endl;
c1->cd(3);
TGraph* graph_sec2 = new TGraph(v_peaks_wix_sec2.size(),x_sec2,y_sec2);
graph_sec2->SetTitle("Alignment, Section 2");
TF1* fit_sec2 = new TF1("fit_sec2","[0]+[1]*sqrt(x)+[2]*x");
graph_sec2->Fit("fit_sec2");
graph_sec2->Draw("AP*");
par_file << fit_sec2->GetParameter(0) << "," << fit_sec2->GetParameter(1) << "," << fit_sec2->GetParameter(2) << endl;
c1->cd(4);
TGraph* graph_sec3 = new TGraph(v_peaks_wix_sec3.size(),x_sec3,y_sec3);
graph_sec3->SetTitle("Alignment, Section 3");
TF1* fit_sec3 = new TF1("fit_sec3","[0]+[1]*sqrt(x)+[2]*x");
graph_sec3->Fit("fit_sec3");
graph_sec3->Draw("AP*");
par_file << fit_sec3->GetParameter(0) << "," << fit_sec3->GetParameter(1) << "," << fit_sec3->GetParameter(2) << endl;



//Double_t x_wix[100], y[100], x_messel[100];
//for (int j = 0; j < min_data_size; j++){
//	x_wix[j] = v_peaks_wix[j];
//	x_messel[j] = v_peaks_messel[j];
//	y[j] = charge_number[j];
//	cout << "Charge number:\t" << y[j] << endl;
//	cout << "Messel points:\t" << x_messel[j] << endl;
//	}	
//TGraph* gr_wix = new TGraph(min_data_size,x_wix,y);
//gr_wix->SetTitle("Wixhausen Alignment");
//TF1 *fit_wix = new TF1("fit_wix","[0]+[1]*sqrt(x)+[2]*x");
//gr_wix->Fit("fit_wix");
//gr_wix->Draw("ap*");
//par_file << fit_wix->GetParameter(0) << "," << fit_wix->GetParameter(1) << "," << fit_wix->GetParameter(2) << endl;
//c1->cd(2);
//TGraph* gr_messel = new TGraph(min_data_size,x_messel,y);
//gr_messel->SetTitle("Messel Alignment");
//TF1 *fit_messel = new TF1("fit_messel","[0]+[1]*sqrt(x)+[2]*x");
//gr_messel->Fit("fit_messel");
//gr_messel->Draw("ap*");
//par_file << fit_messel->GetParameter(0) << "," << fit_messel->GetParameter(1) << "," << fit_messel->GetParameter(2) << endl;
//cout << "first point of messel in x:\t" << x_messel[0] << endl;
//cout << "first y point messel:\t" << y[0] << endl;
//
//




//Double_t x[100], y[100];
//for (int j = 0; j < min_data_size; j++){
//    x[j] = v_peaks_wix[j];
//    y[j] = v_peaks_wix[j];
//	cout << "peaks on wixhausen side:\t" << v_peaks_wix[j] << endl;
//}
//TGraph* gr_wix = new TGraph(min_data_size,x,y);
//TF1 *fit_wix = new TF1("fit_wix","[0]+x*[1]+pow(x,2)*[2]+pow(x,3)*[3]+pow(x,4)*[4]+pow(x,5)*[5]+pow(x,6)*[6]");
//gr_wix->Fit("fit_wix");
//gr_wix->Draw("AC*");
////Double_t slope_wix =  1./(fit_wix->GetParameter(1));
////Double_t offset_wix = -(fit_wix->GetParameter(0))/(fit_wix->GetParameter(1));
////par_file << slope_wix << "," << offset_wix << endl;
//par_file << fit_wix->GetParameter(0) << "," << fit_wix->GetParameter(1) << "," << fit_wix->GetParameter(2) << "," << fit_wix->GetParameter(3) \
// << "," << fit_wix->GetParameter(4) << "," << fit_wix->GetParameter(5) << "," << fit_wix->GetParameter(6) << endl;
//
//c1->cd(2);
//for (int j = 0; j < min_data_size; j++){
//	x[j] = v_peaks_wix[j];
//	y[j] = v_peaks_messel[j];
//	}
//TGraph* gr_messel = new TGraph(min_data_size,x,y);
//TF1 *fit_messel = new TF1("fit_messel","[0]+x*[1]+pow(x,2)*[2]+pow(x,3)*[3]+pow(x,4)*[4]+pow(x,5)*[5]+pow(x,6)*[6]");
//gr_messel->Fit("fit_messel");
//gr_messel->Draw("AC*");
////Double_t slope_messel =  1./(fit_messel->GetParameter(1));
////Double_t offset_messel = -(fit_messel->GetParameter(0))/(fit_messel->GetParameter(1));
////par_file << slope_messel << "," << offset_messel << endl;
//par_file << fit_messel->GetParameter(0) << "," << fit_messel->GetParameter(1) << "," << fit_messel->GetParameter(2) << "," << fit_messel->GetParameter(3) \
// << "," << fit_messel->GetParameter(4) << "," << fit_messel->GetParameter(5) << "," << fit_messel->GetParameter(6) << endl;
par_file.close();
c1->SaveAs("energy_charge_twim_translation.png");
}
