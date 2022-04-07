#ifndef ALL_PARAMS_H
#define ALL_PARAMS_H
#include <iostream>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include <math.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TClonesArray.h"
#include "TChain.h"
#include <sstream>
#include <string>

char hist_name[500];

const Double_t light_c = 299792458.;
const Double_t PI = 3.14159265358979323846;
const Double_t target_shift = 0.;
//these were the old shifts
//const Double_t xMW1_shift = 0.344833+0.7626;
//const Double_t xMW2_shift = 0.451628+0.9026;
const Double_t xMW1_shift = 0.658351;  //offsets from fit on run 36
const Double_t xMW2_shift = 0.784281;  //offsets from fit on run 36
//--------------------------------------------------------------
const Double_t xMW0_shift = 5.72037; //offsets from fit on run 36
//old shift
//const Double_t yMW0_shift = 1.17030e+01;
const Double_t yMW0_shift = 11.6595; //offsets from fit on run 36
//old shift
//const Double_t xMW3_shift = 79.84;
const Double_t xMW3_shift = 79.8035; //offsets from fit on run 36
//old shifts
//const Double_t yMW1_shift = 1.37840e+01;
//const Double_t yMW2_shift = 3.22703e+01;
//const Double_t yMW3_shift = -6.70074e+00;
const Double_t yMW1_shift = 13.7538; //offsets from fit on run 36
const Double_t yMW2_shift = 31.59; //offsets from fit on run 36
const Double_t yMW3_shift = -11.5563; //offsets from fit on run 36
//original value of alpha_g is 14 degrees
const Double_t alpha_G = (14./180.)*PI;
//this is original value of zGm-----
const Double_t zGm = 2577.;
//----------------------------------
const Double_t zMW0 = -2520.;
//probably START and MW0 Detectors are switched, have to check...
//const Double_t zMW0 = -1867.75;
//new, I try to use now also the shift, to get more out of it....
//this is old version of MW0 shift... to be checked which is the right one ..:TODO
//const Double_t xMW0_shift = 5.8959;
const Double_t middle_zM3 = 5937;
const Double_t middle_xM3 = -(middle_zM3-zGm)*tan(PI/10.);
const Double_t middle_zToFW = (6100+760-2777)*cos(PI/10.)+2777;
const Double_t middle_xToFW = -(6100+760-2777)*sin(PI/10.);
const Double_t cutoff_ToFW = middle_xToFW - middle_zToFW*tan(2*PI/5.);
const Double_t zT = -684.5;
const Double_t zM1 = 279. - 81;  //to check again....
const Double_t zM2 = 854.;
const Double_t psi_out_0 = PI/10.;
const Double_t start_to_target = 1183.25;
//const Double_t time_start_target = (1183.25/(0.714549*light_c))*pow(10,6);
const Double_t time_start_target = (1162./(0.714549*light_c))*pow(10,6);
const Double_t current = 1444;
//const Double_t gamma_given = 1.42942;//with ingoing 12C 400 MeV
//const Double_t beta_given = 0.714549;//with ingoing 12C 400 MeV
const Double_t beta_given = 0.709327; //with ingoing 12C 390 MeV
const Double_t gamma_given = 1.41868; //corresponding gamma for 390MeV/u
//const Double_t current = 1498;
const Double_t zGm_cutoff = -tan(PI/2.-alpha_G)*(zGm);
Int_t entries_mw0 = 0.;
Int_t entries_mw1 = 0.;
Int_t entries_mw2 = 0.;
Int_t entries_mw3 = 0.;
Int_t entries_start = 0.;
Int_t entries_tofw = 0.;
Int_t entries_califa = 0.;
Double_t tot_length = 0.;
Double_t charge_val = 0.;
Double_t psi_out;
//Declare more functions here...
Double_t zC;
Double_t xC;
Double_t z_labMW3;
Double_t x_labMW3;
Double_t slope;
Double_t offset_slope;
Double_t psi_out_rec[50];
Double_t z_pos_shift[50];
Double_t x_pos_shift[50];
Double_t zB;
Double_t xB;
Double_t z_D;
Double_t x_D;
//original value of Leff!!----
const Double_t Leff = 2067;
//----------------------------
Double_t l_diff[50];
const Double_t zGLAD1 = (Leff/2.)/(cos(alpha_G));
const Double_t bGLAD_cutoff = -tan(PI/2.-alpha_G)*(zGm-zGLAD1);
const Double_t D_cutoff = -tan(PI/2.-alpha_G)*(zGm+zGLAD1);

//these are the parameters needed for the time of flight calculation
Long64_t nHits =0;
Long64_t nHitsStart = 0;
UShort_t iCh;//Channel number of the TOFD per Detector
UShort_t iDet;//detector number of the TOFD
#define NbDets  28 //detector number of the TOFD
#define  NbChs  2 //number of preams per detector plastic
Double_t iRawTimeNs[NbDets * 2];
Double_t iRawTimeStartNs[2];
UShort_t mult[NbDets * NbChs];
double tof_offs[27] = {};

double offsets[29] = {0,-2.31775e+00,-8.63789e-01,-1.28369e+00,2.02607e+00,6.30893e-01,5.29162e-01,-8.16594e-01,1.69442e-01,-3.52510e+00,-2.16979e+00,1.66105e+00,1.26320e+00,-2.26356e+00,-4.42463e+00,-6.51834e+00,-1.42568e-01,2.20837e+00,6.15364e+00,4.49112e-01,-4.82568e-01,-1.27471e+00,-3.10422e+00,-1.06742e+00,-2.47400e+00,-5.09690e-01,2.38830e-01,-7.83351e-01,0};
void all_params(string runnr){
if (runnr.compare("0082_0001") == 0 || runnr.compare("0083_0001") == 0 || runnr.compare("0078_0001")){
	double tof_offs_1[27] = {
	114.916,
	113.592,
	116.172,
	114.396,
	115.279,
	115.910,
	116.751,
	116.002,
	114.242,
	115.822,
	112.410,
	115.801,
	113.355,
	114.925,
	113.514,
	112.817,
	114.214,
	111.899,
	113.050,
	112.578,
	116.283,
	116.685,
	115.884,
	115.922,
	114.642,
	114.207,
	114.212
	};
	for(int n = 0; n < 27; n++){
		tof_offs[n] = tof_offs_1[n];
		}
	}
else {
	double tof_offs_2[27] = {
	114.4,
	113.588,
	116.132,
	114.358,
	115.261,
	115.886,
	116.734,
	115.982,
	114.222,
	115.801,
	112.402,
	115.786,
	113.347,
	114.917,
	118.516,
	117.823,
	119.217,
	116.913,
	118.063,
	117.591,
	121.295,
	121.690,
	120.896,
	120.931,
	119.649,
	119.225,
	119.274
	};
	for(int n = 0; n < 27; n++){
    tof_offs[n] = tof_offs_2[n];
    }

	}
}


#endif
