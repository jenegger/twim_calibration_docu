#include <iostream>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include <math.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TChain.h"
#include <string>
using namespace std;
//input in mm and output in cm....not so good...
vector<double> path_tof_calc (double psi_in,double xMW1,double xMW3,double t_tof,double t_start){  //xMW1 I  will extrapolate from the calibrated MUSIC
//insert here all detector positions TODO

	//important constants,variables...
	const Double_t light_c = 299792458.;
	const Double_t PI = 3.14159265358979323846;
	//TODO............ all those shifts have to be checked; at the begin may be set to 0....
	const Double_t target_shift = 0.;
	const Double_t xMW1_shift = 0.;
	const Double_t xMW2_shift = 0.;
	//const Double_t xMW0_shift = 0.; for now not interesting...
	//const Double_t yMW0_shift = 0.; for now not interesting...
	const Double_t xMW3_shift = 0.;
	const Double_t yMW1_shift = 0.;
	const Double_t yMW2_shift = 0.;
	const Double_t yMW3_shift = 0.;
	//end of TODO..............
	const Double_t alpha_G = (14./180.)*PI;
	const Double_t zGm = 2577.;  //value in the par file is 2155, I took the one from S444
	const Double_t zMW0 = -2620;
	const Double_t middle_zM3 = 5789.7;
	const Double_t middle_xM3 = -1638.1;
	const Double_t middle_zToFW = 6509.7;
	const Double_t middle_xToFW = -1908.2;
	const Double_t cutoff_ToFW = middle_xToFW - middle_zToFW*tan(2*PI/5.);
	const Double_t zT = -737.5;
	const Double_t zM1 = 300.;
	const Double_t zM2 = 935;
	const Double_t psi_out_0 = PI/10.;
	const Double_t start_to_target = 1130.25;
	const Double_t current = 1444;
	const Double_t beta_given = 0.804555; //value for 637AMeV
	const Double_t gamma_given = 1.68385; //value for 637AMeV
	const Double_t time_start_target = (1130.25/(beta_given*light_c))*pow(10,6);
	const Double_t zGm_cutoff = -tan(PI/2.-alpha_G)*(zGm);
	//declare more constants
	Double_t zC;
	Double_t xC;
	Double_t slope;
	Double_t offset_slope;
	Double_t psi_out_rec[50];
	Double_t z_pos_shift[50];
	Double_t x_pos_shift[50];
	Double_t z_D;
	Double_t x_D;
	//original value of Leff!!----
	const Double_t Leff = 2067;
	Double_t l_diff[50];
	const Double_t zGLAD1 = (Leff/2.)/(cos(alpha_G));
	const Double_t bGLAD_cutoff = -tan(PI/2.-alpha_G)*(zGm-zGLAD1);
	const Double_t D_cutoff = -tan(PI/2.-alpha_G)*(zGm+zGLAD1);
	//-------------------------------

	Double_t T_to_M2 = (zM2-zT)/(cos(psi_in));
	Double_t b1_cutoff = (xMW1+xMW1_shift)-tan(psi_in)*zM1;
	Double_t z_labMW3 = middle_zM3+cos(2*PI/5.)*(xMW3-xMW3_shift);
	Double_t x_labMW3 = middle_xM3+sin(2*PI/5.)*(xMW3-xMW3_shift);
	Double_t zB = (bGLAD_cutoff-b1_cutoff)/(tan(psi_in)-tan(PI/2.-alpha_G));
	Double_t xB = tan(psi_in)*zB+b1_cutoff;
	
        for (Int_t rec = 0; rec < 50; rec++){
            zC = (zGm_cutoff-b1_cutoff)/(tan(psi_in)-tan(PI/2.-alpha_G)) + 2*rec;
            xC = tan(psi_in)*zC+b1_cutoff;
            slope = (x_labMW3-xC)/(z_labMW3-zC);
            offset_slope = x_labMW3-slope*z_labMW3;
            psi_out_rec[rec] = -atan(slope);
            z_D = (D_cutoff-offset_slope)/(slope-tan(PI/2.-alpha_G));
            x_D = tan(PI/2.-alpha_G)*z_D +D_cutoff;
            l_diff[rec] = (sqrt((xC-xB)*(xC-xB)+(zC-zB)*(zC-zB))-sqrt((x_D-xC)*(x_D-xC)+(z_D-zC)*(z_D-zC)));
            z_pos_shift[rec] = zC;
            x_pos_shift[rec] = xC;
        }
        Int_t n = 50;
        TGraph *gr1 = new TGraph(n,psi_out_rec,l_diff);
        gr1->Fit("pol1","Q");
        TF1 *f3 =gr1->GetFunction("pol1");
        Double_t f3_slope = f3->GetParameter(1);
        Double_t f3_offset = f3->GetParameter(0);
        Double_t psi_out = -f3_offset/f3_slope;
        TGraph *gr2 = new TGraph(n,z_pos_shift,l_diff);
        gr2->Fit("pol1", "Q");
        TF1 *f4 = gr2->GetFunction("pol1");
        Double_t f4_slope = f4->GetParameter(1);
        Double_t f4_offset = f4->GetParameter(0);
        Double_t new_z_C = -f4_offset/f4_slope;

        TGraph *gr3 = new TGraph(n,x_pos_shift,l_diff);
        gr3->Fit("pol1","Q");
        TF1 *f5 = gr3->GetFunction("pol1");
        Double_t f5_slope = f5->GetParameter(1);
        Double_t f5_offset = f5->GetParameter(0);
        Double_t new_x_C = -f5_offset/f5_slope;
        Double_t z_diff = new_z_C-z_pos_shift[0];
        Double_t offset_new = new_x_C+tan(psi_out)*new_z_C;
        z_D = (D_cutoff-offset_new)/(-tan(psi_out)-tan(PI/2.-alpha_G));
        x_D = -tan(psi_out)*z_D + offset_new;
        Double_t zToFW = (cutoff_ToFW-offset_new)/(-tan(psi_out)-tan(2*PI/5.));
        Double_t xToFW = -tan(psi_out)*zToFW+offset_new;
        Double_t path_D_to_TOFW = (zToFW-z_D)/(cos(psi_out));
        Double_t BD = sqrt((x_D-xB)*(x_D-xB)+(z_D-zB)*(z_D-zB));

        Double_t m = (cos(psi_out)-cos(psi_in))/(sin(psi_out)+sin(psi_in));
        Double_t rho =  (Leff/(2*sin((psi_in+psi_out)/2)))*sqrt(pow(((m+tan(alpha_G))/(1-m*tan(alpha_G))),2)+1);
        Double_t w = 2*abs(asin(BD/(2*rho)));

        delete gr1;
        delete gr2;
        delete gr3;
	Double_t path_from_target = (abs((zB-zT)/(cos(psi_in)))+rho*w+abs((zToFW-z_D)/cos(psi_out)))*0.1;  //in cm!
	//cout << "path from target\t" << path_from_target << endl;
	//cout << "first part:\t" << abs((zB-zT)/(cos(psi_in))) << endl;
	//cout << "second part:\t" << rho*w+abs((zToFW-z_D)/cos(psi_out)) << endl;
	//cout << "psi in:\t" << psi_in << endl;
	Double_t time_of_flight = t_tof - t_start - time_start_target;
	vector<double> v_return = {path_from_target,time_of_flight,rho};
	return v_return;
}
