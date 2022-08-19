// include headers
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
bool sortcol( const vector<double>& v1,const vector<double>& v2 ) {
                                    return v1[1] < v2[1];
                                    }
static const vector<double>  tof_offset = {
-87.9512, -86.1871, -87.7881, -86.1239, -87.2912, -87.6336, -88.4088, -88.2992, -85.683, -87.1578,-84.9038, -87.8791, -85.3067, -86.7159, -105.254, -104.595, -105.741, -103.607, -104.681, -104.396,
  -106.93, -53.5665, -52.4257, -53.0346, -51.5958, -52.4908, -51.9334, -53.8351
};
static const double tof_lise = 37.0;
////parameters_for TOFW calibration y position
//fstream fin4;
//fin4.open("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/macros/tof_params.csv", ios::in); 
//string line4, word4;
//vector<vector<double> > v_para_tof_y;
//getline(fin4, line4);
//while(fin4.peek()!=EOF) {
//	getline(fin4, line4);
//	stringstream s(line4);
//	vector<double> temp_vec;
//	while (getline(s, word4, ',')) {     temp_vec.push_back(stod(word4));
//	}
//	v_para_tof_y.push_back(temp_vec);
//	
//	temp_vec.clear();
//}
vector<vector<double> >after_glad_pos(TClonesArray* mw3_clone,TClonesArray* tofsingletcal_clone){
R3BMwpcCalData** mwpc3_caldata;
R3BSofTofWSingleTcalData** tofsingletcal_data;
Int_t entries_tof_singletcal = tofsingletcal_clone->GetEntriesFast();
Int_t entries_mw3_cal =  mw3_clone->GetEntriesFast();

mwpc3_caldata = new R3BMwpcCalData*[entries_mw3_cal];
tofsingletcal_data = new R3BSofTofWSingleTcalData*[entries_tof_singletcal];
vector<double> vec_X_mw3(288,0);
vector<double> vec_Y_mw3(120,0);
for (Int_t j = 0; j < entries_mw3_cal;j++){
	mwpc3_caldata[j] = (R3BMwpcCalData*)mw3_clone->At(j);
	Int_t planeId = mwpc3_caldata[j]->GetPlane();
	Int_t padId = mwpc3_caldata[j]->GetPad()-1;
	Double_t charge = mwpc3_caldata[j]->GetQ();
	if (planeId == 1){ //X-PLANE IN MW3
		vec_X_mw3[padId] = charge;
		}
	if (planeId == 3){ //y-plane in MW3
		vec_Y_mw3[padId] = charge;
		}
	}
//get peaks from mw3
vector<double> peaks_X_mw3;
vector<double> peaks_Y_mw3;
for (Int_t i = 0; i < 2;i++){
	Int_t max_element_index_mw3_X = max_element(vec_X_mw3.begin(),vec_X_mw3.end())-vec_X_mw3.begin();
	Double_t max_element_mw3_X = *max_element(vec_X_mw3.begin(),vec_X_mw3.end());
	Int_t max_element_index_mw3_Y = max_element(vec_Y_mw3.begin(),vec_Y_mw3.end())-vec_Y_mw3.begin();
	Double_t max_element_mw3_Y = *max_element(vec_Y_mw3.begin(),vec_Y_mw3.end());
	//MW3 X
	if (max_element_mw3_X > 200){
		Double_t q_left;
        Double_t q_right;
		Double_t q_max = vec_X_mw3[max_element_index_mw3_X];
		if (max_element_index_mw3_X == 0){
		 	q_left=1;
			q_right = vec_X_mw3[max_element_index_mw3_X+1];
			}
		else if (max_element_index_mw3_X == 287){
		 	q_right=1;
			q_left = vec_X_mw3[max_element_index_mw3_X-1];
			}
		else {
			q_left = vec_X_mw3[max_element_index_mw3_X-1];
			q_right = vec_X_mw3[max_element_index_mw3_X+1];
			}
		if (q_left == 0){
			q_left=1;
			}
		if(q_right==0){
			q_right=1;
			}
		peaks_X_mw3.push_back(get_mw3_pos_X(q_max,max_element_index_mw3_X,q_left,q_right));
		if (max_element_index_mw3_X == 0){
			vec_X_mw3[max_element_index_mw3_X] = 0;
			vec_X_mw3[max_element_index_mw3_X+1] = 0;
			}	
		else if (max_element_index_mw3_X == 287){
			vec_X_mw3[max_element_index_mw3_X] = 0;
			vec_X_mw3[max_element_index_mw3_X-1] = 0;
			}	
		else {
			vec_X_mw3[max_element_index_mw3_X] = 0;
			vec_X_mw3[max_element_index_mw3_X+1] = 0;
			vec_X_mw3[max_element_index_mw3_X-1] = 0;
			}
		}
//	//MW3 Y
	if (max_element_mw3_Y > 200){
		Double_t q_left;
        Double_t q_right;
		Double_t q_max = vec_Y_mw3[max_element_index_mw3_Y];
		if (max_element_index_mw3_Y == 0) q_left=1;
		else if (max_element_index_mw3_Y == 287) q_right=1;
		else {
			q_left = vec_Y_mw3[max_element_index_mw3_Y-1];
			q_right = vec_Y_mw3[max_element_index_mw3_Y+1];
			}
		if (q_left == 0) q_left=1;
		if(q_right==0) q_right=1;
		peaks_Y_mw3.push_back(get_mw3_pos_Y(q_max,max_element_index_mw3_Y,q_left,q_right));
		if (max_element_index_mw3_Y == 0){
			vec_Y_mw3[max_element_index_mw3_Y] = 0;
			vec_Y_mw3[max_element_index_mw3_Y+1] = 0;
			}	
		else if (max_element_index_mw3_Y == 287){
			vec_Y_mw3[max_element_index_mw3_Y] = 0;
			vec_Y_mw3[max_element_index_mw3_Y-1] = 0;
			}	
		else {
			vec_Y_mw3[max_element_index_mw3_Y] = 0;
			vec_Y_mw3[max_element_index_mw3_Y+1] = 0;
			vec_Y_mw3[max_element_index_mw3_Y-1] = 0;
			}
		}
	}

//get info about the tof
vector<vector<double> >tof_hits;
//tofsingletcal_data = new R3BSofTofWSingleTcalData*[entries_tof_singletcal];
for (Int_t j = 0; j < entries_tof_singletcal; j++){
	tofsingletcal_data[j] = (R3BSofTofWSingleTcalData*)tofsingletcal_clone->At(j);
    Int_t paddle_ID = tofsingletcal_data[j]->GetDetector();
	Double_t posy_ns = tofsingletcal_data[j]->GetRawPosNs();
	//cout << "posy:\t" << posy_ns << endl;
	Double_t tof_time = tofsingletcal_data[j]->GetRawTimeNs();
	Double_t posx = 840./2. -15. -(Double_t)(paddle_ID-1)*30; //hard coded dimensions (in mm)
	vector<double> temp_vec{posx,posy_ns,tof_time,(Double_t) paddle_ID};
	tof_hits.push_back(temp_vec);
	temp_vec.clear();
	}
//sort tof hits according to the y position
sort(tof_hits.begin(),tof_hits.end(),sortcol);
vector<double> tof_down = tof_hits[0];
vector<double> tof_up = tof_hits[1];
//sort mw3_x
//cout << "peaks_X_mw3.size:\t" << peaks_X_mw3.size() << endl;
vector<vector<double> > positioning_after_glad;
if (peaks_X_mw3.size() == 2){
	sort(peaks_X_mw3.begin(),peaks_X_mw3.end());
if (tof_up[0] > tof_down[0]){
	//cout << "tof_up[0] > tof_down[0]" << endl;
	vector<double> v_temp1{peaks_X_mw3[0],tof_down[0],tof_down[1],tof_down[2],(Double_t) tof_down[3]};
	positioning_after_glad.push_back(v_temp1);
	vector<double> v_temp2{peaks_X_mw3[1],tof_up[0],tof_up[1],tof_up[2],(Double_t) tof_up[3]};
	positioning_after_glad.push_back(v_temp2);
	}

else {
////TJ change again TODO
////	cout << "else" << endl;
	vector<double> v_temp2{peaks_X_mw3[1],tof_down[0],tof_down[1],tof_down[2],(Double_t)tof_down[3]};
	positioning_after_glad.push_back(v_temp2);
	vector<double> v_temp1{peaks_X_mw3[0],tof_up[0],tof_up[1],tof_up[2],(Double_t) tof_up[3]};
	positioning_after_glad.push_back(v_temp1);
////	cout << "end of else" << endl;
//	vector<vector<double> > dummy_vec {
//                {-1000,-1000,-1000,-1000,-1000},
//                {-1000,-1000,-1000,-1000,-1000}
//                };
//    positioning_after_glad = dummy_vec;
	}
//	cout << "positioning vector:\t" << positioning_after_glad[0][0] << positioning_after_glad[0][1] << positioning_after_glad[0][2] << positioning_after_glad[0][3] << endl;
//	cout << "second positioning vector:\t" << positioning_after_glad[1][0] << positioning_after_glad[1][1] << positioning_after_glad[1][2] << positioning_after_glad[1][3] << endl;
	}
else{
	vector<vector<double> > dummy_vec {
				{-1000,-1000,-1000,-1000,-1000},
				{-1000,-1000,-1000,-1000,-1000}
				};
	positioning_after_glad = dummy_vec;
	}
	
	//cout << " I return this" << endl;
	//cout << positioning_after_glad[0][0] << "    " << positioning_after_glad[1][0] << endl;
	delete [] mwpc3_caldata;
	delete [] tofsingletcal_data;
	//cout <<"POS1:\t"<< positioning_after_glad[0][0] << "  " << positioning_after_glad[0][1] << "  " << positioning_after_glad[0][2] << "   " << positioning_after_glad[0][3] << endl;
	//cout <<"POS2:\t"<< positioning_after_glad[1][0] << "  " << positioning_after_glad[1][1] << "  " << positioning_after_glad[1][2] << "   " << positioning_after_glad[1][3] << endl;

	return positioning_after_glad;
}

vector<vector<double> > tof_and_y_calibrated(TClonesArray* mw3_clone,TClonesArray* tofsingletcal_clone){
	fstream fin4;
	fin4.open("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/macros/tof_params.csv", ios::in); 
	string line4, word4;
	vector<vector<double> > v_para_tof_y;
	getline(fin4, line4);
	while(fin4.peek()!=EOF) {
		getline(fin4, line4);
		stringstream s(line4);
		vector<double> temp_vec;
		while (getline(s, word4, ',')) {     temp_vec.push_back(stod(word4));
		}
		v_para_tof_y.push_back(temp_vec);
		
		temp_vec.clear();
	}
	R3BMwpcCalData** mwpc3_caldata;
	Int_t entries_mw3_cal =  mw3_clone->GetEntriesFast();
	mwpc3_caldata = new R3BMwpcCalData*[entries_mw3_cal];	
	vector<double> vec_X_mw3(288,0);
	for (Int_t j = 0; j < entries_mw3_cal;j++){
		mwpc3_caldata[j] = (R3BMwpcCalData*)mw3_clone->At(j);
		Int_t planeId = mwpc3_caldata[j]->GetPlane();
		Int_t padId = mwpc3_caldata[j]->GetPad()-1;
		Double_t charge = mwpc3_caldata[j]->GetQ();
		if (planeId == 1){ //X-PLANE IN MW3
			vec_X_mw3[padId] = charge;
			}
		}
	//get peaks from mw3
	vector<double> peaks_X_mw3;
	
	for (Int_t i = 0; i < 2;i++){
		Int_t max_element_index_mw3_X = max_element(vec_X_mw3.begin(),vec_X_mw3.end())-vec_X_mw3.begin();
		Double_t max_element_mw3_X = *max_element(vec_X_mw3.begin(),vec_X_mw3.end());
		//MW3 X
		if (max_element_mw3_X > 200){
			Double_t q_left;
	        Double_t q_right;
			Double_t q_max = vec_X_mw3[max_element_index_mw3_X];
			if (max_element_index_mw3_X == 0){
			 	q_left=1;
				q_right = vec_X_mw3[max_element_index_mw3_X+1];
				}
			else if (max_element_index_mw3_X == 287){
			 	q_right=1;
				q_left = vec_X_mw3[max_element_index_mw3_X-1];
				}
			else {
				q_left = vec_X_mw3[max_element_index_mw3_X-1];
				q_right = vec_X_mw3[max_element_index_mw3_X+1];
				}
			if (q_left == 0){
				q_left=1;
				}
			if(q_right==0){
				q_right=1;
				}
			peaks_X_mw3.push_back(get_mw3_pos_X(q_max,max_element_index_mw3_X,q_left,q_right));
			if (max_element_index_mw3_X == 0){
				vec_X_mw3[max_element_index_mw3_X] = 0;
				vec_X_mw3[max_element_index_mw3_X+1] = 0;
				}	
			else if (max_element_index_mw3_X == 287){
				vec_X_mw3[max_element_index_mw3_X] = 0;
				vec_X_mw3[max_element_index_mw3_X-1] = 0;
				}	
			else {
				vec_X_mw3[max_element_index_mw3_X] = 0;
				vec_X_mw3[max_element_index_mw3_X+1] = 0;
				vec_X_mw3[max_element_index_mw3_X-1] = 0;
				}
			}
		}

	if (peaks_X_mw3.size() == 2){
	
    sort(peaks_X_mw3.begin(),peaks_X_mw3.end());

	R3BSofTofWSingleTcalData** tofsingletcal_data;
	Int_t entries_tof_singletcal = tofsingletcal_clone->GetEntriesFast();	
	tofsingletcal_data = new R3BSofTofWSingleTcalData*[entries_tof_singletcal];
	vector<vector<double> >tof_hits;
	for (Int_t j = 0; j < entries_tof_singletcal; j++){
		tofsingletcal_data[j] = (R3BSofTofWSingleTcalData*)tofsingletcal_clone->At(j);
		Int_t paddle_ID = tofsingletcal_data[j]->GetDetector();
		//cout << "does this really work:\t" << fTofWHitPar->GetPosOffsetPar(paddle_ID) << endl;

		Double_t posy_ns = tofsingletcal_data[j]->GetRawPosNs();
		Double_t tof_time = tofsingletcal_data[j]->GetRawTofNs();
		tof_time = tof_time - tof_offset[paddle_ID-1] + tof_lise;
		Double_t calibrated_y_pos = posy_ns*v_para_tof_y[paddle_ID][1]+v_para_tof_y[paddle_ID][2];
		//cout << "uncalibrated time:\t" << tofsingletcal_data[j]->GetRawTofNs() << endl;
		//cout << "time:\t" << tof_time << endl;
		vector<double> temp_vec{tof_time,calibrated_y_pos,(Double_t) paddle_ID,(Double_t) paddle_ID};
		tof_hits.push_back(temp_vec);
		temp_vec.clear();
		}
	sort(tof_hits.begin(),tof_hits.end(),sortcol);
	if (tof_hits[0][2] < tof_hits[1][2]){
		tof_hits[0][2] = peaks_X_mw3[1];
		tof_hits[1][2] = peaks_X_mw3[0];

		}
	if (tof_hits[0][2] > tof_hits[1][2]){
		tof_hits[0][2] = peaks_X_mw3[0];
		tof_hits[1][2] = peaks_X_mw3[1];	
		}
	if (tof_hits[0][2] == tof_hits[1][2]){
		tof_hits = {
                {-1000,-1000,-1000,-1000},
                {-1000,-1000,-1000,-1000}
                };
		}
	delete [] tofsingletcal_data;
	return tof_hits; //retrurn two vectors with tof and y, fist entry lower y, second entry higher y
	}
	else {
	vector<vector<double> > dummy_vec {
                {-1000,-1000,-1000,-1000},
                {-1000,-1000,-1000,-1000}
                };
	return dummy_vec;
	}

	}


//single hit 238U reconstruction

vector<double> single_tof_and_y_calibrated(TClonesArray* mw3_clone,TClonesArray* tofsingletcal_clone){
	fstream fin4;
	fin4.open("/scratch8/ge37liw/workingspace/exp_s455/my_macros/twim_calibration_docu/macros/tof_params.csv", ios::in); 
	string line4, word4;
	vector<vector<double> > v_para_tof_y;
	getline(fin4, line4);
	while(fin4.peek()!=EOF) {
		getline(fin4, line4);
		stringstream s(line4);
		vector<double> temp_vec;
		while (getline(s, word4, ',')) {     temp_vec.push_back(stod(word4));
		}
		v_para_tof_y.push_back(temp_vec);
		
		temp_vec.clear();
	}
	R3BMwpcCalData** mwpc3_caldata;
	Int_t entries_mw3_cal =  mw3_clone->GetEntriesFast();
	mwpc3_caldata = new R3BMwpcCalData*[entries_mw3_cal];	
	vector<double> vec_X_mw3(288,0);
	for (Int_t j = 0; j < entries_mw3_cal;j++){
		mwpc3_caldata[j] = (R3BMwpcCalData*)mw3_clone->At(j);
		Int_t planeId = mwpc3_caldata[j]->GetPlane();
		Int_t padId = mwpc3_caldata[j]->GetPad()-1;
		Double_t charge = mwpc3_caldata[j]->GetQ();
		if (planeId == 1){ //X-PLANE IN MW3
			vec_X_mw3[padId] = charge;
			}
		}
	//get peaks from mw3
	vector<double> peaks_X_mw3;
	
	for (Int_t i = 0; i < 2;i++){
		Int_t max_element_index_mw3_X = max_element(vec_X_mw3.begin(),vec_X_mw3.end())-vec_X_mw3.begin();
		Double_t max_element_mw3_X = *max_element(vec_X_mw3.begin(),vec_X_mw3.end());
		//MW3 X
		if (max_element_mw3_X > 200){
			Double_t q_left;
	        Double_t q_right;
			Double_t q_max = vec_X_mw3[max_element_index_mw3_X];
			if (max_element_index_mw3_X == 0){
			 	q_left=1;
				q_right = vec_X_mw3[max_element_index_mw3_X+1];
				}
			else if (max_element_index_mw3_X == 287){
			 	q_right=1;
				q_left = vec_X_mw3[max_element_index_mw3_X-1];
				}
			else {
				q_left = vec_X_mw3[max_element_index_mw3_X-1];
				q_right = vec_X_mw3[max_element_index_mw3_X+1];
				}
			if (q_left == 0){
				q_left=1;
				}
			if(q_right==0){
				q_right=1;
				}
			peaks_X_mw3.push_back(get_mw3_pos_X(q_max,max_element_index_mw3_X,q_left,q_right));
			if (max_element_index_mw3_X == 0){
				vec_X_mw3[max_element_index_mw3_X] = 0;
				vec_X_mw3[max_element_index_mw3_X+1] = 0;
				}	
			else if (max_element_index_mw3_X == 287){
				vec_X_mw3[max_element_index_mw3_X] = 0;
				vec_X_mw3[max_element_index_mw3_X-1] = 0;
				}	
			else {
				vec_X_mw3[max_element_index_mw3_X] = 0;
				vec_X_mw3[max_element_index_mw3_X+1] = 0;
				vec_X_mw3[max_element_index_mw3_X-1] = 0;
				}
			}
		}

	if (peaks_X_mw3.size() == 1 && tofsingletcal_clone->GetEntriesFast() == 1){
		R3BSofTofWSingleTcalData** tofsingletcal_data;
		tofsingletcal_data = new R3BSofTofWSingleTcalData*[1];
		tofsingletcal_data[0] = (R3BSofTofWSingleTcalData*)tofsingletcal_clone->At(0);
		Int_t paddle_ID = tofsingletcal_data[0]->GetDetector();
		Double_t posy_ns = tofsingletcal_data[0]->GetRawPosNs();
		Double_t tof_time = tofsingletcal_data[0]->GetRawTofNs();
		tof_time = tof_time - tof_offset[paddle_ID-1] + tof_lise;
		Double_t calibrated_y_pos = posy_ns*v_para_tof_y[paddle_ID][1]+v_para_tof_y[paddle_ID][2];
		vector<double> tof_hit{tof_time,calibrated_y_pos,peaks_X_mw3[0],(Double_t) paddle_ID};
		delete [] tofsingletcal_data;
		return tof_hit;
		}
	else 
		return {-1000,-1000,-1000,-1000};

	}
