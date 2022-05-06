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
		if (max_element_index_mw3_X == 0) q_left=1;
		else if (max_element_index_mw3_X == 287) q_right=1;
		else {
			q_left = vec_X_mw3[max_element_index_mw3_X-1];
			q_right = vec_X_mw3[max_element_index_mw3_X+1];
			}
		if (q_left == 0) q_left=1;
		if(q_right==0) q_right=1;
		peaks_X_mw3.push_back(get_mw3_pos_X(q_max,max_element_index_mw3_X,q_left,q_right));
		if (max_element_index_mw3_X == 0){
			vec_X_mw3[max_element_index_mw3_X] = 0;
			vec_X_mw3[max_element_index_mw3_X+1] = 0;
			}	
		if (max_element_index_mw3_X == 287){
			vec_X_mw3[max_element_index_mw3_X] = 0;
			vec_X_mw3[max_element_index_mw3_X-1] = 0;
			}	
		else {
			vec_X_mw3[max_element_index_mw3_X] = 0;
			vec_X_mw3[max_element_index_mw3_X+1] = 0;
			vec_X_mw3[max_element_index_mw3_X-1] = 0;
			}
		}
	//MW3 Y
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
		if (max_element_index_mw3_Y == 287){
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
tofsingletcal_data = new R3BSofTofWSingleTcalData*[entries_tof_singletcal];
for (Int_t j = 0; j < entries_tof_singletcal; j++){
	tofsingletcal_data[j] = (R3BSofTofWSingleTcalData*)tofsingletcal_clone->At(j);
    Int_t paddle_ID = tofsingletcal_data[j]->GetDetector();
	Double_t posy_ns = tofsingletcal_data[j]->GetRawPosNs();
	Double_t tof_time = tofsingletcal_data[j]->GetRawTofNs();
	Double_t posx = 840./2. -15. -(Double_t)(paddle_ID-1)*30; //hard coded dimensions (in mm)
	tof_hits.push_back({posx,posy_ns,tof_time});
	}
//sort tof hits according to the y position
sort(tof_hits.begin(),tof_hits.end(),sortcol);
vector<double> tof_down = tof_hits[0];
vector<double> tof_up = tof_hits[1];
//sort mw3_x
if (peaks_X_mw3.size() == 2){
	sort(peaks_X_mw3.begin(),peaks_X_mw3.end());
	}
vector<vector<double> > positioning_after_glad;
if (tof_up[0] > tof_down[0]){
	vector<double> v_temp1 = {peaks_X_mw3[0],tof_down[0],tof_down[1],tof_down[2]};
	positioning_after_glad.push_back(v_temp1);
	vector<double> v_temp2 = {peaks_X_mw3[1],tof_up[0],tof_up[1],tof_up[2]};
	positioning_after_glad.push_back(v_temp2);
	}
else {
	vector<double> v_temp1 = {peaks_X_mw3[0],tof_up[0],tof_up[1],tof_up[2]};
	positioning_after_glad.push_back(v_temp1);
	vector<double> v_temp2 = {peaks_X_mw3[1],tof_down[0],tof_down[1],tof_down[1]};
	positioning_after_glad.push_back(v_temp2);
	}
	return positioning_after_glad; //returns x_pos_mw3,y_pos_tof,tof_time in increasing x order
}
