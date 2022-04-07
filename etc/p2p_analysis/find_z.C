//insert all needed headers
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
#include "all_params.C"

//insert all offsets/constants, geometries etc you need
//
//TH1F* h1_z_all_iso;
//char hist_name[500];

using namespace std;
//char const* par_file = "test_for_file.txt";
double find_z(char const * par_file,char const * root_file){
TH1F* h1_z_all_iso;
sprintf(hist_name, "Z for all isotopes");
h1_z_all_iso = new TH1F(hist_name,hist_name,1000,1,10);
h1_z_all_iso->GetXaxis()->SetTitle("Z (charge)");
h1_z_all_iso->GetYaxis()->SetTitle("# Events ");
h1_z_all_iso->GetXaxis()->CenterTitle(true);
h1_z_all_iso->GetYaxis()->CenterTitle(true);
h1_z_all_iso->GetYaxis()->SetLabelSize(0.045);
h1_z_all_iso->GetYaxis()->SetTitleSize(0.045);

fstream infile(par_file,std::ios::out | std::ios::in | std::ios_base::app);
if (infile.is_open() && !(infile.peek() == std::ifstream::traits_type::eof())){
	cout << "file exists and is not empty " << endl;
	string line, word, temp;
	string check_eof;
	vector<string> temp_vec;
	while (!infile.eof()){
		getline(infile,line);
		stringstream s(line);
		while(getline(s,word,',')){
			temp_vec.push_back(word);
			}
		if (temp_vec[0] == "z_cut"){
			cout << "it will work to read from file...!" << endl;
			if(temp_vec.size() > 1){
				infile.close();
				return stod(temp_vec[1]);
				}
			else { 
				break;
				}
			}	
		temp_vec.clear();
		}
	//infile.close();
	}

TChain* chain = new TChain("evt");
chain->Reset();
chain->Add(root_file);

TClonesArray* SofMwpc3HitData = new TClonesArray("R3BSofMwpcHitData",5);
R3BSofMwpcHitData** sofmwpc3hitdata;
TBranch *branchSofMwpc3HitData = chain->GetBranch("Mwpc3HitData");
branchSofMwpc3HitData->SetAddress(&SofMwpc3HitData);

TClonesArray* SofSciTcalData = new TClonesArray("R3BSofSciTcalData",2);
R3BSofSciTcalData** sofscitcaldata;
TBranch *branchSofSciTcalData = chain->GetBranch("SofSciTcalData");
branchSofSciTcalData->SetAddress(&SofSciTcalData);
//
TClonesArray* SofToFWTcalData = new TClonesArray("R3BSofToFWTcalData",2);
R3BSofToFWTcalData** softofwtcaldata;
TBranch *branchSofToFWTcalData = chain->GetBranch("SofToFWTcalData");
branchSofToFWTcalData->SetAddress(&SofToFWTcalData);
//
TClonesArray* SofMwpc0HitData = new TClonesArray("R3BSofMwpcHitData",5);
R3BSofMwpcHitData** sofmwpc0hitdata;
TBranch *branchSofMwpc0HitData = chain->GetBranch("Mwpc0HitData");
branchSofMwpc0HitData->SetAddress(&SofMwpc0HitData);
//
TClonesArray* SofMwpc1HitData = new TClonesArray("R3BSofMwpcHitData",5);
R3BSofMwpcHitData** sofmwpc1hitdata;
TBranch *branchSofMwpc1HitData = chain->GetBranch("Mwpc1HitData");
branchSofMwpc1HitData->SetAddress(&SofMwpc1HitData);
//
TClonesArray* SofMwpc2HitData = new TClonesArray("R3BSofMwpcHitData",5);
R3BSofMwpcHitData** sofmwpc2hitdata;
TBranch *branchSofMwpc2HitData = chain->GetBranch("Mwpc2HitData");
branchSofMwpc2HitData->SetAddress(&SofMwpc2HitData);
//
TClonesArray* SofToFWMappedData = new TClonesArray("R3BSofToFWMappedData",2);
R3BSofToFWMappedData** softofwmappeddata;
TBranch *branchSofToFWMappedData = chain->GetBranch("SofToFWMappedData");
branchSofToFWMappedData->SetAddress(&SofToFWMappedData);
//
TClonesArray* SofTwimHitData = new TClonesArray("R3BSofTwimHitData",2);
R3BSofTwimHitData** softwimhitdata;
TBranch *branchSofTwimHitData = chain->GetBranch("TwimHitData");
branchSofTwimHitData->SetAddress(&SofTwimHitData);
//
TClonesArray* CalifaHitData = new TClonesArray("R3BCalifaHitData",3);
R3BCalifaHitData** califahitdata;
TBranch *branchCalifaHitData = chain->GetBranch("CalifaHitData");
branchCalifaHitData->SetAddress(&CalifaHitData);

Long64_t nevents = chain->GetEntries();

for(Long64_t i=0;i< nevents;i++){
    Long64_t evtnr = i;
    if (i%100000==0)
        cout<<"Processing event for charge analysis "<<i<<endl;
    chain->GetEvent(i); //event_i_call
    entries_mw0 = SofMwpc0HitData->GetEntries();
    entries_mw1 = SofMwpc1HitData->GetEntries();
    entries_mw2 = SofMwpc2HitData->GetEntries();
    entries_mw3 = SofMwpc3HitData->GetEntries();
    entries_start = SofSciTcalData->GetEntriesFast();
    entries_tofw = SofToFWTcalData->GetEntriesFast();
    entries_califa = CalifaHitData->GetEntries();
    softwimhitdata = new R3BSofTwimHitData*[1];
    softwimhitdata[0] = (R3BSofTwimHitData*)SofTwimHitData->At(0);


	if (entries_califa >= 2 && entries_start ==2 && entries_tofw == 2 && entries_mw0 == 1 && entries_mw1 == 1 && entries_mw2 == 1 && entries_mw3 == 1 && softwimhitdata[0]){
		charge_val = softwimhitdata[0]->GetZcharge();
		sofmwpc0hitdata = new R3BSofMwpcHitData*[1];
        sofmwpc1hitdata = new R3BSofMwpcHitData*[1];
        sofmwpc2hitdata = new R3BSofMwpcHitData*[1];
        sofmwpc3hitdata = new R3BSofMwpcHitData*[1];
        softofwtcaldata = new R3BSofToFWTcalData*[1];
        softofwmappeddata = new R3BSofToFWMappedData*[2];
        califahitdata = new R3BCalifaHitData*[entries_califa];
        sofmwpc0hitdata[0] = (R3BSofMwpcHitData*)SofMwpc0HitData->At(0);
        sofmwpc1hitdata[0] = (R3BSofMwpcHitData*)SofMwpc1HitData->At(0);
        sofmwpc2hitdata[0] = (R3BSofMwpcHitData*)SofMwpc2HitData->At(0);
        sofmwpc3hitdata[0] = (R3BSofMwpcHitData*)SofMwpc3HitData->At(0);
        softofwtcaldata[0] = (R3BSofToFWTcalData*)SofToFWTcalData->At(0);
		califahitdata[0] = (R3BCalifaHitData*)CalifaHitData->At(0);
        Double_t xMW0 = sofmwpc0hitdata[0]->GetX();
        Double_t xMW1 = sofmwpc1hitdata[0]->GetX();
        Double_t xMW2 = sofmwpc2hitdata[0]->GetX();
        Double_t xMW3 = sofmwpc3hitdata[0]->GetX();
        Double_t yMW1 = sofmwpc1hitdata[0]->GetY();
        Double_t yMW2 = sofmwpc2hitdata[0]->GetY();
        Double_t yMW3 = sofmwpc3hitdata[0]->GetY();
		if (xMW1 != -1000 && xMW2 != -1000 && xMW3 != -1000  && xMW0 != -1000 && yMW3 != -1000){
			h1_z_all_iso->Fill(charge_val);
			}
		
		delete [] sofmwpc0hitdata;
        delete [] sofmwpc1hitdata;
        delete [] sofmwpc2hitdata;
        delete [] sofmwpc3hitdata;
        delete [] softofwtcaldata;
        delete [] softofwmappeddata;
        delete [] califahitdata;
		}
		delete [] softwimhitdata;
        }
	Int_t binmax_Z6 = h1_z_all_iso->GetMaximumBin();
	Double_t middle_Z6 = h1_z_all_iso->GetXaxis()->GetBinCenter(binmax_Z6);
	Double_t middle_Z5 = middle_Z6 -1;
	cout << "max bin x position" << middle_Z6 << endl;
	//fit for Z = 6
	TF1 *gauss_6 = new TF1("gauss_6","gaus",middle_Z6-0.5,middle_Z6+0.5);
	h1_z_all_iso->Fit("gauss_6","QR+");
	
	//fit for Z = 5
	TF1 *gauss_5 = new TF1("gauss_5","gaus",middle_Z5-0.5,middle_Z5+0.5);
	h1_z_all_iso->Fit("gauss_5","QR+");
	
	TF1 *fit_g6 = h1_z_all_iso->GetFunction("gauss_6");
	TF1 *fit_g5 = h1_z_all_iso->GetFunction("gauss_5");
	Double_t mean_g6 = fit_g6->GetParameter(1);
	Double_t mean_g5 = fit_g5->GetParameter(1);
	Double_t offs_g6 = fit_g6->GetParameter(0);
	Double_t offs_g5 = fit_g5->GetParameter(0);
	Double_t sigma_g6 = fit_g6->GetParameter(2);
	Double_t sigma_g5 = fit_g5->GetParameter(2);

	TF1 *f1 = new TF1("f1","abs([0]*exp(-0.5*((x-[1])/[2])^2)-([3]*exp(-0.5*((x-[4])/[5])^2)))",mean_g5,mean_g6 );
	f1->SetParameters(offs_g5,mean_g5,sigma_g5,offs_g6,mean_g6,sigma_g6);
	f1->Print("V");
	cout << "new min" << f1->GetMinimumX() << endl;
	const Double_t charge_cut = f1->GetMinimumX();

	delete gauss_6;
	delete gauss_5;
	delete f1;	
	infile.clear();	
	infile << "z_cut" << "," << charge_cut << endl;
	infile.close();
	return charge_cut;
}

