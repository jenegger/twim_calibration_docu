#include "all_params.C"
using namespace std;

//TH1F* h1_a_q_z5;
//TH1F* h1_a_q_z6;
//TH2F* h2_z_vs_a_q_nocut;
vector<double> find_isotopes (char const * par_file,char const * root_file, double const z_cut){
vector<double> dv;
TH1F* h1_a_q_z5;
sprintf(hist_name, "A/q for Z = 5");
h1_a_q_z5 = new TH1F(hist_name,hist_name,400,1,3);
h1_a_q_z5->GetXaxis()->SetTitle("A/q for Z = 5");
h1_a_q_z5->GetYaxis()->SetTitle(" # Events");
h1_a_q_z5->GetXaxis()->CenterTitle(true);
h1_a_q_z5->GetYaxis()->CenterTitle(true);
h1_a_q_z5->GetYaxis()->SetLabelSize(0.045);
h1_a_q_z5->GetYaxis()->SetTitleSize(0.045);
TH1F* h1_a_q_z6;
sprintf(hist_name, "A/q for Z = 6");
h1_a_q_z6 = new TH1F(hist_name,hist_name,400,1,3);
h1_a_q_z6->GetXaxis()->SetTitle("A/q for Z = 6");
h1_a_q_z6->GetYaxis()->SetTitle(" # Events");
h1_a_q_z6->GetXaxis()->CenterTitle(true);
h1_a_q_z6->GetYaxis()->CenterTitle(true);
h1_a_q_z6->GetYaxis()->SetLabelSize(0.045);
h1_a_q_z6->GetYaxis()->SetTitleSize(0.045);
TH2F* h2_z_vs_a_q_nocut;
sprintf(hist_name, "Z versus A/q all isotopes");
h2_z_vs_a_q_nocut = new TH2F(hist_name,hist_name,200,1,3,100,0,10);
h2_z_vs_a_q_nocut->GetXaxis()->SetTitle("A/q");
h2_z_vs_a_q_nocut->GetYaxis()->SetTitle("Z (charge) ");
h2_z_vs_a_q_nocut->GetXaxis()->CenterTitle(true);
h2_z_vs_a_q_nocut->GetYaxis()->CenterTitle(true);
h2_z_vs_a_q_nocut->GetYaxis()->SetLabelSize(0.045);
h2_z_vs_a_q_nocut->GetYaxis()->SetTitleSize(0.045);


fstream infile(par_file,std::ios::out | std::ios::in | std::ios_base::app);
if (infile.is_open() && !(infile.peek() == std::ifstream::traits_type::eof())){
    cout << "file exists " << endl;
    string line, word, temp;
    string check_eof;
    vector<string> temp_vec;
    while (!infile.eof()){
        getline(infile,line);
        stringstream s(line);
        while(getline(s,word,',')){
            temp_vec.push_back(word);
            }
        if (temp_vec[0] == "isotopes_cut"){
            cout << "we can also read isotopes from file, juhuu...!" << endl;
            if(temp_vec.size() == 7){
				temp_vec.erase(temp_vec.begin(),temp_vec.begin()+1);
				transform(temp_vec.begin(), temp_vec.end(), back_inserter(dv), [](const string & astr){ return stod( astr) ; } ) ;
                infile.close();
                return dv;
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
        cout<<"Processing event "<<i<<endl;
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

if (xMW1 != -1000 && xMW2 != -1000 && xMW3 != -1000  && xMW0 != -1000){
        Double_t psi_in = atan(((xMW2+xMW2_shift)-(xMW1+xMW1_shift))/(zM2-zM1));
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


        tot_length = start_to_target+ abs((zB-zT)/(cos(psi_in)))+rho*w+abs((zToFW-z_D)/cos(psi_out));

        delete gr1;
        delete gr2;
        delete gr3;
        if (SofToFWTcalData && SofToFWTcalData->GetEntriesFast())
            {
            for (UShort_t i = 0; i < NbDets; i++)
                {
                for (UShort_t j = 0; j < NbChs; j++)
                    {
                    mult[i * NbChs + j] = 0;
                    }
                }
            for (Int_t ihit = 0; ihit < entries_tofw; ihit++)
                {
                R3BSofToFWMappedData* hitmapped = (R3BSofToFWMappedData*)SofToFWMappedData->At(ihit);
                if (!hitmapped)
                    continue;
                iDet = hitmapped->GetDetector() - 1;
                iCh = hitmapped->GetPmt() - 1;
                mult[iDet * NbChs + iCh]++;
                }
            for (Int_t ihit = 0; ihit < entries_tofw; ihit++)
                {
                R3BSofToFWTcalData* hittcal = (R3BSofToFWTcalData*)SofToFWTcalData->At(ihit);
                if (!hittcal)
                    continue;
                if (hittcal->GetPmt() == 3)
                    continue;
                iDet = hittcal->GetDetector() - 1;
                iCh = hittcal->GetPmt() - 1;
                iRawTimeNs[iDet * 2 + iCh] = hittcal->GetRawTimeNs();
                }
            for (Int_t s_hit = 0 ; s_hit < entries_start; s_hit++)
                {
                R3BSofSciTcalData* start_hits = (R3BSofSciTcalData*)SofSciTcalData->At(s_hit);
                if(!start_hits)
                    continue;
                iRawTimeStartNs[s_hit] = start_hits->GetRawTimeNs();

                }
                sofmwpc3hitdata[0]=(R3BSofMwpcHitData*)SofMwpc3HitData->At(0);
            if ((sofmwpc3hitdata[0]->GetY() > -500.) && (sofmwpc3hitdata[0]->GetX() > -500.))  //cut on the mwpc3
            {
            Double_t raw_mwpc3 = sofmwpc3hitdata[0]->GetY();
            Double_t tofpos_ns = -1000;
            Double_t raw_tofpos_ns = -1000;
            Double_t raw_t_start = -1000000.;
            Double_t raw_time_of_flight = 0.;

            for (UShort_t i = 0; i < NbDets; i++)
                {
                if ((mult[i * NbChs] == 1) && (mult[i * NbChs + 1] == 1))
                    {
                    tofpos_ns = 0.5*(iRawTimeNs[i * NbChs + 1] - iRawTimeNs[i * NbChs]-offsets[i+1]);
                    raw_tofpos_ns = 0.5*(iRawTimeNs[i * NbChs + 1] - iRawTimeNs[i * NbChs]);
                    raw_t_start = 0.5*(iRawTimeStartNs[0]-iRawTimeStartNs[1]);
                    if (raw_t_start !=-1000000. && raw_tofpos_ns != -1000 && tofpos_ns != -1000)
                        {
                        raw_time_of_flight=(0.5*(iRawTimeNs[i*NbChs+1]+iRawTimeNs[i*NbChs]))-0.5*(iRawTimeStartNs[0]+iRawTimeStartNs[1]);
                        Double_t time_target_tof = raw_time_of_flight+tof_offs[i] -time_start_target;
                        Double_t path_from_target = abs((zB-zT)/(cos(psi_in)))+rho*w+abs((zToFW-z_D)/cos(psi_out));
                        Double_t beta = ((path_from_target/time_target_tof)*pow(10,6))/light_c;
                        Double_t mag_field = current*0.0006527728074785267;
                        Double_t a_q = (pow(10,-3)*((((mag_field*rho)/(beta))*1.602176634*pow(10,-19))/(1.66053906660*pow(10,-27)*light_c)))/gamma_given;

                        if (charge_val< z_cut && charge_val> z_cut-1.) //Z = 5
                        {
                        h1_a_q_z5->Fill(a_q);
                        }
                        if (charge_val< z_cut+0.8 && charge_val> z_cut) //Z = 6
                        {
                        h1_a_q_z6->Fill(a_q);
                        }
                        h2_z_vs_a_q_nocut->Fill(a_q,charge_val);
}
}
}
}
}
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

TF1* gauss_10B = new TF1("gauss_10B","gaus",1.85,2.03);
TF1* gauss_11B = new TF1("gauss_11B","gaus",2.1,2.3);
TF1* gauss_11C = new TF1("gauss_11C","gaus",1.7,1.88);
TF1* gauss_12C = new TF1("gauss_12C","gaus",1.9,2.1);


h1_a_q_z6->Fit("gauss_11C","QR+");
h1_a_q_z6->Fit("gauss_12C","QR+");

TF1* fit_11C = h1_a_q_z6->GetFunction("gauss_11C");
TF1* fit_12C = h1_a_q_z6->GetFunction("gauss_12C");
Double_t offs_11C = fit_11C->GetParameter(0);
Double_t offs_12C = fit_12C->GetParameter(0);
Double_t mean_11C = fit_11C->GetParameter(1);
Double_t mean_12C = fit_12C->GetParameter(1);
Double_t sigma_11C = fit_11C->GetParameter(2);
Double_t sigma_12C = fit_12C->GetParameter(2);
TF1 *f2 = new TF1("f2","abs([0]*exp(-0.5*((x-[1])/[2])^2)-([3]*exp(-0.5*((x-[4])/[5])^2)))",mean_11C,mean_12C);
f2->SetParameters(offs_11C,mean_11C,sigma_11C,offs_12C,mean_12C,sigma_12C);
const Double_t aq_cut_z_6 = f2->GetMinimumX();
const Double_t width_11c = abs(mean_11C-aq_cut_z_6);
const Double_t width_12c = abs(mean_12C-aq_cut_z_6);
cout << "cut on z6" << aq_cut_z_6 << endl;
h1_a_q_z5->Fit("gauss_10B","QR+");
h1_a_q_z5->Fit("gauss_11B","QR+");
TF1* fit_10B = h1_a_q_z5->GetFunction("gauss_10B");
TF1* fit_11B = h1_a_q_z5->GetFunction("gauss_11B");
Double_t offs_10B = fit_10B->GetParameter(0);
Double_t offs_11B = fit_11B->GetParameter(0);
Double_t mean_10B = fit_10B->GetParameter(1);
Double_t mean_11B = fit_11B->GetParameter(1);
Double_t sigma_10B = fit_10B->GetParameter(2);
Double_t sigma_11B = fit_11B->GetParameter(2);
TF1 *f3 = new TF1("f3","abs([0]*exp(-0.5*((x-[1])/[2])^2)-([3]*exp(-0.5*((x-[4])/[5])^2)))",mean_10B,mean_11B);
f3->SetParameters(offs_10B,mean_10B,sigma_10B,offs_11B,mean_11B,sigma_11B);
f3->Print("V");
const Double_t aq_cut_z_5 = f3->GetMinimumX();
const Double_t width_10b = abs(mean_10B-aq_cut_z_5);
const Double_t width_11b = abs(mean_11B-aq_cut_z_5);
cout << "cut on z5" << aq_cut_z_5 << endl;


delete gauss_10B;
delete gauss_11B;
delete gauss_11C;
delete gauss_12C;
delete f3;
delete f2;

//I want to return the vector as: cut_5,cut_6,width10B,width11B,width11C,width12C
vector<double> calc_v(6);
calc_v[0] = aq_cut_z_5;
calc_v[1] = aq_cut_z_6;
calc_v[2] = width_10b;
calc_v[3] = width_11b;
calc_v[4] = width_11c;
calc_v[5] = width_12c;

infile.clear();
infile << "isotopes_cut" << "," << aq_cut_z_5 << "," << aq_cut_z_6 << "," <<  width_10b << "," << width_11b << "," << width_11c << "," << width_12c << endl;
infile.close();
return calc_v;
}
