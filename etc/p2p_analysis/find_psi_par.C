#include "all_params.C"
using namespace std;
vector<double> find_psi_par(char const * par_file,char const * root_file, double const z_cut,vector<double> const iso_cuts){
vector<double> dv;
TH2F* h2_theta_out_vs_mw3_11b;
sprintf(hist_name, "Theta_out versus MWPC3.fX for 11B");
h2_theta_out_vs_mw3_11b = new TH2F(hist_name,hist_name,800,-400,400,1200,0,0.6);
h2_theta_out_vs_mw3_11b->GetXaxis()->SetTitle("MWPC3 f.X [mm]");
h2_theta_out_vs_mw3_11b->GetYaxis()->SetTitle("Theta_out [degrees]");
h2_theta_out_vs_mw3_11b->GetXaxis()->CenterTitle(true);
h2_theta_out_vs_mw3_11b->GetYaxis()->CenterTitle(true);
h2_theta_out_vs_mw3_11b->GetYaxis()->SetLabelSize(0.045);
h2_theta_out_vs_mw3_11b->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_theta_out_vs_mw3_10b;
sprintf(hist_name, "Theta_out versus MWPC3.fX for 10B");
h2_theta_out_vs_mw3_10b = new TH2F(hist_name,hist_name,800,-400,400,1200,0,0.6);
h2_theta_out_vs_mw3_10b->GetXaxis()->SetTitle("MWPC3 f.X [mm]");
h2_theta_out_vs_mw3_10b->GetYaxis()->SetTitle("Theta_out [degrees]");
h2_theta_out_vs_mw3_10b->GetXaxis()->CenterTitle(true);
h2_theta_out_vs_mw3_10b->GetYaxis()->CenterTitle(true);
h2_theta_out_vs_mw3_10b->GetYaxis()->SetLabelSize(0.045);
h2_theta_out_vs_mw3_10b->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_theta_out_vs_mw3_12c;
sprintf(hist_name, "Theta_out versus MWPC3.fX for 12C");
h2_theta_out_vs_mw3_12c = new TH2F(hist_name,hist_name,800,-400,400,1200,0,0.6);
h2_theta_out_vs_mw3_12c->GetXaxis()->SetTitle("MWPC3 f.X [mm]");
h2_theta_out_vs_mw3_12c->GetYaxis()->SetTitle("Theta_out [degrees]");
h2_theta_out_vs_mw3_12c->GetXaxis()->CenterTitle(true);
h2_theta_out_vs_mw3_12c->GetYaxis()->CenterTitle(true);
h2_theta_out_vs_mw3_12c->GetYaxis()->SetLabelSize(0.045);
h2_theta_out_vs_mw3_12c->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_theta_out_vs_mw3_11c;
sprintf(hist_name, "Theta_out versus MWPC3.fX for 11C");
h2_theta_out_vs_mw3_11c = new TH2F(hist_name,hist_name,800,-400,400,1200,0,0.6);
h2_theta_out_vs_mw3_11c->GetXaxis()->SetTitle("MWPC3 f.X [mm]");
h2_theta_out_vs_mw3_11c->GetYaxis()->SetTitle("Theta_out [degrees]");
h2_theta_out_vs_mw3_11c->GetXaxis()->CenterTitle(true);
h2_theta_out_vs_mw3_11c->GetYaxis()->CenterTitle(true);
h2_theta_out_vs_mw3_11c->GetYaxis()->SetLabelSize(0.045);
h2_theta_out_vs_mw3_11c->GetYaxis()->SetTitleSize(0.045);


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
        if (temp_vec[0] == "psi_cut"){
            cout << "we can also read psi_cuts from file, juhuu...!" << endl;
            if(temp_vec.size() == 13){
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

Int_t entries_10b = 0;
Int_t entries_11b = 0;
Int_t entries_11c = 0;
Int_t entries_12c = 0;


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

//
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
            if ((sofmwpc3hitdata[0]->GetY() > -500.) && (sofmwpc3hitdata[0]->GetX() > -500.)) 
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
if (charge_val< z_cut && charge_val> z_cut-1.) //Z = 5 before it was charge_val< charge_cut && charge_val> charge_cut-0.8, too straight....
{
                        if (a_q > iso_cuts[0] && a_q < (iso_cuts[0]+2*iso_cuts[3])){  //11B
                            h2_theta_out_vs_mw3_11b->Fill(xMW3,psi_out);
                            entries_11b++;
                        }
                        if (a_q < iso_cuts[0] && a_q > (iso_cuts[0]-0.15)){  //10B
                            h2_theta_out_vs_mw3_10b->Fill(xMW3,psi_out);
                            entries_10b++;
                        }
                        }
                        if (charge_val< z_cut+0.8 && charge_val> z_cut) //Z = 6
                        {
                        if (a_q > iso_cuts[1] && a_q < (iso_cuts[1]+2*iso_cuts[5])){  //12C
                            h2_theta_out_vs_mw3_12c->Fill(xMW3,psi_out);
                            entries_12c++;
                        }
                        if (a_q < iso_cuts[1]-0.06 && a_q > (iso_cuts[1]-1.6*iso_cuts[4])){  //11C before it wasa_q < aq_cut_z_6 && a_q > (aq_cut_z_6-2*width_11c), bad for 11C ...
                            h2_theta_out_vs_mw3_11c->Fill(xMW3,psi_out);
                            entries_11c++;

                        }

                        }
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
//now get ranges for the fits:
////10B
Double_t st_dev_mw3_10b = h2_theta_out_vs_mw3_10b->GetStdDev(1);
Double_t st_dev_psi_out_10b = h2_theta_out_vs_mw3_10b->GetStdDev(2);
Double_t mean_mw3_10b = h2_theta_out_vs_mw3_10b->GetMean(1);
const Double_t mean_psi_out_10b = h2_theta_out_vs_mw3_10b->GetMean(2);

//11B
Double_t st_dev_mw3_11b = h2_theta_out_vs_mw3_11b->GetStdDev(1);
Double_t st_dev_psi_out_11b = h2_theta_out_vs_mw3_11b->GetStdDev(2);
Double_t mean_mw3_11b = h2_theta_out_vs_mw3_11b->GetMean(1);
const Double_t mean_psi_out_11b = h2_theta_out_vs_mw3_11b->GetMean(2);

//12C
Double_t st_dev_mw3_12c = h2_theta_out_vs_mw3_12c->GetStdDev(1);
Double_t st_dev_psi_out_12c = h2_theta_out_vs_mw3_12c->GetStdDev(2);
Double_t mean_mw3_12c = h2_theta_out_vs_mw3_12c->GetMean(1);
const Double_t mean_psi_out_12c = h2_theta_out_vs_mw3_12c->GetMean(2);

//11C
//
Double_t st_dev_mw3_11c = h2_theta_out_vs_mw3_11c->GetStdDev(1);
Double_t st_dev_psi_out_11c = h2_theta_out_vs_mw3_11c->GetStdDev(2);
Double_t mean_mw3_11c = h2_theta_out_vs_mw3_11c->GetMean(1);
const Double_t mean_psi_out_11c = h2_theta_out_vs_mw3_11c->GetMean(2);

TF1 *angle_fit_10b = new TF1("angle_fit_10b","pol1",mean_mw3_10b-2*st_dev_mw3_10b,mean_mw3_10b+2*st_dev_mw3_10b);
TF1 *angle_fit_11b = new TF1("angle_fit_11b","pol1",mean_mw3_11b-2*st_dev_mw3_11b,mean_mw3_11b+2*st_dev_mw3_11b);
TF1 *angle_fit_11c = new TF1("angle_fit_11c","pol1",mean_mw3_11c-2*st_dev_mw3_11c,mean_mw3_11c+2*st_dev_mw3_11c);
TF1 *angle_fit_12c = new TF1("angle_fit_12c","pol1",mean_mw3_12c-2*st_dev_mw3_12c,mean_mw3_12c+2*st_dev_mw3_12c);

//psi parameters for 10B:
h2_theta_out_vs_mw3_10b->Fit("angle_fit_10b", "QR+");
TF1* fit_angle_10b = h2_theta_out_vs_mw3_10b->GetFunction("angle_fit_10b");
const Double_t angle_offs_10b = fit_angle_10b->GetParameter(0);
const Double_t slope_10b = fit_angle_10b->GetParameter(1);

//psi_parameters for 11B:
h2_theta_out_vs_mw3_11b->Fit("angle_fit_11b", "QR+");
TF1* fit_angle_11b = h2_theta_out_vs_mw3_11b->GetFunction("angle_fit_11b");
const Double_t angle_offs_11b = fit_angle_11b->GetParameter(0);
const Double_t slope_11b = fit_angle_11b->GetParameter(1);

//psi_parameters for 11C:
h2_theta_out_vs_mw3_11c->Fit("angle_fit_11c", "QR+");
TF1* fit_angle_11c = h2_theta_out_vs_mw3_11c->GetFunction("angle_fit_11c");
const Double_t angle_offs_11c = fit_angle_11c->GetParameter(0);
const Double_t slope_11c = fit_angle_11c->GetParameter(1);

//psi_parameters for 12C:
h2_theta_out_vs_mw3_12c->Fit("angle_fit_12c", "QR+");
TF1* fit_angle_12c = h2_theta_out_vs_mw3_12c->GetFunction("angle_fit_12c");
const Double_t angle_offs_12c = fit_angle_12c->GetParameter(0);
const Double_t slope_12c = fit_angle_12c->GetParameter(1);

delete angle_fit_10b;
delete angle_fit_11b;
delete angle_fit_11c;
delete angle_fit_12c;

vector<double> psi_v(12);
//for psi_v I want to insert mean_psi_out, slope, angle_offs (for 10b,11b,11c,12c) -> 12 entries
psi_v[0] = mean_psi_out_10b;
psi_v[1] = slope_10b;
psi_v[2] = angle_offs_10b;

psi_v[3] = mean_psi_out_11b;
psi_v[4] = slope_11b;
psi_v[5] = angle_offs_11b;

psi_v[6] = mean_psi_out_11c;
psi_v[7] = slope_11c;
psi_v[8] = angle_offs_11c;

psi_v[9] = mean_psi_out_12c;
psi_v[10] = slope_12c;
psi_v[11] = angle_offs_12c;
infile.clear();
infile << "psi_cut" << "," << mean_psi_out_10b << "," << slope_10b << "," << angle_offs_10b << "," << mean_psi_out_11b << "," << slope_11b << "," << angle_offs_11b << "," << mean_psi_out_11c << "," << slope_11c << "," << angle_offs_11c << "," << mean_psi_out_12c << "," << slope_12c << "," << angle_offs_12c << endl;
return psi_v;
}
