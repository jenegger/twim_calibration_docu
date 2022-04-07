#include "all_params.C"
using namespace std;

void analysis(char const * par_file,char const * root_file, double const z_cut,vector<double> const iso_cuts,vector<double> const psi_par){

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


cout << "finally precise calculation!---------------------------------" << endl;

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
    if (entries_start ==2 && entries_tofw == 2 && entries_mw0 == 1 && entries_mw1 == 1 && entries_mw2 == 1 && entries_mw3 == 1 && softwimhitdata[0]){
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
        Double_t xMW0 = sofmwpc0hitdata[0]->GetX();
        Double_t xMW1 = sofmwpc1hitdata[0]->GetX();
        Double_t xMW2 = sofmwpc2hitdata[0]->GetX();
        Double_t xMW3 = sofmwpc3hitdata[0]->GetX();
        Double_t yMW1 = sofmwpc1hitdata[0]->GetY();
        Double_t yMW2 = sofmwpc2hitdata[0]->GetY();
        Double_t yMW3 = sofmwpc3hitdata[0]->GetY();
if (xMW1 != -1000 && xMW2 != -1000 && xMW3 != -1000  && xMW0 != -1000){
        Double_t psi_in = atan(((xMW2+xMW2_shift)-(xMW1+xMW1_shift))/(zM2-zM1));
        //------------------------------------------------------------------------------
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
                        Double_t tof_diff_up_down = 0.5*(iRawTimeNs[i * NbChs + 1] - iRawTimeNs[i * NbChs]);
                        raw_time_of_flight=(0.5*(iRawTimeNs[i*NbChs+1]+iRawTimeNs[i*NbChs]))-0.5*(iRawTimeStartNs[0]+iRawTimeStartNs[1]);
                        Double_t time_target_tof = raw_time_of_flight+tof_offs[i] -time_start_target;
                        Double_t path_from_target = abs((zB-zT)/(cos(psi_in)))+rho*w+abs((zToFW-z_D)/cos(psi_out));
                        Double_t beta = ((path_from_target/time_target_tof)*pow(10,6))/light_c;
                        Double_t mag_field = current*0.0006527728074785267;
                        Double_t a_q = (pow(10,-3)*((((mag_field*rho)/(beta))*1.602176634*pow(10,-19))/(1.66053906660*pow(10,-27)*light_c)))/gamma_given;
if (charge_val< charge_cut && charge_val> charge_cut-1.) //Z = 5 before it was charge_val< charge_cut && charge_val> charge_cut-0.8, too straight ....
{
if (a_q > aq_cut_z_5 && a_q < (aq_cut_z_5+2*width_11b) && entries_califa >= 2){  //11B
                        psi_in = mean_psi_out_11b -slope_11b*xMW3 - angle_offs_11b;
                        m = (cos(psi_out)-cos(psi_in))/(sin(psi_out)+sin(psi_in));
                        rho =  (Leff/(2*sin((psi_in+psi_out)/2)))*sqrt(pow(((m+tan(alpha_G))/(1-m*tan(alpha_G))),2)+1);
                        w = 2*abs(asin(BD/(2*rho)));
                        Double_t pathlength_precise = abs((zB-zT)/(cos(psi_in)))+rho*w+abs((zToFW-z_D)/cos(psi_out));
                        Double_t beta_precise = ((pathlength_precise/time_target_tof)*pow(10,6))/light_c;
                        Double_t gamma_precise = 1/(sqrt(1-beta_precise*beta_precise));
                        Double_t a_q_precise = (pow(10,-3)*((((mag_field*rho)/(beta_precise))*1.602176634*pow(10,-19))/(1.66053906660*pow(10,-27)*light_c)))/gamma_given;
                        if (charge_val< charge_cut && charge_val> charge_cut-1.) //Z = 5 before it was charge_val< charge_cut && charge_val> charge_cut-0.8, too straight ....
{
                        if (a_q > aq_cut_z_5 && a_q < (aq_cut_z_5+2*width_11b) && entries_califa >= 2){  //11B
                        psi_in = mean_psi_out_11b -slope_11b*xMW3 - angle_offs_11b;
                        m = (cos(psi_out)-cos(psi_in))/(sin(psi_out)+sin(psi_in));
                        rho =  (Leff/(2*sin((psi_in+psi_out)/2)))*sqrt(pow(((m+tan(alpha_G))/(1-m*tan(alpha_G))),2)+1);
                        w = 2*abs(asin(BD/(2*rho)));
                        Double_t pathlength_precise = abs((zB-zT)/(cos(psi_in)))+rho*w+abs((zToFW-z_D)/cos(psi_out));
                        Double_t beta_precise = ((pathlength_precise/time_target_tof)*pow(10,6))/light_c;
                        Double_t gamma_precise = 1/(sqrt(1-beta_precise*beta_precise));
                        Double_t a_q_precise = (pow(10,-3)*((((mag_field*rho)/(beta_precise))*1.602176634*pow(10,-19))/(1.66053906660*pow(10,-27)*light_c)))/gamma_given;
                        Double_t mom_11b = 11*931.5*(beta_precise/(sqrt(1-(beta_precise*beta_precise))));
                        h2_z_vs_a_q->Fill(a_q_precise,charge_val);
                        h2_tof_vs_aq_fix_g_11b->Fill(a_q_precise,time_target_tof);
						
						//CALIFA handling
                        Double_t psi_CALIFA;
                        Double_t theta_CALIFA;
                        Double_t time_CALIFA;
                        Double_t E_Califa;
                        vector<double> vE1;
                        vector<double> vtheta1CALIFA;
                        vector<double> vpsi1CALIFA;
                        vector<double> vE2;
                        vector<double> vtheta2CALIFA;
                        vector<double> vpsi2CALIFA;
                        vector<double> vCALIFA_time_1;
                        vector<double> vCALIFA_time_2;
                        vector<double> v_E_gamma;
                        vector<double> vpsi;
                        vector<double> vtheta;
                        vector<double> vE;
                        vector<double> v_gamma_sorted;
                        vector<double> vtime;
                        vector<double> v_origin_E;
                        Double_t ind_high_energy;
                        Double_t ind_second_high_energy;
                        Double_t ind_first_gamma;
                        Double_t ind_second_gamma;
                        Double_t theta_gamma;
                        Double_t psi1;
                        Double_t psi2;
                        Double_t psi1_2pi;
                        Double_t psi2_2pi;
                        Double_t psi1_2pi_degr;
                        Double_t psi2_2pi_degr;
                        Double_t theta1_degr;
                        Double_t theta2_degr;
                        Double_t complete_gamma_deposited = 0.;
                        for (Int_t j = 0.;j<entries_califa;j++){
                            califahitdata[j] = (R3BCalifaHitData*)CalifaHitData->At(j);
                            psi_CALIFA = califahitdata[j]->GetPhi();
                            theta_CALIFA = califahitdata[j]->GetTheta();
                            E_Califa = califahitdata[j]->GetEnergy();
                            time_CALIFA = califahitdata[j]->GetTime();
                            vpsi.push_back(psi_CALIFA);
                            vtheta.push_back(theta_CALIFA);
                            vE.push_back(E_Califa);
                            v_origin_E.push_back(E_Califa);
                            vtime.push_back(time_CALIFA);
                            }
                        sort(vE.begin(),vE.end());
                        if (vE[vE.size()-1] >30000 && vE[vE.size()-2] > 30000)
                        {
                        vector<vector<double> > energy_theta_v;
                        for (Int_t k =0;k< vE.size();k++)
                        {
                            if (v_origin_E[k] < 10000){
                            vector<double> temp_vec;
                            temp_vec.push_back(vtheta[k]);
                            temp_vec.push_back((v_origin_E[k]/1000)*(1/(sqrt(1-beta_precise*beta_precise)))*(1-beta_precise*cos(vtheta[k])));   //column 1 for energy
                            energy_theta_v.push_back(temp_vec);
                            temp_vec.clear();
                            }
                            if (v_origin_E[k] == vE[vE.size()-1]) {
                            ind_high_energy = k;
                            }
                            else if (v_origin_E[k] == vE[vE.size()-2]){
                            ind_second_high_energy = k;
                            }

                            else if(v_origin_E[k] < 10000)
                            {
                            theta_gamma = vtheta[k];
                            Double_t energy_gamma = (v_origin_E[k]/1000)*(1/(sqrt(1-beta_precise*beta_precise)))*(1-beta_precise*cos(theta_gamma));
                            v_E_gamma.push_back(energy_gamma);
                            h1_gamma_energyE_11B->Fill(energy_gamma);
                            complete_gamma_deposited +=v_origin_E[k]/1000;
                            if (vpsi[k]< 0.){
                            h2_gamma_energyE_psi_11B->Fill(energy_gamma,((2*PI+vpsi[k])/PI)*180);
                            }
                            if (vpsi[k] > 0.){
                            h2_gamma_energyE_psi_11B->Fill(energy_gamma,(vpsi[k]/PI)*180);

                            }
                            }



                        }
                        Double_t E1 = v_origin_E[ind_high_energy]/1000.;
                        Double_t E2 = v_origin_E[ind_second_high_energy]/1000.;
                        Double_t theta1 = vtheta[ind_high_energy];
                        Double_t theta2 = vtheta[ind_second_high_energy];
                        cout << "------------------------------------------------------" << endl;
                        cout << "------------------------------------------------------" << endl;
                        cout << "phi1 = \t" << (vpsi[ind_high_energy]/PI)*180. << endl;
                        cout << "phi2 = \t" << (vpsi[ind_second_high_energy]/PI)*180. << endl;
                        cout << "------------------------------------------------------" << endl;
                        if (vpsi[ind_high_energy] < 0.){
                        psi1 = PI - vpsi[ind_high_energy];
                        psi1_2pi = 2*PI+vpsi[ind_high_energy];
                        }
                        if (vpsi[ind_high_energy] > 0.){
                        psi1_2pi = vpsi[ind_high_energy];
                        }
                        if (vpsi[ind_second_high_energy] < 0.){
                        psi2 = PI - vpsi[ind_second_high_energy];
                        psi2_2pi = 2*PI+vpsi[ind_second_high_energy];
                        }
                        if (vpsi[ind_second_high_energy] > 0.){
                        psi2_2pi = vpsi[ind_second_high_energy];
                        }

                        psi1_2pi_degr = (psi1_2pi/PI)*180.;
                        psi2_2pi_degr = (psi2_2pi/PI)*180.;
                        theta1_degr = (theta1/PI)*180.;
                        theta2_degr = (theta2/PI)*180.;

                        h1_abs_phi_diff_11b->Fill(abs(psi1_2pi_degr-psi2_2pi_degr));


                        Double_t mom_p1 = sqrt(pow((E1+938.272),2)-pow(938.272,2));
                        Double_t mom_p2 = sqrt(pow((E2+938.272),2)-pow(938.272,2));
                        Double_t beta_1 = mom_p1/(sqrt(mom_p1*mom_p1+938.272*938.272));
                        Double_t gamma_1 = 1/(sqrt(1-beta_1*beta_1));
                        Double_t beta_2 = mom_p2/(sqrt(mom_p2*mom_p2+938.272*938.272));               
                        Double_t gamma_2 = 1/(sqrt(1-beta_2*beta_2));
                        Double_t angle_p1_11b = abs(atan(tan(theta1)*cos(psi1_2pi) )-psi_in);
                        Double_t angle_p2_11b = abs(atan(tan(theta2)*cos(psi2_2pi) )-psi_in);
                        Double_t angle_p2_p1 = abs(atan(tan(theta1)*cos(psi1_2pi) )-atan(tan(theta2)*cos(psi2_2pi)));
                        Double_t exc_energy = sqrt(2*938.272*938.272+(11*931.5)*(11*931.5)+2*gamma_1*gamma_2*938.272*938.272*(1-beta_1*beta_2*cos(angle_p2_p1))+gamma_1*gamma_precise*(11*931.5)*938.272*(1-beta_precise*beta_1*cos(angle_p1_11b))+gamma_2*gamma_precise*(11*931.5)*938.272*(1-beta_precise*beta_2*cos(angle_p2_11b)) + 938.272*0.714549*1.42942*938.272*0.714549*1.42942)+complete_gamma_deposited-(12*931.5)-938.272-1.42942*938;
                        cout << "exec energy\t "<< exc_energy << endl;
                        h1_psi_in_11B->Fill(psi_in);
                        TVector2 k_0 =   TVector2(11414.5,0);
                        TVector2 k_1 =   TVector2(mom_p1*cos(theta1),mom_p1*sin(theta1));
                        TVector2 k_2 =   TVector2(mom_p2*cos(theta2),mom_p2*sin(theta2));
                        TVector2 k_a_1 = TVector2(mom_11b*cos(psi_in),mom_11b*sin(psi_in));
                        Double_t cos_distr = ((k_0 - k_1 - k_2)*k_a_1)/(sqrt((k_0 - k_1 - k_2)*(k_0 - k_1 - k_2))*sqrt(k_a_1*k_a_1));
						//CALIFA histos
                        if (v_E_gamma.size()){
                        h1_gamma_energyE_max_val_11B->Fill(*max_element(v_E_gamma.begin(),v_E_gamma.end()));
                        if(v_E_gamma.size() > 1){
                        sort(v_E_gamma.begin(),v_E_gamma.end());
                        h2_gamma_fst_vs_snd_11B->Fill(v_E_gamma[v_E_gamma.size()-2],v_E_gamma[v_E_gamma.size()-1]);
                        }
                        }
                        h1_cos_gamma_11b_p_i->Fill(cos_distr);
                        h1_theta_1_plus_theta_2_CALIFA_11b->Fill(((theta1+theta2)/PI)*180.);
                        h2_E1_vs_E2_CALIFA_11b->Fill(E1,E2);
                        h1_E1_plus_E2_CALIFA_11b->Fill(E1+E2);
                        h2_long_mom_p1p2_long_mom11b->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                        if (v_E_gamma.size() &&  *max_element(v_E_gamma.begin(),v_E_gamma.end()) > 1.){
                        h2_long_mom_p1p2_long_mom11b_large_gamma->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                        }
                        h1_transv_mom_difference_p1_p2_11b->Fill(abs(mom_p1*sin(theta1)-mom_p2*sin(theta2)));
                        h2_E1_plus_E2_CALIFA_vs_full_mom11b->Fill(E1+E2,mom_11b);
                        h1_binding_energy_11b->Fill(11*400-(E1+E2+sqrt(pow(mom_11b,2)+pow(11*938.272,2))-(11*938.272)));
                        h2_theta1_theta2_11b->Fill((theta1/PI)*180.,(theta2/PI)*180.);
                        h2_psi1_psi2_11b->Fill((psi1/PI)*180.,(psi2/PI)*180.);
                        h2_psi1_psi2_11b_2pi->Fill((psi1_2pi/PI)*180.,(psi2_2pi/PI)*180.);
                        h1_theta_1_CALIFA_11b->Fill((theta1/PI)*180.);
                        h1_theta_2_CALIFA_11b->Fill((theta2/PI)*180.);
                        h2_theta_sum_vs_phi_diff_11b->Fill(theta1_degr+theta2_degr,abs(psi1_2pi_degr-psi2_2pi_degr));
                        if (-0.9973*(mom_p1*cos(theta1)+mom_p2*cos(theta2))+10913 > (mom_11b*cos(psi_in)) && (mom_11b*cos(psi_in)) > 9820 && (mom_p1*cos(theta1)+mom_p2*cos(theta2)) > 750) {
                        h2_theta_1_vs_theta_2_CALIFA_11b_cut->Fill((theta1/PI)*180.,(theta2/PI)*180.);
                        h2_E1_vs_E2_CALIFA_11b_cut->Fill(E1,E2);
                        h1_psi_in_11B_cut_triangle->Fill(psi_in);
                        }
                        else {
						}
                        TVector3 p_11b (11*931.5*gamma_precise*beta_precise*sin(psi_in),0.,-gamma_given*0.714549*11*931.5*gamma_precise+11*931.5*gamma_precise*beta_precise*cos(psi_in)*gamma_given);
                        TVector3 p_p1(938.272*gamma_1*beta_1*sin(theta1)*cos(psi1_2pi),0.,-gamma_given*0.714549*938.272*gamma_1+gamma_given*938.272*gamma_1*beta_1*cos(theta1));
                        TVector3 p_p2(938.272*gamma_2*beta_2*sin(theta2)*cos(psi2_2pi),0.,-gamma_given*0.714549*938.272*gamma_2+gamma_given*938.272*gamma_2*beta_2*cos(theta2));
                        TVector3 p_tar(0.,0.,-gamma_given*0.714549*938.272);
                        TVector3 p_initial(p_p1+p_p2-p_tar);
                        cout << "this is momentum of initial projectile\t" << p_initial.Mag() << endl;
                        h1_p_missing_11B->Fill(p_initial.Mag());
                        Double_t my_sweep = 0.048869219056;
                        vector<double> cos_sweep;
                        for (Int_t sweep_1 = -1; sweep_1 < 2; sweep_1++){
                            for (Int_t sweep_2 = -1; sweep_2 < 2; sweep_2++){

                        TVector3 p_11b_sweep (11*931.5*gamma_precise*beta_precise*sin(psi_in),0.,-gamma_given*0.714549*11*931.5*gamma_precise+11*931.5*gamma_precise*beta_precise*cos(psi_in)*gamma_given);
                        TVector3 p_p1_sweep(938.272*gamma_1*beta_1*sin(theta1)*cos(psi1_2pi+sweep_1*my_sweep),0.,-gamma_given*0.714549*938.272*gamma_1+gamma_given*938.272*gamma_1*beta_1*cos(theta1));
                        TVector3 p_p2_sweep(938.272*gamma_2*beta_2*sin(theta2)*cos(psi2_2pi+sweep_2*my_sweep),0.,-gamma_given*0.714549*938.272*gamma_2+gamma_given*938.272*gamma_2*beta_2*cos(theta2));
                        TVector3 p_tar_sweep(0.,0.,-gamma_given*0.714549*938.272);
                        TVector3 p_initial_sweep(p_p1_sweep+p_p2_sweep-p_tar_sweep);
                        Double_t angle_cms_sweep = p_initial_sweep.Angle(p_11b_sweep);
                        cos_sweep.push_back(cos(angle_cms_sweep));
                            }


                        }
                        Double_t min_cos_angle_cms = *min_element(cos_sweep.begin(),cos_sweep.end());
                        cos_sweep.clear();
                        Double_t angle_cms = p_initial.Angle(p_11b);
                        cout << "cos(gamma) questioning....." << endl;
                        if (cos(angle_cms) < -0.6){
                        h2_long_mom_p1p2_long_mom11b_low->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                        }
                        else if (cos(angle_cms)>-0.6){
                        h2_long_mom_p1p2_long_mom11b_high->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                        }
                        if ((theta1_degr+theta2_degr) < 84 ){
                        h2_long_mom_p1p2_long_mom11b_small_opening->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                        }
                        else if ((theta1_degr+theta2_degr) > 84 ){
                        h2_long_mom_p1p2_long_mom11b_large_opening->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                        }
                        cout << "asking for angular restrictions...." << endl;
                        if ((theta1_degr+theta2_degr) < 84 && abs(180-abs(psi1_2pi_degr-psi2_2pi_degr)) < 30){
                            if(energy_theta_v.size() > 1){
                                cout << "AFTER asking!" << endl;
                                sort(energy_theta_v.begin(),energy_theta_v.end(),sortcol);
                                cout << "you sorted it riht..."<< endl;

                                Double_t mev_lab_2_1 = 2.124/((1/(sqrt(1-beta_precise*beta_precise)))*(1-beta_precise*cos(energy_theta_v[0][0])));
                                Double_t mev_lab_4_4 = 4.445/((1/(sqrt(1-beta_precise*beta_precise)))*(1-beta_precise*cos(energy_theta_v[0][0])));
                                Double_t mev_lab_5_0 = 5.020/((1/(sqrt(1-beta_precise*beta_precise)))*(1-beta_precise*cos(energy_theta_v[0][0])));
                                Double_t escape_peak = energy_theta_v[1][1]*(sqrt(1-beta_precise*beta_precise))/(1-beta_precise*cos(energy_theta_v[1][0]));
                                Double_t full_spec_energy = energy_theta_v[0][1];
                                cout << "I reach the lab stuff..." << endl;
                                if ((energy_theta_v[0][1] > (mev_lab_2_1-0.6) && energy_theta_v[0][1] < (mev_lab_2_1-0.4) && escape_peak < 0.6)||( energy_theta_v[0][1] > (mev_lab_2_1-1.1) && energy_theta_v[0][1] < (mev_lab_2_1-0.9) && escape_peak < 1.1)\
                                || (energy_theta_v[0][1] > (mev_lab_4_4-0.6) && energy_theta_v[0][1] < (mev_lab_4_4-0.4) && escape_peak < 0.6)|| (energy_theta_v[0][1] > (mev_lab_4_4-1.1) && energy_theta_v[0][1] < (mev_lab_4_4-0.9) && escape_peak < 1.1)\
                                || (energy_theta_v[0][1] > (mev_lab_5_0-0.6) && energy_theta_v[0][1] < (mev_lab_5_0-0.4) && escape_peak < 0.6)|| (energy_theta_v[0][1] > (mev_lab_5_0-1.1) && energy_theta_v[0][1] < (mev_lab_5_0-0.9) && escape_peak < 1.1)){
                                cout << "this is spec_energy:\t" << full_spec_energy <<endl;
                                for (Int_t spec = 1; spec < energy_theta_v.size(); spec++){
                                    full_spec_energy += (energy_theta_v[spec][1]*(sqrt(1-beta_precise*beta_precise))/(1-beta_precise*cos(energy_theta_v[spec][0])))*((1/(sqrt(1-beta_precise*beta_precise)))*(1-beta_precise*cos(energy_theta_v[0][0]))); //here the energy of 2nd gamma is first transformed back to laboratory then back to the particle system taking the angle of the gamma with the largest engergy

                                    cout << " in the spec loop" <<endl;
                                }
                                cout << "fillig spec histo...." << endl;
                                h1_gamma_spec_escape_angle_cut->Fill(full_spec_energy);
                                h2_gamma_energy_vs_theta_11B->Fill(full_spec_energy,((energy_theta_v[0][0])/PI)*180);
                                }
                                else {
                                    h1_gamma_spec_escape_angle_cut->Fill(full_spec_energy);
                                    h2_gamma_energy_vs_theta_11B->Fill(full_spec_energy,((energy_theta_v[0][0])/PI)*180);
                                }
                            }
                            else if (energy_theta_v.size() == 1){
                                h1_gamma_spec_escape_angle_cut->Fill(energy_theta_v[0][1]);
                                h2_gamma_energy_vs_theta_11B->Fill(energy_theta_v[0][1],((energy_theta_v[0][0])/PI)*180);
                            }
                            h2_long_mom_p1p2_long_mom11b_theta_phi_constr->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                            if (p_initial.Mag() < 250){
                                h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                                if (v_E_gamma.size()){
                                cout << "after questioning vE_gamma_size..." << endl;
                                cout << "this is max element:\t" << *max_element(v_E_gamma.begin(),v_E_gamma.end()) << endl;
                                Double_t max_gamma_test = *max_element(v_E_gamma.begin(),v_E_gamma.end());
                                h1_m_mom_250_gamma_spec->Fill(max_gamma_test);
                                cout << "after filling h1_m_mom_250_gamma_spec ... " << endl;
                                cout << "after h1_m_mom_250_gamma_spec->Fill .." << endl;
                                Int_t gamma_cluster = 0;
                                Int_t proton_cluster = 0;
                                cout << "right before 3rd hello...." << endl;
                                if (vE.size()){
                                    cout << "3rdHELLO" << endl;
                                    for (Int_t cl = 0;cl <vE.size();cl++){
                                        if (vE[cl] > 30000){
                                        cout << "TJ, hello from inside proton cluster" << endl;
                                        proton_cluster+= 1;
                                        }
                                        else{
                                        gamma_cluster+=1;
                                        cout << "TJ, hello from inside gamma cluster" << endl;
                                        }

                                    }
                                }
                                h1_cluster_gamma_11b->Fill(gamma_cluster);
                                h1_cluster_proton_11b->Fill(proton_cluster);
                                }
                                }

                            if (p_initial.Mag() < 200 && v_E_gamma.size()){
                                h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom200->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                                h1_m_mom_200_gamma_spec->Fill(*max_element(v_E_gamma.begin(),v_E_gamma.end()));
                                }
                            if (p_initial.Mag() < 150 &&  v_E_gamma.size()){
                                h2_long_mom_p1p2_long_mom11b_theta_phi_constr_m_mom150->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                                h1_m_mom_150_gamma_spec->Fill(*max_element(v_E_gamma.begin(),v_E_gamma.end()));
                                }
                            if (min_cos_angle_cms < -0.6){
                                h2_long_mom_p1p2_long_mom11b_low_optimized->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                                if (v_E_gamma.size()){
                                    h1_cos_low_gamma->Fill(*max_element(v_E_gamma.begin(),v_E_gamma.end()));
                                }
                            }
                            if (min_cos_angle_cms > 0.6){
                                h2_long_mom_p1p2_long_mom11b_high_optimized->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                            }

                            if (min_cos_angle_cms < -0.8){
                                h2_long_mom_p1p2_long_mom11b_low_optimized_08->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                            }
                            cout << "before triangle cut" <<endl;
                            if (-0.9973*(mom_p1*cos(theta1)+mom_p2*cos(theta2))+10913 > (mom_11b*cos(psi_in)) && (mom_11b*cos(psi_in)) > 9820 && (mom_p1*cos(theta1)+mom_p2*cos(theta2)) > 750) {
                            h1_mw3_fy_11b_cut_triangle->Fill(raw_mwpc3);
                            h1_mw2_fy_11b_cut_triangle->Fill(yMW2);
                            h1_psi_in_11b_cut_triangle->Fill(psi_in);

                            if (v_E_gamma.size()){
                            h1_gamma_energyE_max_val_11B_angle_cut->Fill(v_E_gamma[v_E_gamma.size()-1]);
                            }

                            }
                            else {
                            h1_mw3_fy_11b->Fill(raw_mwpc3);
                            h1_mw2_fy_11b->Fill(yMW2);
                            h1_psi_in_11b->Fill(psi_in);
                            if (v_E_gamma.size()){
                            h1_gamma_energyE_max_val_11B_angle_cut->Fill(v_E_gamma[v_E_gamma.size()-1]);
                            }
                            }
                          if (v_E_gamma.size()){
                          v_E_gamma.clear();
                          }
                        }
                        cout << "before border constraints...." << endl;
                        cout << "psi1 in degrees:\t" << psi1_2pi_degr<< endl;
                        cout << "psi2 in degrees:\t" << psi2_2pi_degr << endl;
                        cout << "theta1 in degrees:\t" << theta1_degr << endl;
                        cout << "theta2 in degrees:\t" << theta2_degr << endl;
                         if( (psi1_2pi_degr > 84 && psi1_2pi_degr < 96) || (psi1_2pi_degr > 264 && psi1_2pi_degr < 276) || (psi2_2pi_degr > 84 && psi2_2pi_degr < 96) || (psi2_2pi_degr > 264 && psi2_2pi_degr)){
                            cout << "hello in first condition" << endl;
                            continue;
                            }
                        if((psi1_2pi_degr> 42 &&  psi1_2pi_degr < 160.5 && theta1_degr < 46) || (psi1_2pi_degr> 199.5 &&  psi1_2pi_degr < 318 && theta1_degr < 46)){
                            cout << "hello in second contidtion" << endl;
                            continue;
                            }
                        if((psi2_2pi_degr> 42 &&  psi2_2pi_degr < 160.5 && theta2_degr < 46) || (psi2_2pi_degr> 199.5 &&  psi2_2pi_degr < 318 && theta2_degr < 46)){
                            cout << "hello in third condition " << endl;
                            continue;
                            }
                        if (theta1_degr < 25 || theta2_degr < 25){
                            cout << "hello in forth condition" << endl;
                            continue;
                            }
                        cout << "After border constraints ...." << endl;
                        h2_long_mom_p1p2_long_mom11b_no_border->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                        h1_cos_gamma_11b_p_i_cms->Fill(cos(angle_cms));
                        h1_cos_gamma_11b_p_i_cms_optimized->Fill(min_cos_angle_cms);
                        cout << "after filling some histograms..." << endl;
                        if (cos(angle_cms) < -0.6){
                        h2_long_mom_p1p2_long_mom11b_no_border_low->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                        }
                        else if (cos(angle_cms)>0.6){
                        h2_long_mom_p1p2_long_mom11b_no_border_high->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_11b*cos(psi_in));
                        }
                        cout << "After asking cos(gamma) <> -0.6...."<< endl;
                        cout << "last calculation....." << endl;
                        for (Int_t k =0;k< vE.size();k++)
                        {
                            cout << "inside the last for loop in 11B" << endl;
                            if(v_origin_E[k] < 10000)
                            {
                            cout << "inside the last if questioning..." << endl;
                            cout << "this is v_origin_E[k]" << v_origin_E[k] << endl;
                            cout << " I shouldn't see this line if v_origin_E[k] is empty..." << endl;
                            Double_t theta_gamma_noborder = vtheta[k];
                            Double_t energy_gamma = (v_origin_E[k]/1000)*(1/(sqrt(1-beta*beta)))*(1-beta*cos(theta_gamma_noborder));
                            h1_gamma_energyE_11B_no_border->Fill(energy_gamma);

                            }



                        }
				}
				}
                        if (a_q < aq_cut_z_5 && a_q > (aq_cut_z_5-0.15) && entries_califa >= 2){  //10B
                        psi_in = mean_psi_out_10b-slope_10b*xMW3 - angle_offs_10b;
                        m = (cos(psi_out)-cos(psi_in))/(sin(psi_out)+sin(psi_in));
                        rho =  (Leff/(2*sin((psi_in+psi_out)/2)))*sqrt(pow(((m+tan(alpha_G))/(1-m*tan(alpha_G))),2)+1);
                        w = 2*abs(asin(BD/(2*rho)));
                        Double_t pathlength_precise = abs((zB-zT)/(cos(psi_in)))+rho*w+abs((zToFW-z_D)/cos(psi_out));
                        Double_t beta_precise = ((pathlength_precise/time_target_tof)*pow(10,6))/light_c;
                        Double_t a_q_precise = (pow(10,-3)*((((mag_field*rho)/(beta_precise))*1.602176634*pow(10,-19))/(1.66053906660*pow(10,-27)*light_c)))/gamma_given;
                        Double_t mom_10b = a_q_precise*5*(beta/(sqrt(1-(beta*beta))))*938.272;

						//CALIFA handling
                        Double_t psi_CALIFA;
                        Double_t theta_CALIFA;
                        Double_t time_CALIFA;
                        Double_t E_Califa;
                        vector<double> vE1;
                        vector<double> vtheta1CALIFA;
                        vector<double> vpsi1CALIFA;
                        vector<double> vE2;
                        vector<double> vtheta2CALIFA;
                        vector<double> vpsi2CALIFA;
                        vector<double> vpsi;
                        vector<double> vtheta;
                        vector<double> vE;
                        vector<double> vtime;
                        vector<double> v_origin_E;
                        Double_t ind_high_energy;
                        Double_t ind_second_high_energy;
                        Double_t theta_gamma;
                        Double_t psi1;
                        Double_t psi2;
                        for (Int_t j = 0.;j<entries_califa;j++){
                            califahitdata[j] = (R3BCalifaHitData*)CalifaHitData->At(j);
                            psi_CALIFA = califahitdata[j]->GetPhi();
                            theta_CALIFA = califahitdata[j]->GetTheta();
                            E_Califa = califahitdata[j]->GetEnergy();
                            time_CALIFA = califahitdata[j]->GetTime();
                            vpsi.push_back(psi_CALIFA);
                            vtheta.push_back(theta_CALIFA);
                            vE.push_back(E_Califa);
                            v_origin_E.push_back(E_Califa);
                            vtime.push_back(time_CALIFA);
                            }
                        sort(vE.begin(),vE.end());
                        if (vE[vE.size()-1] >30000 && vE[vE.size()-2] > 30000)
                        {
                        for (Int_t k =0;k<vE.size();k++)
                        {
                            if (v_origin_E[k] == vE[vE.size()-1]) {
                            ind_high_energy = k;
                            }
                            else if (v_origin_E[k] == vE[vE.size()-2]){
                            ind_second_high_energy = k;
                            }
                            else if(v_origin_E[k] < 10000)
                            {
                            theta_gamma = vtheta[k];
                            Double_t energy_gamma = (v_origin_E[k]/1000)*(1/(sqrt(1-beta*beta)))*(1-beta*cos(theta_gamma));
                            h1_gamma_energyE_10B->Fill(energy_gamma);

                            }



                        }
                        Double_t E1 = v_origin_E[ind_high_energy]/1000.;
                        Double_t E2 = v_origin_E[ind_second_high_energy]/1000.;
                        Double_t theta1 = vtheta[ind_high_energy];
                        Double_t theta2 = vtheta[ind_second_high_energy];

                        if (vpsi[ind_high_energy] < 0.){
                        psi1 = 2*PI + vpsi[ind_high_energy];
                        }
                        if (vpsi[ind_high_energy] > 0.){
                        psi1 = vpsi[ind_high_energy];
                        }
                        if (vpsi[ind_second_high_energy] < 0.){
                        psi2 = 2*PI + vpsi[ind_second_high_energy];
                        }
                        if (vpsi[ind_second_high_energy] > 0.){
                        psi2 = vpsi[ind_second_high_energy];
                        }

                        Double_t mom_p1 = sqrt(pow((E1+938.272),2)-pow(938.272,2));
                        Double_t mom_p2 = sqrt(pow((E2+938.272),2)-pow(938.272,2));
                        h1_theta_1_plus_theta_2_CALIFA_10b->Fill(((theta1+theta2)/PI)*180.);
                        h2_E1_vs_E2_CALIFA_10b->Fill(E1,E2);
                        h1_E1_plus_E2_CALIFA_10b->Fill(E1+E2);
                        h2_long_mom_p1p2_long_mom10b->Fill(mom_p1*cos(theta1)+mom_p2*cos(theta2),mom_10b*cos(psi_in));
                        h1_transv_mom_difference_p1_p2_10b->Fill(abs(mom_p1*sin(theta1)-mom_p2*sin(theta2)));
                        }
                        h2_z_vs_a_q->Fill(a_q_precise,charge_val);
                        h2_tof_vs_aq_fix_g_10b->Fill(a_q_precise,time_target_tof);
                        }
                        }
                        if (charge_val< charge_cut+0.8 && charge_val> charge_cut) //Z = 6
                        {
                        if (a_q > aq_cut_z_6 && a_q < (aq_cut_z_6+2*width_12c)){  //12C
                        psi_in = mean_psi_out_12c-slope_12c*xMW3 - angle_offs_12c;
                        m = (cos(psi_out)-cos(psi_in))/(sin(psi_out)+sin(psi_in));
                        rho =  (Leff/(2*sin((psi_in+psi_out)/2)))*sqrt(pow(((m+tan(alpha_G))/(1-m*tan(alpha_G))),2)+1);
                        w = 2*abs(asin(BD/(2*rho)));
                        Double_t pathlength_precise = abs((zB-zT)/(cos(psi_in)))+rho*w+abs((zToFW-z_D)/cos(psi_out));
                        Double_t beta_precise = ((pathlength_precise/time_target_tof)*pow(10,6))/light_c;
                        Double_t a_q_precise = (pow(10,-3)*((((mag_field*rho)/(beta_precise))*1.602176634*pow(10,-19))/(1.66053906660*pow(10,-27)*light_c)))/gamma_given;
                        if (time_target_tof < 35.13){
                        h1_tof_detnr_strange_12c->Fill(i+1);
                        h1_mw3_fy_strange_12c->Fill(raw_mwpc3);
                        h2_tof_path_strange12c->Fill(pathlength_precise,time_target_tof);
                        h2_radius_path_strange12c->Fill(rho,time_target_tof);
                        h1_tofwall_posy_fast[i]->Fill(tof_diff_up_down);
                        }
                        if (time_target_tof> 35.35 && time_target_tof<35.61){
                        h1_tof_detnr_12c->Fill(i+1);
                        h1_mw3_fy_12c->Fill(raw_mwpc3);
                        h2_tof_path_12c->Fill(pathlength_precise,time_target_tof);
                        h2_radius_path_12c->Fill(rho,time_target_tof);
                        h1_tofwall_posy_medium[i]->Fill(tof_diff_up_down);
                        }
                        if (a_q_precise > 1.997){
                        h1_tof_detnr_strange_12c_large_aq->Fill(i+1);
                        h1_mw3_fy_strange_12c_large_aq->Fill(raw_mwpc3);
                        h2_tof_path_strange_12c_large_aq->Fill(pathlength_precise,time_target_tof);
                        h2_radius_path_strange_12c_large_aq->Fill(rho,time_target_tof);
                        h1_tofwall_posy_slow[i]->Fill(tof_diff_up_down);
                        }
                        if (a_q < aq_cut_z_6-0.06 && a_q > (aq_cut_z_6-1.65*width_11c) && entries_califa >= 2){  //11C before it was a_q < aq_cut_z_6 && a_q > (aq_cut_z_6-2*width_11c), too bad for 11c...
                        psi_in = mean_psi_out_11c-slope_11c*xMW3 - angle_offs_11c;
                        m = (cos(psi_out)-cos(psi_in))/(sin(psi_out)+sin(psi_in));
                        rho =  (Leff/(2*sin((psi_in+psi_out)/2)))*sqrt(pow(((m+tan(alpha_G))/(1-m*tan(alpha_G))),2)+1);
                        w = 2*abs(asin(BD/(2*rho)));
                        Double_t pathlength_precise = abs((zB-zT)/(cos(psi_in)))+rho*w+abs((zToFW-z_D)/cos(psi_out));
                        Double_t beta_precise = ((pathlength_precise/time_target_tof)*pow(10,6))/light_c;
                        Double_t a_q_precise = (pow(10,-3)*((((mag_field*rho)/(beta_precise))*1.602176634*pow(10,-19))/(1.66053906660*pow(10,-27)*light_c)))/gamma_given;
                        h2_z_vs_a_q->Fill(a_q_precise,charge_val);
                        h2_tof_vs_aq_fix_g_11c->Fill(a_q_precise,time_target_tof);
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

//check that I do not have to delete some on heap stuff...





}
