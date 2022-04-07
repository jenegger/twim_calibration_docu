//beginning after isotope has been selected:
//TLorentzVector v_11B_lab(p_x_lab,p_y_lab,p_z_lab,E_lab);
//TLorentzVector v_11B_cms(p_x_cms,p_y_cms,p_z_cms,E_cms);
//fill entries of CALIFA in vectors, no matter if they are protons or gammas..
//variales needed:
vector<double> v_Energy;
vector<double> v_Theta;
vector<double> v_psi;
vector<double> v_CALIFA_evts_lab_energy; //here I insert vectors v = (E_lab,theta_lab,phi_lab)
vector<double> v_CALIFA_evts_cms_energy; //here I inserv vectors v = (E_cms,theta_lab,phi_lab)
vector<TLorentzVector<double> > v_lorentz_CALIFA_evts_lab; //here I insert Lorentz vector v = (p_x_lab,p_y_lab,p_z_lab,E_lab)
vector<TLorentzVector<double> > v_lorentz_CALIFA_evts_cms; //here I insert Lorentz vector v = (p_x_cms,p_y_cms,p_z_cms,E_cms)


for (Int_t j = 0.;j<entries_califa;j++){
	vector<double> v_temp;
	vector<double> v_temp_cms;
	califahitdata[j] = (R3BCalifaHitData*)CalifaHitData->At(j);
	v_temp[0] = (califahitdata[j]->GetEnergy())/1000.;
	v_temp[1] = califahitdata[j]->GetTheta();
	v_temp[2] = califahitdata[j]->GetPhi();
	v_temp_cms[0] = v_temp[0]*(1/(sqrt(1-beta_precise*beta_precise)))*(1-beta_precise*cos(vtheta[k])); //Doppler correction
	v_temp_cms[1] = v_temp[1];
    v_temp_cms[2] = v_temp[2];
	v_CALIFA_evts_lab_energy.push_back(v_temp);
	v_CALIFA_evts_cms_energy.push_back(v_temp_cms);
	v_temp.clear();
	v_temp_cms.clear();
	
}


for (int CALIFA_events = 0; CALIFA_events < entries_califa; CALIFA_events++){

if(v_CALIFA_evts_lab_energy[CALIFA_events][0] < 30) { //gammas
//create lorentz vector in the laboratory system
TLorentzVector temp_lorentz_lab;
temp_lorentz_lab[0] = (v_CALIFA_evts_lab_energy[CALIFA_events][0])*sin(v_CALIFA_evts_lab_energy[CALIFA_events][1])*cos(v_CALIFA_evts_lab_energy[CALIFA_events][2]); //coordinate transform
temp_lorentz_lab[1] = (v_CALIFA_evts_lab_energy[CALIFA_events][0])*sin(v_CALIFA_evts_lab_energy[CALIFA_events][1])*sin(v_CALIFA_evts_lab_energy[CALIFA_events][2]); //coordinate transform
temp_lorentz_lab[2] = (v_CALIFA_evts_lab_energy[CALIFA_events][0])*cos(v_CALIFA_evts_lab_energy[CALIFA_events][1]); //coordinate transform
temp_lorentz_lab[3] = v_CALIFA_evts_lab_energy[CALIFA_events][0];
v_lorentz_CALIFA_evts_lab.push_back(temp_lorentz_lab);
temp_lorentz_lab.clear();

}
else {//protons

TLorentzVector temp_lorentz_lab;
temp_lorentz_lab[0] = sqrt(pow(((v_CALIFA_evts_lab_energy[CALIFA_events][0])+938.272),2)-pow(938.272,2))*sin(v_CALIFA_evts_lab_energy[CALIFA_events][1])*cos(v_CALIFA_evts_lab_energy[CALIFA_events][2]); //coordinate transform
temp_lorentz_lab[1] = sqrt(pow(((v_CALIFA_evts_lab_energy[CALIFA_events][0])+938.272),2)-pow(938.272,2))*sin(v_CALIFA_evts_lab_energy[CALIFA_events][1])*sin(v_CALIFA_evts_lab_energy[CALIFA_events][2]); //coordinate transform
temp_lorentz_lab[2] = sqrt(pow(((v_CALIFA_evts_lab_energy[CALIFA_events][0])+938.272),2)-pow(938.272,2))*cos(v_CALIFA_evts_lab_energy[CALIFA_events][1]); //coordinate transform
temp_lorentz_lab[3] = sqrt(pow((sqrt(pow((v_CALIFA_evts_lab_energy[CALIFA_events][0]+938.272),2)-pow(938.272,2))),2)+pow(938.272,2));
v_lorentz_CALIFA_evts_lab.push_back(temp_lorentz_lab);
temp_lorentz_lab.clear();

}
}
//now sort the vector of lorentz vectors and v_CALIFA_evts_cms_Energy;


