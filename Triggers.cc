#include "Triggers.h"

/******************************************************************************//**
* Default constructor
**********************************************************************************/
Triggers::Triggers(){

	InitializeCounter();
	InitializeFlags();

	  //Jets angular parameter
  R = 0.4;
	Eiso_R = 0.2;
	Miso_R = 0.3;
	Giso_R = 0.4;

  //Lepton isolation cone radius
  LEP_R = 0.5;   //*** Is this correct? ****

  //Calorimeter detector geometry constants
  BARREL_ETA = 1.37, ENDCAP_MIN_ETA = 1.52, ENDCAP_MAX_ETA = 2.37,
  HEC_MIN_ETA = 2.37, HEC_MAX_ETA = 3.1, FCAL_MIN_ETA = 3.1, FCAL_MAX_ETA = 4.9;

  //Tracker geometry constants
	TRACKER_ETA = 2.49;

	//FastJet Clustering Algorithm for General-KT: -1 = anti-kT; 0 = C/A; 1 = kT.
  POWER   = -1;     

  //High level MET and Lepton thresholds 
	MET_HL = 110, ELECTRON_HL = 26, MUON_HL = 26, GAMMA_HL = 140,
					DI_ELECTRON_HL = 17, DI_MUON_HL = 14;

  //High level multijet and HT thresholds 
  HT_HL = 850, THREE_JET_HL = 175, FOUR_JET_HL = 100, FIVE_JET_HL = 70, 
	SIX_JET_HL = 60, ONE_JET_HL_FAT = 420, ONE_JET_HL_THIN = 380;

  //Low level MET and lepton thresholds
  MET_LL = 50, ELECTRON_LL = 22, MUON_LL = 20, GAMMA_LL = 20,
          DI_ELECTRON_LL = 17, DI_MUON_LL = 14;

  //Low level multijet and HT thresholds #3x50 and 4x15
  HT_LL = 0, THREE_JET_LL = 0, FOUR_JET_LL= 50, FIVESIX_JET_LL = 15, ONE_JET_LL = 100;

	//Atlas leption isolation cut
	ISO_CUT_Tight = 0.06;
	ISO_CUT_Loose = 0.15;

}
/******************************************************************************//**
* Destructor
**********************************************************************************/
Triggers::~Triggers(){

}
/******************************************************************************//**
* Initialise trigger counter
**********************************************************************************/
void Triggers::InitializeCounter(){

	// High level counters
	N_JET_ONE							= 0, N_JET_THREE					  = 0,
	N_JET_FOUR 					  = 0, N_JET_FIVE			  			= 0,
	N_JET_SIX							= 0, N_ELECTRON 						= 0,
	N_MUON 							  = 0, N_GAMMA 								= 0,
	N_MET 								= 0, N_HT 									= 0,
	N_JET_ONE_FAT					= 0;
	// Lower level counters
	N_JET_ONE_L						= 0, N_JET_THREE_L				  = 0,
	N_JET_FOUR_L 				  = 0, N_JET_FIVE_L			  		= 0,
	N_JET_SIX_L						= 0, N_ELECTRON_L 					= 0,
	N_MUON_L 						  = 0, N_GAMMA_L 							= 0,
	N_MET_L 							= 0, N_HT_L 								= 0;

}
/******************************************************************************//**
* Initialise trigger call flags
**********************************************************************************/
void Triggers::InitializeFlags(){

	HT_HL_CALL				= false , 	HT_LL_CALL				= false,
	MULTIJET_HL_CALL  = false ,   MULTIJET_LL_CALL  = false,
	EW_HL_CALL				= false , 	EW_LL_CALL				= false,
	MET_HL_CALL 			= false ,	  MET_LL_CALL 			= false;

}
/******************************************************************************//**
* Initialise jet matching in pythia
**********************************************************************************/
void Triggers::InitializeMatching(Pythia *pythia)
{

	UserHooks* matching            = NULL;

  // For jet matching, initialise the respective user hooks code.
  CombineMatchingInput combined;

  matching = combined.getHook( *pythia);
   
  pythia -> setUserHooksPtr(matching);

}
/******************************************************************************//**
* Initialise Pythia8 from a delcared modal/config script
**********************************************************************************/
void Triggers::InitializePythia(Pythia *pythia)
{

	pythia ->readFile("ModA.dat");
	int nEvent = pythia -> mode("Main:numberOfEvents");
	
	pythia -> init();
}
/******************************************************************************//**
* Pythia event for Electro-Weak Initial State Radiation 
**********************************************************************************/
void Triggers::eventElectroWeak(Pythia *pythia, fastjet::JetDefinition jetDef)
{

	electron_hl.resize(0), muon_hl.resize(0), gamma_hl.resize(0),
	electron_ll.resize(0), muon_ll.resize(0), gamma_ll.resize(0),
	fjInputs.resize(0), missingETvec.reset();

	for (int i = 0; i < pythia -> event.size(); ++i) {

		// Selection process for particles to be clustered
		if( pythia -> event[i].isFinal() ){

			if(pythia -> event[i].isVisible() && ( abs( pythia -> event[i].eta() ) 
				< FCAL_MAX_ETA ) ){
				missingETvec += pythia -> event[i].p();
			}

			if( !pythia -> event[i].isVisible() ) continue;

			if(abs( pythia -> event[i].eta() ) > TRACKER_ETA) continue;

			fastjet::PseudoJet particleTemp = pythia -> event[i];

			fjInputs.push_back( particleTemp);

	    if(abs(pythia -> event[i].id()) == 11 && pythia -> event[i].pT() > ELECTRON_HL ){

	    	electron_hl.push_back(i); // fill with electron candidate
	    }
	    if(abs(pythia -> event[i].id()) == 13 && pythia -> event[i].pT() > MUON_HL ){

	      muon_hl.push_back(i); // fill with muon candidate
	    }
	    if(abs(pythia -> event[i].id()) == 22 && pythia -> event[i].pT() > GAMMA_HL ){

	      gamma_hl.push_back(i); // fill with photon candidate
	    }
	    if(abs(pythia -> event[i].id()) == 11 && pythia -> event[i].pT() > ELECTRON_LL ){

	    	electron_ll.push_back(i); // fill with electron candidate
	    }
	    if(abs(pythia -> event[i].id()) == 13 && pythia -> event[i].pT() > MUON_LL ){

	      muon_ll.push_back(i); // fill with muon candidate
	    }
	    if(abs(pythia -> event[i].id()) == 22 && pythia -> event[i].pT() > GAMMA_LL ){

	      gamma_ll.push_back(i); // fill with photon candidate
	    }
		}
	}// end of particle 
	JetClustering( jetDef);

}
/******************************************************************************//**
* Pythia event for Tree level process (No ISR) 
**********************************************************************************/
void Triggers::eventTree(Pythia *pythia)
{

	missingETvec.reset();

	for (int i = 0; i < pythia -> event.size(); ++i) {
		
		// Selection process for particles to be clustered
		if( pythia -> event[i].isFinal() ){

			if(pythia -> event[i].isVisible() && abs( pythia -> event[i].eta() ) 
				< FCAL_MAX_ETA   ){

				missingETvec += pythia -> event[i].p();
			}
		}
	}// end of particle 
}
/******************************************************************************//**
* Pythia event for Tree level process (No ISR) 
**********************************************************************************/
void Triggers::eventTree(Pythia *pythia, fastjet::JetDefinition jetDef)
{

	fjInputs.resize(0);
	missingETvec.reset();

	for (int i = 0; i < pythia -> event.size(); ++i) {
		
		// Selection process for particles to be clustered
		if( pythia -> event[i].isFinal() ){

			if(pythia -> event[i].isVisible() && abs( pythia -> event[i].eta() ) 
				< FCAL_MAX_ETA   ){

				missingETvec += pythia -> event[i].p();
			}

			if( !pythia -> event[i].isVisible() ) continue;

			if(abs( pythia -> event[i].eta() ) > TRACKER_ETA) continue;

			fastjet::PseudoJet particleTemp = pythia -> event[i];

			fjInputs.push_back( particleTemp);

		}
	}// end of particle 
	JetClustering(jetDef);
}
/******************************************************************************//**
* Jet clustering algorithm
**********************************************************************************/
void Triggers::JetClustering( fastjet::JetDefinition jetDef){

	vector <fastjet::PseudoJet> inclusiveJets;

	jets.resize(0);

	fastjet::ClusterSequence clustSeq(fjInputs, jetDef);

	inclusiveJets = clustSeq.inclusive_jets(0);
	jets  = sorted_by_pt(inclusiveJets);

}
/******************************************************************************//**
* MET and MET-Trigger calculation using the TopoCluster style algorithm 
**********************************************************************************/
double Triggers::TopoMET(){

	return missingETvec.pT();
}
/******************************************************************************//**
* MET and MET-Trigger calculation using the Jet style algorithm 
**********************************************************************************/
double Triggers::JetMET(){

	for (unsigned i = 0; i < jets.size(); i++) {

		if(jets[i].perp() > 30){

			// transfering fastjet 4-vector to pythia 4-vector
			missingETvecJet[0] += jets[i].px();
			missingETvecJet[1] += jets[i].py();
			missingETvecJet[2] += jets[i].pz();
			missingETvecJet[3] += jets[i].e();

		} 
	}
	return missingETvecJet.pT();
}
/******************************************************************************//**
* Call to privately stored clustered jets (must initalize jet clustering first)
**********************************************************************************/
vector <fastjet::PseudoJet> Triggers::SortedJets(){ 

	return jets;
}
/******************************************************************************//**
* Call to privately stored leptons (must be found from EW event first)
**********************************************************************************/
vector <double> Triggers::electronsHL(){ 

	return electron_hl;
}
vector <double> Triggers::muonsHL(){ 

	return muon_hl;
}
vector <double> Triggers::gammasHL(){ 

	return gamma_hl;
}
vector <double> Triggers::electronsLL(){ 

	return electron_ll;
}
vector <double> Triggers::muonsLL(){ 

	return muon_ll;
}
vector <double> Triggers::gammasLL(){ 

	return gamma_ll;
}
/******************************************************************************//**
* HT Calulation ( Scalar sum of jet transverse momentum in event)
**********************************************************************************/
double Triggers::HT(){ 

	double Ht = 0;

	for (unsigned i = 0; i < jets.size(); i++) {

		if(jets[i].perp() > 50)  Ht += jets[i].perp();
	}
	return Ht;
}
/******************************************************************************//**
* Lepton isolation algorithm: returns vector of candidate I-variable
**********************************************************************************/
bool Triggers::LeptonISO(Pythia *pythia, vector <double> candidate, int type){

	bool ISO_CUT_FLAG = false;

	//run over lepton candidates
	for (unsigned i = 0; i < candidate.size(); i++){		
	  
	  double Pt_sum = 0;
	  double Et_sum = 0;

		for (unsigned j = 0; j < pythia -> event.size(); ++j){
				
			if( !pythia -> event[j].isFinal() ) continue;

			if( !pythia -> event[j].isVisible() ) continue;

			if(abs( pythia -> event[j].eta() ) > TRACKER_ETA) continue;

			// pseudorapidity difference 
			double dEta = pythia -> event[j].y() - pythia -> event[ candidate[i] ].y();

			// phi coordinate difference
			double dPhi = abs( pythia -> event[j].phi() - pythia -> event[ candidate[i] ].phi() );

			// account for phi coordinate system missmatch
			if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;

			// distance/metric in phi-eta system
			double del_R = sqrt( pow2(dEta) + pow2(dPhi) );

			// sum pt of constituents within cone radius for either electron, muon, gamma
			if(type == 0 && del_R <= Eiso_R ) Pt_sum += pythia -> event[j].pT();
			if(type == 1 && del_R <= Miso_R ) Pt_sum += pythia -> event[j].pT();
			if(type == 2 && del_R <= Giso_R ) Et_sum += pythia -> event[j].eT();
		}
		Et_sum -= pythia -> event[ candidate[i] ].eT();
		double metric = (Pt_sum - pythia -> event[ candidate[i] ].pT() )/pythia -> event[ candidate[i] ].pT();
		double ISO_CUT_Gamma = 0.022 * pythia -> event[ candidate[i] ].eT() + 2.45;

		// Isolation criteria for either lepton/gamma
		if (type == 0 && metric <= ISO_CUT_Loose ) ISO_CUT_FLAG = true;				
		if (type == 1 && metric <= ISO_CUT_Loose ) ISO_CUT_FLAG = true;				
		if (type == 2 && Et_sum <= ISO_CUT_Gamma ) ISO_CUT_FLAG = true;						

	}// end of lepton loop
	return ISO_CUT_FLAG;

}
/******************************************************************************//**
* Efficency callculations
**********************************************************************************/
double Triggers::Efficiency(double accepted, double total){

	return 1.0*(accepted/total);

}
/******************************************************************************//**
* MultiJet Triggers, (1,3,4,5,6) jets 
**********************************************************************************/
void Triggers::MultiJetTriggerHL(){
	
	int count [6] = { };


	for (unsigned i = 0; i < jets.size(); i++) {

		if(jets[i].perp() > ONE_JET_HL_THIN ) 	count[0]++;
		if(jets[i].perp() > THREE_JET_HL    ) 	count[1]++;
		if(jets[i].perp() > FOUR_JET_HL     ) 	count[2]++;
		if(jets[i].perp() > FIVE_JET_HL     ) 	count[3]++;
		if(jets[i].perp() > SIX_JET_HL      ) 	count[4]++;
		if(jets[i].perp() > ONE_JET_HL_FAT  ) 	count[5]++;


	}// end of jet loop

	if( count[0] > 0) 	N_JET_ONE++;
	if( count[5] > 0) 	N_JET_ONE_FAT++;
	if( count[1] > 2) 	N_JET_THREE++;
	if( count[2] > 3) 	N_JET_FOUR++;
	if( count[3] > 4) 	N_JET_FIVE++;
	if( count[4] > 5)	 	N_JET_SIX++;


	MULTIJET_HL_CALL = true;
	
}
void Triggers::MultiJetTriggerLL(){
	
	int count [4] = { };

	for (unsigned i = 0; i < jets.size(); i++) {

		if(jets[i].perp() > ONE_JET_LL      ) 	count[0]++;
		if(jets[i].perp() > THREE_JET_LL    ) 	count[1]++;
		if(jets[i].perp() > FOUR_JET_LL     ) 	count[2]++;
		if(jets[i].perp() > FIVESIX_JET_LL  ) 	count[3]++;

	}// end of jet loop

	if( count[0] > 0) 	N_JET_ONE_L++;
	if( count[2] > 2) 	N_JET_FOUR_L++;
	if( count[3] > 3) 	N_JET_FIVE_L++;
	if( count[3] > 3)	 	N_JET_SIX_L++;


	MULTIJET_LL_CALL = true;
	
}
/******************************************************************************//**
* Missing Transverse Energy (momentum) Trigger 
**********************************************************************************/
void Triggers::METTriggerHL(){

	if( TopoMET() > MET_HL) N_MET++;

	MET_HL_CALL = true;

}
void Triggers::METTriggerLL(){

	if( TopoMET() > MET_LL) N_MET_L++;

	MET_LL_CALL = true;

}
void Triggers::METTriggerHL(int scheme){

	if( scheme == 0){
		if( TopoMET() > MET_HL) N_MET++;
 	}
 	if( scheme == 1){
 		if( JetMET() > MET_HL) N_MET++;
 	}

 	MET_HL_CALL = true;

}
void Triggers::METTriggerLL(int scheme){

	if( scheme == 0){
		if( TopoMET() > MET_LL) N_MET_L++;
 	}
 	if( scheme == 1){
 		if( JetMET() > MET_LL) N_MET_L++;
 	}

 	MET_LL_CALL = true;

}
/******************************************************************************//**
* ElectroWeak Triggers
**********************************************************************************/
void Triggers::EWTriggerHL(Pythia *pythia){

	if( LeptonISO( pythia, electron_hl, 0 ) ) N_ELECTRON++;
	if( LeptonISO( pythia, muon_hl, 1     ) ) N_MUON++;
	if(  gamma_hl.size() > 0 ) N_GAMMA++;

	EW_HL_CALL = true;

}
void Triggers::EWTriggerLL(Pythia *pythia){

	if( LeptonISO( pythia, electron_ll, 0 ) ) N_ELECTRON_L++;
	if( muon_ll.size() > 0  ) N_MUON_L++;
	if( LeptonISO( pythia, gamma_ll, 2 ) ) N_GAMMA_L++;

	EW_LL_CALL = true;

}
/******************************************************************************//**
* HT Trigger using native HT calculation
**********************************************************************************/
void Triggers::HtTriggerHL(){

	if( HT() > HT_HL ) N_HT++;

	HT_HL_CALL = true;

}
void Triggers::HtTriggerLL(){

	if( HT() > HT_LL ) N_HT_L++;

	HT_LL_CALL = true;

}
/******************************************************************************//**
* High level efficiency calls for various triggers
**********************************************************************************/
void Triggers::HLTriggerEfficiencies(char *name, double mass, double total ){

	ofstream file;

	file.open (name, ios::app); 

	file << mass 																<< " ";

	if(HT_HL_CALL){

		file << Efficiency( N_HT , 			   total) << " ";
	}
	if(MULTIJET_HL_CALL){

		file << Efficiency( N_JET_ONE ,    total) << " "
				 << Efficiency( N_JET_ONE_FAT, total) << " "
				 << Efficiency( N_JET_THREE ,  total) << " "
				 << Efficiency( N_JET_FOUR ,   total) << " "
				 << Efficiency( N_JET_FIVE ,   total) << " "
				 << Efficiency( N_JET_SIX ,    total) << " ";
	}
	if(EW_HL_CALL){
		file << Efficiency( N_ELECTRON ,   total) << " "
				 << Efficiency( N_MUON ,       total) << " "
				 << Efficiency( N_GAMMA ,      total) << " ";

	}
	if(MET_HL_CALL){
		file << Efficiency( N_MET ,        total) << " ";
	}

	file << endl;

	file.close();
}
/******************************************************************************//**
* Lower level efficiency calls for various triggers
**********************************************************************************/
void Triggers::LLTriggerEfficiencies(char *name, double mass, double total ){

	ofstream file;

	file.open (name, ios::app); 

	file << mass 																<< " ";


	if(HT_LL_CALL){

		file << Efficiency( N_HT_L , 				total) << " ";
	}
	if(MULTIJET_LL_CALL){

		file << Efficiency( N_JET_ONE_L ,   total) << " "
				 << Efficiency( N_JET_FOUR_L ,  total) << " "
				 << Efficiency( N_JET_FIVE_L ,  total) << " "
				 << Efficiency( N_JET_SIX_L ,   total) << " ";
	}
	if(EW_LL_CALL){
		file << Efficiency( N_ELECTRON_L ,  total) << " "
				 << Efficiency( N_MUON_L ,      total) << " "
				 << Efficiency( N_GAMMA_L ,     total) << " ";

	}
	if(MET_LL_CALL){
		file << Efficiency( N_MET_L ,       total) << " ";
	}

	file << endl;

	file.close();
}