#include "Triggers.h"

/******************************************************************************//**
* Default constructor
**********************************************************************************/
Triggers::Triggers(){

	InitializeCounter();
	InitializeFlags();
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

	N_JET_ONE							= 0;
	N_JET_THREE					  = 0;
	N_JET_FOUR 					  = 0;
	N_JET_FIVE						= 0;
	N_JET_SIX							= 0;
	N_ELECTRON						= 0;
	N_MUON 							  = 0;
	N_GAMMA 							= 0;
	N_MET 								= 0;
	N_HT 									= 0;
}
/******************************************************************************//**
* Initialise trigger call flags
**********************************************************************************/
void Triggers::InitializeFlags(){

	HT_CALL				= kFALSE;
	MULTIJET_CALL = kFALSE;
	EW_CALL				= kFALSE;
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
void Triggers::eventElectroWeak(Pythia *pythia)
{

	electron.resize(0);
	muon.resize(0);
	gamma.resize(0);

	for (int i = 0; i < pythia -> event.size(); ++i) {

		// Selection process for particles to be clustered
		if( pythia -> event[i].isFinal() ){

			if(pythia -> event[i].isVisible() && abs( pythia -> event[i].eta() ) 
				< FCAL_MAX_ETA  ){
				missingETvec += pythia -> event[i].p();
			}

			if( !pythia -> event[i].isVisible() ) continue;

			if(abs( pythia -> event[i].eta() ) > TRACKER_ETA) continue;

	    if(abs(pythia -> event[i].id()) == 11 && pythia -> event[i].pT() > ELECTRON_HL ){

	    	electron.push_back(i); // fill with electron candidate
	    }
	    if(abs(pythia -> event[i].id()) == 13 && pythia -> event[i].pT() > MUON_HL ){

	      muon.push_back(i); // fill with muon candidate
	    }
	    if(abs(pythia -> event[i].id()) == 22 && pythia -> event[i].pT() > GAMMA_HL ){

	      gamma.push_back(i); // fill with photon candidate
	    }
		}
	}// end of particle 
}
/******************************************************************************//**
* Pythia event for Tree level process (No ISR) 
**********************************************************************************/
void Triggers::eventTree(Pythia *pythia)
{

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
	JetClustering(fjInputs, jetDef);
}
/******************************************************************************//**
* Jet clustering algorithm
**********************************************************************************/
void Triggers::JetClustering(vector <fastjet::PseudoJet> fjInputs, 
																		 fastjet::JetDefinition jetDef){

	vector <fastjet::PseudoJet> inclusiveJets;

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
vector <double> Triggers::electrons(){ 

	return electron;
}
vector <double> Triggers::muons(){ 

	return muon;
}
vector <double> Triggers::gammas(){ 

	return gamma;
}
vector <double> Triggers::dielectrons(){ 

	return di_electron;
}
vector <double> Triggers::dimuons(){ 

	return di_muon;
}
vector <double> Triggers::digammas(){ 

	return di_gamma;
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
vector <double> Triggers::LeptonISO(Pythia *pythia, vector <double> candidate){

	vector <double> Iso;

	//run over lepton candidates
	for (unsigned i = 0; i < candidate.size(); i++){		
	  
	  double Pt_sum = 0;

		for (unsigned j = 0; j < pythia -> event.size(); ++j){
				
			// pseudorapidity difference 
			double dEta = pythia -> event[j].y() - pythia -> event[ candidate[i] ].y();

			// phi coordinate difference
			double dPhi = abs( pythia -> event[j].phi() - pythia -> event[ candidate[i] ].phi() );

			// account for phi coordinate system missmatch
			if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;

			// distance/metric in phi-eta system
			double del_R = sqrt( pow2(dEta) + pow2(dPhi) );

			// sum pt of constituents within cone radius
			if(del_R < LEP_R ) Pt_sum += pythia -> event[j].pT();
		}

		Iso.push_back(Pt_sum/pythia -> event[ candidate[i] ].pT());	
	}// end of lepton loop
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
void Triggers::MultiJetTrigger(){
	
	int count [5] = { };


	for (unsigned i = 0; i < jets.size(); i++) {

		if(jets[i].perp() > ONE_JET_HL_THIN ) 	count[0]++;
		if(jets[i].perp() > THREE_JET_HL    ) 	count[1]++;
		if(jets[i].perp() > FOUR_JET_HL     ) 	count[2]++;
		if(jets[i].perp() > FIVE_JET_HL     ) 	count[3]++;
		if(jets[i].perp() > SIX_JET_HL      ) 	count[4]++;

	}// end of jet loop

	if( count[0] > 0) 	N_JET_ONE++;
	if( count[1] > 2) 	N_JET_THREE++;
	if( count[2] > 3) 	N_JET_FOUR++;
	if( count[3] > 4) 	N_JET_FIVE++;
	if( count[4] > 5)	 	N_JET_SIX++;

	MULTIJET_CALL = kTRUE;
	
}
/******************************************************************************//**
* ElectroWeak Triggers
**********************************************************************************/
void Triggers::EWTrigger(Pythia *pythia){


	//if( LeptonISO( pythia, electon ) > ISO_CUT ) N_ELECTRON++;
	//if( LeptonISO( pythia, muon    ) > ISO_CUT ) N_MUON++;
	//if( LeptonISO( pythia, gamma   ) > ISO_CUT ) N_GAMMA++;

}
/******************************************************************************//**
* HT Trigger using native HT calculation
**********************************************************************************/
void Triggers::HtTrigger(){

	if( HT() > HT_HL ) N_HT++;

	HT_CALL = kTRUE;

}
/******************************************************************************//**
* Efficiency calls for various triggers
**********************************************************************************/
void Triggers::TriggerEfficiencies( double total ){

	ofstream file;

	file.open ("Zeff.dat", ios::app); 

	if(HT_CALL){

		file << Efficiency( N_HT , 				total) << " ";
	}
	if(MULTIJET_CALL){

		file << Efficiency( N_JET_ONE ,   total) << " "
				 << Efficiency( N_JET_THREE , total) << " "
				 << Efficiency( N_JET_FOUR ,  total) << " "
				 << Efficiency( N_JET_FIVE ,  total) << " "
				 << Efficiency( N_JET_SIX ,   total) << " ";
	}

	file << endl;

	file.close();
}