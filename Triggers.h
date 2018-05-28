//* 
//* Author: Dylan Lithorne
//*

#ifndef	ATLAS_Triggers_H
#define ATLAS_Triggers_H

// Headers and Namespaces.
#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TFile.h"
#include "Pythia8Plugins/FastJet3.h"
#include "Pythia8Plugins/CombineMatchingInput.h"

using namespace Pythia8;
using namespace std;

class Triggers
{
//**********************************************************************************
// Constants ( Trigger menus & geometry)
//**********************************************************************************
public:
  //Jets angular parameter
  const double R = 0.4;

  //Lepton isolation cone radius
  const double LEP_R = 0.5;   //*** Is this correct? ****

  //Calorimeter detector geometry constants
  const double BARREL_ETA = 1.37, ENDCAP_MIN_ETA = 1.52, ENDCAP_MAX_ETA = 2.37,
          HEC_MIN_ETA = 2.37, HEC_MAX_ETA = 3.1, FCAL_MIN_ETA = 3.1, FCAL_MAX_ETA = 4.9;

  //Tracker geometry constants
	const double TRACKER_ETA = 2.49;

	//FastJet Clustering Algorithm for General-KT: -1 = anti-kT; 0 = C/A; 1 = kT.
  const double POWER   = -1;     

  //High level MET and Lepton thresholds 
	const double MET_HL = 110, ELECTRON_HL = 26, MUON_HL = 26, GAMMA_HL = 140,
					DI_ELECTRON_HL = 17, DI_MUON_HL = 14;

  //High level multijet and HT thresholds 
  const double HT_HL = 850, THREE_JET_HL = 175, FOUR_JET_HL = 100, FIVE_JET_HL      = 70, SIX_JET_HL = 60, ONE_JET_HL_FAT = 420, ONE_JET_HL_THIN = 380;


	//Atlas leption isolation cut
	const double ISO_CUT = 0.1;

 //**********************************************************************************
 // Variables
 //**********************************************************************************
private:
	vector <fastjet::PseudoJet> fjInputs;
  vector <fastjet::PseudoJet> jets;
  vector <int> count;
  vector <double> electron;
  vector <double> muon;
  vector <double> gamma;
  vector <double> di_electron;
  vector <double> di_muon;
  vector <double> di_gamma;

	Vec4 missingETvec;
  Vec4 missingETvecJet; // initialization worries??? //

  int N_JET_ONE, N_JET_THREE, N_JET_FOUR, N_JET_FIVE, N_JET_SIX, N_HT,
    N_ELECTRON, N_MUON, N_GAMMA, N_MET;

  //Flags for OR triggers and calls
  bool MULTIJET_FLAG, MULTIJET_CALL, EW_CALL, HT_CALL;

//**********************************************************************************
// Methods
//**********************************************************************************
public:
	Triggers();
  virtual ~Triggers();

  //Jet Algorithm
  void eventElectroWeak(Pythia *pythia);
  void eventJet(Pythia *pythia);
  void eventTree(Pythia *pythia);
  void eventTree(Pythia *pythia, fastjet::JetDefinition jetDef);

  void InitializeMatching(Pythia *pythia);
  void InitializePythia(Pythia *pythia);
  void InitializeFlags();


  double TopoMET();
  double JetMET();
  double HT();
  double Efficiency(double accepted, double total);
  void   MultiJetTrigger();
  void   HtTrigger();
  void   TriggerEfficiencies(double total);
  void   EWTrigger(Pythia *pythia);

  vector <double> LeptonISO(Pythia *pythia, vector <double> candidate);
  vector <fastjet::PseudoJet> SortedJets();
  vector <double> electrons();
  vector <double> muons();
  vector <double> gammas();
  vector <double> dielectrons();
  vector <double> dimuons();
  vector <double> digammas();

private:
  void JetClustering( vector <fastjet::PseudoJet> fjInputs,  fastjet::JetDefinition jetDef);
  void InitializeCounter();


};

#endif


