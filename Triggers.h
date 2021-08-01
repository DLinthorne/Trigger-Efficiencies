  //* 
//* Author: Dylan Lithorne
//*

#ifndef	ATLAS_Triggers_H
#define ATLAS_Triggers_H

// Headers and Namespaces.
#include "Pythia8/Pythia.h"
//#include "TH1.h"
//#include "TCanvas.h"
//#include "TVirtualPad.h"
//#include "TApplication.h"
//#include "TFile.h"
#include "Pythia8Plugins/FastJet3.h"
#include "Pythia8Plugins/CombineMatchingInput.h"

using namespace Pythia8;
using namespace std;

class Triggers
{
//**********************************************************************************
// Constants ( Trigger menus & geometry) *include constants for scouting*
//**********************************************************************************
public:
  //Jets angular parameter
  double R;

  //Lepton isolation cone radius
  double LEP_R, Eiso_R, Miso_R, Giso_R;

  //Calorimeter detector geometry constants
  double BARREL_ETA, ENDCAP_MIN_ETA, ENDCAP_MAX_ETA,
          HEC_MIN_ETA, HEC_MAX_ETA, FCAL_MIN_ETA, FCAL_MAX_ETA;

  //Tracker geometry constants
  double TRACKER_ETA;

  //FastJet Clustering Algorithm for General-KT: -1 = anti-kT; 0 = C/A; 1 = kT.
  double POWER;     

  //High level MET and Lepton thresholds 
  double MET_HL, ELECTRON_HL, MUON_HL, GAMMA_HL,
          DI_ELECTRON_HL, DI_MUON_HL;

  //High level multijet and HT thresholds 
  double HT_HL , THREE_JET_HL , FOUR_JET_HL , FIVE_JET_HL, SIX_JET_HL , ONE_JET_HL_FAT , ONE_JET_HL_THIN;

  //Low level MET and lepton thresholds
  double MET_LL, ELECTRON_LL , MUON_LL, GAMMA_LL ,
          DI_ELECTRON_LL , DI_MUON_LL;

  //Low level multijet and HT thresholds #3x50 and 4x15
  double HT_LL , THREE_JET_LL, FOUR_JET_LL, FIVESIX_JET_LL , ONE_JET_LL ;

  //Atlas leption isolation cut
  double ISO_CUT_Tight, ISO_CUT_Loose;

 //**********************************************************************************
 // Variables
 //**********************************************************************************
private:
  vector <fastjet::PseudoJet> fjInputs;
  vector <fastjet::PseudoJet> jets;
  vector <int> count;
  vector <double> electron_hl;
  vector <double> muon_hl;
  vector <double> gamma_hl;
  vector <double> electron_ll;
  vector <double> muon_ll;
  vector <double> gamma_ll;

  Vec4 missingETvec;
  Vec4 missingETvecJet; // initialization worries??? //

  int N_JET_ONE, N_JET_THREE, N_JET_FOUR, N_JET_FIVE, N_JET_SIX, N_HT,
    N_ELECTRON, N_MUON, N_GAMMA, N_MET, N_JET_ONE_FAT;
  int N_JET_ONE_L, N_JET_THREE_L, N_JET_FOUR_L, N_JET_FIVE_L, N_JET_SIX_L, N_HT_L, N_ELECTRON_L, N_MUON_L, N_GAMMA_L, N_MET_L;

  //Flags for OR triggers and calls
  bool MULTIJET_HL_FLAG, MULTIJET_HL_CALL, EW_HL_CALL, HT_HL_CALL, MET_HL_CALL;
  bool MULTIJET_LL_FLAG, MULTIJET_LL_CALL, EW_LL_CALL, HT_LL_CALL, MET_LL_CALL;

//**********************************************************************************
// Methods
//**********************************************************************************
public:
  Triggers();
  virtual ~Triggers();

  //Jet Algorithm
  void eventElectroWeak(Pythia *pythia,fastjet::JetDefinition jetDef);
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

  void   MultiJetTriggerHL();
  void   MultiJetTriggerLL();
  void   HtTriggerHL();
  void   HtTriggerLL();
  void   HLTriggerEfficiencies(char* name, double mass, double total);
  void   LLTriggerEfficiencies(char* name, double mass, double total);
  void   EWTriggerHL(Pythia *pythia);
  void   EWTriggerLL(Pythia *pythia);
  void   METTriggerHL();
  void   METTriggerHL(int scheme);
  void   METTriggerLL();
  void   METTriggerLL(int scheme);

  bool LeptonISO(Pythia *pythia, vector <double> candidate, int type);

  vector <fastjet::PseudoJet> SortedJets();
  vector <double> electronsHL();
  vector <double> muonsHL();
  vector <double> gammasHL();
  vector <double> electronsLL();
  vector <double> muonsLL();
  vector <double> gammasLL();

private:
  void JetClustering( fastjet::JetDefinition jetDef);
  void InitializeCounter();


};

#endif


