#include "Triggers.cc"

// Let Pythia8:: be implicit.
int main(int argc, char* argv[]){

	
	Pythia *pythia = new Pythia();

	fastjet::JetDefinition jetDef(fastjet::genkt_algorithm, 0.5, -1);
	
	pythia ->readFile("ModA.dat");
	int nEvent = pythia -> mode("Main:numberOfEvents");

  Triggers *Trig = new Triggers();

  tool -> InitializeMatching(pythia);

  pythia -> init();

	for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

			if (!pythia -> next()) continue;
			Trig -> eventTree(pythia, jetDef);
		  Trig -> MultiJetTrigger();
			Trig -> HtTrigger();
			Trig -> EWTrigger(pythia);
	}

	Trig-> TriggerEfficiencies( pythia -> info.nAccepted() );

	return 0;
}