# Triggers-Efficiencies
Requirements: Pythia8, FASTJET3 and ROOT6

This repository holds classes of CERN lab trigger objects usefull in calculating efficencies for various jet, electroeeak and MET events. efftemp.cc is an example script demonstrates how the Trigger class can be used in efficiencie calculations for a defined trigger menu. The objects are compiled alongside the pythia8 makefile with root6, FASTJET3 plugin flags included. 

***Current Implementation:***

 * Single & Multijet Triggers
    * (working on Fat, R = 1.0 , vs slim, R = 0.4, single jet triggers)
 * Lepton & Photon Triggers
    * Isolation algorithm for tight criteria
 * HT Triggers
 * MET Trigger
    * Topo-cluster style algorithm
    * Jet style algorithm 
 
 ***Work in Progress:***
 
 * All corresponding lower level seed triggers
 * Trigger menu config file
 * Fat (R = 1.0) vs slim (R = 0.4) single jet triggers
    * Multiple jet clustering, recusively using slim jet to seed fat jet

 


  
