# Triggers-Efficiencies                    
**Requirements: Pythia8, FASTJET3 and ROOT6**

[![arXiv](https://img.shields.io/badge/arXiv-2103.08620%20-green.svg)](https://arxiv.org/abs/2103.08620)

This repository contains holds ATLAS triggering objects necessary for making sense of signal efficiencies of emerging jet Monte Carloevents, usually generated from `MadGraph` or `PYTHIA`. See [2103.08620](https://arxiv.org/abs/2103.08620) and [1502.05409
](https://arxiv.org/abs/1502.05409
) for further details on both the triggering stategies and dark showering models. Note that `PYTHIA` (> 8.226) includes the Hidden Valley module used for simulating dark showers. The objects are compiled alongside the `PYTHIA` makefile with `Root6`, `FASTJET3` plugin flags included. 
## Current Implementation:

 * Single & Multijet Triggers: let clustering using `FASTJET3`.
 * Lepton & Photon Triggers: isolation algorithm for tight and loose lepton candidates.
 * HT Triggers: seeded from single/multijet lower level trigger.
 * MET Trigger: both Topo-cluster style & Jet style methods.

## Pythia Example

- With the proper MAKE script
```
make efftemp
```
- Run the example triggering script using model file `model_script_X.dat`
```
./efftemp model_script_X.dat 
```
## References

If you use this code, please cite our paper:

```
@article{Linthorne:2021oiz,
    author = "Linthorne, Dylan and Stolarski, Daniel",
    title = "{Triggering on Emerging Jets}",
    eprint = "2103.08620",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    month = "3",
    year = "2021"
}
```
  
