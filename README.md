                  #################### Quick Start Guide ####################
The easiest way to get this method ("HAPOD Scoring") running and provide you with accurate predictions of the correct ternary protein-PROTAC binding pose, is to start with the few available ternary docking tools, e.g. PRosettaC in its online server version (https://prosettac.weizmann.ac.il/pacb/steps), or its local version (https://github.com/LondonLab/PRosettaC), and also Integrative PROTAC-Model (https://github.com/gaoqiweng/PROTAC-Model). Conjugate usage with PRosettaC outputs are demonstrated in our publication (see also PRosettaC.HAPOD folder here).

To quickly get things running, please download the QuickStartExample.zip file (https://github.com/JL2021MD/PROTACModeling/raw/main/QuickStartExample.zip) from your browser, or if are familiar with downloading specific folders, download the QuickStartExample folder (https://github.com/JL2021MD/PROTACModeling/tree/main/QuickStartExample). It provides you with an example project with 5 candidate poses ready (CRBN-BRD4 system, obtained from the above mentioned local ternary dock tool Integrative PROTAC-Model) for going through the whole HAPOD scoring process. Going through this example project will take only ~15 minutes because all the MD runs can be skipped as outputs are readily available. Once familiar, the same scripts can be used to run your own projects, or rerun the example project without using the already provided results.

If you find our method helpful please cite this paper: In Silico Modeling and Scoring of PROTAC-Mediated Ternary Complex Poses. J. Med. Chem. 2022, 65, 8, 6116–6132, https://pubs.acs.org/doi/full/10.1021/acs.jmedchem.1c02155 , along with any relevant ternary docking method paper(s) if used.

The straightforward idea is that the correct, native pose having the most favorable binding affinity will be most stable and resistant to the short linear heating trials, while non-native poses generally do not. The method is meant to be easy to use involving only classic MD for the best simplicity. The scripts in this repo run natively with AMBER, although the idea and method should work the same with other engines and force fields. Bulk linear heating of the whole system, used in the method, while can be a bit brute, works very effectively due to the large binding interface of protein-proteins, even though it may not be precise enough for small molecule-protein systems. The explicit OPC waters will stabily remain liquid unless the simulated temperature is >650K or so under NPT conditions.

A quick pose refinement before the heating-accelerated pose departure (HAPOD) trials will be important because this refinement will help the interface residues and sidechains fit better with each other, compared to docking outputs, or those from a force field different than the one used in HAPOD.

Notes: 
Please note that in most of the MD runs used here the gamma_ln is set to 0.01 which provides a bit of an increase in conformation change speed compared to gamma_ln=2

For assistance please contact liaojunzhuo at gmail.com
iMessages: +16315900825
