# URCsAnalyticalModel
This repository contains a suite of MATLAB scripts that implement the URCs analytical model for predicting battery discharge performance, which is described in the following paper: "A Fast Analytical Model for Predicting Battery Performance Under Mixed Kinetic Control", Hongxuan Wang, Fan Wang, and Ming Tang

Files in the repository:
1. "URCs_Model_Discharge_ConstEps.m" calculates the normalized galvanostatic discharge capacity and pseudo-steady-state distribution of salt/electrolyte potential/Li surface concentration of battery cells based on input cell parameters and discharge current density. 
2. "DefineParam_NMC_Graphite.m" provides an example of the MATLAB structure array that defines the cell parameters used by "URCs_Model_Discharge_ConstEps.m".
3. "NMC111_OCV.csv" tabulates the open-circuit voltage of the NMC cathode as a function of Li stoichiometry. It defines the cathode OCV profile in "DefineParam_NMC_Graphite.m"
4. "URCs_TrialScript.m" includes two application examples of the URCs model. The first one predicts the rate performance of NMC/Graphite full cells at various cathode thickness and different discharge C rates. The second one makes contour plots of the cell-level specific capacity and energy of NMC half cells at 1C discharge in the parameter space of cathode thickness and porosity. This script calls "URCs_Model_Discharge_ConstEps.m" and "DefineParam_NMC_Graphite.m". 


Contact: mingtang@rice.edu