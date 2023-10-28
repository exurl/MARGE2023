# MARGE2023_2

Version archive of 2023 MARGE computational work, including:
- GVT testing
- wind tunnel testing
- experimental data acquisition, postprocessing, visualization scripts
- aeroservoelastic state-space model generation scripts
- aeroservoelastic state-space model simulation scripts

Note that this repo does not contain most GVT and wind tunnel data as it is too large to upload to github. Wind tunnel data can be found on the Boeing Control Collaboration shared drive.

Scripts for making and analyzing a model
| Name | Description |
|---|---|
| Anthony_ASE_SS_driver.m | defines parameters which define the model. calls functions to generate state-space model. saves state-space model. |
| Anthony_ASE_SS_actuator_generation.m | generates actuator state-space model |
| Anthony_ASE_SS_output_generation.m | generates plant sensing state-space model ([C] & [D]) |
| Anthony_ASE_SS_plant_generation.m | generates plant dynamics state-space model ([A] & [B]) |
| margeResponse.m | compute FRFs of state-space model (and also eigenvalue stability) |
| plotEigenvalues.m | eigenvalue stability plotting |
| margeFreqExpreriment.m | loads FRFs of state-space model (./FRF_ASE_SS.mat) and experiment (./windTunnel/wtData/FRF*) and plots comparison |
| margeTimeExperiment.m | loads experimental time-series data (./windTunnel/wtData/TIME*) and simulates corresponding state-space response and plots comparison. Takes a few minutes to run |

Scripts for optimization of a model
| Name | Description |
|---|---|
| margeOptDriver.m | driver script for model optimization |
| margeObjective.m | objective function which calls other functions below|
| margeComputeFRF.m | computes model FRF |
| margeCompareFRF.m | computes residual between model and experiment |

Miscellaneous scripts
| Name | Description |
|---|---|
| Anthony_ASE_SS_generation_uncorrected.m | old version of Anthony_ASE_SS_driver.m with no model corrections. |
| ASE_Model_Input_info.m | loads Marat's .mat files from ./ASEInputData/ |
| nsStrainStudy.m | quick study on convergence of strain output with number of modes |
| plotMARGEnodes.m | visualize NASTRAN node locations |
| plotMARGEnodesGVT.m | visualize GVT impact points on geometry |
| strainCalibration.m | quick study on strain gauge output compared to linear beam theory |
| czt_FRF.m | John's FRF generation function |
| sysCompareJohn.m | compare my state-space model with John's CIFER state-space model |
