# AvarucciCavicchioliForniZaffaroni_Frequency-Band-Estimation-of-the-Number-of-Factors

The repository contains the Matlab files to reproduce Tables 1-3, Figures 2-4 (main text, section 4) and Tables S1-S9 (Online Supplement) of the paper Frequency-Band Estimation of the Number of Factors, by M. Avarucci, M. Cavicchioli, M. Forni and P. Zaffaroni (2024).
Please note that the random seed has not been fixed in the simulation experiments so the results will differ due to sample variability. 

Folder “Empirics Section 4”: run main.m to reproduce Tables 1,2,3, and Figures 2,4,5

Folder “Simulations OnlineSupplement”:
•	Run Experiment1_HL.m to replicate the First Experiment (Table S.1 and Table S.2). Please select MA or AR loadings (see the preamble in the m-file).
•	Run Experiment2_Onatski.m to replicate the Second Experiment (Table S.3 and Table S.4). Please select MA or AR loadings (see the preamble in the m-file).
•	Run Experiment 3.m to replicate the Third Experiment (Table S.5).  Please select the “size” of the variance of the idiosyncratic component (see the preamble in the m-file).
•	Run Experiment4.m to replicate the Fourth Experiment (Table S.6).
•	Run ExperimentDSGE1.m to replicate the experiment on the DSGE Models of JPT and ACD (Table S.7).
•	Run ExperimentDSGE2.m to replicate the experiment on the DSGE Models of BCR (Table S.8). Please select the sample size (see the preamble in the m-file).
•	Run CalibratingWindowSize.m to replicate the experiment to calibrate the window size (Table S.9).
