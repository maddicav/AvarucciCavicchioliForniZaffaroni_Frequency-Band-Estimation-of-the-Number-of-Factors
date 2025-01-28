# AvarucciCavicchioliForniZaffaroni_Frequency-Band-Estimation-of-the-Number-of-Factors

The repository contains the Matlab files to reproduce Tables 1-3, Figures 2-4 (main text, Section 4) and Tables S1-S9 (Online Supplement) of the paper 
                               “Frequency-Band Estimation of the Number of Factors”
                             by M. Avarucci, M. Cavicchioli, M. Forni and P. Zaffaroni (2025).
                             
In this document we give detailed information regarding the code files necessary to reproduce the methodology and the data used for the empirical analysis.
___________________________________________________________________________________________________________________________
Folder: “Simulations_OnlineSupplement”

This section describes the files, grouped by type, used in the Monte Carlo exercise.

Monte Carlo Simulations

The m.files listed below replicate the Monte Carlo Experiments described in Section S.2  of the Online Supplement to the paper Frequency-Band Estimation of the Number of Factors, by M. Avarucci, M. Cavicchioli, M. Forni and P. Zaffaroni (2025).

Please note that the random seed has not been fixed in the simulation experiments so the results may differ due to sample variability. 

•	Experiment1_HL.m replicates the Experiment described in Section S.2.1, Simulation Design: First Experiment (HL)), Table S.1 and Table S.2.  Users must select MA or AR loadings (additional details are provided in the m-file).

•	Experiment2_Onatski.m replicates the Experiment described in Section S.2.1, Simulation Design: Second Experiment (O)), Table S.3 and S.4. Users must select MA or AR loadings (additional details are provided in the m-file).

•	Experiment3.m replicates the Experiment described in Section S.2.1, Simulation Design: Third Experiment, Table S.5.  Users must select the “size” of the variance of the idiosyncratic component (additional details are provided in the m-file).

•	Experiment4.m replicates the Experiment described in Section 2.2.2, Simulation Design: Fourth Experiment (reduced rank spectral density), Table S.6.

•	ExperimentDSGE1.m replicates the Experiment described in Section 2.2.3, JPT Model and ACD Model, Table S.7.

•	 ExperimentDSGE2.m replicates the Experiment described in Section 2.2.3, BCR Model, Table S.8. Users must select the sample size (additional details are provided in the m-file).

•	CalibratingWindowSize.m replicates the Experiment described in Section 2.3, Table S.9.

Estimators of the number of factors

The m.files  listed below implement the estimators for the number of factors considered in the Monte Carlo Simulations. Additional details are provided in the m-files.

•	ACFZcrit.m returns the estimated number of factors using the DDR, DER and DGR estimators computed as average over all the Fourier frequencies (See Section S.2).  Requires the file DynamicEigenvaluesPERS.m.

•	DDR.m returns the estimated number of factors using the DDR estimator computed as average over the set of Fourier frequencies specified by the user (See Section S.2) .

•	AHcrit.m returns the estimated number of factors using the ER and GR estimators proposed by Ahn and Horenstein (2013).

•	HLcrit.m: returns the estimated number of factors using the log-information criterion proposed by Hallin and Liska (2007).

•	ONcrit.m returns the estimated number of factors using the test proposed by  Onatski (2009), Section 5.3. Requires the files dynamic.m and CVGUE.dat to compute the p-values of Ontaski’s dynamic test.

Data Generating Processes

•	HLmodel.m simulates a (n x T) panel data used in Experiment1_HL.m

•	ONmodel simulates a (n x T) panel data used in Experiment2_Onatski.m

•	DGP3model.m simulates a (n x T) panel data used in Experiment3.m

•	StopBandModel.m simulates a (n x T) panel data used in Experiment4.m

•	TrendCycleModel.m simulates a (n x T) panel data used in Experiment4.m

•	data_generator1.m simulates a (n x T) panel data used in ExperimentDSGE2.m. It requires the functions loadings.m (which loads the files A11.txt, A12.txt, F1.txt, F2.txt, Omega1.txt, Omega2.txt, Omega3.txt, B1.txt, rho.txt, sigma.txt) and parameters.m. The files generate the data as in Onatski and Ruge-Murcia (2013). Additional details are provided in the m-file.

•	solution_acd.mat contains the solution of the Angeletos Collard and Dellas (2020) model. results_jpt.mat contains the results of the replication material for the Justiniano, Primiceri, Tambalotti (2010) model. The files are used to generate the data in ExperimentDSGE1.m. Additional details are provided in the m-file.

Other files

•	standardize.m standardizes the data (the transformed data have mean zero and variance one).

____________________________________________________________________________________________________________________________
Folder : “Empirical Application_Section 4”

This section describes the files, grouped by type, used in Section 4 “Empirical Application: Dissecting the U.S. Economy”

Datasets

•	FREDreduced.xlsx contains QUARTERLY data taken from McCracken and Ng (2020) (FRED-QD) from which we extract the final sample of 216 time series from 1960:Q1 to 2020:Q1. 
The complete list of variables and their transformations is reported in Section S.4 of the Online Supplement.

•	FRED-MD.xlsx contains MONTHLY data for robustness analysis taken from McCracken and Ng (2016) and the final sample consists of 122 time series from 1960:M1 to 2020:M3.

Main program

•	main.m reproduces results reported in Tables 1-3, and in Figures 2-4 and discussed in Section 4 of the main text.

Estimators of the number of factors

•	ACFZcrit.m returns the estimated number of factors using the DDR, DER and DGR estimators computed as average over all the Fourier frequencies (see Section S.3).  Requires the file DynamicEigenvaluesPERS.m.

•	DDR.m returns the estimated number of factors using the DDR estimator computed as average over the set of Fourier frequencies specified by the user (See Section S.3).

•	HLcrit.m: returns the estimated number of factors using the log-information criterion proposed by Hallin and Liska (2007).

•	ONcrit.m returns the estimated number of factors using the test proposed by  Onatski (2009), Section 5.3. Requires the files dynamico.m and CVGUE.dat to compute the p-values of Ontaski’s dynamic test.

Other Subroutines

•	LoadQ.m loads the data from FREDreduced.xlsx and applies transformations to raw quarterly time series.

•	LoadM2.m loads the data from FRED-MD.xlsx and applies transformations to raw monthly time series.

•	standardize.m standardizes the data (the transformed data have mean zero and variance one).

•	DynamicEigenvalues.m obtains the periodogram smoothing estimator.

•	MakeTable.m is an auxiliary file to construct Table 3.
