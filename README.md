The repository contains the Matlab files and the datasets to reproduce Tables 1-3, Figures 2-4 (main text, Section 4) and Tables S1-S9 (Online Supplement) of the manuscript JASA-T&M-2024-0780.R1

In this document, we give detailed information regarding the code files necessary to reproduce the simulation results and the empirical analysis.
_________________________________________________________________________________________
Folder: “Simulations_OnlineSupplement”

This section describes the files used in the Monte Carlo Simulations, grouped by type.

Monte Carlo Simulations

The m.files listed below replicate the Monte Carlo Experiments described in Section S.2, Tables S.1-S.9 of the Online Supplement to the paper Frequency-Band Estimation of the Number of Factors, by M. Avarucci, M. Cavicchioli, M. Forni and P. Zaffaroni (2025).

Please note that the random seed has not been fixed in the simulation experiments so the results may differ due to sample variability.

•Experiment1_HL.m replicates the Experiment described in Section S.2.1, Simulation Design: First Experiment (HL)), Table S.1 and Table S.2. Users must select MA or AR loadings (additional details are provided in the m-file).

•Experiment2_Onatski.m replicates the Experiment described in Section S.2.1, Simulation Design: Second Experiment (O)), Table S.3 and S.4. Users must select MA or AR loadings (additional details are provided in the m-file).

•Experiment3.m replicates the Experiment described in Section S.2.1, Simulation Design: Third Experiment, Table S.5. Users must select the “size” of the variance of the idiosyncratic component (additional details are provided in the m-file).

•Experiment4.m replicates the Experiment described in Section 2.2.2, Simulation Design: Fourth Experiment (reduced rank spectral density), Table S.6.

•ExperimentDSGE1.m replicates the Experiment described in Section 2.2.3, JPT Model and ACD Model, Table S.7.

•ExperimentDSGE2.m replicates the Experiment described in Section 2.2.3, BCR Model, Table S.8. Users must select the sample size (additional details are provided in the m-file).

•CalibratingWindowSize.m replicates the Experiment described in Section 2.3, Table S.9.

Estimators of the number of factors

The m.files listed below implement the estimators for the number of factors considered in the Monte Carlo Simulations. Additional details are provided in the m-files.

•ACFZcrit.m returns the estimated number of factors using the DDR, DER and DGR estimators computed as average over all the Fourier frequencies (See Section S.2). It requires the file DynamicEigenvaluesPERS.m.

•DDR.m returns the estimated number of factors using the DDR estimator computed as average over the set of Fourier frequencies specified by the user (See Section S.2). It requires the file DynamicEigenvaluesPERS.m.

•AHcrit.m returns the estimated number of factors using the ER and GR estimators proposed by Ahn and Horenstein (2013).

•HLcrit.m: returns the estimated number of factors using the log-information criterion proposed by Hallin and Liska (2007).

•ONcrit.m returns the estimated number of factors using the test proposed by Onatski (2009), Section 5.3. It requires the files dynamico.m and CVGUE.dat to compute the p-values of Ontaski’s dynamic test.

Data Generating Processes

•HLmodel.m simulates the panel data used in Experiment1_HL.m

•ONmodel simulates the panel data used in Experiment2_Onatski.m

•DGP3model.m simulates the panel data used in Experiment3.m

•StopBandModel.m simulates the panel data used in Experiment4.m

•TrendCycleModel.m simulates the panel data used in Experiment4.m

•data_generator1.m simulates the panel data used in ExperimentDSGE2.m. It requires the functions loadings.m (which loads the files A11.txt, A12.txt, F1.txt, F2.txt, Omega1.txt, Omega2.txt, Omega3.txt, B1.txt, rho.txt, sigma.txt) and parameters.m. The files generate the data as in Onatski and Ruge-Murcia (2013). Additional details are provided in the m-file

•solution_acd.mat contains the solution of the Angeletos Collard and Dellas (2020) model. results_jpt.mat contains the results of the replication material for the Justiniano, Primiceri, Tambalotti (2010) model. The files are used to generate the data in ExperimentDSGE1.m. Additional details are provided in the m-file.

Other auxiliary files

•DynamicEigenvaluesPERS.m. computes the eigenvalues of the smoothed periodogram.

•standardize.m standardizes the data (the transformed data have mean zero and variance one).

_________________________________________________________________________________
Folder : “Empirical Application_Section 4”

This section describes the files, grouped by type, used in Section 4 “Empirical Application: Dissecting the U.S. Economy”

Datasets and Transformations

•FREDreduced.xlsx contains QUARTERLY data from McCracken and Ng (2020) (FRED-QD_2020m1). The complete list of variables and their transformations is reported in Section S.4 of the Online Supplement.

•FRED-MD.xlsx contains MONTHLY data for robustness analysis taken from McCracken and Ng (2016) – data vintage 2023 – 05.csv

•LoadQ.m loads the data from FREDreduced.xlsx into the Matlab workspace and applies transformations to raw quarterly time series.

•LoadM2.m loads the data from FRED-MD.xlsx into the Matlab workspace and applies transformations to raw monthly time series.

Main program

•main.m reproduces results reported in Tables 1-3, and in Figures 2-4 and discussed in Section 4 of the main text.

Estimators of the number of factors

•ACFZcrit.m returns the estimated number of factors using the DDR, DER and DGR estimators computed as average over all the Fourier frequencies (see Section S.3). It requires the file DynamicEigenvaluesPERS.m. to calculate the eigenvalues of the smoothed periodogram.

•DDR.m returns the estimated number of factors using the DDR estimator computed as average over the set of Fourier frequencies specified by the user (See Section S.3). It requires the file DynamicEigenvaluesPERS.m.

•HLcrit.m: returns the estimated number of factors using the log-information criterion proposed by Hallin and Liska (2007).

•ONcrit.m returns the estimated number of factors using the test proposed by Onatski (2009), Section 5.3. It requires the files dynamico.m and CVGUE.dat to compute the p-values of Ontaski’s dynamic test.


Other auxiliary files

•standardize.m standardizes the data (the transformed data have mean zero and variance one).

•DynamicEigenvalues.m returns the periodogram smoothing estimator, the required numbers of largest eigenvalues and the associated eigenvectors.
•DynamicEigenvaluesPERS.m. which computes the eigenvalues of the smoothed periodogram.

•MakeTable.m is an auxiliary file to construct Table 3.

_________________________________________________________________________________
Folder : “Demo”

While the other folders enable result reproducibility, the large number of files and the length of the codes may make it challenging to follow the methodology. This folder serves an illustrative purpose and includes a demo file demonstrating the implementation of the proposed estimators. 

It includes the files:

•demo.m estimates the number of shocks in simulated dataset using the DER, DGR and DDR criteria considering all frequencies, and the DDR criteria using a use-chosen frequency band. Run the code to obtain the estimated number of shocks in a synthetic dataset using the function ACFZcrit.m and DDR.m

•ACFZcrit.m returns the estimated number of factors using the DDR, DER and DGR estimators computed as average over all the Fourier frequencies). Requires the file DynamicEigenvaluesPERS.m. to calculate the eigenvalues of the smoothed periodogram.

•DDR.m returns the estimated number of factors using the DDR estimator computed as average over the set of Fourier frequencies specified by the user. It requires the file DynamicEigenvaluesPERS.m.

•README_demo is tutorial-style README file describing the file demo.m in details.
