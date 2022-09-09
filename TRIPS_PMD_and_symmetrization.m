clc
clear
close all
%% include the toolbox path 

addpath(genpath('./sspsoct_toolbox'));

%% load data
% load calibration data
load('CalStru.mat');

% load a big frame containing 5 downsampled frames from different locations in a volume
% struct name: bigStokes
load('bigS123.mat');


%% symetrization and PMD correction

PS = PSTriple_Manager_v2(CalStru);

% check the orthogolity of Stokes
PS.checkStokes(bigStokes.S1,bigStokes.S2,bigStokes.S3);

% reconstruct Mueller matrix
PS.load_Stokes(bigStokes.S1,bigStokes.S2,bigStokes.S3);

% calculate symmetric matrix
PS.obtainSymmetrizationCorrection(PS.snrmask(10));

% plot what happens after applying symmetric matrix
PS.checkSymmetrizationCorrection();        

% apply symmetric matrix
PS.applySymmetrizationCorrection();


% search correction vectors of rotation and diattenuation 
PS.obtainPMDCorrection(squeeze(mean(PS.snrmask(10),4)));

% plot what happens after applying correction
PS.checkPMDCorrection();  

% apply PMD correction
PS.applyPMDCorrection();

% save the correction data to CalStru (calibration structure)
CalStru = PS.getCalStru();   
save("CalStru_cor.mat","CalStru");
