% All 'Batch' scripts execute codes for a batch of countries, creating
% outputs that can be read into Excel files for sharing. 

% This code completes model calibration (using MCMC) for all Regions
% listed. The same as Batch_Country_calibrations.m, but for Regions.

% Dependencies:
% =============
% - Setup_model_Regional.m
% - Get_calibrations.m

clear all; 
iso3 = 'region_AFR'
Setup_model_Regional;
Get_calibrations;
close all;

clear all; 
iso3 = 'region_AMR'
Setup_model_Regional;
Get_calibrations;
close all;

clear all; 
iso3 = 'region_EUR'
Setup_model_Regional;
Get_calibrations;
close all;