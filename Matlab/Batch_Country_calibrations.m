% All 'Batch' scripts execute codes for a batch of countries, creating
% outputs that can be read into Excel files for sharing. 

% This code completes model calibration (using MCMC) for all countries
% listed.

% Dependencies:
% =============
% - Setup_model.m
% - Get_calibrations.m

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

tic
clear all; close all;
iso3 = 'THA'
Setup_model;
Get_calibrations;
close all;
toc

tic
clear all; 
iso3 = 'PRK'
Setup_model;
Get_calibrations;
close all;
toc