% All 'Batch' scripts execute codes for a batch of countries, creating
% outputs that can be read into Excel files for sharing. 

% This code completes model calibration (using MCMC) for all countries
% listed.

% Dependencies:
% =============
% - Setup_model_popupdated.m
% - Get_calibrations_popupdated.m

% Iterate over all countries to simulate disruptions

clear all;
ctrys = {'PNG','ZWE','LSO','KEN','PAK','MYS'};

tic
% --- Make calibrations
for ictry = 1:length(ctrys)
    clearvars -except ctrys ictry;
    fprintf('%0.5g: ', ictry);
    ctry = ctrys{ictry};
    iso3 = ctrys{ictry};
    
    Setup_model_popupdated;
    Get_calibrations_popupdated;
    
    fprintf('\n');  close all;
end
toc