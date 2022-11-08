% All 'Batch' scripts execute codes for a batch of countries.

% This code completes alignment of model simulations with WHO estimates,
% ensuring that both have the same uncertainty intervals

% Dependencies:
% =============
% - Extrapolate_Countries.m

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

clear all;
ctrys = {'AGO','AZE','BGD','BRA','CHN','COL','ETH','IDN','KAZ','KEN','KGZ','KHM','LSO','MEX','MMR','MYS','NPL','PAK','PER','PHL','PNG','PRK','ROU','RUS','THA','TLS','TUR','VNM','ZAF','ZWE','region_AFR','region_AMR','region_EUR'};

tic
% --- Extrapolate to fit WHO estimates and write results
for ictry = 1:length(ctrys)
    clearvars -except ctrys ictry;
    fprintf('%0.5g: ', ictry);
    ctry = ctrys{ictry};
    Extrapolate_Countries;
    fprintf('\n');  close all;
end
toc