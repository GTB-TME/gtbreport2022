% All 'Collate' scripts construct final tables of outputs for a batch of
% countries. 

% This script brings together results of forward model projections (from
% 2020 onwards) and organises them into a csv file.

% Dependencies:
% =============
% - Reads in monthly and annual projections (created by Simulate_disruption.m) 
% for incidence and mortality for each of the countries listed

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


clear all;
iso3s = {'AGO','AZE','BGD','BRA','CHN','COL','ETH','IND','IDN','KAZ','KEN','KGZ','KHM','LSO','MEX','MMR','MNG','MYS','NPL','PAK','PER','PHL','PNG','PRK','ROU','RUS','THA','TLS','TUR','VNM','ZAF','ZWE'};

tis = {'mo','yr'};                                                         % Labels for monthly or annual timeseries

% --- Construct the column labels for the output tables
tmp = {};
mos = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
yrs = [2019:2024]; ct = 1;
for yi = 1:length(yrs)
    for mi = 1:length(mos)
        tmp{ct} = [mos{mi}, ' ', num2str(yrs(yi))]; ct = ct+1;
    end
end
colnames_mo = [{'ISO3','Scenario','xx','xx'},tmp];

tmp = {}; ct = 1;
for yi = 1:length(yrs)
    tmp{ct} = [num2str(yrs(yi))]; ct = ct+1;
end
colnames_yr = [{'ISO3','Scenario','xx','xx'},tmp];

colnames = {colnames_mo, colnames_yr};


% --- Then read the files (created by Simulate_disruptions.m) and write into 
% country-specific csv files
for imy = 1:2
    incmat = [];
    mrtmat = [];
    
    for ii = 1:length(iso3s)
        tmp1 = readcell([iso3s{ii},'/Incidence_projections_',tis{imy},'.csv']);
        incmat = [incmat; tmp1];
        
        tmp2 = readcell([iso3s{ii},'/Mortality_projections_',tis{imy},'.csv']);
        mrtmat = [mrtmat; tmp2];
    end
    cell2csv(['Countries_INC',tis{imy},'.csv'], [colnames{imy}; incmat]);
    cell2csv(['Countries_MRT',tis{imy},'.csv'], [colnames{imy}; mrtmat]);
end