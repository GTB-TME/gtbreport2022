% All 'Collate' scripts construct tables of outputs for a batch of
% countries. 

% This script brings together results of forward model projections (from
% 2020 onwards) and organises them into a csv file. Same as
% Batch_country_projections.m, but for Regional outputs.

% Dependencies:
% =============
% - Reads in monthly and annual projections (created by Simulate_disruption.m) 
% for incidence and mortality for each of the Regions listed

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

clear all;
regs = {'region_AFR','region_AMR','region_EUR'};

% --- Construct the column labels for the output tables
tmp = {};
mos = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
yrs = [2019:2024]; ct = 1;
for yi = 1:length(yrs)
    for mi = 1:length(mos)
        tmp{ct} = [mos{mi}, ' ', num2str(yrs(yi))]; ct = ct+1;
    end
end
colnames_mo = [{'Region','ISO3','Scenario','xx','xx'},tmp];

tmp = {}; ct = 1;
for yi = 1:length(yrs)
    tmp{ct} = [num2str(yrs(yi))]; ct = ct+1;
end
colnames_yr = [{'Region','ISO3','Scenario','xx','xx'},tmp];


% --- Then read the files (created by Simulate_disruptions.m) and write into 
% country-specific csv files

incmo = colnames_mo; incyr = colnames_yr; 
mrtmo = colnames_mo; mrtyr = colnames_yr;
for ir = 1:length(regs)
   reg = regs{ir}; 
   
   tmp = readcell([reg,'/Incidence_projections_mo.csv']);
   incmo = [incmo; tmp];
   
   tmp = readcell([reg,'/Incidence_projections_yr.csv']);
   incyr = [incyr; tmp];

   tmp = readcell([reg,'/Mortality_projections_mo.csv']);
   mrtmo = [mrtmo; tmp];

   tmp = readcell([reg,'/Mortality_projections_yr.csv']);
   mrtyr = [mrtyr; tmp];
end

cell2csv('Regions_INCmo.csv',incmo);
cell2csv('Regions_INCyr.csv',incyr);
cell2csv('Regions_MRTmo.csv',mrtmo);
cell2csv('Regions_MRTyr.csv',mrtyr);