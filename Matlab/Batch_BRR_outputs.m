% All 'Batch' scripts execute codes for a batch of countries, creating
% outputs that can be read into Excel files for sharing. This code organises
% organises the outputs required in the BRR framework. 

% Dependencies:
% =============
% BRR_Get_outputs.m

clear all;

% The list of countries
iso3s = {'AGO','BGD','BRA','CHN','IDN','IND','KEN','MMR','PAK','PER','PHL','RUS','UGA','VNM','ZAF'};

computing = 1;
nrows = 17;

if computing
    % Calculate BRR outputs for each country, and save in respective
    % folders as 'BRR_outputs.mat'
    
    for iis = 1:length(iso3s)
        clearvars -except iis iso3s;
        iso3 = iso3s{iis}
        BRR_Get_outputs;
        close all;
    end

else    
    % For each country, read the BRR outputs that have already been
    % calculated, collate across all countries, and save as a csv file named 'BRR_outputs.csv' 
    
    tbl = cell(nrows+1,length(iso3s));
    for iis = 1:length(iso3s)
        iso3 = iso3s{iis};
        load([iso3,'/BRR_outputs']);
        tbl{1,iis} = iso3;
        for irow = 1:nrows
            row = allmat_pct(irow,:);
            if ~isnan(row(2))
               tbl{irow+1, iis} = sprintf('%0.2g (%0.2g - %0.2g)',row(2), row(1), row(3));
            else
               tbl{irow+1, iis} = '--';
            end
        end
    end
    cell2csv('BRR_outputs.csv',tbl);

end