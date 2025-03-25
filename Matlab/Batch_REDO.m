
tic
clear all; close all;
iso3 = 'region_AMR'
Get_calibrations_REDO;
close all;
toc

return;

tic
clear all; close all;
iso3 = 'ZAF'
Setup_model;
Get_calibrations_REDO;
close all;
toc

tic
clear all; close all;
iso3 = 'LSO'
Setup_model;
Get_calibrations_REDO;
close all;
toc

tic
clear all; close all;
iso3 = 'TUR'
Setup_model;
Get_calibrations_REDO;
close all;
toc


% tic
% clear all; close all;
% iso3 = 'ROU'
% Setup_model;
% Get_calibrations_REDO;
% close all;
% toc

% tic
% clear all; 
% iso3 = 'KAZ'
% Setup_model;
% Get_calibrations_REDO;
% close all;
% toc
% 
% tic
% clear all; close all;
% iso3 = 'RUS'
% Setup_model;
% Get_calibrations_REDO;
% close all;
% toc
% 
% tic
% clear all; 
% iso3 = 'ZAF'
% Setup_model;
% Get_calibrations_REDO;
% close all;
% toc