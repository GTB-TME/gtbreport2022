% Set up all aspects of model structure, including definitions of state
% variables; compartment lookups (see Readme file); default parameter
% values; and labels of parameters that are varied in the MCMC

% Dependencies:
% =============
% - get_addresses.m
% - get_distribution_fns.m

% This code is used by:
% =====================
% - Outputs saved in 'Model_setup' serve as inputs for
% Get_calibrations.m

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% Set looping = 0 if this script is not being run by any Batch script, 
% and instead is being used only for a single country
looping = 0;
if ~looping
    clear all; iso3 = 'MNG';
end


% -------------------------------------------------------------------------
% --- Read in information on whether the country has a high HIV burden, or
% has a large private healthcare sector -----------------------------------

[~,~,raw]  = xlsread('Country_modes.xlsx');
col        = find(strcmp(raw(1,:),iso3));
country    = raw{2,col};
opts.hiv   = raw{3,col};
opts.provs = raw{4,col};

gps.hiv   = {'h0','h1','hart'};
gps.provs = {'pu','pr'};

nhiv   = length(gps.hiv);
nprovs = length(gps.provs);

% --- Set up compartment lookups (see 'Readme' file) ----------------------
states1 = {'U','Lf','Ls','I','Rlo','Rhi','R'};
states2 = {'Tx'};
[i, s, d, lim] = get_addresses({states1, gps.hiv}, [], [], [], 0);
[i, s, d, lim] = get_addresses({states2, gps.hiv, gps.provs}, i, s, d, lim);
d = char(d);

s.infectious = s.I;
s.prevalent  = [s.infectious, s.Tx];

% --- Include auxiliary indices: these are for additional integrals used to 
% count things such as incidence or notifications in given time intervals
names = {'inc', 'noti', 'mort'};
lgths =  [nhiv,  nprovs,     2];
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    i.aux.(names{ii}) = inds;
    lim = inds(end);
end
i.nx = lim;


% --- Make aggregators and selectors, used to count incidence-like terms 
% (see 'Readme' file) -----------------------------------------------------

% Selectors for the incidence
tmp = zeros(nhiv,i.nstates);
for ih = 1:nhiv
    tmp(ih,intersect(s.I,s.(gps.hiv{ih}))) = 1;
end
agg.inc = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.I,:) = 1;
if opts.hiv
    tmp(s.h1,s.h0) = 0; tmp(s.hart,s.h1) = 0; tmp(s.h0,s.h1) = 0; tmp(s.h1,s.hart) = 0;
end
sel.inc = tmp - diag(diag(tmp));

% Selectors for notifications
tmp = zeros(nprovs,i.nstates);
for ip = 1:nprovs
    tmp(ip,intersect(s.Tx,s.(gps.provs{ip}))) = 1;
end
agg.noti = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.Tx,:) = 1;
if opts.hiv
    tmp(s.h1,s.h0) = 0; tmp(s.hart,s.h1) = 0; tmp(s.h0,s.h1) = 0; tmp(s.h1,s.hart) = 0;
end
sel.noti = tmp - diag(diag(tmp));

% --- Define variables (parameters varied in the MCMC), and the ranges for
% their prior uniform distributions ---------------------------------------
names = {'r_beta','r_Tx_init','p_HIVlam','r_casefinding','rf_mort_TB','p_Tx_complete','r_ART_init','r_HIV_mort','p_HIV_relrate','rf_self_cure','rf_progression','rf_LTBI_stabil','rf_reactivation','rf_relapse','p_imm'};   
lgths =        [2,          2,         1,              1,           2,              1,           1,           1,              1,             1,               1,               1,                1,           3,      1];

xi = []; lim = 0;
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    xi.(names{ii}) = inds;
    lim = inds(end);
end
xi.calib = xi.p_HIV_relrate;

bds = [];
bds(xi.r_beta,:)          = repmat([0 30],2,1);
bds(xi.r_Tx_init,:)       = repmat([0.1 4],nprovs,1);
bds(xi.p_HIVlam,:)        = [0 10];
bds(xi.r_casefinding,:)   = [0 10];
bds(xi.rf_mort_TB,:)      = [0.5 1.5; 1 10];
bds(xi.p_Tx_complete,:)   = [0.75 0.95];
bds(xi.r_ART_init,:)      = [0 1];
bds(xi.r_HIV_mort,:)      = [0 1];
bds(xi.p_HIV_relrate,:)   = [1 100];
bds(xi.rf_self_cure,:)    = [0.5, 1.5];
bds(xi.rf_progression,:)  = [0.5, 1.5];
bds(xi.rf_LTBI_stabil,:)  = [0.5, 1.5];
bds(xi.rf_reactivation,:) = [0.5, 1.5];
bds(xi.rf_relapse,:)      = repmat([0.5, 1.5],3,1);
bds(xi.p_imm,:)           = [0.5 0.9];
if strcmp(iso3,'MNG')
    bds(xi.r_casefinding,:)   = [0 1e-6];
    bds(xi.r_Tx_init,:)       = repmat([0 4],nprovs,1);
    bds(xi.r_beta,:)          = repmat([0 100],2,1);
    bds(xi.rf_mort_TB,:)      = [0 1.5; 1 10];
end
prm.bounds = bds';


% --- Define baseline parameter values, organising all parameters as
% per-capita rates (r) or proportions (p)

% Natural history
r.progression0  = 0.0826*[1 10 10*0.4];
r.LTBI_stabil0  = 0.872*[1 0 1];
r.reactivation0 = 0.0006*[1 10 10*0.4];
r.self_cure0    = 1/6*[1 0 1];
r.mort_TB0      = 1/6;
r.relapse0      = [0.032 0.14 0.0015];
r.mort          = 1/66;
p.imm           = [0.8 0 0.8];

% Treatment stage
r.Tx          = 2;
p.Tx_complete = 0.8;
r.default     = r.Tx*(1-p.Tx_complete)./p.Tx_complete;

% Bring them all together
prm.p = p; prm.r = r;
ref.i = i; ref.s = s; ref.d = d; ref.xi = xi;
prm.popn_turnover = 1;


% --- Get calibration targets ---------------------------------------------

% Estimated TB incidence and mortality
load('Data/TB Estimates/estim_data_2b.mat');
countryrow = estims(strcmp(estims.iso3,iso3),:);
popn       = countryrow.e_pop_num;
mat        = reshape(countryrow{1,3:end},[3,6])';

% Make adjustments in case upper/lower estimates are same as central (e.g. 
% for MEX)
inds = find(mat(:,1)==mat(:,2));
mat(inds,1) = mat(inds,2)*0.95;
inds = find(mat(:,2)==mat(:,3));
mat(inds,3) = mat(inds,2)*1.05;

data.inc_2019  = mat(1,:);
data.inc_h1    = mat(2,:);
data.mort_H0   = mat(3,:);
data.mort_H1   = mat(4,:);
data.mort_all  = mat(5,:);
data.inc_2014  = mat(6,:);

% In case of DPRK where mortality estimates are unavailable in the csv
% file, incorporate unpublished estimates
if strcmp(iso3,'PRK')
    pt = 71.43; sd = 11.82;                                                % From Philippe's email on 18 Jul 2022
    data.mort_H0 = pt + [-2 0 2]*sd;
end

% Notifications
load('Data/TB Notifications/notif_data.mat');
data.noti_pu = notifs_new{strcmp(notifs_new.iso3,iso3),2}/popn*1e5*[0.9 1 1.1];

% Proportion of ART on TPT
prm.p.covTPT = pTPT{strcmp(pTPT.iso3,iso3),2};

% HIV burden and services
load('Data/HIV data/HIV_estims.mat');

ind = find(strcmp(countries1,country));
if isempty(ind); error('No country match in data'); end
prm.rHIV = HIV_incd(:,2,ind)'/popn;

ind = find(strcmp(countries2,country));
if isempty(ind); error('No country match in data'); end
data.ART_covg = ARTcovg_2019(ind,:)/100;
if data.ART_covg(end) == 1; data.ART_covg(end) = 0.95; end
prm.ART_start = ART_start(ind);

ind = find(strcmp(countries3,country));
if isempty(ind); error('No country match in data'); end
data.HIV_prev = HIVprev_2019(ind,:)/popn;

% --- Set up log-likelihood functions for each of the calibration targets 
show = 0;                                                                  % Select show = 1 to visualise each of the fitted distributions
f1   = get_distribution_fns(data.inc_2019,  'lognorm', show);
f1a  = get_distribution_fns(data.inc_2014,  'lognorm', show);

f2 = get_distribution_fns(data.mort_H0, 'lognorm', show);
f3   = get_distribution_fns(data.noti_pu, 'lognorm', show);
f3pu = get_distribution_fns([0.3 0.5 0.7], 'beta', show);

if opts.hiv
    f4 = get_distribution_fns(data.inc_h1,   'lognorm', show);
    f5 = get_distribution_fns(data.ART_covg, 'lognorm', show);
    f6 = get_distribution_fns(data.HIV_prev, 'lognorm', show);
    f7 = get_distribution_fns(data.mort_H1,  'lognorm', show);
end

% Bring the function handles together to construct an overall
% log-likelihood function
if opts.hiv
    lhd.fn  = @(inc_2019, mort_H0, noti_pu, inc_h1, ART_covg, HIV_prev, mort_H1) f1(inc_2019) + f2(mort_H0) + f3(noti_pu) + f4(inc_h1) + f5(ART_covg) + f6(HIV_prev) + f7(mort_H1);
elseif opts.provs
    lhd.fn  = @(inc_2019, mort_H0, noti_pu) f1(inc_2019) + f2(mort_H0) + f3(noti_pu);
else
    lhd.fn  = @(inc_2019, inc_2014, mort_H0, noti_pu) f1(inc_2019) + f1a(inc_2014) + f2(mort_H0) + f3(noti_pu);
end
lhd.fn_pu = @(prop_pu) f3pu(prop_pu)/10;
lhd.sgn = -Inf;

save([iso3,'/Model_setup']);