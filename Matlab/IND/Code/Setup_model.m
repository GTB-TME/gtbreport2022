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


clear all;

iso3 = 'IND';

% --- Set up compartment lookups (see 'Readme' file) ----------------------
gps.provs = {'pu','pr'};
states1 = {'U','Lf','Ls','I','E','Rlo','Rhi','R'};
states2 = {'Tx'};

[i, s, d, lim] = get_addresses({states1}, [], [], [], 0);
[i, s, d, lim] = get_addresses({states2, gps.provs}, i, s, d, lim);
d = char(d);

s.infectious = [s.I, s.E];
s.prevalent  = [s.infectious, s.Tx];


% --- Include auxiliary indices: these are for additional integrals used to 
% count things such as incidence or notifications in given time intervals
names = {'inc','noti','mort'};
lgths =     [1,     2,     1];
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    i.aux.(names{ii}) = inds;
    lim = inds(end);
end
i.nx = lim;


% --- Make aggregators and selectors, used to count incidence-like terms 
% (see 'Readme' file) -----------------------------------------------------

% Selectors for the incidence
tmp = zeros(1,i.nstates);
tmp(1,s.I) = 1;
agg.inc = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.I,:) = 1;
sel.inc = tmp - diag(diag(tmp));

% Selectors for notifications
tmp = zeros(2,i.nstates);
tmp(1,intersect(s.Tx,s.pu)) = 1;
tmp(2,intersect(s.Tx,s.pr)) = 1;
agg.noti = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.Tx,:) = 1;
sel.noti = tmp - diag(diag(tmp));


% --- Define variables (parameters varied in the MCMC), and the ranges for
% their prior uniform distributions ---------------------------------------
names = {'r_beta','r_Tx_init','r_g','rf_mort_TB','p_Tx_complete','rf_self_cure','rf_progression','rf_LTBI_stabil','rf_reactivation','rf_relapse','p_imm'};
lgths =        [1,          2,    1,           1,              2,             1,               1,               1,                1,           3,      1];
xi = []; lim = 0;
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    xi.(names{ii}) = inds;
    lim = inds(end);
end
% xi.calib = xi.r_mort_TB;
xi.nx = lim;

bds = [];
bds(xi.r_beta,:)          = [0 30];
bds(xi.r_Tx_init,:)       = repmat([0.1 10],2,1);
bds(xi.r_g,:)             = [0 10];
bds(xi.rf_mort_TB,:)      = [0.1 1.9];
bds(xi.p_Tx_complete,:)   = [0.75 0.95; 0.4 0.8];
bds(xi.rf_self_cure,:)    = [0.5, 1.5];
bds(xi.rf_progression,:)  = [0.5, 1.5];
bds(xi.rf_LTBI_stabil,:)  = [0.5, 1.5];
bds(xi.rf_reactivation,:) = [0.5, 1.5];
bds(xi.rf_relapse,:)      = repmat([0.5, 1.5],3,1);
bds(xi.p_imm,:)           = [0.5 0.9];
prm.bounds = bds';


% --- Define baseline parameter values, organising all parameters as
% per-capita rates (r) or proportions (p)

% Natural history
r.progression  = 0.0826;
r.LTBI_stabil  = 0.872;
r.reactivation = 0.0006;
r.self_cure    = 1/6;
r.mort_TB      = 1/6;
r.relapse      = [0.032 0.14 0.0015];
r.mort         = 1/66;
p.imm          = 0.5;

% Treatment stage
r.Tx_init     = [1 1];
r.Tx          = 2;
p.Tx_complete = [.85 0.5];
r.default     = r.Tx*(1-p.Tx_complete)./p.Tx_complete;
p.cure        = [1 1];


% --- Bring them all together
prm.p = p; prm.r = r;
ref.i = i; ref.s = s; ref.d = d; ref.xi = xi;
% prm.bounds = bds';


% --- Get calibration targets ---------------------------------------------

% Estimated TB incidence and mortality
load('../../Data/TB Estimates/estim_data_2b.mat');
countryrow = estims(strcmp(estims.iso3,iso3),:);
popn       = countryrow.e_pop_num;
mat        = reshape(countryrow{1,3:end},[3,6])';

% data.inc_2019  = mat(1,:);
data.inc_2015  = [272 304 347];
data.inc_2019  = [228 259 295];                                            % From Sandip estimates
data.mort_all  = mat(5,:);

popn = countryrow.e_pop_num;                 
data.noti = 1721067/popn*1e5*[0.9 1 1.1];                                  % Taking public only, from 1/1/2019 - 31/12/2019 in Nikshay data

f1a = get_distribution_fns(data.inc_2015,  'lognorm', 0);
f1b = get_distribution_fns(data.inc_2019,  'lognorm', 0);
f2  = get_distribution_fns(data.noti, 'lognorm', 0);
f3  = get_distribution_fns(data.mort_all, 'lognorm', 0);

lhd.fn  = @(inc_2015, inc_2019, noti, mort) f1a(inc_2015) + f1b(inc_2019) + f2(noti) + f3(mort);
lhd.sgn = -Inf;

save Model_setup;