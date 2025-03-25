% Script to calculate the outputs needed for the BRR framework, for any
% given country. For each country, code takes the posterior density
% acquired through model calibration, and performs simple calculations to
% estimate, for example, the duration of untreated TB, together with
% uncertainty intervals. Outputs are saved in each country folder as 'BRR_outputs.csv'

% Dependencies:
% - get_objective2D.m
% - goveqs_basis2.m
% - alloc_parameters.m
% - make_model.m

% This script is used by:
% - BRR_control.m



% Set looping = 0 if this script is not being run by BRR_control, and
% instead is being used only for a single country
looping = 1;
if ~looping
    clear all; iso3 = 'RUS';
end

load([iso3,'/model_fits']);                                                % Posterior sample, created by Get_calibrations_popupdated.m

prm.popn_turnover = 1;

% Take a 'thinned' sample of the posterior density
ix0 = round(size(xsto2,1)/2); nx = 250; dx = round((size(xsto2,1)-ix0)/nx);
xs = xsto2(ix0:dx:end,:,1);

% Make a selector to distinguish recent vs remote sources of incidence (see
% 'Setup_model' for description of selector matrices)
m = zeros(i.nstates);
m(s.I, s.Lf) = 1;
sel.recent = sparse(m);

m = zeros(i.nstates);
m(s.I, s.Ls) = 1;
sel.remote = sparse(m);


% -------------------------------------------------------------------------
% --- Do calculations where simulation is needed --------------------------

mk = round(size(xs,1)/25);
for ii = 1:size(xs,1)
    if mod(ii,mk) == 0; fprintf('%0.5g ',ii/mk); end 

    % Perform a simulation of the TB epidemic leading upto 2019     
    obj = @(x) get_objective2D(x, prm, ref, sel, agg, gps, lhd, opts);
    [out,aux] = obj(xs(ii,:));
    soln      = aux.soln;
    sfin      = soln(end,:);

    % --- Prevalence of LTBI
    ltbi_prev(ii) = sum(sfin([s.Lf, s.Ls]))*100;
    
    % --- Proportion of new incidence arising from recent infection
    allin = aux.allin;
    tmp1 = (allin.*sel.recent)*sfin(1:i.nstates)';
    tmp2 = (allin.*sel.remote)*sfin(1:i.nstates)';
    prop_recent(ii) = sum(tmp1)/sum(tmp1+tmp2)'*100;
    
    % --- Annual Risk of TB Infection
    lam = aux.M0.lambda;
    arti(ii) = sum(lam*sfin(1:i.nstates)')*100;
    
    % --- Incidence in 2019 and 5-year trends
    inc = sum(diff(soln(:,i.aux.inc)),2);
    num = inc(end);
    den = inc(end-5);
    incd_trend(ii) = ((num/den)^0.2-1)*100;
    incd(ii) = num*1e5;

    % --- Mortality in 2019 and 5-year trends
    mrt = sum(diff(soln(:,i.aux.mort)),2);
    num = mrt(end);
    den = mrt(end-5);
    mort_trend(ii) = ((num/den)^0.2-1)*100;
    mort(ii) = num*1e5;
    
    % --- Prevalence in 2019
    prev(ii) = sum(sfin(s.prevalent))*1e5;
    
    % --- HIV indicators: Proportion of TB incidence in 2019 that was HIV
    % coinfected, and prevalence of HIV. For countries where HIV
    % coinfection accounts for <10% of TB incidence, opts.hiv is set to
    % zero.
    if opts.hiv
        tb_incid_hiv(ii) = aux.inc_h1/aux.inc_all_2019;
        hiv_prev(ii) = aux.HIV_prev;
    else
        tb_incid_hiv(ii) = NaN;
        hiv_prev(ii) = NaN;
    end
        
end
fprintf('\n');


% -------------------------------------------------------------------------
% --- Do calculations where only cohort is needed (i.e. no tranmsission) --

r0 = r; p0 = p;

mk = round(size(xs,1)/25);
for ii = 1:size(xs,1)
    if mod(ii,mk) == 0; fprintf('%0.5g ',ii/mk); end 
    
    % Allocate parameter values associated with the 'i'th sample from the
    % posterior density
    [r,p,prm] = alloc_parameters(xs(ii,:),r0,p0,prm,xi,opts);

    % Switch off transmission, population turnover and HIV acquisition, for purpose of
    % cohort modelling
    r.beta    = 0*r.beta;
    r.relapse = 0*r.relapse;
    p.imm     = ones(size(p.imm));
    prm.popn_turnover = 0;
    prm.rHIV = 0*prm.rHIV;

    % Construct model equations
    M = make_model(p, r, i, s, gps, opts);

    
    % ---------------------------------------------------------------------
    % --- 5-year incidence ------------------------------------------------

    % Create a cohort of people with recent infection
    init = zeros(1,i.nx);
    init(i.Lf.h0) = 1;

    % Simulate cohort outcomes
    geq = @(t,in)goveqs_basis2(t, in, M, i, prm, sel, agg);
    [t,soln] = ode15s(geq, [0:100], init, odeset('NonNegative',[1:i.nstates]));

    % Incidence at 5-years, and then over remainder of lifetime
    cinc          = soln(:,i.aux.inc(1));
    cinc_5yr(ii)  = cinc(t==5)*100;
    cinc_life(ii) = (cinc(end) - cinc(t==5))/66*100;
    
    
    % ---------------------------------------------------------------------
    % --- Duration, treated TB --------------------------------------------
    
    % Create a cohort of HIV-ve people having active TB
    init = zeros(1,i.nx);
    init(i.I.h0) = 1;
    
    % Simulate cohort outcomes
    geq = @(t,in)goveqs_basis2(t, in, M, i, prm, sel, agg);
    [t,soln] = ode15s(geq, [0 100], init, odeset('NonNegative',[1:i.nstates],'Refine',64));
    
    % Find the point at which number still having active TB is 50% of
    % intial value - this represents the median duration in this
    % compartment
    ind = find(soln(:,i.I.h0)<0.5,1,'first');
    meandur_treated(ii) = t(ind);
    
    
    % ---------------------------------------------------------------------
    % --- Case fatality and duration, untreated TB ------------------------
    
    % Create a cohort of HIV-ve people having active TB
    init = zeros(1,i.nx);
    init(i.I.h0) = 1;

    % Switch off treatment, and construct updated model equations
    r.Tx_init     = 0*r.Tx_init;
    r.casefinding = 0;
    M             = make_model(p, r, i, s, gps, opts);
    
    % Simulate cohort outcomes 
    geq = @(t,in)goveqs_basis2(t, in, M, i, prm, sel, agg);
    [t,soln] = ode15s(geq, [0 100], init, odeset('NonNegative',[1:i.nstates],'Refine',64));
    
    % Find proportion who have died
    morts = soln(:,i.aux.mort);
    cfr(ii) = morts(end,1)*100;

    % Find the point at which number still having active TB is 50% of
    % intial value - this represents the median duration in this
    % compartment
    ind = find(soln(:,i.I.h0)<0.5,1,'first');
    meandur_untreated(ii) = t(ind);
    
end
fprintf('\n');



% -------------------------------------------------------------------------
% --- Protection given previous exposure ----------------------------------

% Can read these directly as model parameters
priorimm = xs(:,xi.p_imm)'*100;
beta     = xs(:,xi.r_beta(1))';


% -------------------------------------------------------------------------
% --- Bring them all together in same order as BRR form -------------------

mat1 = [cinc_5yr; cinc_life; cfr; meandur_untreated; priorimm];
mat2 = [incd; incd_trend; mort; mort_trend; prev; tb_incid_hiv; hiv_prev];
mat3 = [ltbi_prev; prop_recent; arti; beta; meandur_treated];

allmat = [mat1; mat2; mat3];
allmat_pct = prctile(allmat,[2.5,50,97.5],2);

save([iso3,'/BRR_outputs'],'allmat_pct');
