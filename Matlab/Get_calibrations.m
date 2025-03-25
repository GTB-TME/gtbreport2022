% Code to draw samples from the posterior density for a given country (and
% associated country data). Likelihood functions are as constructed in
% Setup_model. 

% Code has the following steps: 
% (i) Generate 1e4 parameter sets using latin hypercube sampling, and find 
% the 5 best-fitting sets
% (ii) For each set, perform an initial simplex optimisation to increase 
% the likelihood 
% (iii) Using the resulting 5 parameter sets as starting points for 
% independent MCMC chains, perform adaptive Bayesian MCMC to sample from 
% the posterior density. 

% The purpose of the initial optimisations is to assist the MCMC in finding 
% an appropriate starting point, in situations where mixing can be a 
% challenge (e.g. where the posterior density for any given parameter is 
% much narrower than the prior density).

% Dependencies:
% =============
% - get_objective2D.m
% - MCMC_adaptive.m

% This code is used by:
% =====================
% - Batch_Country_calibrations.m
% - Batch_Regional_calibrations.m

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% Set looping = 0 if this script is not being run by any Batch scripts, 
% and instead is being used only for a single country
looping = 1;
if ~looping
    clear all; iso3 = 'TLS';
end
load([iso3,'/Model_setup']);                                               % Created by Setup_Model

% Function to simulate epidemic upto 2019 for any given parameter set 'x'
obj = @(x) get_objective2D(x, prm, ref, sel, agg, gps, lhd, opts);

% --- Generate initial parameter sets using latin hypercube sampling

nsam = 1e4;
% nsam = 100;
nx   = 20;

out = -Inf;
while max(out) == -Inf
    lo = repmat(prm.bounds(1,1:nx),nsam,1);
    hi = repmat(prm.bounds(2,1:nx),nsam,1);
    sams = lo + (hi-lo).*lhsdesign(nsam,nx);
    
    out = zeros(1,nsam);
    inc = zeros(1,nsam);
    mk  = round(nsam/25);
    for ii = 1:nsam
        if mod(ii,mk)==0; fprintf('%0.5g ',ii/mk); end
        [out(ii), aux] = obj(sams(ii,:));
        if out(ii)>-Inf
            inc(ii) = aux.inc_all_2019;
        end
    end
    fprintf('\n');
end

% Order them to find the best-fitting parameter sets
tmp1 = [out;1:nsam]';
tmp2 = tmp1(~isnan(out),:);
mat  = sortrows(tmp2,-1);

ord1  = mat(:,2);
ordx1 = sams(ord1,:);


% --- For the 5 best-fitting parameter sets, perform a simplex optimisation 
% --- to increase their respective likelihoods 

x0 = ordx1(1:5,:);
for ii = 1:size(x0,1)
    [x1(ii,:), fval(ii)] = fminsearch(@(x)-obj(x),x0(ii,:),optimset('PlotFcns',@optimplotfval));    
end

% Order these updated parameter sets according to their likelihoods
tmp1  = [fval; 1:size(x0,1)]';
mat   = sortrows(tmp1,1);
ord2  = mat(:,2);
ordx2 = x1(ord2,:);


% --- For the 3 updated parameter sets having the highest likelihood,
% --- perform an adaptive Bayesian MCMC, running independent chains

xsto   = [];
outsto = [];
niter  = 3;
for ii = 1:niter
    [xsto(:,:,ii), outsto(ii,:)] = MCMC_adaptive(obj, ordx2(ii,:), 3e4, 1, [], [], [], true);
end
xsto1_sto = xsto;

[row,col] = find(outsto==max(outsto));
x0 = xsto(col(1),:,row(1));
mat = xsto(1:end,:,1); cov0 = cov(mat);
[xsto2, outsto2] = MCMC_adaptive(obj, x0, 3e4, 1, [], [], cov0, true);

xsto = xsto2;
outsto = outsto2;

[row,col] = find(outsto==max(outsto));
x0 = xsto(col(1),:,row(1));
mat = xsto(1:end,:,1); cov0 = cov(mat);
[xsto2, outsto2] = MCMC_adaptive(obj, x0, 1.5e5, 1, [], [], cov0, true);

xsto = xsto2;
outsto = outsto2;

save([iso3,'/model_fits']);
