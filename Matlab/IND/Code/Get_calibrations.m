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


clear all; load Model_setup;

tic
obj = @(x) get_objective_wRNTCP(x, prm, ref, sel, agg, gps, lhd);

% --- Generate initial guesses

nsam = 1e3;
nx   = xi.nx;

lo = repmat(prm.bounds(1,1:nx),nsam,1);
hi = repmat(prm.bounds(2,1:nx),nsam,1);
sams = lo + (hi-lo).*lhsdesign(nsam,nx);

out = zeros(1,nsam);
mk  = round(nsam/25);
for ii = 1:nsam
    if mod(ii,mk)==0; fprintf('%0.5g ',ii/mk); end
    out(ii) = obj(sams(ii,:));
end
fprintf('\n');

% Order them to get the best-fitting ones
% out = out(~isnan(out));
% mat = sortrows([out;1:length(out)]',-1);

tmp1 = [out;1:nsam]';
tmp2 = tmp1(~isnan(out),:);
mat  = sortrows(tmp2,-1);

x0   = sams(mat([1:3],2),:);
cov0 = cov(sams(mat([1:10],2),:))/1e3;


% --- Do the sampling -----------------------------------------------------
xsto = [];
niter = 2;
for ii = 1:niter
    [xsto(:,:,ii), outsto(ii,:)] = MCMC_adaptive(obj, x0(ii,:), 1e5, 1, [], [], [], true);
end

% Do a second round with the best-fitting parameter set
[row,col] = find(outsto==max(outsto(:)));
% row1 = row(1); col1 = col(1);

mat   = xsto(:,:,row(1));
cov0  = cov(mat);
xinit = xsto(end,:,row(1));

[xsto2, outsto2] = MCMC_adaptive(obj, xinit, 5e5, 1, [], [], cov0, true);

csvwrite('out.csv',xsto2);

save model_fits2b xsto2 outsto2;
toc