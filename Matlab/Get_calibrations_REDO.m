warning off;

looping = 0;
if ~looping
    clear all; iso3 = 'MNG';
end
% load([iso3,'/Model_setup']);
load([iso3,'/model_fits_popupdated.mat']);

obj = @(x) get_objective2D(x, prm, ref, sel, agg, gps, lhd, opts);

% --- Do the sampling -----------------------------------------------------

[row,col] = find(outsto==max(outsto));
x0 = xsto(col(1),:,row(1));
mat = xsto(1:end,:,1); cov0 = cov(mat);
[xsto2, outsto2] = MCMC_adaptive(obj, x0, 1.5e5, 1, [], [], cov0, true);

xsto = xsto2;
outsto = outsto2;

save([iso3,'/model_fits_popupdated2']);
