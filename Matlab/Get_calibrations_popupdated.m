looping = 0;
if ~looping
    clear all; iso3 = 'MNG';
end
load([iso3,'/Model_setup_popupdated']);

% prm.bounds(2,1) = 200;

obj = @(x) get_objective2D(x, prm, ref, sel, agg, gps, lhd, opts);

% --- Generate initial guesses

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

% Order them to get the best-fitting ones
tmp1 = [out;1:nsam]';
tmp2 = tmp1(~isnan(out),:);
mat  = sortrows(tmp2,-1);

ord1  = mat(:,2);
ordx1 = sams(ord1,:);
% ordi = inc(ord1);
% ordo = out(ord1);
% ordx2 = ordx1(1:5,:);

% % --- Optimise best ones --------------------------------------------------
x0 = ordx1(1:5,:);
% x0 = ordx1(1:2,:);
for ii = 1:size(x0,1)
    [x1(ii,:), fval(ii)] = fminsearch(@(x)-obj(x),x0(ii,:),optimset('PlotFcns',@optimplotfval));    
end

% Order them
tmp1  = [fval; 1:size(x0,1)]';
mat   = sortrows(tmp1,1);
ord2  = mat(:,2);
ordx2 = x1(ord2,:);
% ordo2 = fval(ord2);

% --- Do the sampling -----------------------------------------------------
xsto = [];
outsto = [];
niter = 3;
% niter = 5;
for ii = 1:niter
%     [xsto(:,:,ii), outsto(ii,:)] = MCMC_adaptive(obj, ordx2(ii,:), 1e4, 1, [], [], [], true);
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

save([iso3,'/model_fits_popupdated']);
