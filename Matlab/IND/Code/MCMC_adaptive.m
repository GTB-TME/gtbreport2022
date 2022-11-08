% Adaptive MCMC, using Haario et al:
% https://link.springer.com/article/10.1007/s11222-008-9110-y
% and
% http://probability.ca/jeff/ftpdir/adaptex.pdf

% Arguments:
% ==========
% F:          Function giving log-posterior density for a parameter set x
% x0:         Initial value of parameter set x
% n:          Number of iterations
% cov0:       Initial covariance matrix
% fac:        Scaling factor for covariance matrix. Set fac = 1 for default
% fixinds:    Elements of x that should be held fixed. Set fixinds = [] for full MCMC
% blockinds:  Number of 'epi parameters' (e.g. beta, X2) in x, if we want to vary epi and non-epi parameters as independent 'blocks'. Set blockinds = [] if we want a full covariance matrix
% displ:      Structure with display options. Set displ = true to show progress

% Outputs:
% ========
% xsto:        Matrix giving all parameter values sampled from the posterior density
% outsto:      Vector of values of the posterior density corresponding to xsto
% history:     Record of the outcome of each proposal, i.e. acceptance or rejection
% accept_rate: The proportion of proposed values that were accepted

function [xsto, outsto, history, accept_rate] = MCMC_adaptive(F, x0, n, sigma, fixinds, blockind, cov0, displ)

d = length(x0); b = 0.05; sd = sigma*2.4^2/d;
if ~isempty(fixinds)
    inds = fixinds(1,:); vals = fixinds(2,:);
else
    inds = []; vals = [];
end

% Checks on the initial covariance matrix
if isempty(cov0)
    cov0 = eye(d); cov0(inds,:) = 0; cov0(:,inds) = 0;
    cov0(1:blockind,blockind+1:end) = 0; cov0(blockind+1:end,1:blockind) = 0;
end

% Initiate the output matrices
xsto = zeros(d,n); outsto = zeros(1,n);
history = zeros(d+1,n);                                                    % Rows: 1:d Proposed values 4. Accept or reject

xsto(:,1) = x0(:); xbar = xsto;
FX = F(x0); outsto(1) = FX;
acc = 0;

if displ; figure; end

% --- Start the MCMC loop -------------------------------------------------
for t = 2:n
    
    X = xsto(:,t-1);
    
    % --- Make a proposal from the distribution
    Y0 = mvnrnd(X,0.1^2*cov0*sigma/d);
    if t < 2*d
        Y = max(Y0,0); Y(inds) = vals;
    else
        ind0 = 1; ind1 = t-1;
        covmat = cov(xsto(:,ind0:ind1)');
        covmat(inds,:) = 0; covmat(:,inds) = 0;
        covmat(1:blockind,(blockind+1:end)) = 0; covmat((blockind+1):end,1:blockind) = 0;
        % Need to modify this bit so it goes recursively - faster
        
        covmat = (covmat+covmat')/2+1e-5;
        % Precaution to make sure values don't get negative
        Y = max((1-b)*mvnrnd(X,sd*covmat) + b*Y0,0);
        Y(inds) = vals;
    end
    history(1:d,t) = Y;
    
    % --- Decide whether to accept or not
    FY = F(Y);
    if rand < exp(FY-FX)
        % Accept
        xsel = Y(:);
        FX = FY;
        acc = acc+1;
        history(end,t) = 1;
    else
        % Reject
        xsel = xsto(:,t-1);
    end
    xsto(:,t) = xsel;
    outsto(t) = FX;
    xbar(:,t) = (xbar(:,t-1)*(t-1) + xsel)/t;
    
    % Display options
    if displ && (mod(t,round(n/25))==0); fprintf('%0.5g ', t/n*25); end
    if displ && (mod(t,200)==0)
        plot(xsto(:,1:t-1)'); xlim([0 n]); drawnow;
    end
end

fprintf('\n');
accept_rate = acc/n;
xsto = xsto';
history = history';