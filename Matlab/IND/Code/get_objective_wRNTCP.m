% For a given parameter set x, simulate the epidemic to 2019 and evaluate
% the posterior density.

% Function first simulates a disease-free equilbrium, in the absence or
% ART, before including factors that drive temporal trends: HIV and ART
% (for countries with a substantial burden of HIV); the expansion of public
% sector TB treatment through DOTS; and any increases in case-finding over
% time (for countries without a substantial burden of HIV or a large
% private sector).


% Arguments:
% ==========
% - x: Vector of parameter values being varied in MCMC
% - prm: Any other model parameters
% - sel, agg: Respectively, selector and aggregator matrices used for
% counting incidence-like terms (see Setup_model for explanation)
% - gps: Stratifications for state variables, e.g. HIV status
% - calfn: function to evaluate the posterior density, given inputs for
% incidence, mortality, etc
% - opts: opts.hiv set to 1 for countries with >10% of TB being HIV
% coinfected (0 otherwise) and opts.provs set to 1 for countries with a
% strong role for the private sector (0 otherwise)

% Outputs:
% ========
% - out: log-Posterior density
% - aux: Package of model outputs including incidence, mortality in 2019,
% etc

% Dependencies:
% =============
% - alloc_parameters.m
% - make_model.m
% - goveqs_basis2.m
% - goveqs_scaleup2D.m

% This code is used by:
% =====================
% - Get_calibrations.m

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


function [out, aux] = get_objective_wRNTCP(x, prm, ref, sel, agg, gps, calfn)

r = prm.r; p = prm.p; i = ref.i; s = ref.s; xi = ref.xi;

mat = [prm.bounds(2,1:length(x))-x; x-prm.bounds(1,1:length(x))];
% All values should be positive if x is within the specified ranges

if min(mat(:)) < 0
    % Apply uniform priors: if model parameters are outside uniform ranges
    % specified in prm.bounds, then assign log-likelihood of -Inf
    
    out = Inf*calfn.sgn;
    % +Inf if using sum of squares, -Inf if using likelihoods
    aux = NaN;
    %     keyboard;
else
    
    % Assign all values in parameter vector x to associated model
    % parameters
    [r,p] = alloc_parameters(x,r,p,xi);
    
    % --- Set up the necessary models -------------------------------------
    
    % Final conditions
    p2 = p; r2 = r;
    M2 = make_model(p2, r2, i, s, gps);
    
    % With RNTCP, and without additional case-finding
    p1 = p; r1 = r;
    r1.g = 0;
    M1 = make_model(p1, r1, i, s, gps);
    
    % In absence of RNTCP, and without additional case-finding
    p0 = p; p0.pu = 0;
    r0 = r; r0.Tx_init(1) = 0; r0.g = 0;
    M0 = make_model(p0, r0, i, s, gps);
    
    % Equilibrium model
    init = zeros(1,i.nx); seed = 1e-6;
    init(i.U) = (1-seed); init(i.I) = seed;
    geq = @(t,in)goveqs_basis2(t, in, M0, i, s, p, sel, agg);
    [t0, soln0] = ode15s(geq, [0 2e3], init, odeset('NonNegative',[1:i.nstates]));
    
    
    % --- Solve the models ------------------------------------------------
    
    % RNTCP scale-up
    [t1, soln1] = ode15s(@(t,in) goveqs_scaleup(t, in, M0, M1, [1997 2007], i, s, p, sel, agg), [1997:2015], soln0(end,:), odeset('NonNegative',[1:i.nstates]));
    
    % Increase in case-finding from 2015 onwards
    [t2, soln2] = ode15s(@(t,in) goveqs_scaleup(t, in, M1, M2, [2015 2019], i, s, p, sel, agg), [2015:2019], soln1(end,:), odeset('NonNegative',[1:i.nstates]));
    
    t     = [t1; t2(2:end)];
    soln  = [soln1; soln2(2:end,:)];
    
    sfin  = soln(end,:);
    sdiff = diff(soln,1);
    
    % --- Get the objectives ----------------------------------------------
    
    % Get 2015 incidence
    incd_2015_sim = sdiff(t==2015,i.aux.inc)*1e5;
    incd_2019_sim = sdiff(end,i.aux.inc)*1e5;
    noti_sim = sdiff(end,i.aux.noti)*1e5;
    mort_sim = sdiff(end,i.aux.mort)*1e5;
    
    % Compose the objective function
    out = calfn.fn(incd_2015_sim, incd_2019_sim, noti_sim(1), mort_sim);
    %     out = calfn.fn(incd_sim, noti_sim(1));
    
    % --- Get additional outputs and package ------------------------------
    
    aux.soln      = [soln, t];
    aux.soln0     = [soln0, t0];
    aux.incd_2015 = incd_2015_sim;
    aux.incd_2019 = incd_2019_sim;
    aux.noti      = noti_sim;
    aux.mort      = mort_sim;
    aux.M1        = M1;
end

end
