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


function [out, aux] = get_objective2D(x, prm, ref, sel, agg, gps, calfn, opts)

r = prm.r; p = prm.p; i = ref.i; s = ref.s; xi = ref.xi;
mat = [prm.bounds(2,1:length(x))-x; x-prm.bounds(1,1:length(x))];

if min(mat(:)) < 0
    % Apply uniform priors: if model parameters are outside uniform ranges
    % specified in prm.bounds, then assign log-likelihood of -Inf
    out = Inf*calfn.sgn;
    aux = [];
else
    
    % Assign all values in parameter vector x to associated model
    % parameters
    [r,p,prm] = alloc_parameters(x,r,p,prm,xi,opts);
    
    % --- Set up the necessary models -------------------------------------
    
    if opts.hiv || opts.provs
        
        % Pre-ART, pre-pu conditions
        p0 = p; r0 = r;
        r0.ART_init = 0;
        r0.Tx_init(1) = r0.Tx_init(1)*(1-opts.provs);
        M0 = make_model(p0, r0, i, s, gps, opts);
        
        % Public sector, no ART
        p1 = p; r1 = r;
        r1.ART_init = 0;
        M1 = make_model(p1, r1, i, s, gps, opts);
        
        % ART, no public sector
        p2 = p; r2 = r;
        r2.Tx_init(1) = r2.Tx_init(1)*(1-opts.provs);
        M2 = make_model(p2, r2, i, s, gps, opts);
                
        % --- Simulate the models
        
        % Equilibrium model
        init = zeros(1,i.nx); seed = 1e-6;
        init(i.U.h0) = (1-seed); init(i.I.h0) = seed;
        geq = @(t,in)goveqs_basis2(t, in, M0, i, prm, sel, agg);
        [t0, soln0] = ode15s(geq, [0:2e3], init, odeset('NonNegative',[1:i.nstates]));
        
        % Model programmatic changes after 2000: growth of the public sector
        % from 2000 - 2009 (for public/private models), and implementation of ART
        % (for TB/HIV models)
        init  = soln0(end,:);
        tinit = min([prm.ART_start, 2000]);
        geq = @(t,in) goveqs_scaleup2D(t, in, M0, M1, M2, [2000 2009; prm.ART_start 2019], i, prm, sel, agg);
        [t, soln] = ode15s(geq, [tinit:2019], init, odeset('NonNegative',[1:i.nstates]));
        
        allin = -M0.lin + M1.lin + M2.lin;
    else
        
        % Equilibrium model
        p0 = p; r0 = r;
        r0.casefinding = 0;                                                % Equivalent to parameter 'g' in the model equations
        M0 = make_model(p0, r0, i, s, gps, opts);
        
        % Ramp-up in case-finding
        p1 = p; r1 = r;        
        M1 = make_model(p1, r1, i, s, gps, opts);
        
        % Solve the model
        init = zeros(1,i.nx); seed = 1e-6;
        init(i.U.h0) = (1-seed); init(i.I.h0) = seed;
        geq = @(t,in)goveqs_basis2(t, in, M0, i, prm, sel, agg);
        [t0, soln0] = ode15s(geq, [0:2e3], init, odeset('NonNegative',[1:i.nstates]));
        
        init  = soln0(end,:);
        geq = @(t,in) goveqs_scaleup(t, in, M0, M1, [2014 2019], i, prm, sel, agg);
        [t, soln] = ode15s(geq, [2010:2019], init, odeset('NonNegative',[1:i.nstates]));

        allin = M1.lin;
    end
    
    
    sfin  = soln(end,:);
    sdiff = diff(soln,1);
    
    % --- Get the outcomes used to assess the posterior density -----------
    
    inc_all_2019 = sum(sdiff(end,i.aux.inc))*1e5;                          % Incidence in 2019
    inc_all_2014 = sum(sdiff(end-5,i.aux.inc))*1e5;                        % Incidence in 2014
    inc_h1       = sum(sdiff(end,i.aux.inc([2,3])))*1e5;                   % HIV +ve incidence in 2019
    noti         = sdiff(end,i.aux.noti)*1e5;                              % Total treatment initiations in 2019
    noti_pu      = noti(1);                                                % Treatment initiations amongst notifying providers in 2019
    prop_pu      = noti(1)/sum(noti);                                      % Proportion of treatment initiations amongst notifying providers in 2019     
    
    ART_covg = sum(sfin(s.hart)/sum(sfin([s.h1,s.hart])));                 % ART coverage in 2019
    HIV_prev = sum(sfin([s.h1, s.hart]))/sum(sfin(1:i.nstates));           % HIV prevalence in 2019
    mort_H0  = sdiff(end,i.aux.mort(1))*1e5;                               % HIV-ve TB mortality in 2019
    mort_H1  = sdiff(end,i.aux.mort(2))*1e5;                               % HIV+ve TB mortality in 2019
    

    if inc_all_2019 < 1
        out = Inf*calfn.sgn;
        aux = [];
    else
        
        if opts.hiv
            out = calfn.fn(inc_all_2019, mort_H0, noti_pu, inc_h1, ART_covg, HIV_prev, mort_H1);
        elseif opts.provs
            out = calfn.fn(inc_all_2019, mort_H0, noti_pu);
        else
            out = calfn.fn(inc_all_2019, inc_all_2014, mort_H0, noti_pu);
        end
        if opts.provs 
            out = out + calfn.fn_pu(prop_pu);
        end

        % --- Package all relevant outputs into 'aux' as an output --------
        aux.soln         = [soln, t];
        aux.allsol       = [soln0; soln(2:end,:)];
        aux.inc_all_2019 = inc_all_2019;
        aux.inc_all_2014 = inc_all_2014;
        aux.inct         = sum(sdiff(:,i.aux.inc),2)*1e5;
        aux.inc_h1       = inc_h1;
        aux.noti         = noti;
        aux.noti_pu      = noti_pu;
        aux.prop_pu      = prop_pu;
        aux.ART_covg     = ART_covg;
        aux.HIV_prev     = HIV_prev;
        aux.mort_H0      = mort_H0;
        aux.mort_H1      = mort_H1;
        aux.prev         = sum(sfin(s.prevalent))*1e5;
        aux.allin        = allin;
        aux.M0           = M0;
    end        
    
    if opts.provs || opts.hiv
       aux.M1 = M1;
       aux.M2 = M2;
    else
        aux.M1 = M1;
    end
    
    
end
