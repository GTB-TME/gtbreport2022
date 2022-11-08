% For a given set of parameter values 'x', this code allocates all associated
% model parameters as rates (r) and proportions (p). 

% Arguments:
% ==========
% - x: Vector of parameter values being varied in MCMC
% - r: Values of all rates in the model
% - p: Values of all parameters in the model
% - prm: All remaining parameters
% - xi: records what each element of 'x' means (e.g. typing 'xi.r_beta' 
% will give 1, showing that the first element of 'x' denotes the parameter 
% beta, a rate).
% - opts: specification of whether the country being modelled is a high
% HIV-burden country, has a large private sector, or both

% Outputs:
% ========
% - r: Updated model rates arising from the parameter vector x
% - p: Updated model proportions arising from the parameter vector x
% - prm: Any other model parameters

% Dependencies:
% =============
% None

% This function is used by:
% =========================
% - get_objective2D.m
% - Simulate_disruptions.m
% - BRR_Get_outputs.m

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


function [r,p,prm] = alloc_parameters(x,r,p,prm,xi,opts)

r.beta         = x(xi.r_beta);
r.Tx_init      = x(xi.r_Tx_init);
p.HIVlam       = x(xi.p_HIVlam);
r.mort_TB      = r.mort_TB0*x(xi.rf_mort_TB);
r.default      = r.Tx*(1-x(xi.p_Tx_complete))./x(xi.p_Tx_complete);
r.ART_init     = x(xi.r_ART_init);
r.HIV_mort     = x(xi.r_HIV_mort);
r.casefinding  = x(xi.r_casefinding);
r.casefinding  = 0;

if opts.hiv || opts.provs
    r.casefinding = 0;                                                     % The parameter 'g' in the written equations, only used for non-HIV, non-public/private countries
end

% Construct vector for progression from 'fast' latent infection, each
% element representing a different HIV status: 1.HIV-ve, 2.HIV+ve, 3.On ART
tmp = r.progression0(1)*[1 1 1]*x(xi.rf_progression);
tmp([2,3]) = tmp([2,3])*x(xi.p_HIV_relrate).*[1 0.4-0.24*prm.p.covTPT];    
r.progression = tmp;

% Construct vector for progression from 'slow' latent infection, each
% element representing a different HIV status: 1.HIV-ve, 2.HIV+ve, 3.On ART
tmp = r.reactivation0(1)*[1 1 1]*x(xi.rf_reactivation);
tmp([2,3]) = tmp([2,3])*x(xi.p_HIV_relrate).*[1 0.4-0.24*prm.p.covTPT];
r.reactivation = tmp;

r.LTBI_stabil         = r.LTBI_stabil0*x(xi.rf_LTBI_stabil);
r.self_cure([1,3])    = r.self_cure0([1,3])*x(xi.rf_self_cure);
r.relapse             = r.relapse0.*x(xi.rf_relapse);
p.imm                 = x(xi.p_imm)*[1 0 1];

% For countries where HIV coinfection accounts for <10% of TB incidence,
% opts.hiv is set to zero. Parameters below serve to 'switch off' the HIV
% structure in the model
if ~opts.hiv
    r.beta(2)             = 0;
    p.HIVlam              = 0;
    r.ART_init            = 0;
    r.progression([2,3])  = 0;
    r.reactivation([2,3]) = 0;
    prm.rHIV              = 0*prm.rHIV;
end

% For countries without a large private sector, opts.provs is set to zero.
% The following parameter serves to 'switch off' the public/private
% structure by setting to zero the rate of treatment initiation in the
% non-notifying sector
if ~opts.provs
    r.Tx_init(2) = 0; 
end