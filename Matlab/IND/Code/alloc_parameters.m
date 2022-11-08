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


function [r,p] = alloc_parameters(x,r,p,xi)

r.beta         = x(xi.r_beta);
r.g            = x(xi.r_g);
r.Tx_init      = x(xi.r_Tx_init);
r.mort_TB      = r.mort_TB*x(xi.rf_mort_TB);
r.default      = r.Tx*(1-x(xi.p_Tx_complete))./x(xi.p_Tx_complete);
r.progression  = r.progression*x(xi.rf_progression);
r.LTBI_stabil  = r.LTBI_stabil*x(xi.rf_LTBI_stabil);
r.reactivation = r.reactivation*x(xi.rf_reactivation);
r.relapse      = r.relapse.*x(xi.rf_relapse);
p.imm          = x(xi.p_imm);
