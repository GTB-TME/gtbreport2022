% Function to calculate time derivatives, to be used in the ODE solver.
% Same as goveqs_basis2, but allowing for COVID-related disruptions.

% Arguments: 
% ==========
% - t:  Time 
% - in: State vector at time t
% - M:  Matrix-based model specification (created by make_model)
% - notif_fac: Time-dependent adjustments to the rate of diagnosis,
% during the period of COVID disruptions
% - i: Updated compartment lookup giving index of each model compartment (see 'Read me'
% for more detail on compartment lookups)
% - prm: All model parameters
% - sel, agg: Respectively, selector and aggregator matrices used for 
% counting incidence-like terms (see Setup_model for explanation)

% Outputs:
% ========
% - out: vector of time derivatives
% - lam: Instantaneous force of infection at time t

% Dependencies:
% =============
% None

% This function is used by:
% =========================
% - goveqs_scaleup_disruption.m
% - goveqs_scaleup_disruption2D.m

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


function [out, lam] = goveqs_basis_disruption(t, in, M, notif_fac, i, prm, sel, agg)  

invec = in(1:i.nstates);

rHIV = interp1(0:length(prm.rHIV)-1, prm.rHIV, (t-1980));                  

% Normalise by populations
lam = M.lambda*invec/sum(invec);
fac = interp1(2019+(0:length(notif_fac))/12, [1, notif_fac], t);
if isnan(fac)
    keyboard;
end
% Pull them together
allmat = M.lin + M.Dxlin*fac + rHIV*M.linHIV + lam*M.nlin;
out = allmat*invec;

% Implement deaths
morts = sum(M.mortvec,2).*invec;
out = out - morts;

% Implement births
births = sum(morts);
out(i.U.h0) = out(i.U.h0)+births;

% Get the auxiliaries
out(i.aux.inc)  = agg.inc*(sel.inc.*allmat)*invec;
out(i.aux.noti) = agg.noti*(sel.noti.*allmat)*invec;
out(i.aux.mort(1)) = sum(M.mortvec(:,2).*invec);
out(i.aux.mort(2)) = sum(M.mortvec(:,3).*invec);