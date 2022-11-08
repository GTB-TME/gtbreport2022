% Calculate time derivatives associated with a linear change in model 
% parameters, over a single given time interval. This is essentially a
% wrapper function for goveqs_basis_disruption, allowing for a linear
% change in model parameters. Same as goveqs_scaleup, but allowing for 
% COVID-related disruptions.

% Arguments:
% ==========
% - t:  Time 
% - in: State vector at time t
% - M0: Matrix-based model specification, at the beginning of the time 
% interval (created by make_model)
% - M1: Matrix-based model specification, at the end of the time 
% interval (created by make_model)
% - times: Start and end times for the linear change in model parameters
% - notif_fac: Time-dependent adjustments to the rate of diagnosis,
% during the period of COVID disruptions
% - i:  Compartment lookup giving index of each model compartment (see 'Read me'
% for more detail on compartment lookups)
% - prm: All model parameters
% - sel, agg: Respectively, selector and aggregator matrices used for 
% counting incidence-like terms (see Setup_model for explanation)

% Outputs:
% ========
% - out: Vector of time derivatives

% Dependencies:
% =============
% - goveqs_basis_disruption.m

% This function is used by:
% =========================
% - Simulate_disruptions.m
% - simulate_disruptions_fn.m

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


function out = goveqs_scaleup_disruption(t, in, M0, M1, times, notif_fac, i, prm, sel, agg)

scale = max((t-times(1))/(times(2)-times(1)),0);                           
Mt = M1; Mt.lin = M0.lin + scale*(M1.lin-M0.lin);
out = goveqs_basis_disruption(t, in, Mt, notif_fac, i, prm, sel, agg);