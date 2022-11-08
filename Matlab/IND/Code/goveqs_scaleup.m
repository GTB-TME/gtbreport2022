% Calculate time derivatives associated with a linear change in model 
% parameters, over a single given time interval. This is essentially a
% wrapper function for goveqs_basis2, allowing for a linear
% change in model parameters. 

% Arguments:
% ==========
% - t:  Time 
% - in: State vector at time t
% - M0: Matrix-based model specification, at the beginning of the time 
% interval (created by make_model)
% - M1: Matrix-based model specification, at the end of the time 
% interval (created by make_model)
% - times: Start and end times for the linear change in model parameters
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
% - goveqs_basis2.m

% This function is used by:
% =========================
% - get_objective2D.m
% - Simulate_disruptions.m
% - simulate_disruptions_fn.m

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


function out = goveqs_scaleup(t, in, M0, M1, times, i, s, p, sel, agg)

scale = min((t-times(1))/(times(2)-times(1)),1);
Mt = M1; Mt.lin = M0.lin + scale*(M1.lin-M0.lin);
out = goveqs_basis2(t, in, Mt, i, s, p, sel, agg);
