% Calculate time derivatives associated with a linear change in model 
% parameters. This is essentially a wrapper function for 
% goveqs_basis_disruption, allowing for a linear change in model parameters. 
% Same as goveqs_scaleup_disruption, but allowing for linear change over 
% _two_ independent time intervals.

% Arguments:
% ==========
% - t:  Time 
% - in: State vector at time t
% - M0: Matrix-based model specification, prior to any linear change 
% (created by make_model)
% - M1: Matrix-based model specification, at the end of the first time 
% interval (created by make_model)
% - M2: Matrix-based model specification, at the end of the second time 
% interval (created by make_model)
% - times: Start and end times for the linear change in model parameters.
% Consists of two rows, with each row defining a separate time interval
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


function out = goveqs_scaleup_disruption2D(t, in, M0, M1, M2, times, notif_fac, i, prm, sel, agg)

scale = max((t-times(:,1))./(times(:,2)-times(:,1)),0);                   
scale(1) = min(scale(1),1);

Mt = M1; 

Mt.lin   = M0.lin   + scale(1)*(M1.lin-M0.lin)     + scale(2)*(M2.lin-M0.lin);
Mt.Dxlin = M0.Dxlin + scale(1)*(M1.Dxlin-M0.Dxlin) + scale(2)*(M2.Dxlin-M0.Dxlin);
out = goveqs_basis_disruption(t, in, Mt, notif_fac, i, prm, sel, agg);