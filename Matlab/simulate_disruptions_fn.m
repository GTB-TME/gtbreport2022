function [inct, mort, noti, notq] = simulate_disruptions_fn(xs, fac, betared, prm, ref, agg, sel, tend, lhd, gps, opts)

p  = prm.p; r = prm.r;
i  = ref.i; s = ref.s; xi = ref.xi;

ldown   = 2020 + [4 6.5]/12;                                               % Start and end dates of lockdown

obj     = @(x) get_objective2D(x, prm, ref, sel, agg, gps, lhd, opts);
odeopts = odeset('NonNegative',[1:i.nstates],'Refine',64,'AbsTol',1e-10,'RelTol',1e-10);

r0 = r; p0 = p;
inct = [];
noti = [];

mk = round(size(xs,1)/25);
for ii = 1:size(xs,1)
    if mod(ii,mk) == 0; fprintf('%0.5g ',ii/mk); end
    
    [r,p] = alloc_parameters(xs(ii,:),r0,p0,prm,xi,opts);
    
    % Get the initial conditions
    [out, aux] = obj(xs(ii,:));
    init = aux.soln(end,1:end-1);
    
    % --- Simulate in absence of disruption -------------------------------
    
    if opts.hiv || opts.provs
        M0    = aux.M0;
        M_pu  = aux.M1;
        M_ART = aux.M2;
        geq = @(t,in) goveqs_scaleup2D(t, in, M0, M_pu, M_ART, [2000 2009; prm.ART_start 2019], i, prm, sel, agg);
    else
        M0  = aux.M0;
        M1  = aux.M1;
        geq = @(t,in) goveqs_scaleup(t, in, M0, M1, [2014 2019], i, prm, sel, agg);
    end
    
    [t0, soln0] = ode15s(geq, [2019:1/12:tend], init, odeopts);
    inct(:,ii,1)   = sum(diff(soln0(:,i.aux.inc),1),2);
    mort(:,ii,1)   = sum(diff(soln0(:,i.aux.mort),1),2);
    noti(:,:,ii,1) = sum(diff(soln0(:,i.aux.noti(1)),1),2);
    notq(:,:,ii,1) = sum(diff(soln0(1:3:end,i.aux.noti(1)),1),2);
    
    inct2(:,:,ii,1) = diff(soln0(:,i.aux.inc),1);
    mort2(:,:,ii,1) = diff(soln0(:,i.aux.mort),1);
    noti2(:,:,ii,1) = diff(soln0(:,i.aux.noti(1)),1);
    notq2(:,:,ii,1) = diff(soln0(1:3:end,i.aux.noti(1)),1);
    
    % --- Now simulate disruption -----------------------------------------
    if opts.hiv || opts.provs
        Md0    = M0;     Md0.lambda    = M0.lambda*(1-betared(ii));
        Md_pu  = M_pu;   Md_pu.lambda  = M_pu.lambda*(1-betared(ii));
        Md_ART = M_ART;  Md_ART.lambda = M_ART.lambda*(1-betared(ii));
        geq  = @(t,in) goveqs_scaleup_disruption2D(t, in, M0,  M_pu,  M_ART,  [2000 2009; prm.ART_start 2019], fac, i, prm, sel, agg);
        dgeq = @(t,in) goveqs_scaleup_disruption2D(t, in, Md0, Md_pu, Md_ART, [2000 2009; prm.ART_start 2019], fac, i, prm, sel, agg);
    else
        Md0 = M0;        Md0.lambda = M0.lambda*(1-betared(ii));
        Md1 = M1;        Md1.lambda = M1.lambda*(1-betared(ii));
        geq  = @(t,in) goveqs_scaleup_disruption(t, in, M0,  M1,  [2014 2019], fac, i, prm, sel, agg);
        dgeq = @(t,in) goveqs_scaleup_disruption(t, in, Md0, Md1, [2014 2019], fac, i, prm, sel, agg);
    end
    
    
    % Pre-disruption
    [ta, solna] = ode15s(geq, [2019 ldown(1)], init, odeopts);
    
    % During lockdown
    initb = solna(end,:);
    [tb, solnb] = ode15s(dgeq, [ldown(1) ldown(2)], initb, odeopts);
    
    % After lockdown
    initc = solnb(end,:);
    [tc, solnc] = ode15s(geq, [ldown(2) tend], initc, odeopts);
    
    % Collate solutions
    soln  = [solna; solnb(2:end,:); solnc(2:end,:)];
    t     = [ta; tb(2:end); tc(2:end)];
    soln1 = interp1(t, soln, t0);
    
    inct(:,ii,2)   = sum(diff(soln1(:,i.aux.inc),1),2);
    mort(:,ii,2)   = sum(diff(soln1(:,i.aux.mort),1),2);
    noti(:,:,ii,2) = diff(soln1(:,i.aux.noti(1)),1);
    notq(:,:,ii,2) = diff(soln1(1:3:end,i.aux.noti(1)),1);
    
    inct2(:,:,ii,2) = diff(soln1(:,i.aux.inc),1);
    mort2(:,:,ii,2) = diff(soln1(:,i.aux.mort),1);
    noti2(:,:,ii,2) = diff(soln1(:,i.aux.noti(1)),1);
    notq2(:,:,ii,2) = diff(soln1(1:3:end,i.aux.noti(1)),1);
    
end
fprintf('\n');