% Given calibrations to pre-COVID data, this script simulates the effect of
% disruptions during and after COVID (i.e. from 2020 onwards)

% For every country, the vector vec0 shows the monthly timeseries (or 
% quarterly, depending on what level the notification data is available for) 
% for the time-dependent multiplier on the rate of diagnosis and treatment
% initiation. 

% As an output, the code shows model-simulated treatment initiations, and
% plots these against notification data. The vector vec0 is adjusted so
% that simulations match data as closely as possible.

% Dependencies:
% =============
% - alloc_parameters.m
% - goveqs_basis2.m
% - goveqs_basis_disruption.m

% This code is used by:
% =====================
% - Outputs saved in 'projections_raw' serve as inputs for
% Extrapolate_India.m

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


clear all; load model_fits2; 

popn = 1.38e9;

ix0 = size(xsto,1)/2; nx = 50; dx = round((size(xsto,1)-ix0)/nx);
xs = xsto(ix0:dx:end,:,1);

% --- Get the India notifications
[num, txt, raw] = xlsread('Nikshay.xlsx');
num([1:2],:) = [];
mat = [[num(:,1); num(:,3); num(:,5); num(:,7)], [num(:,2); num(:,4); num(:,6); num(:,8)]];
for ii = 1:2
    vec = mat(:,ii);
    notifs(:,ii) = vec(~isnan(vec))/popn*1e5;
end

tmp = notifs;
tmp(17:27,1)  = smooth(tmp(17:27,1));
tmp(17:27,2)  = smooth(tmp(17:27,2));
tmp(30:end,1) = smooth(tmp(30:end,1));
tmp(30:end,2) = smooth(tmp(30:end,2));
notifs_smth   = tmp;


% --- Lockdown dates 
ldown   = 2020 + [4 7]/12;                                                 % Start and end dates of lockdown
betared = 0.5;                                                             % Reduction in beta during lockdown
tend    = 2025;                                                            % End date for simulation

% % --- Disruption parameters
% vecpu = [1 1 0.63 0.25 0.75 0.5  0.5  0.5  0.5  0.6  0.6  0.6  0.8  0.8  0.8 0.35 0.1 1   0.5 1 0.5 0.85 0.5 0.8 0.4 1 0.7 0.9  0.7  0.75];
% vecpr = [1 1 0.58 0.15 0.7  0.66 0.66 0.65 0.7  0.7  0.75 0.85 0.9  0.9  0.9 0.4  0.5 1.2 0.7 1 0.8 0.8  0.6 0.8 0.5 1 0.7 0.95 0.65 0.7];

vecpu = [1 1 0.63 0.25 1.0 0.4 0.8 0.50 0.70 0.60 0.75 0.75 0.85 0.90 1.0  0.35 0.35 1.0 0.75 1.0 0.80 1.0  0.80 0.95 1.0 1.0  1.1  1.15 1.15 1.1];
vecpr = [1 1 0.58 0.15 1.0 0.4 1.0 0.55 0.95 0.65 1.0  0.80 1.0  1.0  1.15 0.6  0.5  1.2 0.9  1.1 1.0  0.90 1.0  0.90 1.0 1.05 1.05 1.2  1.0  1.1];

fac   = ones(2,(tend-2019)*12); 
fac(1,12+[1:length(vecpu)]) = vecpu;
fac(2,12+[1:length(vecpr)]) = vecpr;
fac   = fac';                                                              % Need to transpose because when there are two columns, interp1 (in goveqs_disruption) only works with column vectors

% --- Do the simulations 
opts = odeset('NonNegative',[1:i.nstates],'Refine',64,'AbsTol',1e-10,'RelTol',1e-10);

r0 = r; p0 = p;
inct = [];
noti = [];

mk = round(size(xs,1)/25);
for ii = 1:size(xs,1)
    if mod(ii,mk) == 0; fprintf('%0.5g ',ii/mk); end 
    
    [r,p] = alloc_parameters(xs(ii,:),prm.r,prm.p,xi);
    
    % Get the initial conditions
    [out, aux] = obj(xs(ii,:));
    init  = aux.soln(end-1,1:end-1);
    dsol0 = diff(aux.soln(:,1:end-1),1);
    dsol0(end,:) = [];
    
    M0 = make_model(p, r, i, s, gps);
    
    % Simulate in absence of disruption
    geq = @(t,in) goveqs_basis2(t, in, M0, i, s, p, sel, agg);
    [t0, soln0] = ode15s(geq, [2019:1/12:tend], init, opts);
    inct(:,ii,1)   = diff(soln0(:,i.aux.inc),1);
    mort(:,ii,1)   = diff(soln0(:,i.aux.mort),1);
    noti(:,:,ii,1) = diff(soln0(:,i.aux.noti),1);
    notq(:,:,ii,1) = diff(soln0(1:3:end,i.aux.noti),1);
    
    % Also aggregate to get overall annual incidence, combined with 1997
    % onwards
    tmp = diff(soln0(1:12:end,i.aux.inc),1);
    inctall(:,ii,1) = [dsol0(:,i.aux.inc); tmp];

    tmp = diff(soln0(1:12:end,i.aux.mort),1);
    mortall(:,ii,1) = [dsol0(:,i.aux.mort); tmp];
    
    % Now simulate disruption
    r1 = r; p1 = p;
    r1.beta = r.beta*(1-betared);
    M1 = make_model(p1, r1, i, s, gps);
    
    geq = @(t,in) goveqs_basis_disruption(t, in, M0, fac, i, sel, agg);
    [ta, solna] = ode15s(geq, [2019 ldown(1)], init, opts);
    
    initb = solna(end,:);
    geq = @(t,in) goveqs_basis_disruption(t, in, M1, fac, i, sel, agg);
    [tb, solnb] = ode15s(geq, [ldown(1) ldown(2)], initb, opts);
    
    initc = solnb(end,:);
    geq = @(t,in) goveqs_basis_disruption(t, in, M0, fac, i, sel, agg);
    [tc, solnc] = ode15s(geq, [ldown(2) tend], initc, opts);
    
    soln  = [solna; solnb(2:end,:); solnc(2:end,:)];
    t     = [ta; tb(2:end); tc(2:end)];
    soln1 = interp1(t, soln, t0);

    inct(:,ii,2)   = diff(soln1(:,i.aux.inc),1);
    mort(:,ii,2)   = diff(soln1(:,i.aux.mort),1);
    noti(:,:,ii,2) = diff(soln1(:,i.aux.noti),1);
    notq(:,:,ii,2) = diff(soln1(1:3:end,i.aux.noti),1);
    
    % Also aggregate to get overall annual incidence, combined with 1997
    % onwards
    tmp = diff(soln1(1:12:end,i.aux.inc),1);
    inctall(:,ii,2) = [dsol0(:,i.aux.inc); tmp];

    tmp = diff(soln1(1:12:end,i.aux.mort),1);
    mortall(:,ii,2) = [dsol0(:,i.aux.mort); tmp];

%     inctall(:,ii,2) = [dsol0(:,i.aux.inc);  inct(:,ii,2)];
%     mortall(:,ii,2) = [dsol0(:,i.aux.mort); mort(:,ii,2)];
end
fprintf('\n');

inc_pct  = permute(prctile(inct,[2.5,50,97.5],2)*1e5,[2,1,3]);
mrt_pct  = permute(prctile(mort,[2.5,50,97.5],2)*1e5,[2,1,3]);
noti_pct = permute(prctile(noti,[2.5,50,97.5],3)*1e5,[3,1,2,4]);           % Dims: 1.Lo/Md/Hi, 2.Month, 3.Pu/Pr, 4.Scenario
notq_pct = permute(prctile(notq,[2.5,50,97.5],3)*1e5,[3,1,2,4]);           % Dims: 1.Lo/Md/Hi, 2.Month, 3.Pu/Pr, 4.Scenario

% save res_forward;
save projections_raw2;
save disruption_vector vecpu vecpr;

% -------------------------------------------------------------------------
% --- Draw figures comparing model simulations with data ------------------

xlbl = []; ct = 1;
mnths = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
years = {'19','20','21','22','23','24','25'};
for iy = 1:length(years)
   for im = 1:length(mnths)
       xlbl{ct} = [mnths{im},' ',years{iy}]; 
       ct = ct+1;
   end
end

ff = figure; fs = 14; lw = 2; ms = 20; ts = 0.1;
skip = 4;

ylbl = {'Monthly notifications per 100k','Monthly notifications, relative to pre-pandemic'};
tis  = {'Public notifications','Private notifications'};

% --- Show public notifications

subplot(1,2,1);
cols = linspecer(2);
jj   = 1;
mat  = squeeze(noti_pct(:,:,jj,:));
for ii = 1:2
    plt = mat(:,:,ii);
    pl(ii,:) = plot(plt(2,:),'Color',cols(ii,:),'linewidth',lw); hold on;
    jbfill(1:size(plt,2),plt(3,:),plt(1,:),cols(ii,:),'None',1,ts); hold on;
end
xlim([1 size(plt,2)]);
pl(3,:) = plot(notifs(:,1),'.-','Color','g','linewidth',lw,'markersize',ms);
yl = ylim; yl(1) = 0; ylim(yl);
xlim([1, size(notifs,1)]);
ylabel('Monthly notifications per 100k');
% legend(pl([3,1,2],:),'Data','Modelled baseline','Modelled disruption','location','SouthEast');
set(gca,'fontsize',fs,'XTick',1:skip:size(plt,2),'XTickLabel',xlbl(1:skip:size(plt,2)));
xtickangle(45);
title(tis{1});
line((12+length(vecpu))*[1 1], ylim, 'linestyle', '--');
ylim([0, max(notifs(:,1))*1.3]);

% --- Show private notifications

subplot(1,2,2);
jj   = 2;
tmp  = squeeze(noti_pct(:,:,jj,:));
mn1  = mean(tmp(2,1:12,1));
mn2  = mean(notifs(1:12,2));
mat  = tmp*mn2/mn1;
for ii = 1:2
    plt = mat(:,:,ii);
    pl2(ii,:) = plot(plt(2,:),'Color',cols(ii,:),'linewidth',lw); hold on;
    jbfill(1:size(plt,2),plt(3,:),plt(1,:),cols(ii,:),'None',1,ts); hold on;
end
% xlim([1 size(plt,2)]);
xlim([1, size(notifs,1)]);
pl2(3,:) = plot(notifs(:,2),'.-','Color','g','linewidth',lw,'markersize',ms);
ylim([0, max(notifs(:,2))*1.3]);

% ylim([0 1.5]);
% ylabel('Monthly notifications per 100k');
legend(pl2([3,1,2],:),'Data','Modelled baseline','Modelled disruption','location','SouthEast','autoupdate','Off');
set(gca,'fontsize',fs,'XTick',1:skip:size(plt,2),'XTickLabel',xlbl(1:skip:size(plt,2)));
xtickangle(45);
title(tis{2});
line((12+length(vecpr))*[1 1], ylim, 'linestyle', '--');

set(ff,'Position',[471   644   940   329]);
set(ff,'Position',[361         404        1000         393]);

t = 1998:2025;
mat  = prctile(inctall,[2.5,50,97.5],2);
mat2 = squeeze(mat(:,2,:));
figure; plot(mat2(t>=2015,:));
yl = ylim; yl(1) = 0; ylim(yl);