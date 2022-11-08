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
% - goveqs_scaleup2D.m
% - goveqs_scaleup.m
% - goveqs_scaleup_disruption2D.m
% - goveqs_scaleup_disruption.m

% This code is used by:
% =====================
% - Outputs saved in 'projections_raw' serve as inputs for
% Extrapolate_Countries.m

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


warning off;

% Set looping_sd = 0 if this script is not being run by a Batch code, and 
% instead is being used only for a single country
looping_sd = 0;
if ~looping_sd
    clear all; ctry = 'MNG';
end
looping_sd = 0;

% Load the country calibrations
load([ctry,'/model_fits']);

% Whether to load previously stored values for vec0, or to use values
% specified below
loading = 1;

% Take a thinned sample of the posterior distribution
ix0 = round(size(xsto,1)/2); nx = 50; dx = round((size(xsto,1)-ix0)/nx);
xs = xsto(ix0:dx:end,:,1);

% --- Read in the monthly or quarterly notifications from 2020 onwards
if contains(ctry,'region')
    load regional_data;
    reg = erase(ctry,'region_');
    notif_rate = regdata.(reg).disruption_notifs;
    notif_smth = notif_rate;
    monthly = 0;
    dt = 1/4;
else
    load('Data/Disruptions/disruption_data.mat');
    ico = find(strcmp(iso3_disrp,ctry));
    monthly = (freq(ico,1)==70);                                           % Assuming here that if 2020 is monthly, then so is 2021, and similarly for quarterly
    
    if ~looping_sd
        freq(ico,:)
    end
    if monthly
        dat = mdata(ico,:,:);
        dt  = 1/12;                                                        
    else
        dat = qdata(ico,:,:);
        dt  = 1/4;
    end
    notif_rate = dat(:)'/popn*1e5;
    notif_smth = notif_rate;
    
    if monthly
        % Create a smoothed timeseries to facilitate matching with
        % simulated values
        tmp = notif_rate;
        tmp(6:end) = smooth(tmp(6:end));
        tmp(isnan(notif_rate)) = nan;
        notif_smth = tmp;
    end
end

tend    = 2025;                                                            % End year for the simulation
ldown   = 2020 + [4 6.5]/12;                                               % Start and end dates of lockdown (for period of reduced TB transmission)
betared = 0.25 + 0.5*rand(1,size(xs,1));                                   % Reduction in beta during lockdown


% --- Disruption parameters
if loading
    load([ctry,'/disruption_vector']);
else
    %     vec0 = [1 0.95 0.9 0.8 1.2 0.8 0.7 0.7 0.65];                                                                                                 % AGO
    %     vec0 = [1 1 1 0.55 0.65 0.9 0.45 0.7 0.45 0.6 0.45 0.45 0.4 0.38 0.4 0.43 0.49 0.54 0.54 0.56 0.6 0.57 0.55 0.55];                            % AZE
    %     vec0 = [1 0.3 0.8 0.85 0.95 1 1.05 1.15];                                                                                                     % BGD
    %     vec0 = [1 0.1 1 0.7 0.9  0.8 0.8 0.8  0.7 0.7 0.55 0.6 0.45 0.55 0.6 0.6 0.55 0.55 0.5 0.5 0.4 0.43 0.35 0.4 0.3 0.4 0.35];                   % CHN
    %     vec0 = [1 1 1 0.5 0.8 0.5 0.6 0.5 0.8 0.6 0.6 0.5 0.3 0.7 0.6 0.5 0.6 0.4 0.7 0.55 0.85 0.5 0.72 0.4 0.6 0.6 0.8];                            % COL
    %     vec0 = [1 0.78 0.9 0.95 0.92 0.95 0.92 0.92 1.1];                                                                                             % ETH
    %     vec0 = [1 0.8 0.9 0.4 0.4 0.7 0.3 0.6 0.25 0.5 0.2 0.5 0.4 0.7 0.45 0.7 0.4 0.55 0.3 0.6 0.45 0.68 0.5 0.7 0.45 0.6 0.35 0.5 0.3];            % IDN
    %     vec0 = [0.7 0.7 0.9 0.4 0.6 0.55 0.55 0.55 0.55 0.5 0.4 0.3 0.45 0.6 0.6 0.5 0.5 0.45 0.45 0.55 0.55];                                        % KEN (smoothed)
    %     vec0 = [1 1 1 0.6 0.88 0.86 0.88 0.88 0.88 0.88 0.65 0.6 0.62 1.1 1 0.95 0.88 0.88 0.83 1.1 1.1];                                             % KEN
    %     vec0 = [1 1 1 0.7 0.35 0.55 0.3 0.55 0.45 0.6 0.5 0.45 0.45 0.45 0.5 0.45 0.45 0.35 0.35 0.3 0.32 0.32 0.34 0.3 0.3 0.33 0.37 0.3 0.33];      % KAZ (smoothed)
    %     vec0 = [1 1 1 0.4 0.6 0.5 0.45 0.39 0.37 0.3 0.3 0.25 0.4 0.4 0.6 0.4 0.4 0.27 0.3 0.25 0.3 0.2 0.25 0.18 0.23 0.18 0.23 0.23];               % KGZ (smoothed)
    %     vec0 = [1 1 0.7 0.7 1 1 1.1 1 1 1 0.9 0.9 0.85 0.89 0.68 0.58 0.5 0.5 0.4 0.45 0.5 0.55 0.65 0.65 0.6 0.65 0.65];                             % KHM
    %     vec0 = [1 0.45 0.48 0.7 0.48 0.53 0.58 0.58 0.58];                                                                                            % LSO
    %     vec0 = [1 1 1 0.4 0.65 0.4 0.65 0.45 0.65 0.6 0.35 0.5 0.4 0.55 0.8 0.55 0.7 0.7 0.55 0.6 0.55 0.55 0.55 0.4 0.55 0.45 0.2];                  % MEX
    %     vec0 = [1 0.6 0.5 0.4 0.25 0.4 0.15 0.3];                                                                                                     % MMR
    %     vec0 = [1 1 1 0.3 1 0.5 1 1 0.5 0.7 1 1 0.4 0.6 0.7 0.5 0.7 0.8 0.2 0.9 0.3 1 0.4 0.8 0.3 0.7 0.3 0.5];                                       % MYS
    %     vec0 = [1 1 1 0.5 0.8 0.8 1 0.8 0.92 0.7 0.73 0.73 0.8 1 1.1 1.1 0.8 0.8 1 1.2 1 1 1];                                                        % NPL
    %     vec0 = [1 0.5 0.9 0.75 0.85 0.9 1.1 1.25 1.45];                                                                                               % PAK
    %     vec0 = [1 1 0.8 0.5 0.3 0.6 0.4 0.6 0.6 0.6 0.6 0.57 0.5 0.5 0.5 0.48 0.46 0.46 0.5 0.5 0.5 0.5 0.45 0.42 0.4 0.45 0.43 0.3 0.3];             % PER
    %     vec0 = [1 1 0.8 0.1 0.6 0.45 0.45 0.4 0.4 0.35 0.35 0.4 0.4 0.5 0.4 0.4 0.35 0.45 0.35 0.43 0.3 0.35 0.25 0.35 0.35 0.45 0.3 0.35];           % PHL
    %     vec0 = [1 1 0.9 0.88 0.86 0.76 0.88 0.33 0.6];                                                                                                % PNG
    %     vec0 = [1 1.3 1.2 1 0.75 0.95 1.1 0.9];                                                                                                       % PRK
    %     vec0 = [1 0.4 0.5 0.3 0.4 0.4 0.42 0.23 0.33];                                                                                                % ROU
    %     vec0 = [1 1 1 1 0.3 0.8 0.5 0.75 0.5 0.75 0.45 0.63 0.52 0.62 0.75 0.68 0.6 0.6 0.57 0.5 0.5 0.55 0.5 0.45 0.45 0.45 0.45];                   % RUS
    %     vec0 = [ones(1,16), 0.8, 0.75 0.7 0.61 0.7 0.67 0.63 0.63 0.63 0.63 0.5 0.5];                                                                 % THA
    %     vec0 = [1 0.55 0.8 0.8 0.77 0.45 0.65 0.85 1.1];                                                                                              % TLS
    %     vec0 = [1 1 1 0.4 0.4 0.9 0.4 0.6 0.4 0.5 0.3 0.5 0.2 0.5 0.4 0.4 0.3 0.4 0.25 0.4 0.2 0.35 0.2 0.35 0.2 0.3 0.25];                           % TUR
    %     vec0 = [1 0.8 1.1 1.1 1.4 1.6 1.8 2];                                                                                                         % UGA
    %     vec0 = [1 1 1 0.5 1 1 0.95 0.9 0.9 0.95 0.94 0.95 0.8 0.8 0.8 0.9 0.6 0.6 0.3 0.3 0.35 0.35 0.4 0.35 0.38 0.4 0.4 0.5];                       % VNM (smoothed)
    %     vec0 = [1 0.75 0.85 1 1.2 0.65 0.8 1];                                                                                                        % ZAF
    %     vec0 = [1 0.65 0.63 0.63 0.63 0.7];                                                                                                           % ZWE
    %     vec0 = [0.9 0.75 0.75 0.77 0.84 0.79 0.82 0.82 0.92];                                                                                         % region AFR
    %     vec0 = [1 0.55 0.65 0.71 0.66 0.59 0.74 0.60 0.62];                                                                                           % region AMR
    %     vec0 = [1 0.67 0.70 0.55 0.58 0.63 0.56 0.45 0.5];                                                                                            % region EUR
    
end

% vec0 starts from January 2020 (or Q1 2020, in the case of quarterly data) and can
% be any length. Here, pad it with ones to convert into a vector ('fac')
% that gives the monthly timeseries from 2019 to 'tend'.
if ~monthly
    tmp = repmat(vec0,3,1);
    vec = tmp(:)';
else
    vec = vec0;
end
fac = ones(1,(tend-2019)*12); fac(12+[1:length(vec)]) = vec;



% -------------------------------------------------------------------------
% --- Do the simulations --------------------------------------------------

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

inc_pct  = permute(prctile(inct,[2.5,50,97.5],2)*1e5,[2,1,3]);
mrt_pct  = permute(prctile(mort,[2.5,50,97.5],2)*1e5,[2,1,3]);
noti_pct = permute(prctile(noti,[2.5,50,97.5],3)*1e5,[3,1,2,4]);           % Dims: 1.Lo/Md/Hi, 2.Month, 3.Pu/Pr, 4.Scenario
notq_pct = permute(prctile(notq,[2.5,50,97.5],3)*1e5,[3,1,2,4]);           % Dims: 1.Lo/Md/Hi, 2.Month, 3.Pu/Pr, 4.Scenario

if ~loading
    save([ctry,'/disruption_vector'],'vec0');
    save([ctry,'/projections_raw']);
end



if ~looping_sd
    % ---------------------------------------------------------------------
    % --- Plot the fit with notifications ---------------------------------
    
    notif_dat = notif_rate;
%     notif_dat = notif_smth;
    
    figure; hold on;
    lw = 1.5; fs = 14; ts = 0.1; ms = 24;
    cols = linspecer(2);
    
    if monthly
        mat = squeeze(noti_pct(:,:,1,:));
        dt = 1/12;
        ylbl = 'Monthly notifications per 100k';
    else
        mat = squeeze(notq_pct(:,:,1,:));
        dt = 1/4;
        ylbl = 'Quarterly notifications per 100k';
    end
    
    for ii = 1:2
        plt = mat(:,:,ii);
        pl2(ii,:) = plot(plt(2,:),'Color',cols(ii,:),'linewidth',lw); hold on;
        jbfill(1:size(plt,2),plt(3,:),plt(1,:),cols(ii,:),'None',1,ts); hold on;
    end
    pl2(3,:) = plot(1/dt + [1:length(notif_dat)], notif_dat,'.-','Color','g','linewidth',lw,'markersize',ms);
    ylabel(ylbl);
    skip = 2; xinds = 1:skip:1/dt+length(notif_dat);
    
    xhi = find(~isnan(notif_dat),1,'last');
    xlim([1, 1/dt+xhi]);
    ym = [0,max(notif_dat)*1.2]; ylim(ym);
    %title([ctry, ' Notifications']);
    legend(pl2([3,1,2],:),'Data','Modelled baseline','Modelled disruption','location','SouthWest','autoupdate','off');
    
    %line((1/dt+length(vec0))*[1 1], ylim, 'linestyle', '--');
    %xlabel('Months since Jan 2019','fontsize',14);
%     mos = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
%     yrs = [2019:2024];
%     xlbl = {}; ct = 1;
%     for iy = 1:length(yrs)
%         for im = 1:length(mos)
%             xlbl{ct} = [mos{im}, ' ', num2str(yrs(iy))];
%             ct = ct+1;
%         end
%     end
%     skip = 3;
%     set(gca,'fontsize',fs,'XTick',1:skip:length(xlbl),'XTIckLabel',xlbl(1:skip:end));
    
end
