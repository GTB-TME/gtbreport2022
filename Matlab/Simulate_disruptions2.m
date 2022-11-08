warning off;

looping_sd = 0;
if ~looping_sd
    clear all; ctry = 'IND';
end

looping_sd = 0;
load([ctry,'/model_fits']);


loading = 1;
tend    = 2025;                                                            % End date for simulation

ix0 = round(size(xsto,1)/2); nx = 50; dx = round((size(xsto,1)-ix0)/nx);
xs = xsto(ix0:dx:end,:,1);

% --- Get the disruption notifications
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
        tmp = notif_rate;
        tmp(6:end) = smooth(tmp(6:end));
        tmp(isnan(notif_rate)) = nan;
        notif_smth = tmp;
    end
end

% --- Lockdown dates
betared = 0.25 + 0.5*rand(1,size(xs,1));                                   % Reduction in beta during lockdown
% betared = ones(size(betared));





vec0 = [1 1 1 0.6];
notif_dat = notif_smth;



% --- Do the simulations

% out = get_disruption_objective(xs, fac, notif_smth, betared, prm, ref, agg, sel, tend, lhd, gps, opts);

if monthly
    vec0 = [1 1 1];
else
    vec0 = [1];
end
tmp = fillmissing(notif_dat,'next');
nsteps = length(tmp) - length(vec0);

for inst = 1:nsteps
    fprintf('Step %0.5g \n',inst);
    vec0 = [vec0, vec0(end)];
    hi = 1.5; lo = 0; prop = 1; pdif = 1;
    while abs(pdif)>0.01
        
        vec0(end) = prop;
        if ~monthly
            tmp = repmat(vec0,3,1);
            vec = tmp(:)';
        else
            vec = vec0;
        end
        fac = ones(1,(tend-2019)*12); fac(12+[1:length(vec)]) = vec;
        [inct, mort, noti, notq] = simulate_disruptions_fn(xs, fac, betared, prm, ref, agg, sel, tend, lhd, gps, opts);
        
        if monthly
            tmp = permute(prctile(noti,[2.5,50,97.5],3)*1e5,[3,1,2,4]);       % Dims: 1.Lo/Md/Hi, 2.Month, 3.Pu/Pr, 4.Scenario
            noti_sim = squeeze(tmp(2,:,:,2));
            dt = 12;
        else
            tmp = permute(prctile(notq,[2.5,50,97.5],3)*1e5,[3,1,2,4]);       % Dims: 1.Lo/Md/Hi, 2.Month, 3.Pu/Pr, 4.Scenario
            noti_sim = squeeze(tmp(2,:,:,2));
            dt = 4;
        end
        
        % Pull out numbers to be compared
        noti_sim(1:dt) = [];
        sim = noti_sim(length(vec0));
        dat = notif_dat(length(vec0));
        out = (1 - sim/dat)^2*100;
        pdif = (1 - sim/dat);
        
        if sim>dat
            hi = prop;
        else
            lo = prop;
        end
        prop = (hi+lo)/2;
    end
    fprintf('\n');
end


return;



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
    
    % notif_dat = notif_rate;
    notif_dat = notif_smth;
    
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
    title([ctry, ' Notifications']);
    legend(pl2([3,1,2],:),'Data','Modelled baseline','Modelled disruption','location','SouthWest','autoupdate','off');
    set(gca,'fontsize',fs);
    line((1/dt+length(vec0))*[1 1], ylim, 'linestyle', '--');
end
