% All countries that did not receive country-specific modelling, yet still
% experienced substantial (>10%) reductions in notifications, are
% aggregated into their respective WHO regions, and the effects of
% disruptions modelled at the Regional level.

% This script pulls together such countries into their respective regions,
% calculating the calibration data for each Region as a population-weighted,
% average across all countries in the Region, as well as the notification
% data during the period of disruption.

% V2: Version to draw population estimates (popn2) from the updated
% estimates sent by Philippe on 4 Aug 2022

clear all;

% -------------------------------------------------------------------------
% --- Get the disruption data, and organise -------------------------------

load('Data/Disruptions/disruption_data.mat');
load('Data/TB Notifications/notif_data_1819.mat');
load('Data/TB Notifications/lookups.mat');
load('Data/Populations/popn_data.mat');
load('country_rankings_global.mat');


% -------------------------------------------------------------------------
% --- Find countries having >10% reduction in notifs at country level -----

pctdrop_indiv = drop./exp;
inds2         = find(pctdrop_indiv>0.1);
tmp           = iso3_disrp(inds2);

% Find countries not represented in global list
for ic = 1:length(tmp)
    ispres(ic) = ~isempty(find(strcmp(ctrs_priority, tmp{ic})));
end
ctrs_regional = tmp(find(1-ispres));

% Now divide the countries into the relevant regions
allregs = {'AFR','AMR','EMR','EUR','SEA','WPR'};
for ir = 1:length(allregs)
    reg = allregs{ir};
    tmp = intersect(reg2iso3s.(reg),ctrs_regional);
    if ~isempty(tmp)
        ctrlist.(reg) = tmp;
    end
end
regs = fieldnames(ctrlist);


% -------------------------------------------------------------------------
% --- Pull together the data on a regional basis

for ir = 1:length(regs)
%     for ir = 1:1
    reg = regs{ir}; list = ctrlist.(reg);
    
    % --- Estimated TB incidence and mortality ----------------------------
    load('Data/TB Estimates/estim_data_2b.mat'); mat2 = [];
    for ic = 1:length(list)
        countryrow   = estims(strcmp(estims.iso3, list{ic}),:);
        mat2(:,:,ic) = reshape(countryrow{1,3:end},[3,6])';
    end
    
    % --- Population ------------------------------------------------------
    popn2 = [];
    for ic = 1:length(list)
        countryrow2   = popns(strcmp(popns.iso3, list{ic}),:);
        popn2(ic)    = countryrow2.e_pop_num;
    end

    % Do a population-weighted average
    if length(list)>1
        tmp  = squeeze(repmat(popn2/sum(popn2),[1,1,size(mat2,1),size(mat2,2)]));
        tmp2 = permute(tmp,[2,3,1]);
        mat3 = sum(mat2.*tmp2,3);
    else
        mat3 = mat2;
    end
    regdata.(reg).inc_2019  = mat3(1,:);
    regdata.(reg).inc_h1    = mat3(2,:);
    regdata.(reg).mort_H0   = mat3(3,:);
    regdata.(reg).mort_H1   = mat3(4,:);
    regdata.(reg).mort_all  = mat3(5,:);
    regdata.(reg).inc_2014  = mat3(6,:);
    
    
    % ---------------------------------------------------------------------
    % --- Notifications ---------------------------------------------------
    load('Data/TB Notifications/notif_data.mat'); noti = []; covTPT = [];
    for ic = 1:length(list)
        countryrow = notifs_new(strcmp(notifs_new.iso3, list{ic}),:);
        noti(ic)   = countryrow.c_newinc;
        
        countryrow = pTPT(strcmp(pTPT.iso3, list{ic}),:);
        covTPT(ic) = countryrow.pTPT;
    end
    
    % Population-weighted average: public notifications
    regdata.(reg).noti_pu = sum(noti)/sum(popn2)*1e5*[0.9 1 1.1];
    
    % Proportion of ART on TPT
    regdata.(reg).pcovTPT = sum(covTPT.*popn2)/sum(popn2);
    
    
    % ---------------------------------------------------------------------
    % --- HIV services ----------------------------------------------------
    load('Data/HIV data/HIV_estims.mat');
    incl         = zeros(1,length(list));
    rHIV_co      = zeros(length(list),size(HIV_incd,1));
    ART_covg_co  = zeros(length(list),3);
    ART_start_co = zeros(1,length(list));
    HIV_prev_co  = zeros(length(list),3);
    
    for ic = 1:length(list)
        country = iso2ctry.(list{ic});
        
        ind = find(strcmp(countries1,country));
        
        if ~isempty(ind)
            incl(ic) = 1;
            rHIV_co(ic,:) = HIV_incd(:,2,ind)';
            
            ind = find(strcmp(countries2,country));
            tmp = ARTcovg_2019(ind,:); tmp(end) = min(tmp(end),95);
            ART_covg_co(ic,:) = tmp/100;
            ART_start_co(ic)  = ART_start(ind);
            
            ind = find(strcmp(countries3,country));
            HIV_prev_co(ic,:) = HIVprev_2019(ind,:);
        end
    end
    
    popden = popn2; popden(incl==0) = 0;
    
    regdata.(reg).rHIV      = sum(rHIV_co.*incl'/sum(popden));
    regdata.(reg).ART_covg  = sum(ART_covg_co.*popden'/sum(popden),1);
    regdata.(reg).ART_start = sum(ART_start_co.*popden/sum(popden));
    regdata.(reg).HIV_prev  = sum(HIV_prev_co.*incl'/sum(popden));
    
    
    % ---------------------------------------------------------------------
    % --- Get aggregated disruption data ----------------------------------
    
    inds = [];
    for ic = 1:length(list)
        inds(ic) = find(strcmp(iso3_disrp,list{ic}));
    end
    num = sum(allqdata(inds,:),1,'omitnan');
    den = sum(popn2);
    regdata.(reg).disruption_notifs = num/den*1e5;
    
end

save regional_data regdata ctrlist






return;

% --- Notifs from a given region
reg = 'EUR';
figure; hold on;
plot(regdata.(reg).disruption_notifs);
plot(xlim, regdata.(reg).noti_pu(1)*[1 1]/4,'linestyle','--');
yl = ylim; yl(1) = 0; ylim(yl);


% --- Visualise notifs from any given country
ctry = 'DEU';
ind  = find(strcmp(iso3_disrp,ctry));
dat  = allqdata(ind,:); div = 4;

% fr   = freq(ind,:);
% if fr(1) == 70
%     dat = squeeze(mdata(ind,:,:));
%     dat = dat(:);
%     div = 12;
% else
%     dat = squeeze(qdata(ind,:,:));
%     dat = dat(:);
%     div = 4;
% end
figure; plot(dat); hold on;
ind = find(strcmp(notifs_19.iso3,ctry));
yy = notifs_19.allnoti(ind)/div;
line(xlim, yy*[1 1], 'linestyle','--');
yl = ylim; yl(1) = 0; ylim(yl);


% --- Grid plot of notifications vs disruptions for all countries in a given
% region

figure;
for ii = 1:length(list)
    subplot(5,4,ii); hold on;
    ynoti = noti(ii)/popn2(ii)*1e5/4;
    yy = allqdata(inds(ii),:)/popn2(ii)*1e5;
    plot(yy);
    plot(xlim, ynoti*[1 1], 'linestyle', '--');
    title(popn2(ii)/sum(popn2));
end
