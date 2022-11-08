% Script to identify the countries that experienced the largest relative
% reductions to their notifications, but were not included in the country
% list created by Rank_countries.m

clear all;

% -------------------------------------------------------------------------
% --- Get the disruption data, and organise -------------------------------

load('Data/Disruptions/disruption_data.mat');
load('Data/TB Notifications/notif_data_1819.mat');


% --- Get the pre-2019 WHO data -------------------------------------------

ctrs = notifs_19.iso3;
% Align countries with WHO disruption data
for ii = 1:length(iso3)
    row  = find(strcmp(ctrs,iso3{ii}));
    el1 = notifs_18{row,end};
    el2 = notifs_19{row,end};
    
    % Pad any nans with equal data
    if isnan(el1) && ~isnan(el2)
        el1 = el2;
    elseif isnan(el2) && ~isnan(el1)
        el2 = el1;
    end
    
    dat2(ii,1) = el1;
    dat2(ii,2) = el2;
    
end

% Extrapolate to get expected 2020 notifications
% tmp1      = dat(:,2).^2./dat(:,1);                                       % Extrapolating assuming exponential decay
% tmp2      = tmp1.^2./dat(:,2);

% tmp1      = 2*dat2(:,2)-dat2(:,1);                                       % Extrapolating assuming linear trend

tmp1      = dat2(:,2);                                                     % Extrapolating assuming constant year trend
tmp2      = 2*tmp1 -dat2(:,2);
notif_est = [tmp1, tmp2];


% --- Use WHO notifications to pad incomplete data in monthly/quarterly data
mdata_aug = mdata;
qdata_aug = qdata;
cinc      = nan(length(iso3),2);
for ii = 1:length(iso3)
    for iy = 1:2
        fr = freq(ii,iy);
        if fr == 70
            vec = mdata(ii,:,iy);
            vec(isnan(vec)) = notif_est(ii,iy)/12;
            mdata_aug(ii,:,iy) = vec;
            cinc(ii,iy) = max(sum(vec),0);
        else
            vec = qdata(ii,:,iy);
            vec(isnan(vec)) = notif_est(ii,iy)/4;
            qdata_aug(ii,:,iy) = vec;
            cinc(ii,iy) = max(sum(vec),0);
        end
    end
end


% -------------------------------------------------------------------------
% --- Find biggest contributors at global level ---------------------------

exp          = sum(notif_est,2);
rep          = sum(cinc,2);
drop         = max(exp-rep,0);
pctdrop_glob = drop/sum(drop);

% Calculate in order of magnitude
mat     = sortrows([pctdrop_glob,[1:length(pctdrop_glob)]'],-1);
nctrs   = 35;
ctrs    = iso3(mat(1:nctrs,2));
props   = log10(mat(1:nctrs,1)*100+1);

figure; hold on;
b1 = bar(props);
set(gca,'XTick',1:nctrs,'XTickLabel',ctrs);
xtickangle(90);
% b1.BarWidth = 0.1;

% Mark existing countries
exis = {'IND','IDN','PHL','BGD','ZAF','CHN','RUS','MMR','AGO','PAK','VNM','UGA','UKR','PER','BRA','TZA'};
exis = {'AGO','BGD','BRA','CHN','IDN','IND','KEN','MMR','PAK','PER','PHL','RUS','UGA','UKR','VNM','ZAF'};

iex  = nan(1,length(exis));
for ii  = 1:length(exis)
    ind = find(strcmp(ctrs,exis{ii}));
    if ~isempty(ind)
        iex(ii) = ind;
    end
end
iex = iex(~isnan(iex));
b2 = bar(iex,props(iex));
b2.BarWidth = 0.8;
set(gca,'fontsize',14);

ylabel('log(Percent contribution to global shortfall + 1)');


% -------------------------------------------------------------------------
% --- Find countries having >10% reduction in notifs at country level -----

tmp  = drop./exp;
inds = find(exp>0);
pctdrop_indiv = tmp(inds);

inds2 = find(pctdrop_indiv>0.1);
ctrs2 = iso3(inds2);

% Find countries not represented in global list
for ic = 1:length(ctrs2)
    ispres(ic) = ~isempty(find(strcmp(ctrs, ctrs2{ic})));
end
ctrs2_mat = ctrs2(find(1-ispres));

reg = {};
% --- Populate the regions for each of these countries --------------------
for ic = 1:length(ctrs2_mat)
    ctr = ctrs2_mat{ic};
    reg{ic} = notifs_18.g_whoregion{find(strcmp(notifs_18.iso3,ctr))};
end

lookup = [ctrs2_mat, reg'];
% Reverse, to get a list of countries in each region
regs = unique(reg);
for ir = 1:length(regs)
    ctrlist.(regs{ir}) = lookup(find(strcmp(lookup(:,2),regs{ir})),1);
end

% --- Get a lookup of each iso3 code to country ---------------------------
for ic = 1:length(ctrs2_mat)
    ctr = ctrs2_mat{ic};
    row = find(strcmp(notifs_18.iso3,ctr));
    iso2ctry.(ctr) = notifs_18.country{row};
end
% Make adjustments according to labelling in Thembisa model
iso2ctry.COD = 'Congo';
iso2ctry.GBR = 'United Kingdom';



for ir = 1:length(regs)
    reg = regs{ir}; list = ctrlist.(reg);
    
    % --- Estimated TB incidence and mortality ----------------------------
    load('Data/TB Estimates/estim_data_2b.mat'); mat2 = [];
    for ic = 1:length(list)
        iso3 = list{ic};
        countryrow  = estims(strcmp(estims.iso3,iso3),:);
        popn2(ic)    = countryrow.e_pop_num;
        mat2(:,:,ic) = reshape(countryrow{1,3:end},[3,6])';
    end
    
    % Do a population-weighted average
    tmp  = squeeze(repmat(popn2/sum(popn2),[1,1,size(mat2,1),size(mat2,2)]));
    tmp2 = permute(tmp,[2,3,1]);
    mat3 = sum(mat2.*tmp2,3);
    
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
        iso3 = list{ic};
        
        countryrow = notifs_new(strcmp(notifs_new.iso3,iso3),:);
        noti(ic)   = countryrow.c_newinc;
        
        countryrow = pTPT(strcmp(pTPT.iso3,iso3),:);
        covTPT(ic) = countryrow.pTPT;
    end
    
    % Population-weighted average: public notifications
    regdata.(reg).noti_pu = sum(noti)/sum(popn2)*1e5*[0.9 1 1.1];
    
    % Proportion of ART on TPT
    regdata.(reg).pcovTPT = sum(covTPT.*popn2)/sum(popn2);    % <---------- Allocate as prm.p.covTPT
    
    
    % ---------------------------------------------------------------------
    % --- HIV services ----------------------------------------------------
    load('Data/HIV data/HIV_estims.mat');
    incl         = zeros(1,length(list));
    rHIV_co      = zeros(length(list),size(HIV_incd,1)); 
    ART_covg_co  = zeros(length(list),3);
    ART_start_co = zeros(1,length(list));
    HIV_prev_co  = zeros(length(list),3);
        
    for ic = 1:length(list)
        iso3 = list{ic};
        country = iso2ctry.(iso3);
        
        ind = find(strcmp(countries1,country));
        
        if ~isempty(ind)
            incl(ic) = 1;
            rHIV_co(ic,:) = HIV_incd(:,2,ind)';
            
            ind = find(strcmp(countries2,country));
            tmp = ARTcovg_2019(ind,:); tmp(end) = min(tmp(end),0.95);
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
end

save regdata regdata



% --- Visualise notifs from any given country
ctry = ctrs2{3};
ctry = 'GBR'
ind  = find(strcmp(iso3,ctry));
fr   = freq(ind,:);
if fr(1) == 70
    dat = squeeze(mdata(ind,:,:));
    dat = dat(:);
else
    dat = squeeze(qdata(ind,:,:));
    dat = dat(:);
end
figure; plot(dat); yl = ylim; yl(1) = 0; ylim(yl);
