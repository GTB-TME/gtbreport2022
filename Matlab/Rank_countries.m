% Script to identify the countries that contributed the most to the global
% drop in notifications between Jan 2020 and Dec 2021

clear all; 

% -------------------------------------------------------------------------
% --- Get the disruption data, and organise -------------------------------

load('Data/Disruptions/disruption_data.mat');
load('Data/TB Notifications/notif_data_1819.mat');


% --- Get the pre-2019 WHO data -------------------------------------------

iso3s = notifs_19.iso3;
% Align countries with WHO disruption data
for ii = 1:length(iso3_disrp)
   row  = find(strcmp(iso3s, iso3_disrp{ii}));
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

% tmp1      = dat(:,2).^2./dat(:,1);                                       % Extrapolating assuming exponential decay
% tmp2      = tmp1.^2./dat(:,2);
% tmp1      = 2*dat2(:,2)-dat2(:,1);                                       % Extrapolating assuming linear trend

tmp1      = dat2(:,2);                                                     % Extrapolating assuming constant year trend
tmp2      = 2*tmp1 -dat2(:,2);
notif_est = [tmp1, tmp2];


% --- Use WHO notifications to pad incomplete data in monthly/quarterly data
mdata_aug = mdata;
qdata_aug = qdata;
cinc      = nan(length(iso3_disrp),2);
for ii = 1:length(iso3_disrp)
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

% --- Finally, compare expected against padded reported -------------------

exp     = sum(notif_est,2);
rep     = sum(cinc,2);
drop    = max(exp-rep,0);
pctdrop = drop/sum(drop);

% Calculate in order of magnitude
mat           = sortrows([pctdrop,[1:length(pctdrop)]'],-1);
nctrs         = 32;
ctrs_priority = iso3_disrp(mat(1:nctrs,2));
props         = log10(mat(1:nctrs,1)*100+1);

figure; hold on;
b1 = bar(props);
set(gca,'XTick',1:nctrs,'XTickLabel',ctrs_priority);
xtickangle(90);

% Mark existing countries
exis = {'IND','IDN','PHL','BGD','ZAF','CHN','RUS','MMR','AGO','PAK','VNM','UGA','UKR','PER','BRA','TZA'};
exis = {'AGO','BGD','BRA','CHN','IDN','IND','KEN','MMR','PAK','PER','PHL','RUS','UGA','UKR','VNM','ZAF'};

iex  = nan(1,length(exis));
for ii  = 1:length(exis)
    ind = find(strcmp(ctrs_priority,exis{ii}));
    if ~isempty(ind)
        iex(ii) = ind; 
    end
end
iex = iex(~isnan(iex));
b2 = bar(iex,props(iex));
b2.BarWidth = 0.8;
set(gca,'fontsize',14);
ylabel('log(Percent contribution to global shortfall + 1)');

save country_rankings_global;