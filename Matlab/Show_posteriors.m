% Show posterior densities for parameter estimates created by
% Get_calibrations.m, for any given country 

clear all; 

iso3 = 'MNG';
load([iso3,'/model_fits_popupdated2']); 

lbls = {'\beta, HIV-ve','\beta, HIV +ve','d_{pu}','d_{pr}',{'Relative exposure,', 'HIV+ve vs -ve'},'Case-finding','\mu_{TB} FACTOR, HIV-ve','\mu_{TB} FACTOR, HIV+ve','Prop Tx completion','Rate ART init','Rate HIV mort',{'Relative progression', 'HIV+ve vs HIV-ve'},'Self cure FACTOR','Progression FACTOR',{'LTBI stabilisation','FACTOR'},'Reactivation FACTOR',{'Relapse (low','risk) FACTOR'},{'Relapse (high','risk) FACTOR'},{'Relapse (long','term) FACTOR'},{'Protection from','prior exposure'}};

ix0 = round(size(xsto2,1)/2);
nx  = 1e4; dx = round(ix0/nx);
xs  = xsto2(ix0:dx:end,:);

inds    = 1:20;
delinds = [];
if ~opts.hiv
    delinds = [delinds, find(ismember(inds,[xi.r_beta(2), xi.p_HIVlam, xi.rf_mort_TB(2), xi.r_ART_init, xi.r_HIV_mort, xi.p_HIV_relrate]))];
    %inds(find(ismember(inds,[xi.r_beta(2), xi.r_casefinding, xi.rf_mort_TB(2), xi.r_ART_init, xi.r_HIV_mort, xi.p_HIV_relrate]))) = [];
end
if ~opts.provs
    delinds = [delinds, find(ismember(inds,[xi.r_Tx_init(2)]))];
    %inds(find(ismember(inds,[xi.r_Tx_init(2), xi.r_casefinding]))) = [];
end
if opts.hiv || opts.provs
    delinds = [delinds, xi.r_casefinding];
end

inds(unique(delinds)) = [];

tmp = sqrt(length(inds));

nc = ceil(tmp);
if nc^2 - length(inds)>nc
    nr = nc-1;
else
    nr = nc;
end

ff=figure;
for ii = 1:length(inds)
   subplot(nr,nc,ii); 
   dat = xs(:,inds(ii));
   hs = histogram(dat,'Normalization','pdf'); 
   hold on;
   vec = prm.bounds(:,inds(ii));
   line(vec, 1/diff(vec)*[1 1], 'Color', 'r');
   xlim(vec);
   
   title(lbls{inds(ii)});   
end

if opts.provs
    subplot(nr,nc,ii+1); hold on;
    mat = xs(:,xi.r_Tx_init);
    rat = mat(:,1)./sum(mat,2);
    histogram(rat,'Normalization','pdf'); 
    xx = linspace(1e-5,1-1e-5,50);
    for ix = 1:length(xx)
        yy(ix) = exp(lhd.fn_pu(xx(ix)));
    end
    plot(xx,yy,'Color','r');    
    title('Proportion Tx public');
end

if opts.hiv
    ind = find(strcmp(lbls(inds),'Rate HIV mort'));
    subplot(nr,nc,ind);
    xlim([0 0.08]);
end

set(ff,'Position',[680   387   738   590]);
saveas(gcf,[iso3,'/Posterior_densities'],'png');
