% Plot figures showing agreement between model projections (created by
% Get_calibrations.m) and data, for given country

looping = 0;
if ~looping
    clear all; iso3 = 'MNG';
end
load([iso3,'/model_fits_popupdated']); 


prm.popn_turnover = 1;

obj = @(x) get_objective2D(x, prm, ref, sel, agg, gps, lhd, opts);

ix0 = round(size(xsto2,1)/2); nx = 200; dx = round((size(xsto2,1)-ix0)/nx);
xs = xsto2(ix0:dx:end,:,1);

mk = round(size(xs,1)/25);
for ii = 1:size(xs,1)
    if mod(ii,mk) == 0; fprintf('%0.5g ',ii/mk); end 
    [out, aux] = obj(xs(ii,:));
    sim(ii,:) = [aux.inc_all_2019, aux.inc_all_2014, aux.noti_pu, aux.inc_h1, aux.mort_H0, aux.mort_H1, aux.ART_covg, aux.HIV_prev];
    prv(ii)    = aux.prev;
    inct(:,ii) = aux.inct;
    ppu(ii)    = aux.prop_pu;
end
fprintf('\n');

sim_pct  = prctile(sim,[2.5,50,97.5],1);
sim_hilo = diff(sim_pct,[],1);


% -------------------------------------------------------------------------
% --- Plot the comparisons between model and data -------------------------

ff=figure; ms = 24; lw = 1.5; fs = 14;

dx = 0.1;
dat = [data.inc_2019; data.inc_2014; data.noti_pu; data.inc_h1; data.mort_H0; data.mort_H1; data.ART_covg; data.HIV_prev]';
dat_hilo = diff(dat,[],1);


allpts  = cat(3, dat(2,:), sim_pct(2,:));
allhilo = cat(3, dat_hilo, sim_hilo);
cols = {'r','b'};

if opts.hiv
    subplot(1,3,1); hold on;

    selinds = [1,4,3]; 
    for iz = 1:2
        xplt = [1:length(selinds)] + (-1)^iz*dx;
        
        mid  = allpts(:,selinds,iz);
        hilo = allhilo(:,selinds,iz);
        
        plot(xplt, mid, '.', 'markersize', ms, 'Color', cols{iz});
        errorbar(xplt, mid, hilo(1,:), hilo(2,:), 'linestyle', 'None', 'linewidth', lw, 'Color', cols{iz});
    end
    xlim([xplt(1)-0.5, xplt(end)+0.5]);
    yl = ylim; yl(1) = 0; ylim(yl);
    set(gca,'XTick',[1:3],'XTicklabel',{'Incidence','Incd HIV+ve','Notifications'},'fontsize',fs);
    ylabel('Rate per 100,000 population');
    
    subplot(1,3,2); hold on;
    selinds = [5,6];
    for iz = 1:2
        xplt = [1:length(selinds)] + (-1)^iz*dx;
        
        mid  = allpts(:,selinds,iz);
        hilo = allhilo(:,selinds,iz);
        
        plot(xplt, mid, '.', 'markersize', ms, 'Color', cols{iz});
        errorbar(xplt, mid, hilo(1,:), hilo(2,:), 'linestyle', 'None', 'linewidth', lw, 'Color', cols{iz});
    end
    xlim([min(xplt)-0.5, max(xplt)+0.5]);
    yl = ylim; yl(1) = 0; ylim(yl);
    set(gca,'XTick',[1,2],'XTicklabel',{'Mortality HIV+ve', 'Mortality HIV-ve'},'fontsize',fs);
    
    subplot(1,3,3); hold on;
    selinds = [7,8];
    for iz = 1:2
        xplt = [1:length(selinds)] + (-1)^iz*dx;
        
        mid  = allpts(:,selinds,iz);
        hilo = allhilo(:,selinds,iz);
        
        lg(iz,:) = plot(xplt, mid, '.', 'markersize', ms, 'Color', cols{iz});
        errorbar(xplt, mid, hilo(1,:), hilo(2,:), 'linestyle', 'None', 'linewidth', lw, 'Color', cols{iz});
    end
    xlim([min(xplt)-0.5, max(xplt)+0.5]);
    yl = ylim; yl(1) = 0; ylim(yl);
    set(gca,'XTick',[1,2],'XTicklabel',{'ART coverage', 'HIV prevalence'},'fontsize',fs);
    ylabel('Proportion');
    
    legend(lg,'Data','Model','Location','NorthEast');
    
elseif opts.provs
    subplot(1,2,1); hold on;

    selinds = [1,3]; cols = {'r','b'};
    for iz = 1:2
        xplt = [1:length(selinds)] + (-1)^iz*dx;
        
        mid  = allpts(:,selinds,iz);
        hilo = allhilo(:,selinds,iz);
        
        plot(xplt, mid, '.', 'markersize', ms, 'Color', cols{iz});
        errorbar(xplt, mid, hilo(1,:), hilo(2,:), 'linestyle', 'None', 'linewidth', lw, 'Color', cols{iz});
    end
    xlim([xplt(1)-0.5, xplt(end)+0.5]);
    yl = ylim; yl(1) = 0; ylim(yl);
    set(gca,'XTick',[1,2],'XTicklabel',{'Incidence','Notifications'},'fontsize',fs);
    ylabel('Rate per 100,000 population');
    
    subplot(1,2,2); hold on;
    selinds = [5];
    for iz = 1:2
        xplt = [1:length(selinds)] + (-1)^iz*dx;
        
        mid  = allpts(:,selinds,iz);
        hilo = allhilo(:,selinds,iz);
        
        lg(iz,:) = plot(xplt, mid, '.', 'markersize', ms, 'Color', cols{iz});
        errorbar(xplt, mid, hilo(1,:), hilo(2,:), 'linestyle', 'None', 'linewidth', lw, 'Color', cols{iz});
    end
    xlim([min(xplt)-0.5, max(xplt)+0.5]);
    yl = ylim; yl(1) = 0; ylim(yl);
    set(gca,'XTick',[1],'XTicklabel',{'Mortality'},'fontsize',fs);
    
    legend(lg,'Data','Model','Location','SouthEast');
    
else    
    
    subplot(1,2,1); hold on;

    selinds = [2,1,3]; cols = {'r','b'};
    for iz = 1:2
        xplt = [1:length(selinds)] + (-1)^iz*dx;
        
        mid  = allpts(:,selinds,iz);
        hilo = allhilo(:,selinds,iz);
        
        plot(xplt, mid, '.', 'markersize', ms, 'Color', cols{iz});
        errorbar(xplt, mid, hilo(1,:), hilo(2,:), 'linestyle', 'None', 'linewidth', lw, 'Color', cols{iz});
    end
    xlim([xplt(1)-0.5, xplt(end)+0.5]);
    yl = ylim; yl(1) = 0; ylim(yl);
    set(gca,'XTick',[1,2,3],'XTicklabel',{'Incidence 2014','Incidence 2019','Notifications'},'fontsize',fs);
    ylabel('Rate per 100,000 population');
    
    subplot(1,2,2); hold on;
    selinds = [5];
    for iz = 1:2
        xplt = [1:length(selinds)] + (-1)^iz*dx;
        
        mid  = allpts(:,selinds,iz);
        hilo = allhilo(:,selinds,iz);
        
        lg(iz,:) = plot(xplt, mid, '.', 'markersize', ms, 'Color', cols{iz});
        errorbar(xplt, mid, hilo(1,:), hilo(2,:), 'linestyle', 'None', 'linewidth', lw, 'Color', cols{iz});
    end
    xlim([min(xplt)-0.5, max(xplt)+0.5]);
    yl = ylim; yl(1) = 0; ylim(yl);
    set(gca,'XTick',[1],'XTicklabel',{'Mortality'},'fontsize',fs);
    
    legend(lg,'Data','Model','Location','SouthEast');
end

set(ff,'Position',[680   638   738   339]);
