clear all;

ff=figure; hold on;
lw = 1.5; fs = 14; ts = 0.1; 

% Construct monthly labels
xlbl_mo = []; ct = 1;
mnths = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
years = {'19','20','21','22','23','24','25'};
for iy = 1:length(years)
   for im = 1:length(mnths)
       xlbl_mo{ct} = [mnths{im},' ',years{iy}];                            
       ct = ct+1;
   end
end

% Construct quarterly labels                                               
xlbl_qt = []; ct = 1;
qrtrs = 1:4;
years = {'19','20','21','22','23','24','25'};
for iy = 1:length(years)
   for iq = 1:length(qrtrs)
       xlbl_qt{ct} = ['Q', num2str(qrtrs(iq)),' ',years{iy}]; 
       ct = ct+1;
   end
end

cols = linspecer(2);

iso3s = {'AGO','BGD','BRA','CHN','IDN','KEN','MMR','PAK','PER','PHL','RUS','UGA','VNM','ZAF'};

for iis = 1:length(iso3s)
    fprintf('%0.5g ',iis);
    subplot(5,3,iis);
    
    if ~ismember(iso3s{iis},{'IND'})
    load([iso3s{iis},'/disruption_vector']);
    
    ms = 10;
    if monthly
        mat  = squeeze(noti_pct(:,:,1,:));
        dt   = 1/12;
        skip = 6;
        xlbl = xlbl_mo;
        ylbl = 'Monthly notifs';
    else
        mat  = squeeze(notq_pct(:,:,1,:));
        dt   = 1/4;
        skip = 2;
        xlbl = xlbl_qt;
        ylbl = 'Quarterly notifs';
    end
    
    for ii = 1:2
        plt = mat(:,:,ii);
        pl2(ii,:) = plot(plt(2,:),'Color',cols(ii,:),'linewidth',lw); hold on;
        jbfill(1:size(plt,2),plt(3,:),plt(1,:),cols(ii,:),'None',1,ts); hold on;
    end
    pl2(3,:) = plot(1/dt + [1:length(notif_rate)], notif_rate,'.-','Color','g','linewidth',lw,'markersize',ms);
    ylabel(ylbl);
    xinds = 1:skip:1/dt+length(notif_rate);
    set(gca,'fontsize',fs,'XTick',xinds,'XTickLabel',xlbl(xinds));
    xlim([1, 1/dt+length(find(~isnan(notif_rate)))]);
    yl = ylim; yl(1) = 0; ylim(yl);
    title(iso3s{iis});
    %legend(pl2([3,1,2],:),'Data','Modelled baseline','Modelled disruption','location','SouthWest');
    end
end
fprintf('\n');
set(ff,'Position',[680   179   786   798]);



