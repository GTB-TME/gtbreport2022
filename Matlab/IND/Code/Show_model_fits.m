clear all; load model_fits2;

ix0 = 2e4; nx = 150; dx = round((size(xsto,1)-ix0)/nx);
xs = xsto(ix0:dx:end,:,1);

% inds = find(outsto==max(outsto));
% xs = xsto(inds(1),:);

mk = round(size(xs,1)/25);
for ii = 1:size(xs,1)
    if mod(ii,mk) == 0; fprintf('%0.5g ',ii/mk); end 
    [out, aux] = obj(xs(ii,:));
    sims(ii,:) = [aux.incd_2015, aux.incd_2019, aux.noti, aux.mort];
    inct(:,ii) = diff(aux.soln(:,i.aux.inc),1)*1e5;
end
fprintf('\n');

inc_pct = prctile(inct,[2.5,50,97.5],2)';


% --- Show all on a plot --------------------------------------------------

ff=figure; fs = 14; lw = 1.5;


% -------------------------------------------------------------------------
% --- Incidence, Notification fits ----------------------------------------

subplot(1,2,1); hold on;
% Plot data
plt  = [data.inc_2015; data.inc_2019; data.noti]';
hilo = diff(plt,1); md = plt(2,:);
xpts = [1:3]-0.1;

plot(xpts, md, 'r.', 'markersize', 24);
errorbar(xpts, md, hilo(1,:), hilo(2,:), 'LineStyle', 'None', 'linewidth', lw, 'Color', 'r');

% Plot simulations
plt = prctile(sims(:,[1:3]),[2.5,50,97.5],1);
hilo = diff(plt,1); md = plt(2,:);
xpts = [1:3]+0.1;

plot(xpts, md, 'b.', 'markersize', 24);
errorbar(xpts, md, hilo(1,:), hilo(2,:), 'LineStyle', 'None', 'linewidth', lw, 'Color', 'b');
xlim([0.5 3.5]);
set(gca,'fontsize',fs,'XTick',1:3,'XTickLabel',{'Incidence 2015','Incidence 2019','Notifications'});
yl = ylim; yl(1) = 0; ylim(yl);


% -------------------------------------------------------------------------
% --- Mortality fits ------------------------------------------------------

subplot(1,2,2); hold on;
plt  = [data.mort_all]';
hilo = diff(plt,1); md = plt(2,:);
xpts = [1]-0.1;

lg(1,:) = plot(xpts, md, 'r.', 'markersize', 24);
errorbar(xpts, md, hilo(1,:), hilo(2,:), 'LineStyle', 'None', 'linewidth', lw, 'Color', 'r');

% Plot simulations
plt = prctile(sims(:,[5]),[2.5,50,97.5],1);
hilo = diff(plt,1); md = plt(2,:);
xpts = [1]+0.1;

lg(2,:) = plot(xpts, md, 'b.', 'markersize', 24);
errorbar(xpts, md, hilo(1,:), hilo(2,:), 'LineStyle', 'None', 'linewidth', lw, 'Color', 'b');
xlim([0.5 1.5]);
set(gca,'fontsize',fs,'XTick',1,'XTickLabel',{'Mortality 2019'});
yl = ylim; yl(1) = 0; ylim(yl);

subplot(1,2,2);
legend(lg,'Calibration targets','Updated model','location','SouthEast');

set(ff,'Position',[680   638   738   339]);
