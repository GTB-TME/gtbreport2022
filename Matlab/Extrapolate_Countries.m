% Script to align model uncertainty with WHO estimates. Because the model
% makes mechanistic links between mortality and incidence, and WHO
% estimates do not, model-based uncertainty is in general narrower than in WHO
% estimates. This code ensures alignment between the two, by (i) estimating
% incidence in year t>2019 as a function F of incidence in 2019, and all
% model parameters. (ii) Drawing samples from an uncertainty
% distribution consistent with WHO estimates in 2019, estimating incidence
% in year t>2019 using the function F. (iii) Doing the same for mortality.

% Dependencies:
% =============
% - get_distribution_fns.m

% This code is used by:
% =====================
% - Batch_Simulate_Extrapolate

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% Set looping = 0 if this script is not being run by Batch_Simulate_Extrapolate, 
% and instead is being used only for a single country
looping  = 0;
if ~looping
    clear all; ctry = 'MNG';
end

% Load forward projections, incorporating COVID disruptions, created by 
% 'Simulate_disruptions' or 'Simulate_disruptions_AUTOMATIC'
load([ctry,'/projections_raw']);

% Whether the code will plot the outcomes of extrapolation,
% comparing with WHO estimates
plotting = 1;

% Whether the code will save results of the extrapolation, as annual
% timeseries for incidence and mortality
printing = 1;


% --- Get the data
alldat = [data.inc_2019; data.inc_h1; data.mort_all; data.mort_H1]/12;

% --- Collate the model simulations
tmp1 = squeeze(sum(inct2(:,[2,3],:,:),2));
tmp2 = squeeze(mort2(:,2,:,:));
allsim = cat(4,inct,tmp1,mort,tmp2);

% For countries where HIV doesn't play a substantial role, take only the
% HIV-ve outputs
if ~opts.hiv
    alldat = alldat([1,3],:);
    allsim = allsim(:,:,:,[1,3]);
end

% For DPRK, exclude HIV+ve data, which is unavailable in WHO estimates
if strcmp(ctry,'PRK')
    alldat = [data.inc_2019; data.mort_H0]/12;
end


% -------------------------------------------------------------------------
% --- Fit a log-normal distribution to WHO 2019 estimates for incidence and
% --- mortality, and sample from that distribution

for iz = 1:size(alldat,1)
    dat = alldat(iz,:);
    [~, out] = get_distribution_fns(dat, 'lognorm', 0);
    cdif = 1;
    while cdif>1e-3
        sam = lognrnd(out(1),out(2),size(xs,1),1);
        vec = prctile(sam,[2.5,50,97.5])./dat-1;
        cdif = sum(vec.^2);
    end
    fprintf('Agreement with WHO estimates:\n'); prctile(sam,[2.5,50,97.5])./dat
    sam_who(:,iz) = sam;
end


% -------------------------------------------------------------------------
% --- Using regression, estimate the function relating incidence (or
% --- mortality) in year t>2019 to model parameters, and its value in 2019

for iz = 1:size(allsim,4)                                                  % Incidence and mortality
    for ii = 1:size(allsim,3)                                              % Baseline vs disruption scenario
        mat = allsim(:,:,ii,iz)'*1e5;
        X = [mat(:,1), xs];
        
        beta = [];
        for ic = 2:size(mat,2)
            Y = mat(:,ic);
            beta(:,ic) = mvregress(X,Y);
        end
        beta(1,1) = 1;
        
        % Make extrapolations
        X = [sam_who(:,iz), xs];
        
        mat1 = X*beta;
        timeser(:,:,ii,iz) = mat1;
        all_pct_mo(:,:,ii,iz) = prctile(mat1,[2.5,50,97.5],1);             % Dims: 1.Lo/Md/Hi, 2.Month, 3.Baseline/disruption, 4.Inc/Mort/HIV data
    
        % Also get annualised rates
        dims = size(mat1);
        tmp  = reshape(mat1,[dims(1), 12, dims(2)/12]);
        mat2 = squeeze(sum(tmp,2));
        all_pct_yr(:,:,ii,iz) = prctile(mat2,[2.5,50,97.5],1);             % Dims: 1.Lo/Md/Hi, 2.Year, 3.Baseline/disruption, 4.Inc/Mort/HIV data
    end
    betasto(:,:,iz) = beta;
end


% -------------------------------------------------------------------------
% --- Plot the results to compare alignment with WHO uncertainty intervals 

if plotting
    
    figure; hold on;
    tp = 0.1; fs = 14; lw = 1.5;
    cols = linspecer(3);
    xpts = (1:size(all_pct_mo,2))+1;
    
    if ~opts.hiv
        nr = 1; nc = 2; inds = [1,2];
    else
        nr = 2; nc = 2; inds = [1,3,2,4];
    end
    
    for iz = 1:size(all_pct_mo,4)
        subplot(nr,nc,inds(iz));
        for ii = 1:size(all_pct_mo,3)
            plt = all_pct_mo(:,:,ii,iz);
            lg(ii,:) = plot(xpts, plt(2,:),'Color',cols(ii,:),'linewidth',lw); hold on;
            jbfill(xpts,plt(3,:),plt(1,:),cols(ii,:),'None',1,tp); hold on;
        end
        dat = alldat(iz,:);
        difdat = diff(dat);
        lg(ii+1,:) = errorbar(1,dat(2),difdat(1),difdat(2),'linewidth',lw,'Color',cols(ii+1,:));
        plot(1,dat(2),'.','markersize',24,'Color',cols(ii+1,:))
        
        yl = ylim; yl(1) = 0; ylim(yl);
        xl = xlim; xl(2) = xpts(end)+10; xlim(xl);
        
        xlabel('Months from Jan 2019');
        ylabel('Monthly rate per 100k population');
        set(gca,'fontsize',fs);
        if iz == 1
            legend(lg,'Model posterior-based outputs','Statistically extrapolated','WHO uncertainty ranges, 2019','Location','SouthEast');
        end
    end
end


% -------------------------------------------------------------------------
% --- Print the results to csv files in the appropriate country folder ----

if printing
    
    sources = {all_pct_mo, all_pct_yr};
    tis = {'mo','yr'};
    
    for si = 1:length(sources)
        
        source = sources{si};
        
        % Stack baseline and disruption scenarios
        incout = [];
        mrtout = [];
        if ~opts.hiv
            incout = [source(:,:,1,1); source(:,:,2,1)];
            mrtout = [source(:,:,1,2); source(:,:,2,2)];
            
            % Construct column names
            cols = repmat({'xx'},6,4);
            cols(1,1) = {ctry};
            
            cols(1,2) = {'Baseline'};
            cols(4,2) = {'COVID'};
            
            cols([1,4],4) = {'lo'};
            cols([2,5],4) = {'md'};
            cols([3,6],4) = {'hi'};
            
            inctbl = [cols, num2cell(incout)]; 
            mrttbl = [cols, num2cell(mrtout)]; 
            
            
        else
            incout = [source(:,:,1,1); source(:,:,1,2); source(:,:,2,1); source(:,:,2,2)];
            mrtout = [source(:,:,1,3); source(:,:,1,4); source(:,:,2,3); source(:,:,2,4)];
            
            % Construct column names
            cols = repmat({'xx'},12,4);
            cols(1,1) = {ctry};
            
            cols(1,2) = {'Baseline'};
            cols(7,2) = {'COVID'};
            
            cols([1,7],3) = {'All'};
            cols([4,10],3) = {'HIV +ve'};
            
            cols([1,4,7,10],4) = {'lo'};
            cols([2,5,8,11],4) = {'md'};
            cols([3,6,9,12],4) = {'hi'};
            
            inctbl = [cols, num2cell(incout)]; 
            mrttbl = [cols, num2cell(mrtout)];
        end
        cell2csv([iso3,'/Incidence_projections_',tis{si},'.csv'],inctbl);
        cell2csv([iso3,'/Mortality_projections_',tis{si},'.csv'],mrttbl);
        
    end
    
    save([ctry,'/projections_extrapolated.mat'], 'all_pct_yr', 'opts');
end
