% Script to align model uncertainty with WHO estimates. Because the model
% makes mechanistic links between mortality and incidence, and WHO
% estimates do not, model-based uncertainty is in general narrower than in WHO
% estimates. This code ensures alignment between the two, by (i) estimating
% incidence in year t>2019 as a function F of incidence in 2019, and all
% model parameters. (ii) Drawing samples from an uncertainty
% distribution consistent with WHO estimates in 2019, estimating incidence
% in year t>2019 using the function F. (iii) Doing the same for mortality.

% Script is essentially the same as Extrapolate_Countries, but applied to
% Regional projections instead

% Dependencies:
% =============
% - get_distribution_fns.m

% This code is used by:
% =====================
% - Batch_Simulate_Extrapolate

clear all;

reg = 'EUR'; dir = ['region_',reg,'/'];
load([dir,'res_forward']);


% -------------------------------------------------------------------------
% --- Get the incidence and mortality rate ratios -------------------------

% IRR for total incidence
IRR_all = (inct(:,:,1)./inct(1,:,1))';
% IRR for HIV+ve
tmp = squeeze(sum(inct2(:,[2,3],:,:),2));
IRR_hiv = (tmp(:,:,1)./tmp(1,:,1))';
% MRR for total mortality
MRR_all = (mort(:,:,1)./mort(1,:,1))';
% MRR for HIV+ve
tmp = squeeze(sum(mort2(:,2,:,:),2));
MRR_hiv = (tmp(:,:,1)./tmp(1,:,1))';
% Bring them all together
if opts.hiv
    RRs_bsline = cat(3, IRR_all, IRR_hiv, MRR_all, MRR_hiv);
else
    RRs_bsline = cat(3, IRR_all, MRR_all);
end

% IRR for total incidence
IRR_all = (inct(:,:,2)./inct(1,:,2))';
% IRR for HIV+ve
tmp = squeeze(sum(inct2(:,[2,3],:,:),2));
IRR_hiv = (tmp(:,:,2)./tmp(1,:,2))';
% MRR for total mortality
MRR_all = (mort(:,:,2)./mort(1,:,2))';
% MRR for HIV+ve
tmp = squeeze(sum(mort2(:,2,:,:),2));
MRR_hiv = (tmp(:,:,2)./tmp(1,:,2))';
% Bring them all together
if opts.hiv
    RRs_disrupt = cat(3, IRR_all, IRR_hiv, MRR_all, MRR_hiv);
else
    RRs_disrupt = cat(3, IRR_all, MRR_all);
end


load('Data/TB Estimates/estim_data_2b.mat');

if opts.hiv
    names = {'inc_2019','inc_h1','mort_all','mort_H1'};
else
    names = {'inc_2019','mort_all'};
end

tblinc = [];
tblmrt = [];


for ic = 1:length(ctrlist.(reg))
    iso3 = ctrlist.(reg){ic};
    
    % Get the data for incidence and mortality
    countryrow = estims(strcmp(estims.iso3,iso3),:);
    mat        = reshape(countryrow{1,3:end},[3,6])';
    
    % Make adjustments in case boundary estimates are same as central (e.g. for
    % MEX)
    inds = find(mat(:,1)==mat(:,2));
    mat(inds,1) = mat(inds,2)*0.95;
    
    inds = find(mat(:,2)==mat(:,3));
    mat(inds,3) = mat(inds,2)*1.05;
    
    data.inc_2019  = mat(1,:)/12;
    data.inc_h1    = mat(2,:)/12;
    data.mort_H1   = mat(4,:)/12;
    data.mort_all  = mat(5,:)/12;
    
    incmrt = [];
    
    for ni = 1:length(names)
        name = names{ni};
        
        [~, out] = get_distribution_fns(data.(name), 'lognorm', 0);
        cdif = 1; count = 0;
        while cdif>1e-3 && count<1e4   % <--- TEMPORARY: fix when finalising
            sam = lognrnd(out(1),out(2),size(RRs_bsline,1),1);
            vec = prctile(sam,[2.5,50,97.5])./data.(name)-1;
            cdif = sum(vec.^2); count = count+1;
            if count>1e4
                error(sprintf('Not providing good samples, iso3: %s',iso3));
            end
        end
        
        tmp1 = sam.*RRs_bsline(:,:,ni);
        tmp2 = sam.*RRs_disrupt(:,:,ni);
        tmp3 = cat(3,tmp1,tmp2);
        % Get annual aggregations
        dims = size(tmp3);
        tmp4 = reshape(tmp3,[dims(1), 12, dims(2)/12, dims(3)]);
        tmp5 = squeeze(sum(tmp4,2));
        % Get percentiles
        outmat_mo(:,:,:,ni,ic)  = prctile(tmp3,[2.5,50,97.5],1);           % 1.Lo/Md/Hi 2.Month 3.Baseline/COVID 4.Inc/Mrt by HIV 5.Country
        outmat_yr(:,:,:,ni,ic)  = prctile(tmp5,[2.5,50,97.5],1);           % 1.Lo/Md/Hi 2.Year 3.Baseline/COVID 4.Inc/Mrt by HIV 5.Country
                
    end
    
    % --- Save the projections under country names, to allow later
    % plotting: e.g. with Plot_projections_Regional
    incmrt = outmat_yr(:,:,:,:,ic);
    dimnames = {'1. Lo/Md/hi'; '2. Year'; '3. Baseline/COVID'; '4. Inc/Mort by HIV'};
    save([dir,'country projections/inc_mrt_',iso3,'.mat'],'incmrt','dimnames');
end

if opts.hiv
   incmat_mo = squeeze(cat(1,outmat_mo(:,:,:,1,:),outmat_mo(:,:,:,2,:)));  % Dims: 1.Lo/Md/Hi, then HIV 2.Year 3.Baseline/COVID 4.Country
   mrtmat_mo = squeeze(cat(1,outmat_mo(:,:,:,3,:),outmat_mo(:,:,:,4,:)));

   incmat_yr = squeeze(cat(1,outmat_yr(:,:,:,1,:),outmat_yr(:,:,:,2,:)));  % Dims: 1.Lo/Md/Hi, then HIV 2.Year 3.Baseline/COVID 4.Country
   mrtmat_yr = squeeze(cat(1,outmat_yr(:,:,:,3,:),outmat_yr(:,:,:,4,:)));
else
   incmat_mo = squeeze(outmat_mo(:,:,:,1,:));
   mrtmat_mo = squeeze(outmat_mo(:,:,:,2,:));
   
   incmat_yr = squeeze(outmat_yr(:,:,:,1,:));
   mrtmat_yr = squeeze(outmat_yr(:,:,:,2,:));
end


incout_mo = []; mrtout_mo = [];
for iz1 = 1:size(incmat_mo,4)
    for iz2 = 1:size(incmat_mo,3)
        incout_mo = [incout_mo; incmat_mo(:,:,iz2,iz1)];
        mrtout_mo = [mrtout_mo; mrtmat_mo(:,:,iz2,iz1)];
    end
end

incout_yr = []; mrtout_yr = [];
for iz1 = 1:size(incmat_yr,4)
    for iz2 = 1:size(incmat_yr,3)
        incout_yr = [incout_yr; incmat_yr(:,:,iz2,iz1)];
        mrtout_yr = [mrtout_yr; mrtmat_yr(:,:,iz2,iz1)];
    end
end

% --- Construct the row labels

mat = repmat({'xx'}, 3*2*length(names)/2,length(ctrlist.(reg)));
mat{1,1} = reg;
col0 = mat(:);

mat = repmat({'xx'}, 3*2*length(names)/2,length(ctrlist.(reg)));
mat(1,:) = ctrlist.(reg)';
col1 = mat(:);

mat = repmat({'xx'},3*2*length(names)/2,length(ctrlist.(reg)));
mat(1,:) = {'Baseline'};
mat(3*2*length(names)/4+1,:) = {'With COVID'};
col2 = mat(:);

mat = repmat({'xx'},3*2*length(names)/2,length(ctrlist.(reg)));
if opts.hiv
    rows1 = [1,7];
    rows2 = [4,10];
    mat(rows1,:) = {'All'};
    mat(rows2,:) = {'HIV +ve'};
end
col3 = mat(:);

vec = {'lo','md','hi'}';
col4 = repmat(vec,length(col3)/3,1);

inctbl_mo = [col0, col1, col2, col3, col4, num2cell(incout_mo)];
mrttbl_mo = [col0, col1, col2, col3, col4, num2cell(mrtout_mo)];
inctbl_yr = [col0, col1, col2, col3, col4, num2cell(incout_yr)];
mrttbl_yr = [col0, col1, col2, col3, col4, num2cell(mrtout_yr)];

cell2csv([dir,'/Incidence_projections_mo.csv'],inctbl_mo);
cell2csv([dir,'/Mortality_projections_mo.csv'],mrttbl_mo);
cell2csv([dir,'/Incidence_projections_yr.csv'],inctbl_yr);
cell2csv([dir,'/Mortality_projections_yr.csv'],mrttbl_yr);
