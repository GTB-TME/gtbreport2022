clear all; load res_forward.mat;

% --- Get excess cases and deaths
allmat = cat(4,inct,mort);
outcms = {'cases','deaths'};
for ia = 1:2
    tmp1         = squeeze(sum(allmat(13:24,:,:,ia),1));                                  % For 2020 only
    tmp2         = squeeze(sum(allmat(13:end,:,:,ia),1));                                 % For 2020 - 2025
    mat          = cat(3,tmp1,tmp2);
    ex_inc       = squeeze(mat(:,2,:)./mat(:,1,:) - 1);
    mat1(:,:,ia) = prctile(ex_inc, [2.5, 50. 97.5], 1)*100;    
end

% --- Incidence and mortality timeseries
tmp  = cat(4,inc_pct,mrt_pct);
% Drop 2019 from projections
tmp(:,1:12,:,:) = [];
% Append baseline and disruption scenarios
mat2 = squeeze(cat(1,tmp(:,:,1,:),tmp(:,:,2,:)));

% Bring them together
tmp1 = nan(size(mat1)); 
tmp2 = cat(1,tmp1,mat1);
out  = [tmp2, mat2];

writematrix(out(:,:,1), 'Incd_results.csv');
writematrix(out(:,:,2), 'Mort_results.csv');