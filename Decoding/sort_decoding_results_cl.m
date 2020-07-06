function[mean_out, std_out, pval_out] = sort_decoding_results_cl(which_stim,which_alg)
if strcmp(which_alg,'knn')
    load('results_stimulus_predict_knn_review.mat');
else
    load('results_stimulus_predict_rf_review.mat');
end
Nres = numel(res); Ntrial = 344;
%
correct = [];
ind_correct = [];
count = 0;
for n = 1:Nres
    if isequal(which_stim, res(n).indS)
        count = count+1;
        correct{count} = res(n).correct;
        %correct{count} = correct{count}(:,1:7,:);
    end
end
N = numel(correct);
%identify the best combination of parameters
for n = 1:N
    mean_correct = squeeze(mean(correct{n},1));
    [ind_best_K,ind_best_Npca] = find(mean_correct==max(mean_correct(:)),1);
    correct{n} = 100*correct{n}(:,ind_best_K,ind_best_Npca);
end
%output means & stds
for n = 1:N
    mean_out(n) = mean(correct{n});
    std_out(n) = std(correct{n});
end



