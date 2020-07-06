function[pval,es] = test_mean_difference_cl(X)
Nstim = numel(X);
stim_ind = [];
for n = 1:Nstim
    Ntrial(n) = size(X{n},1);
    stim_ind = [stim_ind n*ones(1,Ntrial(n))];
end
X = vertcat(X{:});
%calculate distance between mean responses
D = test_mean_difference_cl_sub(X,stim_ind,Nstim);
%shuffle
rng(1);
Nshuf = 1000; 
for n = 1:Nshuf
    stim_ind = stim_ind(randperm(numel(stim_ind)));
    Dshuf(:,:,n) = test_mean_difference_cl_sub(X,stim_ind,Nstim);
end
%calculate p-value
pval = 1 - sum((D>Dshuf),3)/Nshuf;
pval(1,1) = NaN; pval(2,2) = NaN; pval(2,1) = NaN; 
es = (D-mean(Dshuf,3))./mean(Dshuf,3);
es(1,1) = NaN; es(2,2) = NaN; es(2,1) = NaN; 


function[D] = test_mean_difference_cl_sub(X,stim_ind,Nstim)
for n = 1:Nstim
    Xm(n,:) = nanmean(X(stim_ind==n,:));
end
for n = 1:Nstim
    for m = n+1:Nstim
        D(n,m) = norm(Xm(n,:)-Xm(m,:));
    end
end