function [group,Ngroup,C,Y] = cluster_movements_multi_timestep_review_v3(indX, Ngroup)
%here use empirical cost criterion
%e.g.:
%indX = {137:143, 144:150, 151:157, 158:164, 165:171, 172:178, 179:185};
%Ngroup = 5; 
%[LL,group,C,Y] = cluster_movements_multi_timestep_review_v1(indX, Ngroup);
%C are centroids?
%load data
rng(2);
which_ind = 1:9; 
Nts = numel(indX);
load('movements_quantified_final.mat');
%take good trials
X = X(flag_val(:,1) == 1,:,:);
flag_val = flag_val(flag_val(:,1) == 1,:);
N = size(X,1);
%reduce stimuli to flash/loom/sound
ind_stim_vec0 = flag_val(:,5);
ind_stim_vec = ind_stim_vec0;
ind_stim_vec((ind_stim_vec0==1)|(ind_stim_vec0==2)) = 1;
ind_stim_vec((ind_stim_vec0==3)|(ind_stim_vec0==4)) = 2;
Nstim = 3;
which_stim = unique(ind_stim_vec);
for n = 1:Nstim
    ind_stim{n} = find(ind_stim_vec==which_stim(n));
end
%generate merge matrix
merge = []; merge2 = [];
for n = 1:Nts
    temp = cluster_movements_multi_timestep_review_v3_sub1(X,indX{n},which_ind);
    merge(:,:,n) = temp;
    merge2 = [merge2; temp];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CLUSTER RESP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, Y2, vars] = pca(merge2);
THpca = 0.8;
Npca = min(find((cumsum(vars)/sum(vars))>THpca));
%Npca = 2;
Y2 = Y2(:,1:Npca);


%%%%%%%%%%%%%%%%%%%%CLUSTER%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nrep = 200;  
[group2,C] = kmeans(Y2, Ngroup, 'Replicates',Nrep,'Distance','sqeuclidean');

%%%%%%%%%%%%%%%%%%%%RESHAPE GROUP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
group = zeros(N,Nts);
Y = zeros(N,Npca,Nts);
for n = 1:Nts
    group(:,n) = group2((n-1)*N+1:N*n);
    Y(:,:,n) = Y2((n-1)*N+1:N*n,:);
end

%%%%%%%%%%%%%%%%%%%SAVE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('clustered for %s clusters',num2str(Ngroup)));
%save('clustered_response','Y','group','C','D','merge','ind_stim_vec');


%%%%%%%%%%%%%%%%%%%SUBFUN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[merge] = cluster_movements_multi_timestep_review_v3_sub1(X,indX,which_ind)
[N,Nt,~] = size(X);
Ndim = numel(which_ind);
Nind = numel(indX);
merge = zeros(N, Nind*Ndim);
for n = 1:N
    for m = 1:Ndim
        merge(n,(m-1)*Nind+1:m*Nind) = squeeze(X(n,indX,which_ind(m)));
    end
end

