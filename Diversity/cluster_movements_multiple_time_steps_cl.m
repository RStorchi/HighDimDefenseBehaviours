function[group_ts,Y,Nts_pre,Nts_post] = cluster_movements_multiple_time_steps_cl(Ngroup,Ntime_steps)
%This function provides clustering of responses across different epochs of
%the response and across trials. The clustering results is then used to estimate Variable-order Markov Models. 
%INPUT: 
%Ngroup: to set the number of clusters
%OUTPUT:
%group_ts: matrix of cluster ID; in rows are trials, in columns are
%different epochs within a trial
%Y: PC reduction of each cluster (explains >80% variance). dim1 = trials;
%dim2 = scores on the principal components; dim3 = epochs within a trial
%Nts_pre: number of epochs before stimulus onset
%Nts_post: number of epochs after stimulus onset

if Ntime_steps == 2
    %2 time steps: 
    indXpre = {136:150};
    indXpost = {151:165, 166:180};
elseif Ntime_steps == 3
    %3 time steps: 
    indXpre = {141:150};
    indXpost = {151:160, 161:170, 171:180};
elseif Ntime_steps == 5
    %5 time steps: 
    indXpre = {139:144, 145:150};
    indXpost = {151:156, 157:162, 163:168, 169:174, 175:180};
elseif Ntime_steps == 6
    %6 time steps: 
    indXpre = {141:145, 146:150};
    indXpost = {151:155, 156:160, 161:165, 166:170, 171:175, 176:180};
elseif Ntime_steps == 10
    %10 time steps: 
    indXpre = {139:141,142:144, 145:147, 148:150};
    indXpost = {151:153, 154:156, 157:159, 160:162, 163:165, 166:168, 169:171, 172:174, 175:177, 178:180};
elseif Ntime_steps == 15
    %15 time steps: 
    indXpre = {141:142, 143:144, 145:146, 147:148, 149:150};
    indXpost = {151:152, 153:154, 155:156, 157:158, 159:160, 161:162, 163:164, 165:166, 167:168, 169:170, 171:172, 173:174, 175:176, 177:178, 179:180};
else
    disp('allowed Ntime_steps are: 2, 3, 5, 6, 10, 15');
    group_ts = [];
    Y = [];
    return;
end
%
Nts_post = numel(indXpost);
Nts_pre = numel(indXpre);
Nts = Nts_pre+Nts_post;
%load data
which_ind = 1:9; 
load('../3Dreconstruction/data_3D_all','flag_val'); 
load('../3Dreconstruction/measures_postures_movements.mat');
%select good trials
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
%generate merge matrix for pre-stimulus
merge = []; merge2 = [];
count = 0;
for n = 1:Nts_pre
    count = count+1;
    temp = cluster_movements_multi_timestep_review_v3_sub1(X,indXpre{n},which_ind);
    merge(:,:,count) = temp;
    merge2 = [merge2; temp];
end
%generate a merge matrix for response
for n = 1:Nts_post
    count = count+1;
    temp = cluster_movements_multi_timestep_review_v3_sub1(X,indXpost{n},which_ind);
    merge(:,:,count) = temp;
    merge2 = [merge2; temp];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CLUSTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, Y2, vars] = pca(merge2);
THpca = 0.8;
Npca = min(find((cumsum(vars)/sum(vars))>THpca));
Y2 = Y2(:,1:Npca);
Nrep = 100;  
rng(1);
group2_ts = kmeans(Y2, Ngroup, 'Replicates',Nrep,'Distance','sqeuclidean');

%%%%%%%%%%%%%%%%%%%%RE-ORDER CLUSTERS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
group_ts = zeros(N,Nts);
Y = zeros(N,Npca,Nts);
for n = 1:Nts
    group_ts(:,n) = group2_ts((n-1)*N+1:N*n);
    Y(:,:,n) = Y2((n-1)*N+1:N*n,:);
end

%%%%%%%%%%%%%%%%%%%SAVE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('clustered for %s clusters',num2str(Ngroup)));
% filename = ['review_antonio_multistep_cluster' num2str(Ngroup)]; 
% save(filename,'Y','group_ts','Nts_pre','Nts_post','merge','ind_stim_vec');


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

