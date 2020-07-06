function[SI,SI_loc] = calculate_SI_cl(which_stim,Kneigh,epoch,make_fig)
%calculate SI 
%INPUT: 
%which_stim: for flash vs loom = [0 1]; for flash vs sound = [0 2]; for loom vs sound = [1 2];
%Kneigh: number of neighbours; 
%epoch: for 0-1s type 'early'; for 1-2s type 'intermediate';for 2-3s type 'late'; 0-2s type
%'all;
%OUTPUT: 

%load data
load('../3Dreconstruction/measures_postures_movements_10bin.mat','X'); 
load('../3Dreconstruction/data_3D_all.mat','flag_val');

%epoch selection
if strcmp(epoch,'early')
    indX = 151:165;   
elseif strcmp(epoch,'intermediate')
    indX = 166:180; 
elseif strcmp(epoch,'late') 
    indX = 181:195; 
elseif strcmp(epoch,'all')
    indX = 151:180; 
else
   disp('epoch not recognized, see help');      
end

%select good trials
ind_good = find(flag_val(:,1) == 1);
X = X(ind_good,:,:);
flag_val = flag_val(ind_good,:);
%
ind_stim_vec0 = flag_val(:,5);
ind_stim_vec = ind_stim_vec0;
ind_stim_vec((ind_stim_vec0==1)|(ind_stim_vec0==2)) = 1;
ind_stim_vec((ind_stim_vec0==3)|(ind_stim_vec0==4)) = 2;
%
[N,~,~] = size(X);
Nind = numel(indX);
which_dim = [1:9];
Ndim = numel(which_dim);
%
merge = zeros(N, Nind*Ndim);
merge_locomotion = zeros(N, Nind);
for n = 1:N
    merge_locomotion(n,:) = squeeze(X(n,indX,4));
    for m = 1:Ndim
        merge(n,(m-1)*Nind+1:m*Nind) = squeeze(X(n,indX,m));
    end
end
%
if numel(which_stim) == 2
    ind_stim = find((ind_stim_vec==which_stim(1))|(ind_stim_vec==which_stim(2)));
    merge = merge(ind_stim,:);
    X = X(ind_stim,:,:);
    merge_locomotion = merge_locomotion(ind_stim,:);
    ind_stim_vec = ind_stim_vec(ind_stim);
end
[N,Nt] = size(merge);

%PCA reduction
[~, score,vars] = pca(merge);
[~, score_loc,vars_loc] = pca(merge_locomotion);
if Nind == 15
    Npca = 15; Npca_loc = 15;
else
    Npca = 45; Npca_loc = 30;
end

D = score(:,1:Npca); 
Dlocomotion = score_loc(:,1:Npca_loc);

%calculate SI
SI = calculate_SI_cl_sub(D,ind_stim_vec, Kneigh);
SI_loc = calculate_SI_cl_sub(Dlocomotion,ind_stim_vec, Kneigh);

if make_fig
    figure; 
    subplot(1,1,1); hold on;
    plot(1:Npca, nanmean(SI),'.-','LineWidth',2,'MarkerSize',18)
    plot(1:Npca_loc, nanmean(SI_loc),'k.-','LineWidth',2,'MarkerSize',18)
    xlabel('#PC'); ylabel('SI'); 
    legend('Full Set','Locomtion');
end

function[spec] = calculate_SI_cl_sub(D,ind_stim_vec, Kneigh)
[N, Npca] = size(D);
spec = zeros(N,Npca);
Dval = pdist(D); Dval = Dval(Dval>0); Dmin = min(Dval);
for n = 1:N
    for m = 1:Npca
        ind_neigh = knnsearch(D(:,1:m), D(n,1:m),'K',Kneigh+1,'NSMethod','kdtree');
        ind_neigh = ind_neigh(2:end);
        same_stim = ind_stim_vec(ind_neigh)==ind_stim_vec(n);
        [~, dist] = knnsearch(D(n,1:m), D(ind_neigh,1:m),'K',Kneigh,'NSMethod','kdtree');
        for nd = 1:numel(dist)
            if dist(nd) == 0
               dist(nd) = Dmin; 
            end
        end
        spec(n,m) = sum(squeeze(same_stim)./dist)/sum(1./dist);
    end
end
