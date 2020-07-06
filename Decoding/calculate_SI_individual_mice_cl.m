function[same_mouse,same_mouse_loc,p,ploc] = calculate_SI_individual_mice_cl()
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

Kneigh = 1;
epoch = 'all';
which_stim = [0 1 2];

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
which_mouse = flag_val(:,3); 
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

%calculate same mouse
Nsh = 1000; 
same_mouse = calculate_SI_individual_mice_cl_sub(D,which_mouse);
same_mouse_loc = calculate_SI_individual_mice_cl_sub(Dlocomotion,which_mouse);
for n = 1:Nsh
    same_mouse_sh(n) = calculate_SI_individual_mice_cl_sub(D,which_mouse(randperm(N)));
    same_mouse_loc_sh(n) = calculate_SI_individual_mice_cl_sub(Dlocomotion,which_mouse(randperm(N)));
end
p = 1-sum(same_mouse> same_mouse_sh)/Nsh;
ploc = 1-sum(same_mouse_loc> same_mouse_loc_sh)/Nsh;
same_mouse  = 100*same_mouse/N;
same_mouse_loc  = 100*same_mouse_loc/N;




function[same_mouse] = calculate_SI_individual_mice_cl_sub(D,which_mouse)
ind = knnsearch(D,D,'k',2);
ind = ind(:,2); 
same_mouse = sum(which_mouse==which_mouse(ind));
