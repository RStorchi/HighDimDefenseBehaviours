function[group,Y] = cluster_movements_cl(is_prestimulus,Ngroup)
%kmeans++ clustering for ehavioural responses and a pre-stimulus epoch.
%INPUT: 
%is_prestimulus: type false to cluster 0-2s behvioural response; if true the clustering will be performed on 0.5s epoch preceding the stimulus onset.
%Ngroup: to set the number of clusters to extract from the kmeans analysis
%OUTPUT:
%group: returns cluster ID for each trial
%Y: returns the PC reduction for each trial. The number of components is set to explain >80% variance 



if is_prestimulus
    indX = 144:150;
else
    indX = 151:180;
end
%load data
which_ind = 1:9; Ngroup = 7;
load('../3Dreconstruction/data_3D_all','flag_val'); 
load('../3Dreconstruction/measures_postures_movements.mat');
%take good trials
X = X(flag_val(:,1) == 1,:,:);
flag_val = flag_val(flag_val(:,1) == 1,:);
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
[N,Nt,~] = size(X);
Ndim = numel(which_ind);
Nind = numel(indX);
merge = zeros(N, Nind*Ndim);
for n = 1:N
    for m = 1:Ndim
        merge(n,(m-1)*Nind+1:m*Nind) = squeeze(X(n,indX,which_ind(m)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CLUSTER RESP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, Y, vars] = pca(merge);
THpca = 0.80;
Npca = min(find((cumsum(vars)/sum(vars))>THpca));
Y = Y(:,1:Npca);


%%%%%%%%%%%%%%%%%%%%CLUSTER%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ngroup_max = 10; 
Nrep = 100;  
rng(1);
C = []; D = [];
for n = 1:Ngroup_max
    [group(:,n),C{n},~,temp] = kmeans(Y, n, 'Replicates',Nrep,'Distance','sqeuclidean');
    D{n} = min(temp');
end
group = group(:,Ngroup);
C = C{Ngroup}; D = D{Ngroup};

%%%%%%%%%%%%%%%%%%%%REORDER%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gval = unique(group);
group_size = hist(group,gval);
[~,ind_size] = sort(group_size,'descend');
group_ord = group;
for n = 1:Ngroup
    group_ord(group==gval(ind_size(n))) = n;
end
gval = unique(group_ord);


%%%%%%%%%%%%%%%%%%%FIGURES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
merge_ordered = [];
group_ordered = [];
Yordered = [];
ind_stim_vec_ord = [];
ind_ordered = [];
for n = 1:numel(gval)
    ind_temp = find(group_ord == gval(n));
    ind_ordered = [ind_ordered; ind_temp];
    merge_ordered = [merge_ordered; merge(ind_temp,:)];
    group_ordered = [group_ordered; group_ord(ind_temp)];
    ind_stim_vec_ord = [ind_stim_vec_ord; ind_stim_vec(ind_temp)]; 
    Yordered = [Yordered; Y(ind_temp,:)];
end


make_pcolor_fig(merge_ordered,group_ordered,gval,true,'Clustered',Nind,Ndim);
prob = make_hist_fig(group_ordered,gval,ind_stim_vec_ord);

%save('single_clustering_results_final','flag_val','group_ord');
%save('clustered_response','Y','group','C','D','merge','ind_stim_vec');

%%%%%%%%%%%%%%%%%%FIG FUNCS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[hg] = make_hist_fig(group,gval,ind_stim_vec)
which_stim = unique(ind_stim_vec);
Nstim = numel(unique(which_stim));
Ng = numel(gval);
for n = 1:Nstim
   hg(n,:) = hist(group(ind_stim_vec==which_stim(n)),gval);
end
for n = 1:Ng
    hg(:,n) = hg(:,n)/sum(hg(:,n));
end

hg_p = [[hg; zeros(1,Ng)] zeros(Nstim+1,1)];
%
fig1 = figure; 
set(fig1,'Position',[150 50 400 900]);
tt = {'Flash','Loom','Sound'};
for n = 1:3
    h1(n) = subplot(3,1,n); hold on;
    for m = 1:Ng
        bar(m,hg(n,m),'BarWidth',0.5,'FaceColor',ones(1,3)*0.66);
    end
    xlabel('#group','FontSize',14);
    title(tt{n});
    set(h1(n),'XTick',[1:Ng],'XTickLabel',[1:Ng],'FontSize',14);
    set(h1(n),'YTick',[0:0.25:1],'FontSize',14);
    ylabel('\DeltaP(S|group)','FontSize',14);
    xlim([0 Ng+1]); ylim([0 1]);
end


function[] = make_pcolor_fig(merge,group,gval,ordered,tt,Nind,Ndim)
N = size(merge,1);
fig = figure;
set(fig,'Position',[50 50 500 600]);
h = subplot(1,1,1); hold on;
if ordered
    p = pcolor(merge); set(p,'LineStyle','none');
    count = 0;
    for n = 1:numel(gval)
        count = count+sum(group==gval(n));
        line([1 Ndim*Nind],[count count],'Color','r','LineWidth',2);
        text(-20,count,num2str(n),'FontSize',12) 
    end
    for n = 1:Ndim
        line(Nind*[n n],[0 N],'Color','g','LineWidth',2);
    end
else
    p = pcolor(merge); set(p,'LineStyle','none');
end
xlim([1 size(merge,2)]);
ylim([1 size(merge,1)]);
title(tt);
set(h,'XTick',[],'YTick',[]);
