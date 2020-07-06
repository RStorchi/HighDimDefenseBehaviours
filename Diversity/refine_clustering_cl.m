function[group, ind_ok, ind_stim_vec] = refine_clustering_cl(pre_stimulus,graph,ind_ok1)
%This script improves goodness of clustering by removing trials whose cluster membership is most ambiguous.
%In practice we remove the trials that are furthest from their centroid.
%Goodness of clustering is measured with Davies-Bouldin (DB) and Silhouette
%(SL) indexes
%NOTE:
%DB: Lower value is better clustering
%Sil: Higher value if better clustering
%INPUT: 
%pre_stimulus: if true the clusters of pre-stimulus epoch are analysed, if
%false the analysis is performed on clusters associated with 0-2s epoch
%after stimulus onset
%graph: if true a graph of the clustering results are generated
%ind_ok1: if provided allows to restrict analysis to the trials that are flagged as true and removes those flagged as false
%OUTPUT: 
%group: provides the cluster ID of a restricted set of trials (50%, this can be changed in line 34)
%ind_ok: flags 1 when a trial is preserved, 0 when it is removed
%ind_stim_vec: vector with stimulus identity fr each trial; 0 = flash, 1 = loom; 2 = sound  

if pre_stimulus
    load('clustered_prestimulus.mat');
else
    load('clustered_response.mat');
end

Ndim = 9; 
Nts = size(merge,2)/Ndim;

X = Y;
N0 = size(X,1); K = size(C,1);
Nd = size(X,2);

%refine clustering
frac = [1:-0.05:0.5]; 
Nfrac = numel(frac);
ind_ok = true(1,N0);
N = N0;
for n = 1:Nfrac
    %trim dataset
    for m = 1:K
        ind = find(group==m);
        TH = quantile(D(ind),frac(n));
        ind_ok(ind) = D(ind)<=TH;
    end
    %evaluate clustering
    temp = evalclusters(X(ind_ok,:),group(ind_ok),'DaviesBouldin');
    DB(n) = temp.CriterionValues;   
    S = silhouette(X(ind_ok,:),group(ind_ok));
    Sq(n,:) = quantile(S,[0:0.1:1]);
    N(n) = sum(ind_ok);
end

%re-order group labels from large to small
gval = unique(group);
group_size = hist(group,gval);
[~,ind_size] = sort(group_size,'descend');
group_ord = group;
for n = 1:K
    group_ord(group==gval(ind_size(n))) = n;
end
gval = unique(group_ord);

%further restrain
if nargin == 3
    ind_ok = ind_ok & ind_ok1;
end

%re-order rows
group_ordered = [];
merge_ordered = [];
ind_ok_ordered = [];
for n = 1:numel(gval)
    ind_temp = find(group_ord == gval(n));
    group_ordered = [group_ordered; group_ord(ind_temp)];
    merge_ordered = [merge_ordered; merge(ind_temp,:)];
    ind_ok_ordered = [ind_ok_ordered ind_ok(ind_temp)];
end
count_all = hist(group_ordered,gval);
count_red = hist(group_ordered(ind_ok_ordered==true),gval);


%results
if graph
    %
    fig = figure; 
    set(fig,'Position',[100 100 600 250]);
    h1 = subplot(1,2,1); hold on;
    plot(N,DB,'.-','LineWidth',2,'MarkerSize',15); 
    ylabel('DB index','FontSize',16);
    xlabel('#Trials','FontSize',16); set(h1,'FontSize',14); 
    h2 = subplot(1,2,2); hold on;
    plot(N,Sq(:,6),'.-','LineWidth',2,'MarkerSize',15); 
    ylabel('SL index','FontSize',16);
    xlabel('#Trials','FontSize',16); set(h2,'FontSize',14); 
    %%%%
    figure; 
    subplot(1,2,1); hold on;
    p = pcolor(merge_ordered); set(p,'LineStyle','none');
    caxis([0 1]);
    for n = 1:Ndim
        line(n*Nts*[1 1]+1,[1 N(1)],'LineWidth',2,'Color','g');
    end
    for n = 1:numel(gval)
        line(Ndim*Nts*[0 1],sum(count_all(1:n))*[1 1]+1,'LineWidth',2,'Color','r');
    end
    xlim([1 Ndim*Nts]); ylim([1 N(1)]); 
    %%%%
    subplot(1,2,2); hold on;
    merge1 = merge_ordered(find(ind_ok_ordered),:); %merge1(~ind_ok_ordered,:) = 0;
    p = pcolor(merge1); set(p,'LineStyle','none');
    caxis([0 1]);
    for n = 1:Ndim
        line(n*Nts*[1 1]+1,[1 N(end)],'LineWidth',2,'Color','g');
    end
    for n = 1:numel(gval)
        line(Ndim*Nts*[0 1],sum(count_red(1:n))*[1 1]+1,'LineWidth',2,'Color','r');
    end
    xlim([1 Ndim*Nts]); ylim([1 sum(ind_ok)]); 
end